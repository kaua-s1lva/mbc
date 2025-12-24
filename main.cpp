#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include "header.h"

#define DBG

int main(int argc, char* argv[]) {
	double arg_TI = 798.27;
	double arg_TC = 0.12;
	double arg_TR = 0.8656;
	int arg_SAMAX = 151; // Se 0, vamos definir baseado no num_nos depois
	double arg_TEM_MAX = 300;
	int seed = 123456;
	std::string instance_path = "instancias/SW-1000-4-0d1-trial1.txt";

	// 2. Parsing da Linha de Comando
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--ti" && i + 1 < argc) arg_TI = std::atof(argv[++i]);
		else if (arg == "--tc" && i + 1 < argc) arg_TC = std::atof(argv[++i]);
		else if (arg == "--tr" && i + 1 < argc) arg_TR = std::atof(argv[++i]);
		else if (arg == "--samax" && i + 1 < argc) arg_SAMAX = std::atoi(argv[++i]);
		else if (arg == "--seed" && i + 1 < argc) seed = std::atoi(argv[++i]);
		else if ((arg == "--instancia" || arg == "-i") && i + 1 < argc) instance_path = argv[++i];
	}

	// 3. Configurar Semente e Leitura
	srand(seed); // Importante para o irace!

	// Converter string para char* para sua funcao ler_arquivo
	ler_arquivo(instance_path.c_str());

	// Se SAMAX não foi passado (ou é 0), usa a lógica original (num_nos)
	if (arg_SAMAX == 0) arg_SAMAX = num_nos;

	Solucao sol;

	// Construção Inicial
	heu_const_ale(sol);
	calcular_fo(sol);

	// Variáveis de Saída (dummy)
	double TEM_TOT, TEM_MEL;
	int NUM_SOL;

	// 4. Executar o SA com os parametros lidos
	simulated_annealing(sol, arg_TI, arg_TC, arg_TR, arg_SAMAX, arg_TEM_MAX, TEM_TOT, TEM_MEL, NUM_SOL);

	// 5. Saída para o irace (APENAS O NÚMERO DA FO)
	// O irace minimiza por padrão. Como sua FO parece ser minimizar, está ok.
	// Use printf ou cout, mas garanta que seja a única saída numérica na última linha.
	escrever_solucao(sol);
	std::cout << sol.fo << " " << TEM_TOT << std::endl;

	return 0;
}

void simulated_annealing(Solucao& s, const double& TI, const double& TC,
	const double& TR, const int& SAMAX, const double& TEM_MAX,
	double& TEM_TOT, double& TEM_MEL, int& NUM_SOL)
{
	Solucao s_viz, s_atu;
	clock_t h = clock();
	TEM_TOT = TEM_MEL = 0.0;
	NUM_SOL = 1;
	memcpy(&s_atu, &s, sizeof(Solucao));
	while (true)
	{
		double temperatura = TI;
		while (temperatura > TC)
		{
			for (int i = 0; i < SAMAX; i++)
			{
				memcpy(&s_viz, &s_atu, sizeof(Solucao));
				gerar_vizinho(s_viz);
				//heu_BL_rand(s_viz, 1 * (num_moc + 1) * num_obj);
				//calcular_fo_solucao(s_viz);
				calcular_fo(s_viz);
				NUM_SOL++;
				double delta = s_atu.fo - s_viz.fo;
				if (delta > 0)
				{
					memcpy(&s_atu, &s_viz, sizeof(Solucao));
					if (s_viz.fo < s.fo)
					{
						memcpy(&s, &s_viz, sizeof(Solucao));
						TEM_MEL = (double)(clock() - h) / CLOCKS_PER_SEC;
#ifdef DBG
						printf("FO: %d\tTempo: %.2f\n", s.fo, TEM_MEL);
#endif
					}
				}
				else
				{
					double x = rand() % 1001;
					x /= 1000.0;
					if (x < exp(-delta / temperatura))
						memcpy(&s_atu, &s_viz, sizeof(Solucao));
				}
				TEM_TOT = (double)(clock() - h) / CLOCKS_PER_SEC;
				if (TEM_TOT > TEM_MAX)
					goto FIM;
			}
			temperatura *= TR;
		}
	}
FIM:;
}

void heu_const_ale_gul(Solucao& sol) {
	int vet_aux[MAX_NOS];

	memcpy(&vet_aux, &vet_qtd_rel, sizeof(vet_qtd_rel));

	int pivo = 0;
	for (int i = 0; i < NOS_FONTE; i++) {
		for (int j = 0; j <= num_nos; j++) {
			if (vet_aux[j] > vet_aux[pivo]) {
				pivo = j;
			}
		}
		vet_aux[pivo] = -1;
		sol.vet_sol[i] = pivo;
	}

	for (int i = 0; i < num_nos; i++) {
		vet_aux[i] = i + 1;
	}

	for (int i = 0; i < NOS_FONTE; i++) {
		for (int j = sol.vet_sol[i] - 1; j < num_nos; j++) {
			vet_aux[j] = vet_aux[j+1];
		}
	}

	for (int i = NOS_FONTE; i < num_nos; i++) {
		int pos = i + rand() % (num_nos - i);
		sol.vet_sol[i] = vet_aux[pos];
		int aux = vet_aux[i];
		vet_aux[i] = vet_aux[pos];
		vet_aux[pos] = aux;
	}
}

void heu_const_gul(Solucao& sol) {
	//construir o vetor solução baseado na ordem do vetor_aux, colocar os nós com maior número de arestas primeiro
	int vet_aux[MAX_NOS];

	memcpy(&vet_aux, &vet_qtd_rel, sizeof(vet_qtd_rel));

	int pivo = 0;
	for (int i = 0; i < num_nos; i++) {
		for (int j = 0; j <= num_nos; j++) {
			if (vet_aux[j] > vet_aux[pivo]) {
				pivo = j;
			}
		}
		vet_aux[pivo] = -1;
		sol.vet_sol[i] = pivo;
	}
}

void heu_const_ale(Solucao& sol) {
	int vet_aux[MAX_NOS];

	for (int i = 0; i < num_nos; i++) {
		vet_aux[i] = i+1;
	}

	for (int i = 0; i < num_nos; i++) {
		int pos = i + rand() % (num_nos - i);
		sol.vet_sol[i] = vet_aux[pos];
		int aux = vet_aux[i];
		vet_aux[i] = vet_aux[pos];
		vet_aux[pos] = aux;
	}
}

//trocar NÓ FONTE por NÓ COMUM
void gerar_vizinho(Solucao& sol) {
	int pos_origem = rand() % NOS_FONTE;
	int pos_destino = (rand() % num_nos - NOS_FONTE) + NOS_FONTE;

	int aux = sol.vet_sol[pos_origem];
	sol.vet_sol[pos_origem] = sol.vet_sol[pos_destino];
	sol.vet_sol[pos_destino] = aux;
}

//trocar QUALQUER NÓ
void gerar_vizinho2(Solucao& sol) {
	int pos_origem = rand() % num_nos;
	int pos_destino = rand() % num_nos;

	while (pos_destino == pos_origem) {
		pos_destino = rand() % num_nos;
	}

	int aux = sol.vet_sol[pos_origem];
	sol.vet_sol[pos_origem] = sol.vet_sol[pos_destino];
	sol.vet_sol[pos_destino] = aux;
}

void calcular_fo(Solucao& sol) {
	sol.fo = 0;

	int order[MAX_NOS];
	int vet_aux[MAX_NOS];

	for (int i = 0; i < NOS_FONTE; i++) {
		order[i] = sol.vet_sol[i];
	}

	for (int i = NOS_FONTE; i < num_nos; i++) {
		vet_aux[i - NOS_FONTE] = sol.vet_sol[i];
	}

	int tam_order = NOS_FONTE;
	int aux;

	while (tam_order < num_nos) {
		sol.fo++;
		aux = tam_order;
		for (int i = 0; i < aux; i++) {
			for (int j = 0; j < num_nos - NOS_FONTE; j++) {
				if (mat_bin_rel[order[i]][vet_aux[j]]) {
					order[tam_order] = vet_aux[j];
					vet_aux[j] = 0;
					tam_order++;
					break;
				}
			}
		}
	}
}

void ler_arquivo(const char* path) {
	int a, b;
	FILE* f = fopen(path, "r");

	if (f == NULL) {
		fprintf(stderr, "Erro: Nao foi possivel abrir o arquivo %s\n", path);
		exit(1);
	}

	memset(&vet_qtd_rel, 0, sizeof(vet_qtd_rel));
	memset(&mat_bin_rel, 0, sizeof(mat_bin_rel));

	fscanf(f, "%d %d %d", &num_nos, &num_arestas, &num_arestas);

	for (int i = 0; i < num_arestas; i++) {
		fscanf(f, "%d %d", &a, &b);
		mat_bin_rel[a][b] = 1;
		vet_qtd_rel[a]++;
		mat_bin_rel[b][a] = 1;
		vet_qtd_rel[b]++;
	}

	fclose(f);
}

void escrever_dados(const char* arq) {
	FILE* f;

	if (!strcmp(arq, " ")) f = stdout;
	else f = fopen(arq, "w");

	fprintf(f, "%d %d %d\n", num_nos, NOS_FONTE, num_arestas);

	for (int i = 1; i <= num_nos; i++) {
		for (int j = 1; j <= num_nos; j++) {
			fprintf(f, "%d ", mat_bin_rel[i][j]);
		}
		fprintf(f, "\n");
	}
}

void escrever_solucao(Solucao& sol) {
	printf("\nValor da fo: %d", sol.fo);

	printf("\nNos fonte:");
	for (int i = 0; i < NOS_FONTE; i++) {
		printf(" %d ", sol.vet_sol[i]);
	}

	printf("\nNos:");
	for (int i = NOS_FONTE; i < num_nos; i++) {
		printf(" %d ", sol.vet_sol[i]);
	}
}