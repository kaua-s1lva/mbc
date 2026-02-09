#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include "header.h"

#define DBG
//#define IRACE_MODE
//#define DEBUG_MODE
#define RELEASE_MODE

int main(int argc, char* argv[]) {

#ifdef IRACE_MODE

	double arg_TI = 1418.49;
	double arg_TC = 0.79;
	double arg_TR = 0.8976;
	int arg_SAMAX = 348;
	double arg_TEM_MAX = 300;

	double TEM_TOT, TEM_MEL;
	int NUM_SOL;

	int seed = 123456;
	std::string instance_path = "instancias/SW-100-6-0d2-trial1.txt";

	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--ti" && i + 1 < argc) arg_TI = std::atof(argv[++i]);
		else if (arg == "--tc" && i + 1 < argc) arg_TC = std::atof(argv[++i]);
		else if (arg == "--tr" && i + 1 < argc) arg_TR = std::atof(argv[++i]);
		else if (arg == "--samax" && i + 1 < argc) arg_SAMAX = std::atoi(argv[++i]);
		else if (arg == "--seed" && i + 1 < argc) seed = std::atoi(argv[++i]);
		else if ((arg == "--instancia" || arg == "-i") && i + 1 < argc) instance_path = argv[++i];
	}

	srand(seed);

	ler_arquivo(instance_path.c_str());

	if (arg_SAMAX == 0) arg_SAMAX = num_nos;

	Solucao sol;

	heu_const_ale(sol);

	calcular_fo(sol);

	simulated_annealing(sol, arg_TI, arg_TC, arg_TR, arg_SAMAX, arg_TEM_MAX, TEM_TOT, TEM_MEL, NUM_SOL);

	std::cout << sol.fo << " " << TEM_TOT << std::endl;
#endif


#ifdef DEBUG_MODE
	double arg_TI = 1418.49;
	double arg_TC = 0.79;
	double arg_TR = 0.8976;
	int arg_SAMAX = 348;
	double arg_TEM_MAX = 300;
	double TEM_TOT, TEM_MEL;
	int NUM_SOL;
	char instance[] = "instancias/SW-100-6-0d2-trial1.txt";

	srand(time(NULL));

	ler_arquivo(instance);

	Solucao sol;

	heu_const_ale(sol);

	calcular_fo(sol);

	simulated_annealing(sol, arg_TI, arg_TC, arg_TR, arg_SAMAX, arg_TEM_MAX, TEM_TOT, TEM_MEL, NUM_SOL);
#endif

#ifdef RELEASE_MODE
	double arg_TI = 1418.49;
	double arg_TC = 0.79;
	double arg_TR = 0.8976;
	int arg_SAMAX = 348;
	double arg_TEM_MAX = 300;
	double TEM_TOT, TEM_MEL;
	int NUM_SOL;
	int total_exec = 3;
	char instances[][50] = { "H10_30", "H9_30", "SW-100-4-0d1-trial1", "SW-100-5-0d1-trial1", "SW-100-6-0d1-trial3", "SW-100-6-0d2-trial1" };

	FILE* f = fopen("resultados.csv", "w");

	fprintf(f, "Instância;Melhor FO;FO Média;Desvio (%%);Tempo Médio (seg.);T. Melhor (seg.)\n");

	srand(time(NULL));

	for (int i = 0; i < (sizeof(instances) / sizeof(instances[0])); i++) {

#ifdef DBG
		printf("\nExecutando instancia: %s\n", instances[i]);
#endif 

		char instance[100] = "instancias/";
		strcat(instance, instances[i]);
		strcat(instance, ".txt");
		ler_arquivo(instance);

		Solucao sol, mel_sol;
		mel_sol.fo = INT_MAX;

		int total_fo=0;
		double total_tempo=0, total_mel_tempo=0;

		for (int j = 0; j < total_exec; j++) {
#ifdef DBG
			printf("\nEtapa: %d\n", j+1);
#endif 
			heu_const_ale(sol);

			calcular_fo(sol);

			simulated_annealing(sol, arg_TI, arg_TC, arg_TR, arg_SAMAX, arg_TEM_MAX, TEM_TOT, TEM_MEL, NUM_SOL);

			total_fo += sol.fo;
			total_tempo += TEM_TOT;
			total_mel_tempo += TEM_MEL;

			if (sol.fo < mel_sol.fo) {
				memcpy(&mel_sol, &sol, sizeof(Solucao));
			}
		}
		double fo_media = total_fo / total_exec;
		double desvio = (fabs((total_fo / total_exec) - mel_sol.fo) / mel_sol.fo);
		double tempo_medio = (total_tempo / total_exec);
		double mel_tempo_medio = (total_mel_tempo / total_exec);

		fprintf(f, "%s;%d;%.3f;%.3f%%;%.3f;%.3f\n", 
			instances[i], //instancia
			mel_sol.fo, //Melhor FO
			fo_media, //FO Média
			desvio, //Desvio (%)
			tempo_medio, //Tempo Médio (seg.)
			mel_tempo_medio //T. Melhor (seg.)
		);

		char arq_solucao[100] = "solucoes/";
		strcat(arq_solucao, instances[i]);
		strcat(arq_solucao, ".txt");

		FILE* s = fopen(arq_solucao, "w");

		fprintf(s, "--- Relatorio da Solucao Otimizada ---\n\n");

		fprintf(s, "Tempo de transmissao: %d\n\n", mel_sol.fo);

		fprintf(s, "Nos Escolhidos:\n");

		for (int i = 0; i < NOS_FONTE; i++) {
			fprintf(s, "  - %d\n", mel_sol.vet_sol[i]);
		}

		fprintf(s, "\nArestas Usadas no Caminho:\n");
		fprintf(s, "tempo, origem, destino\n");

		int order[MAX_NOS];
		int vet_aux[MAX_NOS];

		for (int i = 0; i < NOS_FONTE; i++) {
			order[i] = mel_sol.vet_sol[i];
		}

		for (int i = NOS_FONTE; i < num_nos; i++) {
			vet_aux[i - NOS_FONTE] = mel_sol.vet_sol[i];
		}

		int tam_order = NOS_FONTE;
		int aux;

		for(int t=1; t<=mel_sol.fo; t++) {
			aux = tam_order;
			for (int i = 0; i < aux; i++) {
				for (int j = 0; j < num_nos - NOS_FONTE; j++) {
					if (mat_bin_rel[order[i]][vet_aux[j]]) {
						fprintf(s, "%d, %d, %d\n", t, order[i], vet_aux[j]);
						order[tam_order] = vet_aux[j];
						vet_aux[j] = 0;
						tam_order++;
						break;
					}
				}
			}
		}

		fclose(s);
	}

	fclose(f);

#endif

	return 0;
}

void grasp(Solucao& s, const double& LRC, const double& TEM_MAX,
           double& TEM_TOT, double& TEM_MEL, int& NUM_SOL)
{
	Solucao s_viz;
	
	clock_t h = clock();
	TEM_TOT = TEM_MEL = 0.0;
	NUM_SOL = 1;
	s.fo = INT_MAX;
	while (TEM_TOT < TEM_MAX)
	{
		heu_const_ale_gul(s_viz);
		calcular_fo(s_viz);
		//heu_BL_MM(s_viz);
		
		if (s_viz.fo < s.fo)
		{
			memcpy(&s, &s_viz, sizeof(Solucao));
			TEM_MEL = (double)(clock() - h) / CLOCKS_PER_SEC;

            #ifdef DBG
			    printf("FO: %d\tTempo: %.2f\n", s.fo, TEM_MEL);
            #endif
		}
		
		NUM_SOL++;
		printf("\n%d", NUM_SOL);

		TEM_TOT = (double)(clock() - h) / CLOCKS_PER_SEC;
	}
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
				static void (*vizinhancas[])(Solucao&) = { gerar_vizinho3, gerar_vizinho };
    			vizinhancas[!(NUM_SOL % 100)](s_viz);
				//gerar_vizinho3(s_viz);
				//heu_BL_rand(s_viz, 1 * (num_moc + 1) * num_obj);
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

void heu_BL_rand(Solucao& s, const int& iter) {
    int mel_fo = s.fo;
    while (true)
    {
        int flag = 1;
        for (int i = 0; i < iter; i++)
        {
            int no = rand() % num_nos;
            int no_origem = s.vet_sol[no];
            int no_destino;
            do
                no_destino = rand() % num_nos;
            while (no_destino == no_origem);
            int fo_ori = s.fo;
            s.vet_sol[no] = no_destino;
            calcular_fo(s);
            if (s.fo < mel_fo)
            {
                mel_fo = s.fo;
                flag = 0;
            }
            else
            {
                s.vet_sol[no] = no_origem;
                s.fo = fo_ori;
            }
        }
        if (flag)
            break;
    }
}

void heu_BL_MM(Solucao& s) {
    int mel_fo = s.fo;
    while (true)
    {
        int mel_no_i, mel_no_j;
        int flag = 0;
        for (int j = NOS_FONTE; j < num_nos; j++)
        {
            //int no_ori = s.vet_sol[j];
            for (int i = j+1; i < num_nos; i++)
            {
				int aux = s.vet_sol[j];
				s.vet_sol[j] = s.vet_sol[i];
				s.vet_sol[i] = aux;
                calcular_fo(s);
                if (s.fo < mel_fo)
                {
                    mel_fo = s.fo;
                    mel_no_j = j;
                    mel_no_i = i;
                    flag = 1;
					printf("\nFO: %d", s.fo);
                }
				aux = s.vet_sol[j];
				s.vet_sol[j] = s.vet_sol[i];
				s.vet_sol[i] = aux;
				printf("\n[%d][%d]", j, i);
            }
        }
        s.fo = mel_fo;
        if (flag) {
			int aux = s.vet_sol[mel_no_j];
			s.vet_sol[mel_no_j] = s.vet_sol[mel_no_i];
			s.vet_sol[mel_no_i] = aux;
		}
        else
            break;
    }
}

void heu_BL_PM(Solucao& s)
{
    int mel_fo = s.fo;
    INICIO : ;
    for (int j = NOS_FONTE; j < num_nos; j++)
    {
        int fo_ori = s.fo;
        for (int i = j+1; i < num_nos; i++)
        {
			int aux = s.vet_sol[j];
			s.vet_sol[j] = s.vet_sol[i];
			s.vet_sol[i] = aux;
			calcular_fo(s);
            if (s.fo < mel_fo)
            {
                mel_fo = s.fo;
				printf("\nFO: %d", s.fo);
                goto INICIO;
            }
            else
            {
				aux = s.vet_sol[j];
				s.vet_sol[j] = s.vet_sol[i];
				s.vet_sol[i] = aux;
                s.fo = fo_ori;
            }
        }
    }
    s.fo = mel_fo;
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
	//construir o vetor solu��o baseado na ordem do vetor_aux, colocar os n�s com maior n�mero de arestas primeiro
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

//trocar N� FONTE por N� COMUM
void gerar_vizinho(Solucao& sol) {
	int pos_origem = rand() % NOS_FONTE;
	int pos_destino = (rand() % num_nos - NOS_FONTE) + NOS_FONTE;

	int aux = sol.vet_sol[pos_origem];
	sol.vet_sol[pos_origem] = sol.vet_sol[pos_destino];
	sol.vet_sol[pos_destino] = aux;
}

//trocar QUALQUER N�
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

void gerar_vizinho3(Solucao& sol) {
	int pos_no_fonte, no_fonte, pos_vizinhos[MAX_NOS], num_vizinhos=0, pos_origem, pos_destino, aux;
	
	pos_no_fonte = rand() % NOS_FONTE;
	no_fonte = sol.vet_sol[pos_no_fonte];

	if (vet_qtd_rel[no_fonte] > 1) {

		for (int i=NOS_FONTE; i<num_nos; i++) {
			if (mat_bin_rel[no_fonte][sol.vet_sol[i]]) {
				pos_vizinhos[num_vizinhos] = i;
				num_vizinhos++;
			}
		}

		pos_origem = rand() % num_vizinhos;
		pos_destino = rand() % num_vizinhos;

		while (pos_destino == pos_origem) {
			pos_destino = rand() % num_vizinhos;
		}

		aux = sol.vet_sol[pos_vizinhos[pos_origem]];
		sol.vet_sol[pos_vizinhos[pos_origem]] = sol.vet_sol[pos_vizinhos[pos_destino]];
		sol.vet_sol[pos_vizinhos[pos_destino]] = aux;
	} else {
		gerar_vizinho(sol);
	}
}

void calcular_fo2(Solucao& sol) {
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

void calcular_fo(Solucao& sol) {
	sol.fo = 0;

    // Estruturas estáticas para evitar alocação de memória repetitiva
    static int prioridade[MAX_NOS];
    static bool infectado[MAX_NOS];
    static int fila_transmissores[MAX_NOS]; // Substitui o 'order'
    
    // Inicialização O(N)
    // Se o MAX_NOS for muito grande (ex: > 10.000), use memset. Para < 1000, for loop é ok.
    for (int i = 0; i <= num_nos; i++) {
        infectado[i] = false;
        // Inicializa com um valor maior que qualquer índice possível
        prioridade[i] = num_nos + 1; 
    }

    // 1. Configurar Fontes e Prioridades
    // Os nós fonte entram na fila e são marcados
    for (int i = 0; i < NOS_FONTE; i++) {
        int no = sol.vet_sol[i];
        infectado[no] = true;
        fila_transmissores[i] = no;
    }

    // Mapeia a prioridade dos candidatos.
    // O nó que está em sol.vet_sol[k] tem prioridade 'k'.
    // Isso permite verificar em O(1) qual vizinho aparece antes no cromossomo.
    for (int i = NOS_FONTE; i < num_nos; i++) {
        prioridade[sol.vet_sol[i]] = i;
    }

    int total_infectados = NOS_FONTE;
    int tam_fila = NOS_FONTE;

    // Vetor auxiliar para armazenar quem recebeu mensagem neste turno
    // para adicioná-los à fila apenas no final do pulso de clock.
    static int novos_neste_turno[MAX_NOS]; 

    // 2. Loop de Tempo (Simulação)
    while (total_infectados < num_nos) {
        sol.fo++;
        int qtd_novos = 0;
        int snapshot_tam_fila = tam_fila; // Apenas quem JÁ tinha a mensagem transmite

        // Loop sobre os transmissores ativos (O(N_infectados))
        for (int i = 0; i < snapshot_tam_fila; i++) {
            int transmissor = fila_transmissores[i];
            
            int melhor_vizinho = -1;
            int melhor_prio = num_nos + 1;

            // Loop sobre os vizinhos (O(Grau_do_no)) - Muito mais rápido que O(N)
            for (int k = 0; k < vet_qtd_rel[transmissor]; k++) {
                int vizinho = mat_rel[transmissor][k];

                if (!infectado[vizinho]) {
                    // Verifica se esse vizinho tem prioridade maior (menor índice)
                    // do que o melhor encontrado até agora para este transmissor
                    if (prioridade[vizinho] < melhor_prio) {
                        melhor_prio = prioridade[vizinho];
                        melhor_vizinho = vizinho;
                    }
                }
            }

            // Se o transmissor conseguiu encontrar um alvo livre
            if (melhor_vizinho != -1) {
                infectado[melhor_vizinho] = true; // Marca como usado para que outro transmissor NESTE MESMO turno não o pegue
                novos_neste_turno[qtd_novos++] = melhor_vizinho;
                total_infectados++;
            }
        }

        // Adiciona os novos infectados à fila principal de transmissores para os próximos turnos
        for (int i = 0; i < qtd_novos; i++) {
            fila_transmissores[tam_fila++] = novos_neste_turno[i];
        }
        
        // Se em um turno ninguém for infectado e ainda restam nós, temos um grafo desconexo
        // ou isolado (break de segurança para evitar loop infinito, opcional)
        if (qtd_novos == 0 && total_infectados < num_nos) break; 
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

		mat_rel[a][ vet_qtd_rel[a] ] = b;
		mat_rel[b][ vet_qtd_rel[b] ] = a;

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