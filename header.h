#pragma once

#define MAX_NOS 1001
#define NOS_FONTE 1

typedef struct tSolucao {
	int vet_sol[MAX_NOS];
	int fo;
} Solucao;

int num_nos;
int num_arestas;
int mat_rel[MAX_NOS][MAX_NOS];
int mat_bin_rel[MAX_NOS][MAX_NOS];
int vet_qtd_rel[MAX_NOS];

void ler_arquivo(const char* path);
void escrever_dados(const char* arq);
void calcular_fo(Solucao& sol);
void escrever_solucao(Solucao& sol);
void gerar_vizinho(Solucao& sol);
void gerar_vizinho2(Solucao& sol);
void heu_const_ale(Solucao& sol);
void heu_const_gul(Solucao& sol);
void heu_const_ale_gul(Solucao& sol);
void simulated_annealing(Solucao& s, const double& TI, const double& TC,
	const double& TR, const int& SAMAX, const double& TEM_MAX,
	double& TEM_TOT, double& TEM_MEL, int& NUM_SOL);