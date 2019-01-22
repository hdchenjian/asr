#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "nrutil.h"
#include "hmm.h"

void FreeHMM(HMM *hmm)
{
	free_dmatrix(hmm->A, 1, hmm->N, 1, hmm->N);
	free_dmatrix(hmm->B, 1, hmm->N, 1, hmm->M);
	free_dvector(hmm->pi, 1, hmm->N);
}

void save_hmm(char nome_file[], HMM *hmm)
{
    FILE *fp = fopen(nome_file,"wb");
    if(fp == NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(double);
    size_t elementi = hmm->N;
    size_t ritorno;

    int i;
    for(i=1; i<hmm->N; i++) {
        ritorno = fwrite(hmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = hmm->M;
    for(i=1; i<hmm->N; i++) {
        ritorno = fwrite(hmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = hmm->N;
    ritorno = fwrite(hmm->pi, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}

/*  Function that loads a binary file on hidden Markov model */
void load_hmm(char nome_file[], HMM *hmm, int N, int M)
{
    int i;
    hmm->N = N;
    hmm->M = M;
    hmm->A = (double **) dmatrix(1, hmm->N, 1, hmm->N);
    hmm->B = (double **) dmatrix(1, hmm->N, 1, hmm->M);
    hmm->pi = (double *) dvector(1, hmm->N);

    FILE *fp = fopen(nome_file,"rb");
    if(fp==NULL) print_error_exit("Errore di apertura del file.");

    size_t dimensione = sizeof(double);
    size_t elementi = hmm->N;
    size_t ritorno;

    for(i=1; i<hmm->N; i++) {
        ritorno = fread(hmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = hmm->M;
    for(i=1; i<hmm->N; i++) {
        ritorno = fread(hmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = hmm->N;
    ritorno = fread(hmm->pi, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}

void InitHMM(HMM *hmm, int N, int M, int seed)
{
	int i, j, k;
	double sum;
	srand(seed);
    hmm->M = M;
    hmm->N = N;
    hmm->A = (double **) dmatrix(1, hmm->N, 1, hmm->N);
    for(i = 1; i <= hmm->N; i++) {
		sum = 0.0;
        for(j = 1; j <= hmm->N; j++) {
            if(j<i || j>=i+3) hmm->A[i][j] = 0.0;
            else hmm->A[i][j] = (double)rand() / RAND_MAX;
			sum += hmm->A[i][j];
		}
        for(j = 1; j <= hmm->N; j++) hmm->A[i][j] /= sum;
	}

    hmm->B = (double **) dmatrix(1, hmm->N, 1, hmm->M);
    for(j = 1; j <= hmm->N; j++) {
        sum = 0.0;
        for(k = 1; k <= hmm->M; k++) {
            hmm->B[j][k] = (double)rand() / RAND_MAX;
			sum += hmm->B[j][k];
		}
        for(k = 1; k <= hmm->M; k++) hmm->B[j][k] /= sum;
	}

    hmm->pi = (double *) dvector(1, hmm->N);
    hmm->pi[1] = 1;
    for(i = 2; i <= hmm->N; i++) {
        hmm->pi[i] = 0;
	}
}
