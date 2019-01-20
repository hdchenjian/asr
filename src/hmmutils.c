#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "nrutil.h"
#include "hmm.h"

void FreeHMM(HMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->B, 1, phmm->N, 1, phmm->M);
	free_dvector(phmm->pi, 1, phmm->N);
}

void save_hmm(char nome_file[], HMM *phmm)
{
    FILE *fp = fopen(nome_file,"wb");
    if(fp==NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(double);
    size_t elementi = phmm->N;
    size_t ritorno;

    int i;
    for(i=1; i<phmm->N; i++) {
        ritorno = fwrite(phmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            print_error_exit("Errore scrittura su file.");
    }

    elementi = phmm->M;
    for(i=1; i<phmm->N; i++) {
        ritorno = fwrite(phmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = phmm->N;
    ritorno = fwrite(phmm->pi, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}

/*  Function that loads a binary file on hidden Markov model */
void load_hmm(char nome_file[], HMM *phmm, int N, int M)
{
    int i;
    phmm->N = N;
    phmm->M = M;
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    phmm->pi = (double *) dvector(1, phmm->N);

    FILE *fp = fopen(nome_file,"rb");
    if(fp==NULL) print_error_exit("Errore di apertura del file.");

    size_t dimensione = sizeof(double);
    size_t elementi = phmm->N;
    size_t ritorno;

    for(i=1; i<phmm->N; i++) {
        ritorno = fread(phmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = phmm->M;
    for(i=1; i<phmm->N; i++) {
        ritorno = fread(phmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }

    elementi = phmm->N;
    ritorno = fread(phmm->pi, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}

void InitHMM_SC1(HMM *phmm, int N, int M, int seed)
{
	int i, j, k;
	double sum;
	hmmsetseed(seed);
    phmm->M = M;
    phmm->N = N;
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    for(i = 1; i <= phmm->N; i++) {
		sum = 0.0;
        for(j = 1; j <= phmm->N; j++) {
            if(j<i || j>=i+3) phmm->A[i][j] = 0.0;
            else phmm->A[i][j] = hmmgetrand();
			sum += phmm->A[i][j];
		}
        for(j = 1; j <= phmm->N; j++) phmm->A[i][j] /= sum;
	}

    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    for(j = 1; j <= phmm->N; j++) {
        sum = 0.0;
        for(k = 1; k <= phmm->M; k++) {
            phmm->B[j][k] = hmmgetrand();
			sum += phmm->B[j][k];
		}
        for(k = 1; k <= phmm->M; k++) phmm->B[j][k] /= sum;
	}

    phmm->pi = (double *) dvector(1, phmm->N);
    phmm->pi[1] = 1;
    for(i = 2; i <= phmm->N; i++) {
        phmm->pi[i] = 0;
	}
}
