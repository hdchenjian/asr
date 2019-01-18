#include <stdio.h>
#include <stdlib.h>

void print_error_exit(char *testo);
double *array_double(int dimensione);
int *array_int(int dimensione);
double **matrix_double(int righe, int colonne);
int **matrice_i(int righe, int colonne);
double **trasposta_d(double **x, int righe, int colonne);
double **gemm(double **A, int righe_A, int colonne_A, double **B, int colonne_B);
void save_matrix_double(char nome_file[], double **dati, int righe, int colonne);
void load_matrix_double(char nome_file[], double **dati, int righe, int colonne);
void save_array_int(char nome_file[], int *dati, int elementi);
void load_array_int(char nome_file[], int *dati, int elementi);
