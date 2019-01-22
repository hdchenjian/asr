#include <stdio.h>
#include <stdlib.h>

void print_error_exit(char *testo);
double **matrix_double(int row, int column);
void free_matrix_double(double **m, int row);
int **matrice_i(int row, int column);
double **trasposta_d(double **x, int row, int column);
double **gemm(double **A, int row_A, int column_A, double **B, int column_B);
void save_matrix_double(char nome_file[], double **dati, int row, int column);
void load_matrix_double(char nome_file[], double **dati, int row, int column);
void save_array_int(char nome_file[], int *dati, int elementi);
void load_array_int(char nome_file[], int *dati, int elementi);
