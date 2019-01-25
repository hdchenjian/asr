#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

void print_error_exit(char *info)
{
    fprintf(stderr,"%s\n",info);
	exit(1);
}

void free_matrix_double(double **m, int row)
{
	for(int i = 0; i< row; i++) free(m[i]);
	free(m);
}

double **matrix_double(int row, int column)
{
	double **m = (double **) malloc(sizeof(double*) * row);
	for(int i = 0; i < row; i++) {
		m[i] = (double *)calloc(column, sizeof(double));
	}
	return m;
}

double **dmatrix(int row, int col)
{
	double **m=(double **) calloc(row + 1, sizeof(double*));
	for(int i = 0; i <= row; i++) {
		m[i] = (double *)calloc(col + 1, sizeof(double));
	}
	return m;
}

void free_dmatrix(double **m, int row)
{
	for(int i = 0; i <= row; i++) free(m[i]);
	free(m);
}

int **imatrix(int row, int col)
{
	int **m=(int **) calloc(row + 1, sizeof(int*));
	for(int i = 0; i <= row; i++) {
		m[i] = (int *)calloc(col + 1, sizeof(int));
	}
	return m;
}

void free_imatrix(int **m, int row)
{
	for(int i = 0; i <= row; i++) free(m[i]);
	free(m);
}
