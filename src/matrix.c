#include "matrix.h"

void print_error_exit(char *testo)
{
    fprintf(stderr,"%s\n",testo);
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
	for(int i = 0; i< row; i++) {
		m[i] = (double *)calloc(column, sizeof(double));
	}
	return m;
}

int **matrice_i(int row, int column)
{
	int **m = (int **) malloc(sizeof(int*) * row);
	for(int i=0; i<row; i++) {
		m[i] = (int *)calloc(column, sizeof(int));
	}
	return m;
}

/*  Function that calculates the transpose of a matrix of type double */
double **trasposta_d(double **x, int row, int column)
{
    double **y = matrix_double(column, row);
    int n,m;
    for(n=0; n<column; n++)
        for(m=0; m<row; m++) y[n][m]=x[m][n];
    return y;
}

double **gemm(double **A, int row_A, int column_A, double **B, int column_B)
{
    /*  If A is MxN and B id NxP then C is MxP */
    double temp;
    int n,m,i;
    double **C = matrix_double(row_A, column_B);

    for(n=0; n<row_A; n++){
        for(m=0; m<column_B; m++) {
            temp = 0.0;
            for(i=0; i<column_A; i++) temp += A[n][i]*B[i][m];
            C[n][m]=temp;
        }
    }
    return C;
}

