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
	if (m==NULL) print_error_exit("malloc failed");
	for(int i = 0; i< row; i++) {
		m[i] = (double *)calloc(column, sizeof(double));
		if (m[i]==NULL) print_error_exit("malloc failed");
	}
	return m;
}

int **matrice_i(int row, int column)
{
	int **m = (int **) malloc(sizeof(int*) * row);
	if (m==NULL) print_error_exit("malloc failed");
	for(int i=0; i<row; i++) {
		m[i] = (int *)calloc(column, sizeof(int));
		if (m[i]==NULL) print_error_exit("malloc failed");
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

void save_matrix_double(char nome_file[], double **dati, int row, int column)
{
    FILE *fp = fopen(nome_file,"wb") ;
    if(fp==NULL) print_error_exit("open file failed.");
    size_t dimensione = sizeof(double);
    size_t elementi = column ;
    size_t ritorno;
    for(int i=0; i<row; i++) {
        ritorno = fwrite(dati[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }
    fclose(fp);
}

void load_matrix_double(char nome_file[], double **dati, int row, int column)
{
    FILE *fp = fopen(nome_file,"rb") ;
    if(fp==NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(double);
    size_t elementi = column ;
    size_t ritorno;
    for(int i=0; i < row; i++) {
        ritorno = fread(dati[i], dimensione, elementi, fp);
        if (ritorno != column) print_error_exit("Errore scrittura su file.");
    }
    fclose(fp);
}

void save_array_int(char nome_file[], int *dati, int elementi)
{
    FILE *fp =fopen(nome_file,"wb") ;
    if(fp==NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(int);
    size_t ritorno;
    ritorno = fwrite(dati, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}

void load_array_int(char nome_file[], int *dati, int elementi)
{
    FILE *fp = fopen(nome_file,"rb") ;
    if(fp==NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(int);
    size_t ritorno;
    ritorno = fread(dati, dimensione, elementi, fp);
    if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    fclose(fp);
}
