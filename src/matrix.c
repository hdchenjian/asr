#include "matrix.h"

void print_error_exit(char *testo)
{
    fprintf(stderr,"%s\n",testo);
	exit(1);
}

/*  Function that allocates space in memory to store a double type vector */
double *array_double(int dimensione)
{
	double *v = (double *)calloc(dimensione, sizeof(double));
	if (v == NULL) print_error_exit("array_double failed");
    return v;
}

/*  Function that allocates space in memory to store a int type vector */
int *array_int(int dimensione)
{
	int *v = (int *)calloc(dimensione, sizeof(int));
	if (v == NULL) print_error_exit("array_int failed");
	return v;
}

double **matrix_double(int righe, int colonne)
{
	double **m = (double **) malloc(sizeof(double*) * righe);
	if (m==NULL) print_error_exit("malloc failed");
	for(int i=0; i<righe; i++) {
		m[i] = (double *)calloc(colonne, sizeof(double));
		if (m[i]==NULL) print_error_exit("malloc failed");
	}
	return m;
}

int **matrice_i(int righe, int colonne)
{
	int **m = (int **) malloc(sizeof(int*) * righe);
	if (m==NULL) print_error_exit("malloc failed");
	for(int i=0; i<righe; i++) {
		m[i] = (int *)calloc(colonne, sizeof(int));
		if (m[i]==NULL) print_error_exit("malloc failed");
	}
	return m;
}

/*  Function that calculates the transpose of a matrix of type double */
double **trasposta_d(double **x, int righe, int colonne)
{
    double **y = matrix_double(colonne, righe);
    int n,m;

    for(n=0; n<colonne; n++)
        for(m=0; m<righe; m++)
            y[n][m]=x[m][n];

    return y;
}

double **gemm(double **A, int righe_A, int colonne_A, double **B, int colonne_B)
{
    /*  If A is MxN and B id NxP then C is MxP */
    double temp;
    int n,m,i;
    double **C = matrix_double(righe_A, colonne_B);

    for(n=0; n<righe_A; n++){
        for(m=0; m<colonne_B; m++) {
            temp = 0.0;
            for(i=0; i<colonne_A; i++) temp += A[n][i]*B[i][m];
            C[n][m]=temp;
        }
    }
    return C;
}

void save_matrix_double(char nome_file[], double **dati, int righe, int colonne)
{
    FILE *fp = fopen(nome_file,"wb") ;
    if(fp==NULL) print_error_exit("open file failed.");
    size_t dimensione = sizeof(double);
    size_t elementi = colonne ;
    size_t ritorno;
    for(int i=0; i<righe; i++) {
        ritorno = fwrite(dati[i], dimensione, elementi, fp);
        if (ritorno != elementi) print_error_exit("Errore scrittura su file.");
    }
    fclose(fp);
}

void load_matrix_double(char nome_file[], double **dati, int righe, int colonne)
{
    FILE *fp = fopen(nome_file,"rb") ;
    if(fp==NULL) print_error_exit("Errore di apertura del file.");
    size_t dimensione = sizeof(double);
    size_t elementi = colonne ;
    size_t ritorno;
    for(int i=0; i<righe; i++) {
        ritorno = fread(dati[i], dimensione, elementi, fp);
        if (ritorno != colonne) print_error_exit("Errore scrittura su file.");
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
