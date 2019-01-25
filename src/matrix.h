void print_error_exit(char *info);
double **matrix_double(int row, int column);
void free_matrix_double(double **m, int row);

double **dmatrix(int row, int col);
void free_dmatrix(double **m, int row);

int **imatrix(int row, int col);
void free_imatrix(int **m, int row);
