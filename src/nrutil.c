#include <malloc.h>
#include <stdio.h>

/*  Prepare a memory location that will contain a vector of type double */
double *dvector(int nl, int nh)
{
	double *v = (double *)calloc((unsigned) (nh-nl+1),sizeof(double));
	return v-nl;
}

/*  Prepare a memory location that will contain a matrix of type double */
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
	{
		m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
		m[i] -= ncl;
	}
	return m;
}

/*  Prepare a memory location that will contain a matrix of type int */
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)calloc((unsigned) (nrh-nrl+1),sizeof(int*));
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)calloc((unsigned) (nch-ncl+1),sizeof(int));
		m[i] -= ncl;
	}
	return m;
}

void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for(i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for(i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
