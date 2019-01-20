/*  Given the sequence of observations O=O1,O2,...,OT and the hidden Markov model
    lambda=(A, B, pi), how to can be calculated efficiently P(O|lambda) */

#include "hmm.h"

/*  Where alpha is the auxiliary variable forward, defined as the probability of
    partial sequence of observations up to time t and the state Si,
    alpha: frame_num * N, alpha_t(i)=P(O1 O2...Ot, q_t=S_i|lambda) */
void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob)
{
    /*  alpha_1(i) = pi_i * b_i(O1), with 1 <= i <= N */
    for(int i = 1; i <= phmm->N; i++) alpha[1][i] = phmm->pi[i]* phmm->B[i][O[0]];

    /*  for t=1,2...T-1, with 1 <= j <= N
        alpha_t+1(j) = {sum alpha_t(i)*a_ij}*b_j(Ot+1) */
    for(int t = 1; t <= T - 1; t++) {
        for(int j = 1; j <= phmm->N; j++) {
            double sum = 0;
            for(int i = 1; i <= phmm->N; i++) sum += alpha[t][i]* (phmm->A[i][j]);
            alpha[t+1][j] = sum*(phmm->B[j][O[t]]);
        }
    }

    *pprob = 0.0;
    for(int i = 1; i <= phmm->N; i++) *pprob += alpha[T][i];
}

void ForwardWithScale(HMM *phmm, int T, int *O, double **alpha, double *scale, double *pprob)
/*  pprob is the probability value relating to the logarithm */
{
	int	i, j;
	int	t;

	double sum;

	/* 1. Initialization */
	scale[1] = 0.0;
	for(i = 1; i <= phmm->N; i++)
	{
		alpha[1][i] = phmm->pi[i]* (phmm->B[i][O[1]]);
		scale[1] += alpha[1][i];
	}
	for(i = 1; i <= phmm->N; i++)
		alpha[1][i] /= scale[1];

	/* 2. Induction */
	for(t = 1; t <= T - 1; t++)
	{
		scale[t+1] = 0.0;
		for(j = 1; j <= phmm->N; j++)
		{
		    sum = 0.0;
			for(i = 1; i <= phmm->N; i++)
				sum += alpha[t][i]* (phmm->A[i][j]);

			alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
			scale[t+1] += alpha[t+1][j];
		}
		for(j = 1; j <= phmm->N; j++)
			alpha[t+1][j] /= scale[t+1];
	}

	/* 3. Termination */
	*pprob = 0.0;
	for(t = 1; t <= T; t++)
		*pprob += log(scale[t]);
}
