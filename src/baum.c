/* **	Purpose: Changed the convergence criterion from ratio to absolute value.
   Solution to the problem of learning: how to adjust the parameters of the hidden Markov model
   lambda = (A, B, pi), so as to maximize P(O|lambda)
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrutil.h"
#include "hmm.h"

#define DELTA 0.001
#define EPSILON 0.0001

/*  Function that calculates the gamma as: gamma_t(i) = { alpha_t(i) * beta_t(i) } / P(O|lambda) */
void ComputeGamma(HMM *hmm, int T, double **alpha, double **beta, double **gamma)
{
    for(int t = 1; t <= T; t++) {
    	double sum = 0.0;
        for(int j = 1; j <= hmm->N; j++) {
            gamma[t][j] = alpha[t][j] * beta[t][j];
            sum += gamma[t][j];
        }
        for(int i = 1; i <= hmm->N; i++) gamma[t][i] /= sum;
    }
}

void BaumWelch(HMM *hmm, int T, int *O, double **alpha, double **beta,	double **gamma, int *pniter,
		double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, logprobprev;

	xi = alloc_3d_matrix(T, hmm->N);
	scale = dvector(1, T);

	ForwardWithScale(hmm, T, O, alpha, scale, &logprobf);
	*plogprobinit = logprobf; /* log P(O |intial model) */
	BackwardWithScale(hmm, T, O, beta, scale, &logprobb);
	ComputeGamma(hmm, T, alpha, beta, gamma);
	compute_epsilon(hmm, T, O, alpha, beta, xi);
	logprobprev = logprobf;

	do
	{
        /*  Re-estimates the frequency of the state i over the time t=1 */
        /*  The robability of being in state S_i at time t */
		for(i = 1; i <= hmm->N; i++)
			hmm->pi[i] = .001 + .999*gamma[1][i];

        /*  Re-estimating the transition probability matrix and the different
            observation symbols per state */

		for(i = 1; i <= hmm->N; i++)
		{   /*  a^_ij = sum xi_t(i,j) / sum gamma_t(i) */
			denominatorA = 0.0;
			for(t = 1; t <= T - 1; t++)
				denominatorA += gamma[t][i];

			for(j = 1; j <= hmm->N; j++)
			{
				numeratorA = 0.0;
				for(t = 1; t <= T - 1; t++)
					numeratorA += xi[t][i][j];
				hmm->A[i][j] = .001 + .999*numeratorA/denominatorA;
			}

            /*  b^_j(k) = sum gamma_t(j) / sum gamma_t(j) */
			denominatorB = denominatorA + gamma[T][i];
			for(k = 1; k <= hmm->M; k++)
			{
				numeratorB = 0.0;
				for(t = 1; t <= T; t++)
				{
					if (O[t] == k)
						numeratorB += gamma[t][i];
				}
				hmm->B[i][k] = .001 + .999*numeratorB/denominatorB;
			}
		}

		ForwardWithScale(hmm, T, O, alpha, scale, &logprobf);
		BackwardWithScale(hmm, T, O, beta, scale, &logprobb);
		ComputeGamma(hmm, T, alpha, beta, gamma);
		compute_epsilon(hmm, T, O, alpha, beta, xi);

        /*  Calculates the difference between the logarithmic probability
            of two iterations, the current and the previous */
		delta = logprobf - logprobprev;
		logprobprev = logprobf;
        /*  Increments the counter that tracks the number of iterations */
		l++;
	}while(delta > DELTA);
	/*  Exits from the do while if the logarithmic probability does not change much */

	*pniter = l;
	/*  log P(O|estimated models) */
	*plogprobfinal = logprobf;
	free_3d_matrix(xi, T, hmm->N);
	free_dvector(scale, 1, T);
}

/*  Function that calculates the xi as: xi_t(i,j) = { alpha_t(i) * beta_t+1(j) * a_ij * b_j } / P(O|lambda) */
void compute_epsilon(HMM* hmm, int T, int *O, double **alpha, double **beta, double ***xi)
{
	for(int t = 1; t <= T - 1; t++) {
		double sum = 0.0;
		for(int i = 1; i <= hmm->N; i++){
			for(int j = 1; j <= hmm->N; j++) {
				xi[t][i][j] = alpha[t][i] * beta[t+1][j] * hmm->A[i][j] * hmm->B[j][O[t]];
				sum += xi[t][i][j];
			}
		}
		for(int i = 1; i <= hmm->N; i++) {
			for(int j = 1; j <= hmm->N; j++) {
				xi[t][i][j] /= sum;
			}
		}
	}
}

/*  Allocates the array of 3 dimensions x: C * N * N. */
double ***alloc_3d_matrix(int C, int N)
{
	double ***x = (double ***)malloc(C * sizeof(double**));
	for(int t = 0; t < C; t++) x[t] = dmatrix(1, N, 1, N);
	return x;
}

void free_3d_matrix(double ***x, int C, int N)
{
	for(int t = 0; t < C; t++) free_dmatrix(x[t], 1, N, 1, N);
	free(x);
}

void BaumWelch_C(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma)
{
    double observe_prob, logprobb;
    Forward(hmm, frame_num, O, alpha, &observe_prob);
    Backward(hmm, frame_num, O, beta, &logprobb);
    ComputeGamma(hmm, frame_num, alpha, beta, gamma);
    double ***epsilon = alloc_3d_matrix(frame_num, hmm->N);
    compute_epsilon(hmm, frame_num, O, alpha, beta, epsilon);

    double observe_prob_previous = observe_prob;
    double delta;
    do {
        for(int i = 1; i <= hmm->N; i++) hmm->pi[i] = .001 + .999*gamma[1][i];;

        for(int i = 1; i <= hmm->N; i++) {   /*  a^_ij = sum epsilon_t(i,j) / sum gamma_t(i) */
            double sum = 0.0;
            for(int t = 1; t <= frame_num - 1; t++) sum += gamma[t][i];

            for(int j = 1; j <= hmm->N; j++) {
            	double tmp = 0.0;
                for(int t = 1; t <= frame_num - 1; t++) tmp += epsilon[t][i][j];
                hmm->A[i][j] = 0.001f + 0.999f * tmp / sum;
            }

            /*  b^_j(k) = sum gamma_t(j) / sum gamma_t(j) */
            double sum_gamma = sum + gamma[frame_num][i];
            for(int k = 1; k <= hmm->M; k++) {
                sum = 0.0;
                for(int t = 1; t <= frame_num; t++) {
                    if (O[t] == k) sum += gamma[t][i];
                }
                hmm->B[i][k] = 0.001f + 0.999f * sum / sum_gamma;
            }
        }

        Forward(hmm, frame_num, O, alpha, &observe_prob);
        Backward(hmm, frame_num, O, beta, &logprobb);
        ComputeGamma(hmm, frame_num, alpha, beta, gamma);
        compute_epsilon(hmm, frame_num, O, alpha, beta, epsilon);

        /* Calculates the difference between the logarithmic probability of two iterations */
        delta = log(observe_prob) - log(observe_prob_previous);
        observe_prob_previous = observe_prob;
    } while(delta > DELTA);
    free_3d_matrix(epsilon, frame_num, hmm->N);
}
