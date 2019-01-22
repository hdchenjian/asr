/** Solution to the problem of decoding:
    given the sequence of observations O = O1O2 ... OT and the Hidden Markovian Model lambda = (A, B, pi),
    how to choose a corresponding state sequence Q = q1 q2 ... qT that is optimal, to justify the better observations? */
#include <math.h>
#include "hmm.h"
#include "nrutil.h"

#define VITHUGE  100000000000.0
#define VITHUGE2  100000000.0
//#define VITHUGE  -100000000000.0
/*  Want to find Q (sequences of states) that maximizes P(O,Q|lambda).
    It defines the quantity:
    delta_t(i) = max P(q1 q2 ... qt = i, O2 O1...Ot|lambda)
    which is the highest probability of a single path, at time t, which justifies
    the first t observations and ends in the state Si.
    By induction we have:
    delta_t+1(j) = {max delta_t(s)*a_ij}+b_j(Ot+1)
    It keeps track of the argument that maximizes the above formula for each t e j
    in psi_t (j)
*/

/* delta psi: frame_num * status_num, q:  1 * frame_num. */
void Viterbi(HMM *hmm, int T, int *O, int *path, double *pprob)
{
	double **delta = dmatrix(1, T, 1, hmm->N);
	int	**psi = imatrix(1, T, 1, hmm->N);
	for(int i = 1; i <= hmm->N; i++) {
		delta[1][i] = hmm->pi[i] * (hmm->B[i][O[0]]);
		psi[1][i] = 0;
	}

	/*  2. Recursion */
	/*  delta_t(j) = max{delta_t-1(i)*a_ij}*b_j(Ot), with 2 <= t <= T e 1 <= j <= N
        psi_t(j) = arg max{delta_t-1(i)*a_ij}, with 2 <= t <= T e 1 <= j <= N */
	for(int t = 2; t <= T; t++) {
		for(int j = 1; j <= hmm->N; j++) {
			double maxval = 0.0;
			int maxvalind = 1;
			for(int i = 1; i <= hmm->N; i++) {
				double val = delta[t-1][i]*(hmm->A[i][j]);
				if (val > maxval) {
				    maxval = val;
					maxvalind = i;
				}
			}
			delta[t][j] = maxval*(hmm->B[j][O[t - 1]]);
			psi[t][j] = maxvalind;
		}
	}

	*pprob = 0.0;
	path[T - 1] = 1;
	for(int i = 1; i <= hmm->N; i++) {
	    if (delta[T][i] > *pprob) {
			*pprob = delta[T][i];
			path[T - 1] = i;
		}
	}

    /*  4. Backtracking on the path (sequence of states) q*_t = psi_t+1 (q*_t+1), for t=T-1, T-2, ..., 1*/
	for(int t = T - 1; t >= 1; t--) path[t - 1] = psi[t+1][path[t]];
	free_dmatrix(delta, 1, T, 1, hmm->N);
	free_imatrix(psi, 1, T, 1, hmm->N);
}

/*  Function with the Viterbi algorithm, with operations with logarithms. */
void ViterbiLog(HMM *hmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
    int i, j;
    int t;

    int     maxvalind;
    double  maxval, val;

	double  **biot;

	/*  0. Preprocessing */

	for(i = 1; i <= hmm->N; i++)
		hmm->pi[i] = log(hmm->pi[i]);
	for(i = 1; i <= hmm->N; i++)
		for(j = 1; j <= hmm->N; j++)
		{
			hmm->A[i][j] = log(hmm->A[i][j]);
		}

	biot = dmatrix(1, hmm->N, 1, T);
	for(i = 1; i <= hmm->N; i++)
		for(t = 1; t <= T; t++)
		{
			biot[i][t] = log(hmm->B[i][O[t]]);
		}

    /*  1. Initialisation  */

    for(i = 1; i <= hmm->N; i++)
    {
        delta[1][i] = hmm->pi[i] + biot[i][1];
        psi[1][i] = 0;
    }

    /*  2. Recursion */
    for(t = 2; t <= T; t++)
    {
        for(j = 1; j <= hmm->N; j++)
        {
            maxval = -VITHUGE;
            maxvalind = 1;
            for(i = 1; i <= hmm->N; i++)
            {
                val = delta[t-1][i] + (hmm->A[i][j]);
                if (val > maxval)
                {
                    maxval = val;
                    maxvalind = i;
                }
            }
            delta[t][j] = maxval + biot[j][t];
            psi[t][j] = maxvalind;
        }
    }

    /*  3. Termination */
    *pprob = -VITHUGE;
    q[T] = 1;
    for(i = 1; i <= hmm->N; i++)
    {
        if (delta[T][i] > *pprob)
        {
            *pprob = delta[T][i];
            q[T] = i;
        }
    }

    /* 4. Backtracking on the path (sequence of states) */
    for(t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
}

void ViterbiLog_C(HMM *hmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
    int     maxvalind;
    double  maxval, val;
   	for(int i = 1; i <= hmm->N; i++) {
	    if(hmm->pi[i]*hmm->B[i][O[1]] != 0) delta[1][i] = log (hmm->pi[i] * (hmm->B[i][O[1]]));
        else delta[1][i] = -VITHUGE2;
        psi[1][i] = 0;
	}

    /*  2. Recursion */
    for(int t = 2; t <= T; t++) {
        for(int j = 1; j <= hmm->N; j++) {
            maxval = -VITHUGE;
            maxvalind = 1;
            for(int i = 1; i <= hmm->N; i++) {
                if(hmm->A[i][j]*(hmm->B[j][O[t]]) != 0)
                    val = delta[t-1][i] + log(hmm->A[i][j]*(hmm->B[j][O[t]]));
                else
                    val = delta[t-1][i];
                if (val > maxval)
                {
                    maxval = val;
                    maxvalind = i;
                }
            }
            delta[t][j] = maxval;
            psi[t][j] = maxvalind;
        }
    }

    /*  3. Termination */
    *pprob = -VITHUGE;
    q[T] = 1;
    for(int i = 1; i <= hmm->N; i++) {
        if (delta[T][i] > *pprob) {
            *pprob = delta[T][i];
            q[T] = i;
        }
    }

    /* 4. Backtracking on the path (sequence of states) */
    for(int t = T - 1; t >= 1; t--) q[t] = psi[t+1][q[t+1]];
}
