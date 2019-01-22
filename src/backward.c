#include "hmm.h"

/*  Where beta is the auxiliary variable backward, defined as the probability of
    partial sequence of observations from t+1 to the end, given the state Si at time t,
    and the model lambda=(A, B, pi). beta: frame_num * N, beta_t(i)=P(Ot+1 Ot+2...OT, q_t=S_i|lambda) */

void Backward(HMM *hmm, int T, int *O, double **beta, double *pprob)
{
    /*  beta_T(i) = 1, with 1 <= i <= N */
    for(int i = 1; i <= hmm->N; i++) beta[T][i] = 1.0;

    /*  2. Induction */
    /*  for t=T-1,T-2...1, with 1 <= i <= N
        beta_t+1(j) = {sum a_ij*b_j(Ot+1)}*beta_t+1(j) */
    for(int t = T - 1; t >= 1; t--) {
        for(int i = 1; i <= hmm->N; i++) {
            double sum = 0.0;
            for(int j = 1; j <= hmm->N; j++) sum += hmm->A[i][j] * (hmm->B[j][O[t]]) * beta[t+1][j];
            beta[t][i] = sum;
        }
    }

    //*pprob = 0.0;
    //for(int i = 1; i <= hmm->N; i++) *pprob += beta[1][i];
}

void BackwardWithScale(HMM *hmm, int T, int *O, double **beta,	double *scale, double *pprob)
{
    int i, j;
    int t;
	double sum;

    /* 1. Initialization */
    for(i = 1; i <= hmm->N; i++)
        beta[T][i] = 1.0/scale[T];

    /* 2. Induction */
    for(t = T - 1; t >= 1; t--)
    {
        for(i = 1; i <= hmm->N; i++)
        {
			sum = 0.0;
            for(j = 1; j <= hmm->N; j++)
                sum += hmm->A[i][j] * (hmm->B[j][O[t+1]])*beta[t+1][j];
            beta[t][i] = sum/scale[t];

        }
    }

    /* 3. Termination
        This phase in the Backward there is not.  */
}
