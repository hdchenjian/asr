/* **    Purpose: Changed the convergence criterion from ratio to absolute value.
   Solution to the problem of learning: how to adjust the parameters of the hidden Markov model
   lambda = (A, B, pi), so as to maximize P(O|lambda)
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "hmm.h"
#include "matrix.h"

#define DELTA 0.001

void FreeHMM(HMM *hmm)
{
    free_dmatrix(hmm->A, hmm->N);
    free_dmatrix(hmm->B, hmm->N);
    free(hmm->pi);
}

void InitHMM(HMM *hmm, int N, int M, int seed)
{
    int i, j, k;
    double sum;
    srand(seed);
    hmm->M = M;
    hmm->N = N;
    hmm->A = (double **) dmatrix(hmm->N, hmm->N);
    for(i = 1; i <= hmm->N; i++) {
        sum = 0.0;
        for(j = 1; j <= hmm->N; j++) {
            //if(j < i || j >= i+3) hmm->A[i][j] = 0.0;
            hmm->A[i][j] = (double)rand() / RAND_MAX;
            sum += hmm->A[i][j];
        }
        for(j = 1; j <= hmm->N; j++) hmm->A[i][j] /= sum;
    }

    hmm->B = (double **) dmatrix(hmm->N, hmm->M);
    for(j = 1; j <= hmm->N; j++) {
        sum = 0.0;
        for(k = 1; k <= hmm->M; k++) {
            hmm->B[j][k] = (double)rand() / RAND_MAX;
            sum += hmm->B[j][k];
        }
        for(k = 1; k <= hmm->M; k++) hmm->B[j][k] /= sum;
    }

    hmm->pi = calloc(hmm->N + 1, sizeof(double));
    sum = 0.0;
    for(i = 1; i <= hmm->N; i++) {
        hmm->pi[i] = (double)rand() / RAND_MAX;
        sum += hmm->pi[i];
    }
    for(i = 1; i <= hmm->N; i++) hmm->pi[i] /= sum;
    //for(i = 1; i <= hmm->N; i++) printf("pi %d %f\n", i, hmm->pi[i]);
}

/** Solution to the problem of decoding:
    given the sequence of observations O = O1O2 ... OT and the Hidden Markovian Model lambda = (A, B, pi),
    how to choose a corresponding state sequence Q = q1 q2 ... qT that is optimal, to justify the better observations */
/* delta psi: frame_num * status_num, q:  1 * frame_num. */
void Viterbi(HMM *hmm, int T, int *O, int *path, double *pprob)
{
    double **delta = dmatrix(T, hmm->N);
    int    **psi = imatrix(T, hmm->N);
    for(int i = 1; i <= hmm->N; i++) {
        delta[1][i] = hmm->pi[i] * (hmm->B[i][O[0]]);
        psi[1][i] = 0;
        printf("pi %d %f\n", i, hmm->pi[i]);
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
            printf("max delta %d %d %f %d %f\n", t, j, maxval, O[t - 1], hmm->B[j][O[t - 1]]);

            delta[t][j] = maxval*(hmm->B[j][O[t - 1]]);
            psi[t][j] = maxvalind;
        }
    }

    *pprob = 0.0;
    path[T - 1] = 1;
    for(int i = 1; i <= hmm->N; i++) {
        printf("delta %d %f\n", i, delta[T][i]);
        if (delta[T][i] > *pprob) {
            *pprob = delta[T][i];
            path[T - 1] = i;
        }
    }

    /*  4. Backtracking on the path (sequence of states) q*_t = psi_t+1 (q*_t+1), for t=T-1, T-2, ..., 1*/
    for(int t = T - 1; t >= 1; t--) path[t - 1] = psi[t+1][path[t]];
    free_dmatrix(delta, T);
    free_imatrix(psi, T);
}

void Backward_with_scale(HMM *hmm, int T, int *O, double **beta, double *scale)
{
    for(int i = 1; i <= hmm->N; i++) beta[T][i] = 1.0 / scale[T];

    /*  for t=T-1,T-2...1, with 1 <= i <= N
        beta_t(i) = {sum a_ij*b_j(Ot+1)}*beta_t+1(j) */
    for(int t = T - 1; t >= 1; t--) {
        for(int i = 1; i <= hmm->N; i++) {
            double sum = 0.0;
            for(int j = 1; j <= hmm->N; j++) sum += hmm->A[i][j] * (hmm->B[j][O[t]]) * beta[t+1][j];
            beta[t][i] = sum / scale[t];
        }
    }
}

/*  Where beta is the auxiliary variable backward, defined as the probability of
    partial sequence of observations from t+1 to the end, given the state Si at time t,
    and the model lambda=(A, B, pi). beta: frame_num * N, beta_t(i)=P(Ot+1 Ot+2...OT, q_t=S_i|lambda) */
void Backward(HMM *hmm, int T, int *O, double **beta)
{
    for(int i = 1; i <= hmm->N; i++) beta[T][i] = 1.0;

    /*  for t=T-1,T-2...1, with 1 <= i <= N
        beta_t(i) = {sum a_ij*b_j(Ot+1)}*beta_t+1(j) */
    for(int t = T - 1; t >= 1; t--) {
        for(int i = 1; i <= hmm->N; i++) {
            double sum = 0.0;
            for(int j = 1; j <= hmm->N; j++) sum += hmm->A[i][j] * (hmm->B[j][O[t]]) * beta[t+1][j];
            beta[t][i] = sum;
        }
    }
}

void Forward_with_scale(HMM *hmm, int T, int *O, double **alpha, double *scale, double *pprob)
{
    /*  alpha_1(i) = pi_i * b_i(O1), with 1 <= i <= N */
    scale[1] = 0.0;
    for(int i = 1; i <= hmm->N; i++){
        alpha[1][i] = hmm->pi[i]* hmm->B[i][O[0]];
        scale[1] += alpha[1][i];
    }
    for(int i = 1; i <= hmm->N; i++) alpha[1][i] /= scale[1];

    /*  for t=1,2...T-1, with 1 <= j <= N
        alpha_t+1(j) = {sum alpha_t(i)*a_ij}*b_j(Ot+1) */
    for(int t = 1; t <= T - 1; t++) {
        scale[t+1] = 0.0;
        for(int j = 1; j <= hmm->N; j++) {
            double sum = 0;
            for(int i = 1; i <= hmm->N; i++) sum += alpha[t][i]* (hmm->A[i][j]);
            alpha[t+1][j] = sum*(hmm->B[j][O[t]]);
            scale[t+1] += alpha[t+1][j];
        }
        for(int j = 1; j <= hmm->N; j++) alpha[t+1][j] /= scale[t+1];
    }

    *pprob = 0.0;
    for (int t = 1; t <= T; t++) *pprob += log(scale[t]);
}

/*  Given the sequence of observations O=O1,O2,...,OT and the hidden Markov model
    lambda=(A, B, pi), how to can be calculated efficiently P(O|lambda)
    Where alpha is the auxiliary variable forward, defined as the probability of
    partial sequence of observations up to time t and the state Si,
    alpha: frame_num * N, alpha_t(i)=P(O1 O2...Ot, q_t=S_i|lambda) */
void Forward(HMM *hmm, int T, int *O, double **alpha, double *pprob)
{
    /*  alpha_1(i) = pi_i * b_i(O1), with 1 <= i <= N */
    for(int i = 1; i <= hmm->N; i++) alpha[1][i] = hmm->pi[i]* hmm->B[i][O[0]];


    /*  for t=1,2...T-1, with 1 <= j <= N
        alpha_t+1(j) = {sum alpha_t(i)*a_ij}*b_j(Ot+1) */
    for(int t = 1; t <= T - 1; t++) {
        for(int j = 1; j <= hmm->N; j++) {
            double sum = 0;
            for(int i = 1; i <= hmm->N; i++) sum += alpha[t][i] * (hmm->A[i][j]);
            alpha[t+1][j] = sum*(hmm->B[j][O[t]]);
        }
    }

    *pprob = 0.0;
    for(int i = 1; i <= hmm->N; i++) *pprob += alpha[T][i];
}

/* gamma_t(i) = { alpha_t(i) * beta_t(i) } / P(O|lambda) */
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

/* xi_t(i,j) = { alpha_t(i) * beta_t+1(j) * a_ij * b_j } / P(O|lambda) */
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
    for(int t = 0; t < C; t++) x[t] = dmatrix(N, N);
    return x;
}

void free_3d_matrix(double ***x, int C, int N)
{
    for(int t = 0; t < C; t++) free_dmatrix(x[t], N);
    free(x);
}

void BaumWelch(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma)
{
    double observe_prob;
    //double *scale = malloc((frame_num + 1) * sizeof(double));
    Forward(hmm, frame_num, O, alpha, &observe_prob);
    Backward(hmm, frame_num, O, beta);
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
                double tmp = 0.0;
                for(int t = 1; t <= frame_num; t++) {
                    if (O[t - 1] == k) tmp += gamma[t][i];
                }
                hmm->B[i][k] = 0.001f + 0.999f * tmp / sum_gamma;
            }
        }

        Forward(hmm, frame_num, O, alpha, &observe_prob);
        Backward(hmm, frame_num, O, beta);
        ComputeGamma(hmm, frame_num, alpha, beta, gamma);
        compute_epsilon(hmm, frame_num, O, alpha, beta, epsilon);

        /* Calculates the difference between the logarithmic probability of two iterations */
        delta = log(observe_prob) - log(observe_prob_previous);
        observe_prob_previous = observe_prob;
    } while(delta > DELTA);
    free_3d_matrix(epsilon, frame_num, hmm->N);
    //free(scale);
}

void BaumWelch_with_scale(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma)
{
    double observe_prob;
    double *scale = malloc((frame_num + 1) * sizeof(double));
    Forward_with_scale(hmm, frame_num, O, alpha, scale, &observe_prob);
    Backward_with_scale(hmm, frame_num, O, beta, scale);
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
                double tmp = 0.0;
                for(int t = 1; t <= frame_num; t++) {
                    if (O[t - 1] == k) tmp += gamma[t][i];
                }
                hmm->B[i][k] = 0.001f + 0.999f * tmp / sum_gamma;
            }
        }

        Forward_with_scale(hmm, frame_num, O, alpha, scale, &observe_prob);
        Backward_with_scale(hmm, frame_num, O, beta, scale);
        ComputeGamma(hmm, frame_num, alpha, beta, gamma);
        compute_epsilon(hmm, frame_num, O, alpha, beta, epsilon);

        /* Calculates the difference between the logarithmic probability of two iterations */
        delta = log(observe_prob) - log(observe_prob_previous);
        observe_prob_previous = observe_prob;
    } while(delta > DELTA);
    free_3d_matrix(epsilon, frame_num, hmm->N);
    free(scale);
}
