#include <math.h>

typedef struct {
    /* N number of states;  Q={1,2,...,N} */
    int N;
    /* M number of symbols of the observation; V={1,2,...,M}*/
    int M;
    /*  A[1..N][1..N]. a[i][j] is the distribution of the transition probability of change
        from state i at time t to state j at time t + 1 */
    double	**A;
    /*  B[1..N][1..M]. b[j][k] is the probability distribution of observation k symbols in state j */
    double	**B;
    /*  pi[1..N] pi[i] is the initial distribution of states */
    double	*pi;
} HMM;

void FreeHMM(HMM *phmm);
void InitHMM_SC1(HMM *phmm, int N, int M, int seed);
void save_hmm(char nome_file[], HMM *phmm);
void load_hmm(char nome_file[], HMM *phmm, int N, int M);

/*  Functions belonging to the random number generator */
unsigned int hmmgetseed(void);
void hmmsetseed(unsigned int seed);
double hmmgetrand(void);

void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob);
void ForwardWithScale(HMM *phmm, int T, int *O, double **alpha, double *scale, double *pprob);
void Backward(HMM *phmm, int T, int *O, double **beta, double *pprob);
void BackwardWithScale(HMM *phmm, int T, int *O, double **beta, double *scale, double *pprob);
void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,
               double **gamma, int *niter, double *plogprobinit, double *plogprobfinal);

void BaumWelch_C(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma);

double *** alloc_3d_matrix(int T, int N);
void free_3d_matrix(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta, double **gamma);
void compute_epsilon(HMM* phmm, int T, int *O, double **alpha, double **beta, double ***xi);
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);
void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);

void ViterbiLog_C(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);
