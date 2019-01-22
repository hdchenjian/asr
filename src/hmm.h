typedef struct {
    int N; // number of states
    int M; // number of symbols of the observation
    /*  A[1..N][1..N]. a[i][j] distribution of the transition probability of change
        from state i at time t to state j at time t + 1 */
    double	**A;
    /*  B[1..N][1..M]. b[j][k] is the probability distribution of observation k symbols in state j */
    double	**B; // probability distribution of observation k symbols in state j
    /*  pi[1..N] pi[i] is the initial distribution of states */
    double	*pi;  // initial distribution of states
} HMM;

void FreeHMM(HMM *hmm);
void InitHMM(HMM *hmm, int N, int M, int seed);
void save_hmm(char nome_file[], HMM *hmm);
void load_hmm(char nome_file[], HMM *hmm, int N, int M);

void Forward(HMM *hmm, int T, int *O, double **alpha, double *pprob);
void ForwardWithScale(HMM *hmm, int T, int *O, double **alpha, double *scale, double *pprob);
void Backward(HMM *hmm, int T, int *O, double **beta, double *pprob);
void BackwardWithScale(HMM *hmm, int T, int *O, double **beta, double *scale, double *pprob);
void BaumWelch(HMM *hmm, int T, int *O, double **alpha, double **beta,
               double **gamma, int *niter, double *plogprobinit, double *plogprobfinal);
void BaumWelch_C(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma);

double *** alloc_3d_matrix(int T, int N);
void free_3d_matrix(double *** xi, int T, int N);
void ComputeGamma(HMM *hmm, int T, double **alpha, double **beta, double **gamma);
void compute_epsilon(HMM* hmm, int T, int *O, double **alpha, double **beta, double ***xi);
void Viterbi(HMM *hmm, int T, int *O, int *path, double *pprob);
void ViterbiLog(HMM *hmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);
void ViterbiLog_C(HMM *hmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);
