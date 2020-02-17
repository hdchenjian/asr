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
void BaumWelch(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma);
void BaumWelch_with_scale(HMM *hmm, int frame_num, int *O, double **alpha, double **beta, double **gamma);
void Viterbi(HMM *hmm, int T, int *O, int *path, double *pprob);
