#include <string.h>

#include "mfcc.h"
#include "hmm.h"
#include "nrutil.h"
#include "lbg_vq.h"

/*  Number of codevectors, allowed only power of 2 */
//#define NVQ 32
#define NVQ 64
//#define NVQ 128

/*
  If SET is 0 it compares with words not used in the training
  If SET is 1 it compares with words used in the training
  If SET is 2 it compares with words not used in the training contaminated with Gaussian noise
  If SET is 3 it compares with words used in the training contaminated with Gaussian noise */

/*  If SC is 1 it uses the model with serial constraint with single and double transitions
    If SC is 2 it uses the model with serial constraint with single transitions */
#define SC 1

int main()
{
    list_audio *lista = NULL;
    lista = add_audio("audio/00_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/00_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/00_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/00_01_pietro.r.wav", lista);
    lista = add_audio("audio/00_02_pietro.r.wav", lista);
    lista = add_audio("audio/00_03_pietro.r.wav", lista);
    lista = add_audio("audio/01_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/01_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/01_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/01_01_pietro.r.wav", lista);
    lista = add_audio("audio/01_02_pietro.r.wav", lista);
    lista = add_audio("audio/01_03_pietro.r.wav", lista);
    lista = add_audio("audio/02_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/02_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/02_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/02_01_pietro.r.wav", lista);
    lista = add_audio("audio/02_02_pietro.r.wav", lista);
    lista = add_audio("audio/02_03_pietro.r.wav", lista);
    lista = add_audio("audio/03_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/03_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/03_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/03_01_pietro.r.wav", lista);
    lista = add_audio("audio/03_02_pietro.r.wav", lista);
    lista = add_audio("audio/03_03_pietro.r.wav", lista);
    lista = add_audio("audio/04_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/04_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/04_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/04_01_pietro.r.wav", lista);
    lista = add_audio("audio/04_02_pietro.r.wav", lista);
    lista = add_audio("audio/04_03_pietro.r.wav", lista);
    lista = add_audio("audio/05_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/05_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/05_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/05_01_pietro.r.wav", lista);
    lista = add_audio("audio/05_02_pietro.r.wav", lista);
    lista = add_audio("audio/05_03_pietro.r.wav", lista);
    lista = add_audio("audio/06_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/06_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/06_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/06_01_pietro.r.wav", lista);
    lista = add_audio("audio/06_02_pietro.r.wav", lista);
    lista = add_audio("audio/06_03_pietro.r.wav", lista);
    lista = add_audio("audio/07_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/07_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/07_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/07_01_pietro.r.wav", lista);
    lista = add_audio("audio/07_02_pietro.r.wav", lista);
    lista = add_audio("audio/07_03_pietro.r.wav", lista);
    lista = add_audio("audio/08_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/08_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/08_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/08_01_pietro.r.wav", lista);
    lista = add_audio("audio/08_02_pietro.r.wav", lista);
    lista = add_audio("audio/08_03_pietro.r.wav", lista);
    lista = add_audio("audio/09_01_giuseppe.s.wav", lista);
    lista = add_audio("audio/09_02_giuseppe.s.wav", lista);
    lista = add_audio("audio/09_03_giuseppe.s.wav", lista);
    lista = add_audio("audio/09_01_pietro.r.wav", lista);
    lista = add_audio("audio/09_02_pietro.r.wav", lista);
    lista = add_audio("audio/09_03_pietro.r.wav", lista);

    int frame_num_total = 0;
    list_audio *current_audio =lista;
    int train_num = 0;
    while(current_audio != NULL) {
        frame_num_total += current_audio->frame_num;
        current_audio = current_audio->next;
        train_num += 1;
    }

    int lettere = 128;
    char **set = (char **)malloc(sizeof(char*) * (train_num + 1));
    if (set == NULL) print_error_exit("malloc failed");
    for(int i = 0; i < train_num + 1; i++) {
        set[i] = (char *)malloc(sizeof(char) * lettere);
        if (set[i] == NULL) print_error_exit("malloc failed");
    }
    int mfcc_coeff_num = N_COEFF_CEP*3;

    double **training = matrix_double(frame_num_total, mfcc_coeff_num);
    double **mfcc_feature = NULL;
    current_audio = lista;
	int train_set_index = 0, feature_index = 0;
    while(current_audio != NULL) {
        mfcc_feature = matrix_double(current_audio->frame_num, mfcc_coeff_num);
        load_matrix_double(current_audio->filename_car, mfcc_feature, current_audio->frame_num, mfcc_coeff_num);
        for(int n = 0; n < current_audio->frame_num; n++){
            for(int k=0; k< mfcc_coeff_num; k++) {
                training[feature_index + n][k] = mfcc_feature[n][k];
            }
        }
        strcpy(set[train_set_index],current_audio->filename);
        //printf("%s, frame_num: %d\n", current_audio->filename, current_audio->frame_num);
        train_set_index += 1;
        feature_index += current_audio->frame_num;
        current_audio = current_audio->next;
    }
    printf("train_num %d, mfcc_coeff_num %d, frame_num_total %d\n", train_num, mfcc_coeff_num, frame_num_total);


    double **codebook = lbg(training, frame_num_total, mfcc_coeff_num, NVQ);
    exit(-1);

    FILE *codebook_file = fopen("codebook.txt","w");
    for(int n = 0; n < NVQ; n++)
    {
        fprintf(codebook_file,"%d\n", n);
        for(int k = 0; k<mfcc_coeff_num; k++) fprintf(codebook_file,"%f\t", codebook[n][k]);
        fprintf(codebook_file,"\n");
    }

    current_audio = lista;
    char *nome = NULL;
    int *sequenza = NULL;
    char *pos = NULL;

    FILE *tutte_sequenze = NULL;
    tutte_sequenze = fopen("tutte_sequenze.txt","w");

    int N_vq = NVQ + 1;
    /*  After the codebook, it determines the vector quantization for every vector of the training set */
    while(current_audio != NULL)
    {
        mfcc_feature = matrix_double(current_audio->frame_num, mfcc_coeff_num-1);
        load_matrix_double(current_audio->filename_car, mfcc_feature, current_audio->frame_num, mfcc_coeff_num-1);

        /*  The quantized vector is the vector of sequences that will be used in the Hidden Markov Models */
        sequenza = vq(codebook, mfcc_feature, current_audio->frame_num, mfcc_coeff_num, N_vq);
        nome = strdup(current_audio->filename);
        /*  Find the last character '.' */
        pos = strrchr(nome, '.' );
        int i = pos-nome;
        /*  After '.' changes the extension */
        nome[i+1] = 's';
        nome[i+2] = 'e';
        nome[i+3] = 'q';

        fprintf(tutte_sequenze,"%s\n", current_audio->filename);
        for(int n = 1; n<=current_audio->frame_num; n++)
        {
            fprintf(tutte_sequenze, "%d\t", sequenza[n]);
        }
        fprintf(tutte_sequenze, "\n\n");


        strcpy(current_audio->filename_seq, nome);
        save_array_int(current_audio->filename_seq, sequenza, current_audio->frame_num +1);

        current_audio = current_audio->next;
    }

    /*  Set the seed of the random. */
    unsigned int seed = 333;

    current_audio = lista;

    HMM hmm;
    /*  In the Hidden Markov Models:
        - M is the number of the observed symbols
        - N is the number of states */
    int M = NVQ;
    int N = 5;

    double	logprobinit, logprobfinal;
    double 	**alpha;
    double	**beta;
    double	**gamma;

    int T, niter, T_max = 0;

    while(current_audio != NULL) {
        nome = strdup(current_audio->filename);
        pos = strrchr(nome, '.' );
        int i = pos-nome;
        nome[i+1] = 'h';
        nome[i+2] = 'm';
        nome[i+3] = 'm';

        strcpy(current_audio->filename_hmm, nome);
        int *O = array_int(current_audio->frame_num +1);
        load_array_int(current_audio->filename_seq, O, current_audio->frame_num +1);
        T = current_audio->frame_num;
        if(T_max<T) T_max = T;

#if SC==1
        InitHMM_SC1(&hmm, N, M, seed);
#elif SC==2
        InitHMM_SC2(&hmm, N, M, seed);
#endif

        alpha = dmatrix(1, T, 1, hmm.N);
        beta = dmatrix(1, T, 1, hmm.N);
        gamma = dmatrix(1, T, 1, hmm.N);

        BaumWelch_C(&hmm, T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal);

        Salva_HMM(current_audio->filename_hmm, &hmm);

        current_audio = current_audio->next;
    }

    double 	proba[train_num];
    int	*q;
    double **delta;
    int	**psi;
    int frequenza_campionamento, campioni_segnale, n_finestre;
    double *segnale = NULL;
    char **analizza = NULL;
    int parole  =  40;
    analizza = (char **) malloc(sizeof(char*) * parole);
    if (analizza == NULL) print_error_exit("malloc failed");
    for(int i = 0; i<parole; i++) {
        analizza[i]  =  (char *) malloc(sizeof(char) * lettere);
        if (analizza[i]==NULL)
            print_error_exit("malloc failed");
    }

    analizza[0] = "audio/00_04_giuseppe.s.wav";
    analizza[1] = "audio/00_05_giuseppe.s.wav";
    analizza[2] = "audio/00_04_pietro.r.wav";
    analizza[3] = "audio/00_05_pietro.r.wav";
    analizza[4] = "audio/01_04_giuseppe.s.wav";
    analizza[5] = "audio/01_05_giuseppe.s.wav";
    analizza[6] = "audio/01_04_pietro.r.wav";
    analizza[7] = "audio/01_05_pietro.r.wav";
    analizza[8] = "audio/02_04_giuseppe.s.wav";
    analizza[9] = "audio/02_05_giuseppe.s.wav";
    analizza[10] = "audio/02_04_pietro.r.wav";
    analizza[11] = "audio/02_05_pietro.r.wav";
    analizza[12] = "audio/03_04_giuseppe.s.wav";
    analizza[13] = "audio/03_05_giuseppe.s.wav";
    analizza[14] = "audio/03_04_pietro.r.wav";
    analizza[15] = "audio/03_05_pietro.r.wav";
    analizza[16] = "audio/04_04_giuseppe.s.wav";
    analizza[17] = "audio/04_05_giuseppe.s.wav";
    analizza[18] = "audio/04_04_pietro.r.wav";
    analizza[19] = "audio/04_05_pietro.r.wav";
    analizza[20] = "audio/05_04_giuseppe.s.wav";
    analizza[21] = "audio/05_05_giuseppe.s.wav";
    analizza[22] = "audio/05_04_pietro.r.wav";
    analizza[23] = "audio/05_05_pietro.r.wav";
    analizza[24] = "audio/06_04_giuseppe.s.wav";
    analizza[25] = "audio/06_05_giuseppe.s.wav";
    analizza[26] = "audio/06_04_pietro.r.wav";
    analizza[27] = "audio/06_05_pietro.r.wav";
    analizza[28] = "audio/07_04_giuseppe.s.wav";
    analizza[29] = "audio/07_05_giuseppe.s.wav";
    analizza[30] = "audio/07_04_pietro.r.wav";
    analizza[31] = "audio/07_05_pietro.r.wav";
    analizza[32] = "audio/08_04_giuseppe.s.wav";
    analizza[33] = "audio/08_05_giuseppe.s.wav";
    analizza[34] = "audio/08_04_pietro.r.wav";
    analizza[35] = "audio/08_05_pietro.r.wav";
    analizza[36] = "audio/09_04_giuseppe.s.wav";
    analizza[37] = "audio/09_05_giuseppe.s.wav";
    analizza[38] = "audio/09_04_pietro.r.wav";
    analizza[39] = "audio/09_05_pietro.r.wav";


    FILE *risultati = fopen("risultati.txt","w");
    for(int j = 0; j<parole; j++)
    {
        current_audio = lista;

        frequenza_campionamento = campioni_segnale = n_finestre = 0;

        segnale = open_wav(analizza[j], &frequenza_campionamento, &campioni_segnale);
        mfcc_feature = extract_mfcc(segnale, frequenza_campionamento, campioni_segnale,
        		FREQUENZA_MAX, DIM_FINESTRA, DIM_PASSO, N_COEFF_CEP, &n_finestre);

        sequenza = vq(codebook, mfcc_feature, n_finestre, mfcc_coeff_num, N_vq);

        T = n_finestre;

        q = ivector(1,T);
        delta = dmatrix(1, T, 1, N);
        psi = imatrix(1, T, 1, N);

        int i = 0;
        while(current_audio != NULL)
        {
            proba[i] = 0;
            Carica_HMM(current_audio->filename_hmm, &hmm, N, M);

            ViterbiLog_C(&hmm, T, sequenza, delta, psi, q, &proba[i]);

            i++;
            current_audio = current_audio->next;
        }

        double alt_proba;
        double basso_proba;
        int val_alto, val_basso;

        alt_proba = basso_proba = proba[0];
        val_alto = val_basso = 0;

        for(i = 0; i<train_num; i++) {
            if(proba[i]>alt_proba) {
                alt_proba = proba[i];
                val_alto = i;
            }
            if(proba[i]<basso_proba) {
                basso_proba = proba[i];
                val_basso = i;
            }
        }
        fprintf(risultati,"%s\tlo associo a\t%s\n", analizza[j], set[val_alto]);
    }
    return 0;
}
