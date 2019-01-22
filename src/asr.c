#include <string.h>
#include <time.h>

#include "lista.h"
#include "mfcc.h"
#include "hmm.h"
#include "nrutil.h"
#include "lbg_vq.h"

/*  Number of codevectors, allowed only power of 2 */
//#define NVQ 32
#define VQ_NUM 64
//#define NVQ 128

/*
  If SET is 0 it compares with words not used in the training
  If SET is 1 it compares with words used in the training
  If SET is 2 it compares with words not used in the training contaminated with Gaussian noise
  If SET is 3 it compares with words used in the training contaminated with Gaussian noise */

int main()
{
    list_audio *train_list = 0;
    add_audio("audio/00_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/00_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/00_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/00_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/00_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/00_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/01_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/01_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/01_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/01_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/01_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/01_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/02_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/02_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/02_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/02_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/02_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/02_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/03_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/03_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/03_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/03_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/03_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/03_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/04_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/04_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/04_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/04_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/04_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/04_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/05_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/05_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/05_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/05_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/05_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/05_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/06_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/06_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/06_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/06_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/06_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/06_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/07_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/07_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/07_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/07_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/07_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/07_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/08_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/08_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/08_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/08_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/08_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/08_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/09_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/09_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/09_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/09_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/09_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/09_03_pietro.r.wav", train_list);

    int frame_num_total = 0;
    list_audio *current_audio =train_list;
    int train_num = 0;
    while(current_audio != NULL) {
        frame_num_total += current_audio->frame_num;
        current_audio = current_audio->next;
        train_num += 1;
    }

    int mfcc_coeff_num = MFCC_COEFF_NUM*3;
    double **training = matrix_double(frame_num_total, mfcc_coeff_num);
    current_audio = train_list;
	int train_set_index = 0, feature_index = 0;
    while(current_audio != NULL) {
        for(int n = 0; n < current_audio->frame_num; n++){
            for(int k=0; k< mfcc_coeff_num; k++) {
                training[feature_index + n][k] = current_audio->mfcc[n][k];
            }
        }
        //printf("%s, frame_num: %d\n", current_audio->filename, current_audio->frame_num);
        train_set_index += 1;
        feature_index += current_audio->frame_num;
        current_audio = current_audio->next;
    }
    printf("train_num %d, mfcc_coeff_num %d, frame_num_total %d\n", train_num, mfcc_coeff_num, frame_num_total);

    double **codebook = lbg_get_codebook(training, frame_num_total, mfcc_coeff_num, VQ_NUM);
    free_matrix_double(training, frame_num_total);
    printf("lbg_get_codebook, M %d, K %d, N %d\n", frame_num_total, mfcc_coeff_num, VQ_NUM);

    current_audio = train_list;
    while(current_audio != NULL) {
        int *vq_result = lbg_encode(codebook, current_audio->mfcc, current_audio->frame_num, mfcc_coeff_num, VQ_NUM);
        current_audio->vq_result = vq_result;
        current_audio = current_audio->next;
    }

    HMM hmm;
    int observed_num = VQ_NUM;
    int status_num = 5;
    current_audio = train_list;
    while(current_audio != NULL) {
        strcpy(current_audio->filename_hmm, current_audio->filename);
        char *pos = strrchr(current_audio->filename_hmm, '.' );
        int i = pos - current_audio->filename_hmm;
        current_audio->filename_hmm[i+1] = 'h';
        current_audio->filename_hmm[i+2] = 'm';
        current_audio->filename_hmm[i+3] = 'm';

        InitHMM(&hmm, status_num, observed_num, time(0));
        double **alpha = dmatrix(1, current_audio->frame_num, 1, hmm.N);
        double **beta = dmatrix(1, current_audio->frame_num, 1, hmm.N);
        double **gamma = dmatrix(1, current_audio->frame_num, 1, hmm.N);
        BaumWelch_C(&hmm, current_audio->frame_num, current_audio->vq_result, alpha, beta, gamma);
        save_hmm(current_audio->filename_hmm, &hmm);
        free_dmatrix(alpha, 1, current_audio->frame_num, 1, hmm.N);
        free_dmatrix(beta, 1, current_audio->frame_num, 1, hmm.N);
        free_dmatrix(gamma, 1, current_audio->frame_num, 1, hmm.N);
        FreeHMM(&hmm);
        current_audio = current_audio->next;
    }
    printf("train hmm over\n");

    int test_num  =  40;
    char **test_list = (char **)malloc(sizeof(char*) * test_num);
    for(int i = 0; i<test_num; i++) test_list[i]  =  (char *)malloc(sizeof(char) * 128);

    test_list[0] = "audio/00_04_giuseppe.s.wav";
    test_list[1] = "audio/00_05_giuseppe.s.wav";
    test_list[2] = "audio/00_04_pietro.r.wav";
    test_list[3] = "audio/00_05_pietro.r.wav";
    test_list[4] = "audio/01_04_giuseppe.s.wav";
    test_list[5] = "audio/01_05_giuseppe.s.wav";
    test_list[6] = "audio/01_04_pietro.r.wav";
    test_list[7] = "audio/01_05_pietro.r.wav";
    test_list[8] = "audio/02_04_giuseppe.s.wav";
    test_list[9] = "audio/02_05_giuseppe.s.wav";
    test_list[10] = "audio/02_04_pietro.r.wav";
    test_list[11] = "audio/02_05_pietro.r.wav";
    test_list[12] = "audio/03_04_giuseppe.s.wav";
    test_list[13] = "audio/03_05_giuseppe.s.wav";
    test_list[14] = "audio/03_04_pietro.r.wav";
    test_list[15] = "audio/03_05_pietro.r.wav";
    test_list[16] = "audio/04_04_giuseppe.s.wav";
    test_list[17] = "audio/04_05_giuseppe.s.wav";
    test_list[18] = "audio/04_04_pietro.r.wav";
    test_list[19] = "audio/04_05_pietro.r.wav";
    test_list[20] = "audio/05_04_giuseppe.s.wav";
    test_list[21] = "audio/05_05_giuseppe.s.wav";
    test_list[22] = "audio/05_04_pietro.r.wav";
    test_list[23] = "audio/05_05_pietro.r.wav";
    test_list[24] = "audio/06_04_giuseppe.s.wav";
    test_list[25] = "audio/06_05_giuseppe.s.wav";
    test_list[26] = "audio/06_04_pietro.r.wav";
    test_list[27] = "audio/06_05_pietro.r.wav";
    test_list[28] = "audio/07_04_giuseppe.s.wav";
    test_list[29] = "audio/07_05_giuseppe.s.wav";
    test_list[30] = "audio/07_04_pietro.r.wav";
    test_list[31] = "audio/07_05_pietro.r.wav";
    test_list[32] = "audio/08_04_giuseppe.s.wav";
    test_list[33] = "audio/08_05_giuseppe.s.wav";
    test_list[34] = "audio/08_04_pietro.r.wav";
    test_list[35] = "audio/08_05_pietro.r.wav";
    test_list[36] = "audio/09_04_giuseppe.s.wav";
    test_list[37] = "audio/09_05_giuseppe.s.wav";
    test_list[38] = "audio/09_04_pietro.r.wav";
    test_list[39] = "audio/09_05_pietro.r.wav";

    double probability[train_num];
    for(int j = 0; j < test_num; j++) {
        current_audio = train_list;
        int sample_rate = 0;
        int audio_data_num = 0;
        int frame_num = 0;
        double *signal = open_wav(test_list[j], &sample_rate, &audio_data_num);
        double **mfcc_feature = extract_mfcc(signal, sample_rate, audio_data_num,
        		FREQUENZA_MAX, FRAME_LENGTH, FRAME_SHIFT, MFCC_COEFF_NUM, &frame_num);
        free(signal);
        int *audio_vq = lbg_encode(codebook, mfcc_feature, frame_num, mfcc_coeff_num, VQ_NUM);
        //printf("audio_vq %p\n", audio_vq);

        free_matrix_double(mfcc_feature, frame_num);
        printf("test %s\n", test_list[j]);

        int index = 0;
        int	*path = (int *)calloc(frame_num, sizeof(int));
        while(current_audio != NULL) {
            load_hmm(current_audio->filename_hmm, &hmm, status_num, observed_num);
            Viterbi(&hmm, frame_num, audio_vq, path, &probability[index]);
            index++;
            current_audio = current_audio->next;
        }
        free(path);
        //printf("audio_vq %p\n", audio_vq);
        free(audio_vq);


        double max_prob = -1;
        double min_prob = -1;
        int max_index, min_index;
        max_index = - 1;
        min_index = -1;
        for(int i = 0; i < train_num; i++) {
            if(probability[i] > max_prob) {
                max_prob = probability[i];
                max_index = i;
            }
            if(probability[i] < min_prob) {
                min_prob = probability[i];
                min_index = i;
            }
        }
        printf("%s %d %d\n", test_list[j], max_index, min_index);
    }
    for(int i = 0; i<test_num; i++) free(test_list[i]);
    free(test_list);
    free_matrix_double(codebook, VQ_NUM);
    return 0;
}
