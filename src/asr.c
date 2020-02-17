#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "mfcc.h"
#include "hmm.h"
#include "lbg_vq.h"
#include "matrix.h"
#include "data.h"
#include "utils.h"

/*  Number of codevectors, allowed only power of 2 */
#define VQ_NUM 64
#define CLASS_NUM 10

int main()
{
    list_audio *train_list = 0;
    train_list = add_audio("audio/train/00_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/00_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/00_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/00_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/00_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/00_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/01_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/01_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/01_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/01_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/01_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/01_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/02_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/02_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/02_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/02_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/02_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/02_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/03_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/03_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/03_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/03_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/03_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/03_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/04_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/04_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/04_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/04_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/04_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/04_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/05_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/05_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/05_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/05_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/05_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/05_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/06_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/06_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/06_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/06_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/06_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/06_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/07_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/07_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/07_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/07_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/07_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/07_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/08_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/08_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/08_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/08_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/08_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/08_03_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/09_01_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/09_02_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/09_03_giuseppe.s.wav", train_list);
    train_list = add_audio("audio/train/09_01_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/09_02_pietro.r.wav", train_list);
    train_list = add_audio("audio/train/09_03_pietro.r.wav", train_list);
    printf("load train data over\n");

    int frame_num_total = 0;
    list_audio *current_audio = train_list;
    int train_num = 0;
    while(current_audio != NULL) {
        frame_num_total += current_audio->frame_num;
        current_audio = current_audio->next;
        train_num += 1;
    }
    double **training = matrix_double(frame_num_total, MFCC_COEFF_ALL);
    current_audio = train_list;
    int feature_index = 0;
    while(current_audio != NULL) {
        for(int n = 0; n < current_audio->frame_num; n++){
            for(int k=0; k< MFCC_COEFF_ALL; k++) {
                training[feature_index + n][k] = current_audio->feature[n][k];
            }
        }
        feature_index += current_audio->frame_num;
        current_audio = current_audio->next;
    }
    double **codebook = lbg_get_codebook(training, frame_num_total, MFCC_COEFF_ALL, VQ_NUM);
    free_matrix_double(training, frame_num_total);
    printf("lbg_get_codebook, M %d, K %d, N %d, train_num %d\n", frame_num_total, MFCC_COEFF_ALL, VQ_NUM, train_num);

    current_audio = train_list;
    int max_frame_num = 0;
    char **labels = get_labels("audio/train_label.txt");
    while(current_audio != NULL) {
        if(current_audio->frame_num > max_frame_num) max_frame_num = current_audio->frame_num;
        int *vq_result = lbg_encode(codebook, current_audio->feature, current_audio->frame_num, MFCC_COEFF_ALL, VQ_NUM);
        current_audio->vq_result = vq_result;
        for(int i = 0; i < CLASS_NUM; ++i){
            if(strstr(current_audio->filename, labels[i])){
                current_audio->class = i;
                break;
            }
        }
        //printf("class %s %d\n", current_audio->filename, current_audio->class);
        current_audio = current_audio->next;
    }

    int observed_num = VQ_NUM;
    int status_num = 5;
    HMM hmm[CLASS_NUM];
    for(int k=0; k< CLASS_NUM; k++) InitHMM(&hmm[k], status_num, observed_num, time(0) + k*77);
    current_audio = train_list;
    double **alpha = dmatrix(max_frame_num, status_num);
    double **beta = dmatrix(max_frame_num, status_num);
    double **gamma = dmatrix(max_frame_num, status_num);
    current_audio = train_list;

    while(current_audio != NULL) {
    	BaumWelch_with_scale(&hmm[current_audio->class], current_audio->frame_num, current_audio->vq_result, alpha, beta, gamma);
        break;
        current_audio = current_audio->next;
    }
    free_dmatrix(alpha, max_frame_num);
    free_dmatrix(beta, max_frame_num);
    free_dmatrix(gamma, max_frame_num);
    printf("train hmm over\n");

    int test_num  = 40;
    char **test_list = (char **)malloc(sizeof(char*) * test_num);
    //test_list[0] = "audio/test/00_04_giuseppe.s.wav";
    test_list[0] = "audio/train/00_01_giuseppe.s.wav";
    test_list[1] = "audio/test/00_05_giuseppe.s.wav";
    test_list[2] = "audio/test/00_04_pietro.r.wav";
    test_list[3] = "audio/test/00_05_pietro.r.wav";
    test_list[4] = "audio/test/01_04_giuseppe.s.wav";
    test_list[5] = "audio/test/01_05_giuseppe.s.wav";
    test_list[6] = "audio/test/01_04_pietro.r.wav";
    test_list[7] = "audio/test/01_05_pietro.r.wav";
    test_list[8] = "audio/test/02_04_giuseppe.s.wav";
    test_list[9] = "audio/test/02_05_giuseppe.s.wav";
    test_list[10] = "audio/test/02_04_pietro.r.wav";
    test_list[11] = "audio/test/02_05_pietro.r.wav";
    test_list[12] = "audio/test/03_04_giuseppe.s.wav";
    test_list[13] = "audio/test/03_05_giuseppe.s.wav";
    test_list[14] = "audio/test/03_04_pietro.r.wav";
    test_list[15] = "audio/test/03_05_pietro.r.wav";
    test_list[16] = "audio/test/04_04_giuseppe.s.wav";
    test_list[17] = "audio/test/04_05_giuseppe.s.wav";
    test_list[18] = "audio/test/04_04_pietro.r.wav";
    test_list[19] = "audio/test/04_05_pietro.r.wav";
    test_list[20] = "audio/test/05_04_giuseppe.s.wav";
    test_list[21] = "audio/test/05_05_giuseppe.s.wav";
    test_list[22] = "audio/test/05_04_pietro.r.wav";
    test_list[23] = "audio/test/05_05_pietro.r.wav";
    test_list[24] = "audio/test/06_04_giuseppe.s.wav";
    test_list[25] = "audio/test/06_05_giuseppe.s.wav";
    test_list[26] = "audio/test/06_04_pietro.r.wav";
    test_list[27] = "audio/test/06_05_pietro.r.wav";
    test_list[28] = "audio/test/07_04_giuseppe.s.wav";
    test_list[29] = "audio/test/07_05_giuseppe.s.wav";
    test_list[30] = "audio/test/07_04_pietro.r.wav";
    test_list[31] = "audio/test/07_05_pietro.r.wav";
    test_list[32] = "audio/test/08_04_giuseppe.s.wav";
    test_list[33] = "audio/test/08_05_giuseppe.s.wav";
    test_list[34] = "audio/test/08_04_pietro.r.wav";
    test_list[35] = "audio/test/08_05_pietro.r.wav";
    test_list[36] = "audio/test/09_04_giuseppe.s.wav";
    test_list[37] = "audio/test/09_05_giuseppe.s.wav";
    test_list[38] = "audio/test/09_04_pietro.r.wav";
    test_list[39] = "audio/test/09_05_pietro.r.wav";

    char **test_labels = get_labels("audio/test_label.txt");
    for(int j = 0; j < test_num; j++) {
        int test_class = -1;
        for(int i = 0; i < CLASS_NUM; ++i){
            if(strstr(test_list[j], test_labels[i])){
                test_class = i;
                break;
            }
        }
        int frame_num = 0;
        WaveData data = wavRead(test_list[j], strlen(test_list[j]));
        double **feature = extract_mfcc(data.data, data.sampleRate, data.size, &frame_num);
        free(data.data);
        int *audio_vq = lbg_encode(codebook, feature, frame_num, MFCC_COEFF_ALL, VQ_NUM);
        free_matrix_double(feature, frame_num);

        double probability;
        int *path = (int *)calloc(frame_num, sizeof(int));
        for(int i = 0; i < CLASS_NUM; ++i){
            Viterbi(&hmm[i], frame_num, audio_vq, path, &probability);
            printf("%d %s class %d %f\n", i, test_list[j], test_class, probability);
            break;
        }
        free(path);
        free(audio_vq);
        break;
    }
    for(int i = 0; i < CLASS_NUM; ++i) FreeHMM(&hmm[i]);
    free_matrix_double(codebook, VQ_NUM);
    free(test_list);
    free_ptrs((void**)labels, CLASS_NUM);
    free_ptrs((void**)test_labels, CLASS_NUM);
    return 0;
}
