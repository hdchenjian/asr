#include "matrix.h"

/*  The windows is defined in milliseconds */
#define FRAME_LENGTH 25
#define FRAME_SHIFT 5
/*  The frequency is defined in Hertz */
#define FREQUENZA_MAX 5000
/*  Number of considered cepstral coefficients */
#define MFCC_COEFF_NUM 13

/*  Data structure for the list of words */
typedef struct list_audio
{
    char filename[128];
    int frame_num;
    int sample_rate;
    int audio_frame;
    char filename_hmm[128];
    double **mfcc;
    int *vq_result;
    struct list_audio *next;
}list_audio;

list_audio *add_audio(char *filename, list_audio *lista);
