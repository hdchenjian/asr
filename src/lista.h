#include "matrix.h"

typedef struct list_audio
{
    char filename[128];
    int frame_num;
    int sample_rate;
    int audio_frame;
    char filename_hmm[128];
    double **feature;
    int *vq_result;
    struct list_audio *next;
}list_audio;

list_audio *add_audio(char *filename, list_audio *lista);
