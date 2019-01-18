#include "matrix.h"

/*  The windows is defined in milliseconds */
#define DIM_FINESTRA 25
#define DIM_PASSO 5
/*  The frequency is defined in Hertz */
#define FREQUENZA_MAX 5000
/*  Number of considered cepstral coefficients */
#define N_COEFF_CEP 13

/*  Data structure for the list of words */
typedef struct list_audio
{
    char filename[128];
    int frame_num;
    int sample_rate;
    int audio_frame;
    char filename_car[128];
    char filename_seq[128];
    char filename_hmm[128];
    struct list_audio *next;
}list_audio;

list_audio *add_audio(char *filename, list_audio *lista);
void print_list_audio(list_audio *lista);
