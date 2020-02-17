#include <math.h>
#include <stdint.h>

#define MFCC_COEFF_NUM 13
#define MFCC_COEFF_ALL 40
#define FRAME_LENGTH 23
#define FRAME_SHIFT 10
#define FFT_NUM 1024


typedef struct {
    int16_t *data;
    int size;
    int sampleRate;
} WaveData;

typedef struct list_audio
{
	int class;
    char filename[128];
    int frame_num;
    int sample_rate;
    double **feature;
    int *vq_result;
    struct list_audio *next;
} list_audio;

list_audio *add_audio(char *filename, list_audio *lista);

double **extract_mfcc(int16_t *audio_data, int sample_rate, int audio_data_num, int *frame_num);
WaveData wavRead(char fileName[], size_t fileNameSize);
