#include "matrix.h"
#include <math.h>
#include <complex.h>

#define MFCC_COEFF_NUM 13
#define FRAME_LENGTH 25
#define FRAME_SHIFT 10

typedef struct {
    int16_t *data;
    int size;
    int sampleRate;
} WaveData;

double **extract_mfcc(int16_t *audio_data, int sample_rate, int audio_data_num, int *frame_num);
WaveData wavRead(char fileName[], size_t fileNameSize);
