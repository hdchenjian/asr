#include "matrix.h"
#include <math.h>
#include "sndfile.h"
#include "fftw3.h"
#include <complex.h>

typedef struct {
    int16_t *data;
    int size;
    int sampleRate;
} WaveData;

double **fft(double **x, int frame_num, int campioni_finestra, int frequenza_campionamento, int freq_max, int *n_freq, double *freq_passo);
double freq_mel(double freq);
double mel_freq(double mel);
double **filtri_mel(int n_filtri, int n_frequenze, int freq_max, double frequenza_passo);
double **mfcc(double **absX, int frame_num, int n_frequenze, int campioni_finestra, int n_filtri, double **banco_filtri, int n_coeff_cepstrali, double *energia_log);
double *en_log(double **x, int frame_num, int campioni_finestra);
double **regressione_lineare(double **x, int frame_num, int n_coeff);
double **extract_mfcc(int16_t *segnale, int frequenza_campionamento, int campioni_segnale, int freq_max, int dim_finestra, int dim_passo, int n_coeff_cepstrali, int *n_fin);
double *open_wav(char *file_name, int *sample_rate, int *audio_frame);

WaveData wavRead(char fileName[], size_t fileNameSize);
