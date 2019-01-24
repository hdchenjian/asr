#include <string.h>
#include "stdint.h"

#include "mfcc.h"

#define FFT_NUM    512
#define PI   3.14159265358979

typedef struct {
    char chunkID[4];            // 1-4      "RIFF"
    int32_t chunkSize;          // 5-8
    char format[4];             // 9-12     "WAVE"
    char subchunkID[4];         // 13-16    "fmt\0"
    int32_t subchunkSize;       // 17-20
    uint16_t audioFormat;       // 21-22    PCM = 1
    uint16_t numChannels;       // 23-24
    int32_t sampleRate;         // 25-28
    int32_t bytesPerSecond;     // 29-32
    uint16_t blockAlign;        // 33-34
    uint16_t bitDepth;          // 35-36    16bit support only
    char dataID[4];             // 37-40    "data"
    int32_t dataSize;           // 41-44
} WaveHeader;

WaveData wavRead(char fileName[], size_t fileNameSize)
{
    if (fileName[fileNameSize] != '\0') {
        fprintf(stderr, "wavRead: Invalid fileName: %s.\n", fileName);
        exit(-1);
    }
    FILE* filePtr = fopen(fileName, "r");
    if (filePtr == NULL) {
        fprintf(stderr, "unable to open file: %s.\n", fileName);
        exit(-1);
    }
    WaveHeader header;
    fread(&header, sizeof(header), 1, filePtr);
    if(strncmp(header.chunkID, "RIFF", 4) || strncmp(header.format,"WAVE", 4)||
       strncmp(header.subchunkID, "fmt", 3) || strncmp(header.dataID,"data", 4) ||
       header.audioFormat != 1 || header.bitDepth != 16) {
        fprintf(stderr, "Unsupported file type.\n");
        fclose(filePtr);
        exit(-1);
    }
    WaveData data;
    data.sampleRate = header.sampleRate;
    data.size = header.dataSize / sizeof(int16_t);
    data.data = (int16_t*)malloc(data.size * sizeof(int16_t));
    fread(data.data, sizeof(int16_t), data.size, filePtr);
    fclose(filePtr);
    return data;
}

/*  Function that calculates the pre-emphasis, Filter:H(z) = 1-a*z^(-1) */
void pre_emphasize(int16_t *input, double *output, int length, double alpha)
{
    /*  From the second, it is possible to apply the pre-emphasis */
    output[0] = input [0];
    for(int i = 1; i < length; i++) output[i] = (double)input[i] - alpha * (double)input[i - 1];
}

double *split_frame(double *x, int audio_data_num, int sample_rate, int audio_frame_num,
                    int frame_shift, int *sample_per_window_ret, int *frame_num_ret)
{
    int frame_length = floor((float)sample_rate * audio_frame_num / 1000 + 0.5);
    *sample_per_window_ret = frame_length;
    int frame_step = floor((float)sample_rate * frame_shift / 1000 + 0.5);
    int frame_num = 0;
    if(audio_data_num <= frame_length) frame_num = 1;
    else frame_num = 1 + (int)ceil((float)(audio_data_num - frame_length) / frame_step);
    *frame_num_ret = frame_num;
    printf("audio_data_num %d frame_num %d  frame_length %d \n", audio_data_num, frame_num, frame_length);

    double *output = malloc(frame_num * frame_length * sizeof(double));
    for(int i = 0; i < frame_num; i++) {
        for(int j = 0; j < frame_length; j++) {
            int out_index = i * frame_length + j;
            int in_index = (i == 0) ? j : i * frame_step +j;
            //printf("i, j %d %d in_index index %d %d\n", i, j, in_index, out_index);
            if(in_index < audio_data_num){
                output[out_index] = x[in_index];
            }
            else output[out_index] = 0;
        }
    }
    return output;
}

unsigned int get_bit_one_num(uint32_t n)
{
    unsigned int c =0 ;
    for (c = 0; n; ++c) {
        n &= (n - 1) ; // remove last 1
    }
    return c ;
}

uint32_t floor_log2(uint32_t x)
{
    x |= (x>>1);
    x |= (x>>2);
    x |= (x>>4);
    x |= (x>>8);
    x |= (x>>16);
    return (get_bit_one_num(x>>1));
}

typedef struct {
    double real;
    double image;
} COMPLEX;

int fft_complex(COMPLEX *x, uint32_t N)
{
    int i,j,l,k,ip;
    static uint32_t M = 0;
    static int le,le2;
    static double sR,sI,tR,tI,uR,uI;
    M = floor_log2(N);

    /* bit reversal sorting */
    l = N >> 1;
    j = l;
    ip = N - 2;
    for (i = 1; i <= ip; i++) {
        if (i < j) {
            tR = x[j].real;
            tI = x[j].image;
            x[j].real = x[i].real;
            x[j].image = x[i].image;
            x[i].real = tR;
            x[i].image = tI;
        }
        k = l;
        while (k <= j) {
            j = j - k;
            k = k >> 1;
        }
        j = j + k;
    }

    for (l=1; l<=M; l++) {   /* loop for ceil{log2(N)} */
        le  = (int)(1 << l);
        le2 = (int)(le >> 1);
        uR = 1;
        uI = 0;
        k = floor_log2(le2);
        sR = cos(PI / le2);
        sI = -sin(PI / le2);
        for (j=1; j<=le2; j++) {   /* loop for each sub DFT */
            for (i=j-1; i<N; i+=le) {  /* loop for each butterfly */
                ip = i + le2;
                tR = x[ip].real * uR - x[ip].image * uI;
                tI = x[ip].real * uI + x[ip].image * uR;
                x[ip].real = x[i].real - tR;
                x[ip].image = x[i].image - tI;
                x[i].real += tR;
                x[i].image += tI;
            }
            tR = uR;
            uR = tR * sR - uI * sI;
            uI = tR * sI + uI *sR;
        }
    }
    return 0;
}

int fft_real(COMPLEX *x, int N)
{
    int i,j,l,k;
    int M = 0;
    int ND4 = 0;
    double sR,sI,tR,tI,uR,uI;
    M = N / 2;
    /* N/2 points FFT */
    fft_complex(x, M);

    /* Even/Odd frequency domain decomposition */
    ND4 = N >> 2;
    for (i=1; i<ND4; i++) {
        j = M - i;
        k = i + M;
        l = j + M;
        x[k].real = (x[i].image + x[j].image) / 2;
        x[l].real = x[k].real;
        x[k].image = -(x[i].real - x[j].real) / 2;
        x[l].image = -x[k].image;
        x[i].real = (x[i].real + x[j].real) / 2;
        x[j].real = x[i].real;
        x[i].image = (x[i].image - x[j].image) / 2;
        x[j].image = -x[i].image;
    }
    x[N-ND4].real = x[ND4].image;
    x[M].real = x[0].image;
    x[N-ND4].image = 0;
    x[M].image = 0;
    x[ND4].image = 0;
    x[0].image = 0;

    /* Complete last stage FFT */
    uR = 1;
    uI = 0;
    k = floor_log2(M);
    sR = cos(PI / M);
    sI = -sin(PI / M);
    for (i=0; i<M; i++) {   /* loop for each sub DFT */
        k = i + M;
        tR = x[k].real * uR - x[k].image * uI;
        tI = x[k].real * uI + x[k].image * uR;
        x[k].real = x[i].real - tR;
        x[k].image = x[i].image - tI;
        x[i].real += tR;
        x[i].image += tI;

        tR = uR;
        uR = tR * sR - uI * sI;
        uI = tR * sI + uI *sR;
    }
    return 0;
}

double hz2mel(double hz)
{
    return 2595 * log10(1 + hz / 700.0f);
}

double mel2hz(double mel)
{
    return 700 * (pow(10, (mel / 2595.0)) - 1);
}

void get_filterbanks(int filter_num, int fft_num, int sample_rate, double *fbank)
{
    int high_freq = sample_rate / 2;
    int low_freq = 0;
    double low_mel = hz2mel(low_freq);
    double high_mel = hz2mel(high_freq);
    int mel_point_num = filter_num + 2;
    double *mel_points = malloc(mel_point_num * sizeof(double));
    int *bin = malloc(mel_point_num * sizeof(int));
    double mel_point_interval = (high_mel - low_mel) / (mel_point_num - 1);
    for(int i = 0; i < mel_point_num; i++){
        mel_points[i] = low_mel + i * mel_point_interval;
        bin[i] = floor(mel2hz(mel_points[i]) * (fft_num + 1) / sample_rate);
    }

    int fft_num_local = (fft_num / 2 + 1);
    for(int j = 0; j < filter_num; j++){
        for(int i = bin[j]; i < bin[j + 1]; i++){
            fbank[j * fft_num_local + i] = (double)(i - bin[j]) / (double)(bin[j + 1] - bin[j]);
        }
        for(int i = bin[j + 1]; i < bin[j + 2]; i++){
            fbank[j * fft_num_local + i] = (double)(bin[j + 2] - i) / (double)(bin[j + 2] - bin[j + 1]);
        }
    }
    free(bin);
    free(mel_points);
}

void gemm_nt_log(int M, int N, int K, double *A, double *B, double *C)
{
    //printf("gemm M %d N %d K %d\n", M, N, K);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            double sum = 0;
            for(int m = 0; m < K; m++){
                sum += A[i * K + m] * B[j * K + m];
            }
            C[i * N + j] = log(sum);
        }
    }
}

void gemm_nt_lift(int M, int N, int K, double *A, double *B, double **C)
{
    //printf("gemm M %d N %d K %d\n", M, N, K);
    int L = 22;
    double *lift = malloc(MFCC_COEFF_NUM * sizeof(double));
    for(int i = 0; i < MFCC_COEFF_NUM; i++){
    	lift[i] = 1 + (L / 2.0f) * sin(PI * i / L);
    }

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            double sum = 0;
            for(int m = 0; m < K; m++){
                sum += A[i * K + m] * B[j * K + m];
            }
            C[i][j] = sum * lift[j];
        }
    }
    free(lift);
}

void get_dct_coeff(int dct_num, int n, double *coeff)
{
    double scale_0 = sqrt(1.0f / n);
    double scale = sqrt(2.0f / n);
    for(int i = 0; i < dct_num; i++) {
        for(int j = 0; j < n; j++) {
            if(i == 0){
                coeff[i * n + j] = scale_0 * cos(PI * (double)i / (double)n * ((double)j + 0.5));
            } else {
                coeff[i * n + j] = scale * cos(PI * (double)i / (double)n * ((double)j + 0.5));
            }
        }
    }
}

double **extract_mfcc(int16_t *audio_data, int sample_rate, int audio_data_num, int *frame_num)
{
    double *data = (double *)malloc(audio_data_num * sizeof(double));
    pre_emphasize(audio_data, data, audio_data_num, 0.97);
    int sample_per_window;
    double *frame_data = split_frame(data, audio_data_num, sample_rate, FRAME_LENGTH,
                                     FRAME_SHIFT, &sample_per_window, frame_num);
    free(data);
    if(sample_per_window > FFT_NUM){
        printf("frame length %d, is greater than fft_complex size %d, frame will be truncated. Increase FFT_NUM to avoid\n",
               sample_per_window, FFT_NUM);
    }

    int fft_num_local = (FFT_NUM / 2 + 1);
    double *magnitude_spectrum = malloc(*frame_num * fft_num_local * sizeof(double));
    double *energy = calloc(*frame_num, sizeof(double));
    for (int m = 0; m < *frame_num; m++ ) {
        COMPLEX data_fft[FFT_NUM / 2 + 1];
        int fft_num_half = FFT_NUM / 2;
        for (int i = 0; i < fft_num_half && i * 2 < sample_per_window; i++ ){
            data_fft[i].real = frame_data[m * sample_per_window + 2 * i];
            data_fft[i].image = frame_data[m * sample_per_window + 2 * i + 1];
        }
        fft_real(data_fft, FFT_NUM);
        double scale = 1.0f / FFT_NUM;
        for(int i = 0; i < fft_num_local; i++){
            int index = i + m * fft_num_local;
            magnitude_spectrum[index] = scale * (data_fft[i].real * data_fft[i].real + data_fft[i].image * data_fft[i].image);
            energy[m] += magnitude_spectrum[index];
        }
    }
    free(frame_data);

    int filter_num = 26;
    double *fbank = calloc(filter_num * fft_num_local, sizeof(double));
    get_filterbanks(filter_num, FFT_NUM, sample_rate, fbank);
    double *fbak_feature = (double *)malloc(*frame_num * filter_num * sizeof(double));
    gemm_nt_log(*frame_num, filter_num, fft_num_local, magnitude_spectrum, fbank, fbak_feature);


    double *dct_coeff = malloc(MFCC_COEFF_NUM * filter_num * sizeof(double));
    get_dct_coeff(MFCC_COEFF_NUM, filter_num, dct_coeff);
    double **feature = matrix_double(*frame_num, MFCC_COEFF_NUM);
    gemm_nt_lift(*frame_num, MFCC_COEFF_NUM, filter_num, fbak_feature, dct_coeff, feature);
    for (int m = 0; m < MFCC_COEFF_NUM; m++ ) {
        printf("feature %d %f %f %f \n", m, fbak_feature[m], dct_coeff[m], feature[0][m]);
    }

    for (int m = 0; m < *frame_num; m++ ) {
    	feature[m][0] = log(energy[m]);
    }
    for (int m = 0; m < MFCC_COEFF_NUM; m++ ) {
    	printf("feature %d %f \n", m, feature[0][m]);
    }
    free(fbank);
    free(magnitude_spectrum);
    free(energy);
    free(fbak_feature);
    free(dct_coeff);
    return feature;
}
