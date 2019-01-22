#include <string.h>

#include "mfcc.h"

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
    data.size = header.dataSize;
    data.data = (int16_t*)malloc(header.dataSize * sizeof(int16_t));
    fread(data.data, sizeof(float), header.dataSize, filePtr);
    fclose(filePtr);
    return data;
}

typedef struct
{
  SNDFILE *audio;
  SF_INFO info;
  char *file_name;
}data_audio;

double *open_wav(char *file_name, int *sample_rate, int *audio_frame)
{
    data_audio *audio_data = (data_audio*)malloc(sizeof(data_audio));
    audio_data->file_name = file_name;
    audio_data->audio = sf_open(audio_data->file_name, SFM_READ, &audio_data->info); /*Open WAV file read-only*/
    if(audio_data->audio == NULL) {
        fprintf(stderr, "audio file not exist %s\n", audio_data->file_name);
        exit(-1);
    }

    *audio_frame = audio_data->info.frames * audio_data->info.channels;
    int *audio_data_point = (int *)calloc(*audio_frame, sizeof(int));
    /*  Read the integer type samples by placing them in the array, */
    sf_readf_int(audio_data->audio, audio_data_point, *audio_frame);
    for(int i = 0; i < 120; i++) printf("%d %s %d\n", i, file_name, audio_data_point[i]);
    /*  Maximum values (signed or unsigned) representable by an int (4 bytes) */
    double *val_n = (double *)calloc(*audio_frame, sizeof(double));
    *sample_rate = audio_data->info.samplerate;
    sf_close(audio_data->audio);

    for(int i = 0; i < *audio_frame; i++) {
        //if(audio_data_point[i] > val_max_n) val_max_n = audio_data_point[i];
        //if(audio_data_point[i] < val_min_n) val_min_n = audio_data_point[i];
        val_n[i] = audio_data_point[i]; // / int_max;
    }
    free(audio_data);
    free(audio_data_point);
    return val_n;
}

double **merge_feature(int frame_num, int n_coeff, double **mfcc, double **delta, double **delta_delta)
{
    double **y = matrix_double(frame_num, 3*n_coeff);
    for(int n=0; n<frame_num; n++) {
        for(int m=0; m<n_coeff; m++) {
            y[n][m] = mfcc[n][m];
            y[n][m+n_coeff] = delta[n][m];
            y[n][m+2*n_coeff] = delta_delta[n][m];
        }
    }
    return y;
}

/*  Function that calculates the pre-emphasis, Filter:H(z) = 1-a*z^(-1) */
void pre_emphasize(int16_t *input, double *output, int length, double alpha)
{
    /*  From the second, it is possible to apply the pre-emphasis */
	output[0] = input [0];
    for(int i = 1; i < length; i++) output[i] = input[i] - alpha * input[i - 1];
}

double **make_frames_hamming(double *x, int audio_data_num, int sample_rate, int frame_length,
		int frame_shift, int *sample_per_window_ret, int *frame_num_ret)
{
    int i,n,m;
    /* how many samples are in a window */
    int sample_per_window = floor( (float)sample_rate * frame_length / 1000);
    *sample_per_window_ret = sample_per_window;
    int campioni_passo = floor( (float)sample_rate * frame_shift / 1000 );

    /*  Determine the number of windows that will be in the signal */
    int frame_num = floor( (float)(audio_data_num - sample_per_window) / campioni_passo );
    *frame_num_ret = frame_num;

    /*  Determine the displacement vector that contains the sample number from which start the window */
    int *spiazzamento = calloc(frame_num, sizeof(int));

    for(i=0; i<frame_num; i++)
        spiazzamento[i] = i*campioni_passo;

    double *hann = (double *)calloc(sample_per_window, sizeof(double));
    for(i=0; i<sample_per_window; i++)
        hann[i] = 0.5*(1-cos( (2*M_PI*i) / (sample_per_window-1) ));

    /*  Prepare the matrix that contains the data for each window */
    double **y = matrix_double(frame_num, sample_per_window);

    /*  Multiply the data for the used window */
    for(n=0; n<frame_num ; n++ )
        for(m=0; m<sample_per_window; m++)
            y[n][m] = x[spiazzamento[n] + m ] * hann[m];
    return y;
}

double **extract_mfcc(int16_t *audio_data, int sample_rate, int audio_data_num, int freq_max,
		int frame_length, int frame_shift, int coeff_num, int *frame_num)
{
	double *data = (double *)malloc(audio_data_num * sizeof(double));
    pre_emphasize(audio_data, data, audio_data_num, 0.97);
    for(int i = 0; i < 3; i++) printf("%d %d %f\n", audio_data_num, i, data[i]);

    /*  DIVISION SIGNAL IN WINDOWS  */
    double **dati_finestra = NULL;
    int sample_per_window;

    dati_finestra = make_frames_hamming(data, audio_data_num, sample_rate, frame_length,
    		frame_shift, &sample_per_window, frame_num);

    /*  SHORT-TIME AVERAGE ENERGY  */
    /*  Determine the Short-Time Average Energy, that is, the average energy relative to a window of a certain number of seconds. */
    /*  Calculates the logarithmic energy of the audio signal */
    double *energia_log=NULL;
    energia_log = en_log(dati_finestra, *frame_num, sample_per_window);

    double **dati_fft = NULL;
    int n_frequenze;
    double frequenza_passo;
    dati_fft = fft(dati_finestra, *frame_num, sample_per_window, sample_rate, freq_max, &n_frequenze, &frequenza_passo);

    /* ******************************* */
    /*  BENCH OF FILTERS ON SCALE MEL  */
    /* ******************************* */

    double **banco_filtri_mel = NULL;
    int n_filtri=24;

    banco_filtri_mel = filtri_mel(n_filtri, n_frequenze, freq_max, frequenza_passo);

    /* ********************************************* */
    /*  MEL-FREQUENCY CEPSTRAL COEFFICIENTS (MFCCs)  */
    /* ********************************************* */

    double **cepstrali = NULL;

    cepstrali=mfcc(dati_fft, *frame_num, n_frequenze, sample_per_window, n_filtri, banco_filtri_mel ,coeff_num, energia_log);

    /* ******************** */
    /*  DELTA COEFFICIENTS  */
    /* ******************** */

    double **delta=NULL;
    /*  The Delta coefficients represent the features changes over time.
       As coefficient 0 has been inserted the logarithmic energy. */
    delta=regressione_lineare(cepstrali, *frame_num, coeff_num);

    /* *************************** */
    /*  DOUBLE-DELTA COEFFICIENTS  */
    /* *************************** */

    double **delta_delta=NULL;
    /*  The Double-Delta coefficients represent the features of the variation speed in time. */
    delta_delta=regressione_lineare(delta, *frame_num, coeff_num);

    /* ******************** */
    /*  EXTRACTED FEATURES  */
    /* ******************** */

    double **caratteristiche=NULL;

    caratteristiche = merge_feature(*frame_num, coeff_num, cepstrali, delta, delta_delta);

    return caratteristiche;
}



/*  Function that computes the Fast Fourier Transform */
double **fft(double **x, int frame_num, int campioni_finestra, int frequenza_campionamento, int freq_max, int *n_freq, double *freq_passo)
{
    int i,n,m;

    /*  Declare the variables used to calculate the FFT */
    /*  The input data are real numbers */
    double *in;

    /*  The output of the FFT will have real and imaginary values
       fftw_complex as default identifies a two-dimensional array of double
       where the first column identifies the real values and the second ones
       imaginary values. */
    fftw_complex *out;

    /*  This variable will contain the plan, with all the info that the FFTW algorithm needs to compute the FFT.*/
    fftw_plan p;

    /*  Number of samples processed in the FFT */
    /*  Because the algorithm works faster with a number of samples equal to 2^n,
        set the value as follows */
    int n_fft = pow(2, ceil(log2(campioni_finestra)));
    /*  Only the first N/2 samples are significant and have the information, others are redundant. */
    int nc = (n_fft/2) + 1;

    double frequenza_passo = (double)frequenza_campionamento / n_fft;
    *freq_passo = frequenza_passo;

    int n_frequenze = ceil( freq_max/frequenza_passo );
    *n_freq = n_frequenze;

    in = fftw_malloc ( sizeof ( double ) * n_fft * frame_num);
    out = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * nc * frame_num);

    /*  Prepare the vector to be processed
       Because the algorithm accepts a vector instead of a matrix,
       it is necessary to put the rows side by side. */
    for(i=0; i<(frame_num*n_fft); )
        for(n=0; n<frame_num; n++)
            for(m=0; m<n_fft; m++)
            {
                if(m < campioni_finestra)
                {
                    in[i]=x[n][m];
                    i++;
                }
                else
                {
                    /*  Zero padding */
                    /*  Fill with zero positions of excess samples. */
                    in[i]=0;
                    i++;
                }
            }

    /*  The plan that needs to be elaborated needs of this function:
       fftw_plan_dft_r2c_2d(int n0, n1 int, double *in, fftw_complex *out, unsigned flags);
       where n0 is the number of rows and n1 that of the columns of the matrix,
       the variable flags can be FFTW_ESTIMATE or FFTW_MEASURE */
    p = fftw_plan_dft_r2c_2d( frame_num, n_fft, in, out, FFTW_ESTIMATE );
    /*  After creating the plan this function is performed to compute the FFT.
       It will get the vector output with the values of the transformed. */
    fftw_execute ( p );

    /*  Fetch the real and imaginary values of the Fourier transform, and then
       inserts in a matrix that facilitates the understanding of the algorithm. */
    double complex *Y[frame_num];
    double absY_max=0;

    double **absY = (double **)malloc( sizeof(double *) * frame_num );

    for(n=0; n<frame_num; n++)
    {
        Y[n] = malloc ( sizeof (double complex) * n_frequenze );
        absY[n] = (double *)malloc ( sizeof (double) * n_frequenze );
    }

    /*  The values of the FFT have as unit of measurement Pa/Hz */
    for(n=0; n<frame_num; n++)
    {
        /*  Fetch up to a frequency that interests,
            to pick up all the calculated frequency it is necessary to put nc instead n_frequency. */
        for(i=0; i<n_frequenze; i++)
        {
            Y[n][i] = out[i + nc*n][0] + I*out[i + nc*n][1] ;
            absY[n][i] = cabs(Y[n][i]);
            /*  Find the maximum amplitude value */
            if(absY[n][i] > absY_max)
                absY_max = absY[n][i];
        }
    }

    /*  After the FFT calculation the plan is deallocated.*/
    fftw_destroy_plan ( p);
    /*  The array are deallocated. */
    fftw_free ( in );
    fftw_free ( out );

    return absY;
}

/*  Function that converts the scaled frequency mel */
double freq_mel(double freq)
{
    /*  Calculates the value of mel scale by frequency */
    double mel;
    mel = 2595 * log10(1 + freq/700);
    return mel;
}

/*  Function that converts the value in the mel-scale frequency */
double mel_freq(double mel)
{
    /*  Calculates the frequency value given a value in the scale mel */
    double freq;
    freq = 700*(pow(10, (mel/2595)) -1);
    return freq;
}

/*  Function that creates the filter bank with mel scale */
double **filtri_mel(int n_filtri, int n_frequenze, int freq_max, double frequenza_passo)
{
    int i;

    /*  By imposing the maximum frequency, it determines the value in the mel scale */
    double mel_max = freq_mel(freq_max);

    /*  The increases in the mel scale will be determined according to the maximum value and the number of filters. */
    double mel_passo = mel_max / (n_filtri + 1);

    /*  The centre frequencies will have a value in the mel scale determined by: */
    double mel_centrali[n_filtri];

    for(i=0; i<n_filtri; i++)
        mel_centrali[i]=(i+1)*mel_passo;

    /*  The previous values are determined in frequency */
    double freq_centrali[n_filtri];

    for(i=0; i<n_filtri; i++)
        freq_centrali[i] = mel_freq(mel_centrali[i]);

    /*  Quantizes the central indices in function of the FFT */
    int indice_centrale[n_filtri];

    for(i=0; i<n_filtri; i++)
    {
        indice_centrale[i] = (int) round(freq_centrali[i] / frequenza_passo );
    }

    /*  Determine the initial indices */
    int indice_iniziale[n_filtri];

    for(i=0; i<n_filtri; i++)
    {
        if(i==0)
            indice_iniziale[i]=1;
        else
            indice_iniziale[i]=indice_centrale[i-1];
    }

    /* Determine the final indices */
    int indice_finale[n_filtri];

    for(i=0; i<n_filtri; i++)
    {
        if(i<(n_filtri-1))
            indice_finale[i]=indice_centrale[i+1];
        else
            indice_finale[i]=n_frequenze;
    }

    double **y = matrix_double(n_filtri, n_frequenze);
    int n,m;

    double incremento, decremento;

    for(n=0; n<n_filtri; n++)
    {
        /*  Left ramp */
        incremento = 1.0/(indice_centrale[n] - indice_iniziale[n]);
        for(m=indice_iniziale[n]; m<indice_centrale[n]; m++)
            y[n][m] = (m - indice_iniziale[n])*incremento;

        /*  Right ramp */
        decremento = 1.0/(indice_finale[n] - indice_centrale[n]);
        for(m=indice_centrale[n]; m<indice_finale[n]; m++)
            y[n][m] = 1 - (m - indice_centrale[n])*decremento;
    }

    y=trasposta_d(y, n_filtri, n_frequenze);

    return y;
}

/*  Function that determines the mel-frequency cepstral coefficients */
double **mfcc(double **absX, int frame_num, int n_frequenze, int campioni_finestra, int n_filtri, double **banco_filtri, int n_coeff_cepstrali, double *energia_log)
{
    int n,m,i;
    double **y_mel = matrix_double(frame_num, n_filtri);
    double **y_log = matrix_double(frame_num, n_filtri);
    double **y = matrix_double(frame_num, n_coeff_cepstrali);

    /*  Apply the triangular filter */
    y_mel = gemm(absX, frame_num, n_frequenze, banco_filtri, n_filtri);

    /*  Calculate the spectral energy density */
    for(n=0; n<frame_num; n++)
        for(m=0; m<n_filtri; m++)
            y_log[n][m] = log10(pow(y_mel[n][m],2));

    /*  Calculate the discrete cosine transform */
    double *in;
    double *out;
    fftw_plan p;

    int n_dct = pow(2, ceil(log2(n_filtri)));

    in = fftw_malloc ( sizeof ( double ) * frame_num * n_dct);
    out = fftw_malloc ( sizeof ( double ) * frame_num * n_dct);

    for(i=0; i<(frame_num*n_dct); )
        for(n=0; n<frame_num; n++)
            for(m=0; m<n_dct; m++)
            {
                if(m < n_filtri)
                {
                    in[i]=y_log[n][m];
                    i++;
                }
                else
                {
                    /*  Zero padding */
                    in[i]=0;
                    i++;
                }
            }

    p = fftw_plan_r2r_2d(frame_num, n_dct, in, out, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);

    fftw_execute ( p );

    for(n=0; n<frame_num; n++)
        for(m=0; m<n_coeff_cepstrali; m++)
        {
            if(m==0)
            {
                /*  As a first coefficient insert the logarithmic energy of the signal */
                y[n][m]=energia_log[n];
            }
            else
                y[n][m] = out[m + n_dct*n];
        }

    fftw_destroy_plan ( p);
    fftw_free ( in );
    fftw_free ( out );

    return y;
}

/*  Function that calculates the logarithmic energy */
double *en_log(double **x, int frame_num, int campioni_finestra)
{
    double *E_log = (double *)calloc(frame_num, sizeof(double));
    int n,m;

    for(n=0; n<frame_num; n++)
    {
        /*  It is possible to add directly because the vettore_d function,
            in addition to allocating the vector, initializes the values at zero */
        for(m=0; m<campioni_finestra; m++)
            E_log[n] += pow( x[n][m], 2 );

        E_log[n]=log10(E_log[n]);
    }

    return E_log;
}

/*  Function defining the linear regression */
double **regressione_lineare(double **x, int frame_num, int n_coeff)
{
    double **y = matrix_double(frame_num, n_coeff);
    int n,m;

    for(n=0; n<frame_num; n++)
        for(m=0; m<n_coeff; m++)
        {
            if(m!=0 && m!=n_coeff-1)
                y[n][m] = (x[n][m+1]-x[n][m-1])/2;
            else if(m==0)
                /*  Usa la differenza semplice di primo ordine */
                y[n][m] = x[n][m+1] - x[n][m];
            else
                y[n][m] = x[n][m]-x[n][m-1];
        }

    return y;
}
