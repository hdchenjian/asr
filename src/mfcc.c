#include "mfcc.h"

double *open_wav(char *file_name, int *sample_rate, int *frame_num)
{
    data_audio *audio_data = (data_audio*)malloc(sizeof(data_audio));
    audio_data->file_name = file_name;
    audio_data->audio = sf_open(audio_data->file_name, SFM_READ, &audio_data->info); /*Open WAV file read-only*/
    if(audio_data->audio == NULL) {
        fprintf(stderr, "audio file not exist %s\n", audio_data->file_name);
        exit(-1);
    }

    *frame_num = audio_data->info.frames * audio_data->info.channels;
    int *audio_data_point = array_int(*frame_num);
    /*  Read the integer type samples by placing them in the array, */
    sf_readf_int(audio_data->audio, audio_data_point, *frame_num);
    /*  Maximum values (signed or unsigned) representable by an int (4 bytes) */
    double *val_n = array_double(*frame_num);
    *sample_rate = audio_data->info.samplerate;
    sf_close(audio_data->audio);

    /*  CALCULATION OF VALUES OF NORMALISED SIGNAL  */
    //int val_max_n = 0, val_min_n = 0;
    double int_max = pow(2, 31) - 1;
    for(int i = 0; i < *frame_num; i++) {
        //if(audio_data_point[i] > val_max_n) val_max_n = audio_data_point[i];
        //if(audio_data_point[i] < val_min_n) val_min_n = audio_data_point[i];
        val_n[i] = audio_data_point[i] / int_max;
    }
    //printf("val_max_n, val_min_n %d %d\n", val_max_n, val_min_n);
    return val_n;
}

double **extract_mfcc(double *segnale, int frequenza_campionamento, int campioni_segnale, int freq_max,
		int dim_finestra, int dim_passo, int n_coeff_cepstrali, int *n_fin)
{
    /*  CALCULATION OF PRE-EMPHASIS  */
    double *val_n_p=NULL;

    val_n_p=pre_enfasi(segnale, campioni_segnale, 0.97);

    /*  DIVISION SIGNAL IN WINDOWS  */
    double **dati_finestra = NULL;
    int campioni_finestra, campioni_passo, n_finestre;

    dati_finestra = preleva_dati_finestra(val_n_p, campioni_segnale, frequenza_campionamento, dim_finestra,
    		dim_passo, &campioni_finestra, &campioni_passo, &n_finestre);
    *n_fin = n_finestre;

    /* *************************** */
    /*  SHORT-TIME AVERAGE ENERGY  */
    /* *************************** */

    /*  Determine the Short-Time Average Energy, that is, the average energy relative to a window of a certain number of seconds. */
    /*  Calculates the logarithmic energy of the audio signal */
    double *energia_log=NULL;
    energia_log = en_log(dati_finestra, n_finestre, campioni_finestra);

    /* ************************ */
    /*  FAST FOURIER TRANSFORM  */
    /* ************************ */

    double **dati_fft = NULL;
    int n_frequenze;
    double frequenza_passo;
    dati_fft = fft(dati_finestra, n_finestre, campioni_finestra, frequenza_campionamento, freq_max, &n_frequenze, &frequenza_passo);

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

    cepstrali=mfcc(dati_fft, n_finestre, n_frequenze, campioni_finestra, n_filtri, banco_filtri_mel ,n_coeff_cepstrali, energia_log);

    /* ******************** */
    /*  DELTA COEFFICIENTS  */
    /* ******************** */

    double **delta=NULL;
    /*  The Delta coefficients represent the features changes over time.
       As coefficient 0 has been inserted the logarithmic energy. */
    delta=regressione_lineare(cepstrali, n_finestre, n_coeff_cepstrali);

    /* *************************** */
    /*  DOUBLE-DELTA COEFFICIENTS  */
    /* *************************** */

    double **doppio_delta=NULL;
    /*  The Double-Delta coefficients represent the features of the variation speed in time. */
    doppio_delta=regressione_lineare(delta, n_finestre, n_coeff_cepstrali);

    /* ******************** */
    /*  EXTRACTED FEATURES  */
    /* ******************** */

    double **caratteristiche=NULL;

    caratteristiche = estrai_caratteristiche(n_finestre, n_coeff_cepstrali, cepstrali, delta, doppio_delta);

    return caratteristiche;
}

/*  Function that prints on the screen the audio file information */
void print_audio_info(data_audio *dato)
{
    printf("print_audio_info: file name %s\n", dato->file_name);
    printf("frames: %ld samplerate: %d, frames / samplerate: %f, channels: %d, formato : %d, sections: %d, Seekable: %d\n",
           dato->info.frames, dato->info.samplerate, (float)dato->info.frames / dato->info.samplerate,
		   dato->info.channels, dato->info.format, dato->info.sections, dato->info.seekable);
}

/*  Function that calculates the pre-emphasis  */
double *pre_enfasi(double *x, int campioni_segnale, double alpha)
{
    int i;
    double *y = array_double(campioni_segnale);

    /*  The first value is equal to that of the signal */
    y[0]=x[0];

    /*  From the second, it is possible to apply the pre-emphasis */
    for(i=1; i<campioni_segnale; i++)
        y[i]=x[i]-alpha*x[i-1];

    return y;
}

/*  Function that selects the window containing the signal data */
double **preleva_dati_finestra(double x[], int campioni_segnale, int freq_camp, int dim_finestra, int dim_passo, int *camp_fin, int *camp_passo, int *n_fin)
{
    int i,n,m;

    /*  Determine how many samples are in a window */
    int campioni_finestra = floor( (float)freq_camp * dim_finestra/1000);
    *camp_fin = campioni_finestra;

    /*  Determine how many samples are in a step */
        /*  ceil approximates by excess
            floor approximates by defect */
    int campioni_passo = floor( (float)freq_camp * dim_passo/1000 );
    *camp_passo = campioni_passo;

    /*  Determine the number of windows that will be in the signal */
    int n_finestre = floor( (float)(campioni_segnale - campioni_finestra)/campioni_passo );
    *n_fin = n_finestre;

    /*  Determine the displacement vector that contains the sample number from which start the window */
    int *spiazzamento=array_int(n_finestre);

    for(i=0; i<n_finestre; i++)
        spiazzamento[i] = i*campioni_passo;

    /*  Prepare the Hann window */
    double *hann=array_double(campioni_finestra);
    for(i=0; i<campioni_finestra; i++)
        hann[i] = 0.5*(1-cos( (2*M_PI*i) / (campioni_finestra-1) ));

    /*  Prepare the matrix that contains the data for each window */
    double **y = matrix_double(n_finestre, campioni_finestra);

    /*  Multiply the data for the used window */
    for(n=0; n<n_finestre ; n++ )
        for(m=0; m<campioni_finestra; m++)
            y[n][m]=x[ spiazzamento[n] + m ] * hann[m];

    /*  Return the double pointer memory address */
    return y;
}

/*  Function that computes the Fast Fourier Transform */
double **fft(double **x, int n_finestre, int campioni_finestra, int frequenza_campionamento, int freq_max, int *n_freq, double *freq_passo)
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

    in = fftw_malloc ( sizeof ( double ) * n_fft * n_finestre);
    out = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * nc * n_finestre);

    /*  Prepare the vector to be processed
       Because the algorithm accepts a vector instead of a matrix,
       it is necessary to put the rows side by side. */
    for(i=0; i<(n_finestre*n_fft); )
        for(n=0; n<n_finestre; n++)
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
    p = fftw_plan_dft_r2c_2d( n_finestre, n_fft, in, out, FFTW_ESTIMATE );
    /*  After creating the plan this function is performed to compute the FFT.
       It will get the vector output with the values of the transformed. */
    fftw_execute ( p );

    /*  Fetch the real and imaginary values of the Fourier transform, and then
       inserts in a matrix that facilitates the understanding of the algorithm. */
    double complex *Y[n_finestre];
    double absY_max=0;

    double **absY = (double **)malloc( sizeof(double *) * n_finestre );

    for(n=0; n<n_finestre; n++)
    {
        Y[n] = malloc ( sizeof (double complex) * n_frequenze );
        absY[n] = (double *)malloc ( sizeof (double) * n_frequenze );
    }

    /*  The values of the FFT have as unit of measurement Pa/Hz */
    for(n=0; n<n_finestre; n++)
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
double **mfcc(double **absX, int n_finestre, int n_frequenze, int campioni_finestra, int n_filtri, double **banco_filtri, int n_coeff_cepstrali, double *energia_log)
{
    int n,m,i;
    double **y_mel = matrix_double(n_finestre, n_filtri);
    double **y_log = matrix_double(n_finestre, n_filtri);
    double **y = matrix_double(n_finestre, n_coeff_cepstrali);

    /*  Apply the triangular filter */
    y_mel = gemm(absX, n_finestre, n_frequenze, banco_filtri, n_filtri);

    /*  Calculate the spectral energy density */
    for(n=0; n<n_finestre; n++)
        for(m=0; m<n_filtri; m++)
            y_log[n][m] = log10(pow(y_mel[n][m],2));

    /*  Calculate the discrete cosine transform */
    double *in;
    double *out;
    fftw_plan p;

    int n_dct = pow(2, ceil(log2(n_filtri)));

    in = fftw_malloc ( sizeof ( double ) * n_finestre * n_dct);
    out = fftw_malloc ( sizeof ( double ) * n_finestre * n_dct);

    for(i=0; i<(n_finestre*n_dct); )
        for(n=0; n<n_finestre; n++)
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

    p = fftw_plan_r2r_2d(n_finestre, n_dct, in, out, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);

    fftw_execute ( p );

    for(n=0; n<n_finestre; n++)
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
double *en_log(double **x, int n_finestre, int campioni_finestra)
{
    double *E_log = array_double(n_finestre);
    int n,m;

    for(n=0; n<n_finestre; n++)
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
double **regressione_lineare(double **x, int n_finestre, int n_coeff)
{
    double **y = matrix_double(n_finestre, n_coeff);
    int n,m;

    for(n=0; n<n_finestre; n++)
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

/*  Function that takes the features and insert one after the other */
double **estrai_caratteristiche(int n_finestre, int n_coeff, double **mfcc, double **delta, double **doppio_delta)
{
    int n,m;
    double **y = matrix_double(n_finestre, 3*n_coeff);

    for(n=0; n<n_finestre; n++)
        for(m=0; m<n_coeff; m++)
        {
            y[n][m] = mfcc[n][m];
            y[n][m+n_coeff] = delta[n][m];
            y[n][m+2*n_coeff] = doppio_delta[n][m];
        }

    return y;
}
