#include  <math.h>
#include "lbg_vq.h"

/* training M * k, output: N * K, N is the desired number of codevectors, which together form
    the codebook C = {c_1, c_2, ..., c_N} where every codevector has the dimension K
    c_n = (c_n,1, c_n,2, ..., c_n,K)    with n = 1, 2, ..., N
    The relationship M > =  1000*N should normally be valid */
double **lbg(double **data, int M, int K, int N)
{
    printf("lbg, M %d, K %d, N %d\n", M, K, N);
    double epsilon =  0.001;
    int k,m,n;
    double **codebook = matrix_double(N, K);
    double **codebook_star = matrix_double(N, K);
    int quanti[N], indice;
    int flag_epsilon = 0;
    double distance_min, distance_temp;
    double distort, distort_star, rate_distort;
    int codevector_num = 1;

    for (k = 0; k < K; k++) {
        for (m = 0; m < M; m++){
            codebook_star[0][k] += data[m][k];
        }
        codebook_star[0][k] /=  (float)M;
    }

    /*  Until is not reached the desired number of codevectors */
    while (codevector_num < N-1)
    {
        /* DIVISION */
        for(n = 1; n <= codevector_num; n++)
            for(k = 1; k < K; k++)
            {
                codebook[n][k] = (1+epsilon)*codebook_star[n][k];
                codebook[n+codevector_num][k] = (1-epsilon)*codebook_star[n][k];
            }

        codevector_num *= 2;

        /*  Set a high value to distort_star so the cycle is repeated more than once */
        distort_star = 10.0E15;

        flag_epsilon = 0;
        /*  Until rate_distort is greater than epsilon */
        while (flag_epsilon == 0)
        {
            /*  ITERATION */
            distort = 0.0;

            /*  Preparation of vectors */
            for(n = 1; n <= codevector_num; n++)
            {
                for(k = 1; k < K; k++)
                    codebook_star[n][k] = 0.0;
                quanti[n] = 0;
            }

            for(m = 1; m < M; m++)
            {
                distance_min = 0.0;

                /*  Nearest Neighbor Condition
                    Associate the closest codevector */

                /*  Because there is a comparison later, it imposes the minimum distance values for n = 1 */
                for(k = 1; k < K; k++)
                {
                    /*  Take in consideration codebook[n][k] with n = 1 */
                    distance_min +=  pow( (data[m][k] - codebook[1][k]), 2 );
                }
                /*  Store the index n = 1 */
                indice = 1;

                /*  Start from n = 2, because n = 1 it has already processed before */
                for(n = 2; n <= codevector_num; n++)
                {
                    distance_temp = 0.0;
                    for(k = 1; k < K; k++)
                        distance_temp +=  pow( (data[m][k] - codebook[n][k]), 2 );

                    if(distance_temp < distance_min)
                    {
                        distance_min = distance_temp;
                        indice = n;
                    }
                }
                /*  Increase of one the n-th position of the codevector that satisfied the previous criteria */
                ++quanti[indice];

                /*  Start to calculate the codebook that will be updated */
                for(k = 1; k < K; k++)
                    codebook_star[indice][k] +=  data[m][k];

                distort = distort+distance_min;
            }

            /*  Centroid Condition
                The codevector c_n should be the average of all vectors of the training that are in the codified region S_n */
            for(n = 1; n <= codevector_num; n++)
            {
                if(quanti[n] > 0)
                {
                    for(k = 1; k < K; k++)
                    {
                        codebook_star[n][k] /=  (float) quanti[n];
                        codebook[n][k] = codebook_star[n][k];
                    }
                }
            }

            distort /= (float) M*K;

            rate_distort = (distort_star-distort)/distort_star;

            flag_epsilon = 1;

            if(rate_distort > 0.001) {
                flag_epsilon = 0;
                distort_star = distort;
            }
        }
    }

    return codebook;
}

/*  Function that quantises the data passed on the function of codebook */
int *vq(double **codebook, double **dati, int n_fin, int K, int N)
{
    double distance_min, distance_temp;
    int m, k, n, indice;
    int *quali = array_int(n_fin+1);

    for(m = 0; m < n_fin; m++)
    {
        /*  Nearest Neighbor Condition
            Associate the closest codevector */
        distance_min = 0.0;

        for(k = 1; k < K; k++)
        {
            distance_min +=  pow( (dati[m][k-1]-codebook[1][k]), 2);
        }
        indice = 1;

        for(n = 2; n < N; n++)
        {
            distance_temp = 0.0;
            for(k = 1; k < K; k++)
                distance_temp +=  pow( (dati[m][k-1]-codebook[n][k]), 2);

            if(distance_temp < distance_min)
            {
                distance_min = distance_temp;
                indice = n;
            }
        }
        quali[m+1] = indice;
    }
    return quali;
}

