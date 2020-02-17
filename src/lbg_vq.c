#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"
#include "lbg_vq.h"

/* training M * k, output: N * K, N is the desired number of codevectors that has the dimension K  */
double **lbg_get_codebook(double **data, int M, int K, int N)
{
    double epsilon =  0.001;
    double **codebook = matrix_double(N, K);
    double **codebook_previous = matrix_double(N, K);
    int *category_num = calloc(N, sizeof(int));;

    for(int k = 0; k < K; k++) {
        for(int m = 0; m < M; m++){
            codebook_previous[0][k] += data[m][k];
        }
        codebook_previous[0][k] /=  (float)M;
    }

    int codevector_num = 1;
    while(codevector_num < N) {
        for(int n = 0; n < codevector_num; n++){         /* split */
            for(int k = 0; k < K; k++) {
                codebook[n][k] = (1 + epsilon) * codebook_previous[n][k];
                codebook[n + codevector_num][k] = (1 - epsilon) * codebook_previous[n][k];
            }
        }
        codevector_num *= 2;
        double error_previous = DBL_MAX;

        int flag_epsilon = 0;
        while(flag_epsilon == 0) {
            double error_current = 0.0;
            for(int n = 0; n < codevector_num; n++) {
                for(int k = 0; k < K; k++){
                    codebook_previous[n][k] = 0;
                }
                category_num[n] = 0;
            }

            for(int m = 0; m < M; m++) {
                int index = -1;
                double distance_min = DBL_MAX;  // todo: replease DBL_MAX to FLT_MAX
                /*  Nearest Neighbor Condition Associate the closest codevector */
                for(int n = 0; n < codevector_num; n++) {
                    double distance_temp = 0;
                    for(int k = 0; k < K; k++){
                        double tmp = data[m][k] - codebook[n][k];
                        distance_temp += (tmp * tmp);
                    }
                    if(distance_temp < distance_min) {
                        distance_min = distance_temp;
                        index = n;
                    }
                }
                category_num[index] += 1;
                for(int k = 0; k < K; k++) codebook_previous[index][k] +=  data[m][k];
                error_current = error_current + distance_min;
            }

            for(int n = 0; n < codevector_num; n++) {
                if(category_num[n] > 0) {
                    for(int k = 0; k < K; k++) {
                        codebook_previous[n][k] /= (float)category_num[n];
                        codebook[n][k] = codebook_previous[n][k];
                    }
                }
            }

            error_current /= (float)(M * K);
            double rate_distort = (error_previous - error_current) / error_previous;
            if(rate_distort > 0.001) {
                flag_epsilon = 0;
                error_previous = error_current;
            } else {
                flag_epsilon = 1;
                printf("iterate end: error diff %f, codevector_num %d\n", rate_distort, codevector_num);
            }
        }
    }
    free(category_num);
    free_matrix_double(codebook_previous, N);
    return codebook;
}

int *lbg_encode(double **codebook, double **data, int M, int K, int N)
{
    int *result_vq = calloc(M, sizeof(int));
    for(int m = 0; m < M; m++) {
        double distance_min = DBL_MAX;
        int index = -1;

        for(int n = 0; n < N; n++) {
            double distance_temp = 0.0;
            for(int k = 0; k < K; k++){
                double tmp = data[m][k] - codebook[n][k];
                distance_temp += (tmp * tmp);
            }
            if(distance_temp < distance_min) {
                distance_min = distance_temp;
                index = n;
            }
        }
        result_vq[m] = index + 1;
    }
    return result_vq;
}
