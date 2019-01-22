#include <string.h>

#include "lista.h"
#include "mfcc.h"

list_audio *add_audio(char *filename, list_audio *lista)
{
    list_audio *new = (list_audio *)calloc(1, sizeof(list_audio));
    strcpy(new->filename , filename);
    strcpy(new->filename_hmm , "");
    new->next = NULL;
    WaveData data = wavRead(new->filename, strlen(new->filename));
    //for(int i = 0; i < 120; i++) printf("%d %d\n", i, data.data[i]);
    //double *signal = open_wav(new->filename, &new->sample_rate, &new->audio_frame);
    double **mfcc = extract_mfcc(data.data, new->sample_rate, new->audio_frame, FREQUENZA_MAX, FRAME_LENGTH,
                               FRAME_SHIFT, MFCC_COEFF_NUM, &new->frame_num);
    printf("audio mfcc: %s\n", new->filename);
    for(int i = 0; i < MFCC_COEFF_NUM * 3; i++) printf("%f, ", mfcc[0][i]);
    exit(0);
    new->mfcc = mfcc;
    free(data.data);

    /*  If it is the first element, it adds on the top */
    if(lista == NULL) {
        lista = new;
    } else {
    	list_audio *previous = lista;
        list_audio *next_node = lista->next;
        while(next_node != NULL) {
        	previous = next_node;
            next_node = next_node->next;
        }
        previous->next = new;
    }
    return lista;
}

