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
    new->sample_rate = data.sampleRate;
    new->audio_frame = data.size;
    new->feature = extract_mfcc(data.data, new->sample_rate, new->audio_frame, &new->frame_num);
    free(data.data);

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

