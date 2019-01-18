#include <string.h>

#include "lista.h"
#include "mfcc.h"

list_audio *add_audio(char *filename, list_audio *lista)
{
    list_audio *new = (list_audio *)calloc(1, sizeof(list_audio));
    if(new == 0){
        printf("malloc failed.\n");
        exit(1);
    }

    /*  Copy the file name */
    strcpy(new->filename , filename);
    strcpy(new->filename_car , "");
    strcpy(new->filename_seq , "");
    strcpy(new->filename_hmm , "");
    new->next = NULL;

    double *signal = open_wav(new->filename, &new->sample_rate, &new->audio_frame);
    double **mfcc=extract_mfcc(signal, new->sample_rate, new->audio_frame, FREQUENZA_MAX, DIM_FINESTRA,
                               DIM_PASSO, N_COEFF_CEP, &new->frame_num);

    char *filename_bak = strdup(new->filename);
    /*  Search for the last occurrence of the character '.' */
    char *pos = strrchr(filename_bak, '.');

    int i = pos - filename_bak;
    /*  After the '.', it changes the extension */
    filename_bak[i+1] = 'c';
    filename_bak[i+2] = 'a';
    filename_bak[i+3] = 'r';
    strcpy(new->filename_car, filename_bak);
    save_matrix_double(new->filename_car, mfcc, new->frame_num, 3*N_COEFF_CEP);

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

/*  Function that prints a list on the screen */
void print_list_audio(list_audio *lista)
{
    list_audio *corrente;
    int i=1;

    corrente = lista;
    /*  Start from the top and scan the list */
    while(corrente != NULL)
    {
        printf("Parola numero:%3i\nNome File: %s\n", i++, corrente->filename);
        printf("N.finestre: %d\nFreq. camp.: %d\n", corrente->frame_num, corrente->sample_rate);
        printf("N.campioni: %d\n", corrente->audio_frame);
        printf("Nome Carat: %s\n", corrente->filename_car);
        printf("Nome Seq: %s\n", corrente->filename_seq);
        printf("Nome HMM: %s\n", corrente->filename_hmm);
        printf("\n\n");
        corrente = corrente->next;
    }
}
