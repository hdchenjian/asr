#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils.h"

void file_error(const char *s)
{
    fprintf(stderr, "Couldn't open file: %s\n", s);
    exit(0);
}

void free_ptrs(void **ptrs, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(ptrs[i]);
        ptrs[i] = 0;
    }
    free(ptrs);
}

char *fgetl(FILE *fp)
{
    if(feof(fp)) return 0;
    size_t size = 512;
    char *line = malloc(size*sizeof(char));
    if(!fgets(line, size, fp)){
        free(line);
        return 0;
    }

    size_t curr = strlen(line);

    while((line[curr-1] != '\n') && !feof(fp)){
        if(curr == size-1){
            size *= 2;
            line = realloc(line, size*sizeof(char));
            if(!line) {
                printf("%ld\n", size);
                fprintf(stderr, "Malloc error\n");
                exit(-1);
            }
        }
        size_t readsize = size-curr;
        if(readsize > INT_MAX) readsize = INT_MAX-1;
        char *return_point = fgets(&line[curr], readsize, fp);
        if(return_point == 0){
            fprintf(stderr, "fread read 0 object");
        }
        curr = strlen(line);
    }
    if(line[curr-1] == '\n') line[curr-1] = '\0';

    return line;
}
