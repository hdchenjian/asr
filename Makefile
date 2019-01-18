DEBUG=1

VPATH=./src/
EXEC=asr
OBJDIR=./obj/

CC=gcc
LDFLAGS= -lm -pthread -Lthirdpart/fftw-3.3.8/.libs/ -lfftw3 -Lthirdpart/libsndfile-1.0.28/src/.libs/ -lsndfile
COMMON= -Isrc/ -Ithirdpart/libsndfile-1.0.28/src/ -Ithirdpart/fftw-3.3.8/api/
CFLAGS=-Wall -Wno-unknown-pragmas -Wfatal-errors -fPIC --std=gnu11 -Wunused-but-set-variable -Wno-unused-result
OPTS=-Ofast
ifeq ($(DEBUG), 1) 
OPTS=-O0 -g
endif
CFLAGS+=$(OPTS)

OBJ=mfcc.o backward.o baum.o forward.o hmmrand.o hmmutils.o lbg_vq.o lista.o matrix.o nrutil.o sequence.o viterbi.o
EXECOBJA=asr.o
EXECOBJ = $(addprefix $(OBJDIR), $(EXECOBJA))
OBJS = $(addprefix $(OBJDIR), $(OBJ))

all: obj $(EXEC)

$(EXEC): $(EXECOBJ) $(OBJS)
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

obj:
	mkdir -p obj

.PHONY: clean

clean:
	rm -rf $(OBJS) $(EXEC) $(EXECOBJ)
