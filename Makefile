DEBUG=1

VPATH=./src/
EXEC=asr
OBJDIR=./obj/

CC=gcc
LDFLAGS= -lm -pthread
COMMON= -Isrc/
CFLAGS=-Wall -Wno-unknown-pragmas -Wfatal-errors -fPIC --std=gnu11 -Wunused-but-set-variable -Wno-unused-result
OPTS=-Ofast
ifeq ($(DEBUG), 1) 
OPTS=-O0 -g
endif
CFLAGS+=$(OPTS)

OBJ=mfcc.o lbg_vq.o matrix.o hmm.o data.o utils.o list.o
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
