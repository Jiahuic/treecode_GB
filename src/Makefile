
# linux
# LIBS = -llapack -lblas -lm
LIBS = -lm

#CFLAGS = -O2
CFLAGS  = -g
CC = gcc -c
LD = gcc


HEADERS = global.h tree.h
# OBJ = gkDisplay.o gkSetup.o input.o gkGlobal.o expan.o moments.o numQuad.o fmm.o kernel.o solver.o direct.o
OBJ = input.o global.o tree_display.o tree.o

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) -o $@ $<

main: main.o $(OBJ)
	$(LD) -o main main.o $(OBJ) $(LIBS) $(LIBSLA)

clean:
	\rm -f *.o
