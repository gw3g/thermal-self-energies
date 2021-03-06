CC	= gcc
VPATH	=thy aux
CFLAGS= -I. -lgsl -lgslcblas -lm
DEPS	= core.h
ODIR	= obj
_OBJ	= main.o htl.o thermal.o gluon.o photon.o quark.o funcs.o disp.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tl: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

