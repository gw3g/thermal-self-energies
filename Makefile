CC	= gcc
VPATH	=thy
CFLAGS= -I. -lgsl -lgslcblas -lm
DEPS	= core.h
ODIR	= obj
_OBJ	= main.o htl.o thermal.o QCD_gluon.o QED_photon.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

HTL: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

