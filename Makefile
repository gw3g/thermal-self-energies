CC	= gcc
#VPATH	=intgrand:x-sec
CFLAGS= -I. -lgsl -lgslcblas -lm
DEPS	= core.h
ODIR	= obj
_OBJ	= main.o htl.o integrand.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

HTL: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

