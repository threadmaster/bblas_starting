include ../Makefile.inc

# Timing Library targets 

all: lib

vvm.o : vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mmm.o : mmm.c
	$(CC) $(CFLAGS)  -c mmm.c

dls.o : dls.c
	$(CC) $(CFLAGS)  -c dls.c

lib: vvm.o mmm.o dls.o 
	ar -rc libserial_bblas.a vvm.o mmm.o dls.o 
	ranlib libserial_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
