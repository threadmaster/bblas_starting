# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -O2
CFLAGS = -O2
LDLIBS = -L./serial -L./lbstime -lserial_blas -llbs_time
