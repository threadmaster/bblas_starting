# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -O1
CFLAGS = -O1
LDLIBS = -L./serial -lserial_bblas -L./lbstime -llbstime
