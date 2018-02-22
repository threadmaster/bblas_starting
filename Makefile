# Top level makefile
#
include Makefile.inc

all : driver serial

driver: driver.o
	$(F90) driver.o -o driver $(LDLIBS)  

serial: 
	cd serial && $(MAKE)


