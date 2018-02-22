# Top level makefile
#
include Makefile.inc

all : driver serial lbstime

driver: driver.o serial lbstime
	$(F90) driver.o -o driver $(LDLIBS)  

driver.o: driver.f90
	$(F90) $(FFLAGS) driver.f90 -c  

serial: 
	cd serial && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

clean:
	cd serial && $(MAKE) clean
	cd lbstime && $(MAKE) clean

pristine:
	cd serial && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine

#This next target get "made" every time
.PHONY: serial lbstime
