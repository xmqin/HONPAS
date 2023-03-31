#
#  Makefile for FDF library and FDF example.
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: check_flags
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Add items to the INCFLAGS as set by the arch.make file
#
INCFLAGS:=-I$(VPATH) $(INCFLAGS)    # For VPATH operation
#
# This is needed to pick up top-level module files, such as mpi_siesta
#
LOCAL_INCFLAGS=-I../
INCFLAGS:=$(INCFLAGS) $(LOCAL_INCFLAGS)
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libfdf.a
default: module
module: libfdf.a
	cp libfdf.a ..
	@cp -p *.*d ..     
#
libfdf.a: check_flags iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o 
	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a \
                iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o 
	-$(RANLIB) libfdf.a

check_flags:
	@echo "In fdf, INCFLAGS is: $(INCFLAGS)"
#
fdf.o: io_fdf.o parse.o prec.o utils.o
io_fdf.o: utils.o prec.o
parse.o: utils.o prec.o
utils.o: prec.o

#
# Tests for FDF library (Sequential, MPI, IO)
#
sample: sample.F90 libfdf.a
	$(FC) $(LDFLAGS) -o $@ $^ $(MPIFLAGS)

#
# Cleaning
#
clean:
	rm -f *.o *.mod libfdf.a
	rm -f sample sample_mpi sample_cluster sample_blocking
	rm -f fdf.debug* sample.out* sample.fdf.*
