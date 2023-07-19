# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
# Makefile for MPI tests
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: MPI clean
#
default: pi3 blacs_prb pblas_prb # mpi_test
#
OBJDIR=Obj
VPATH=../../Src
#
include ../../$(OBJDIR)/arch.make
#
INCFLAGS=-I$(MPI_INCLUDE)
INCFLAGS:=$(INCFLAGS) -I../../$(OBJDIR)/MPI # To pick up mod files
#
MPI_INTERFACE_LIB= ../../$(OBJDIR)/MPI/libmpi_f90.a 
#
pi3:    MPI pi3.o
	$(FC) -o pi3 $(LDFLAGS) pi3.o $(MPI_INTERFACE_LIB) $(MPI_LIBS) 
#
# Blacs test
#
blacs_prb:  blacs_prb.o
	$(FC) -o blacs_prb $(LDFLAGS) blacs_prb.o $(LIBS) 
#
# Pblas test. You might need to uncomment the following line
# in those systems where pdlaprnt is not in the standard library.
#
EXTRA_ROUTINE=pdlaprnt.o
#
pblas_prb:   pblas_prb.o $(EXTRA_ROUTINE)
	@echo
	@echo --Attempting to build Pblas test. 
	@echo --You might need to enable compilation of pdlaprnt.f in the Makefile
	@echo --if that routine is not in the standard library on your system...
	@echo --..or if its interface is not compatible, i.e., segfaults
	@echo
	$(FC) -o pblas_prb $(LDFLAGS) pblas_prb.o $(EXTRA_ROUTINE) $(LIBS) 
#
# Siesta-specific MPI test --- iffy on some systems
#----------------------------------------------------------------------------
MT_OBJS=precision.o m_mpi_utils.o sys.o parallel.o pxf.o
mpi_test: MPI $(MT_OBJS) mpi_test.o
	$(FC) -o mpi_test \
	       $(LDFLAGS) mpi_test.o $(MT_OBJS) $(MPI_INTERFACE_LIB) $(LIBS) 
#----------------------------------------------------------------------------
#
clean:
	rm -f *.o *.*d  
	rm -f pi3 blacs_prb pblas_prb
#
#
MPI:
	@echo "Making sure that the MPI interface is compiled..."
	(cd ../../$(OBJDIR)/MPI ; make "VPATH=$(VPATH)/MPI")
#
# DO NOT DELETE THIS LINE - used by make depend
#
m_mpi_utils.o: sys.o
sys.o: parallel.o
