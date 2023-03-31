# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for MPI tests
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90

.PHONY: MPI clean

default: pi3 blacs_prb pblas_prb # mpi_test

OBJDIR=Obj
VPATH=$(shell pwd)/../../Src
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

# Remove timings
FFLAGS := $(filter-out -DMPI_TIMING, $(FFLAGS))
FPPFLAGS := $(filter-out -DMPI_TIMING, $(FPPFLAGS))

MPI_INTERFACE= libmpi_f90.a 

pi3: $(MPI_INTERFACE) pi3.o timer_stub.o
	$(FC) -o pi3 $(LDFLAGS) pi3.o $(MPI_INTERFACE) timer_stub.o

newcomm:    MPI newcomm.o 
	$(FC) -o newcomm $(LDFLAGS) newcomm.o $(MPI_INTERFACE_LIB) $(MPI_LIBS) 
DISTOBJS=class_BlockCyclicDist.o class_PEXSIDist.o class_Dist.o \
         m_matrix.o m_transfers.o m_redist.o m_uuid.o \
         dielocal.o object_debug.o pxf.o dist.o
dist:    MPI $(DISTOBJS)
	$(FC) -o $@ $(LDFLAGS) $(DISTOBJS) $(MPI_INTERFACE_LIB) $(MPI_LIBS) 
#
# Blacs test
#
blacs_prb: blacs_prb.o
	$(FC) -o blacs_prb $(LDFLAGS) blacs_prb.o $(LIBS)

#
# Pblas test. You might need to uncomment the following line
# in those systems where pdlaprnt is not in the standard library.
#
EXTRA_ROUTINE=pdlaprnt.o
#
pblas_prb: pblas_prb.o $(EXTRA_ROUTINE)
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
mpi_test: $(MPI_INTERFACE) $(MT_OBJS) mpi_test.o
	$(FC) -o mpi_test \
	       $(LDFLAGS) mpi_test.o $(MT_OBJS) $(MPI_INTERFACE) $(LIBS)


# Note that -DMPI_TIMING might not be removed here in some systems,
# hence the timer_stub routine.
#
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
libmpi_f90.a: 
	-mkdir -p MPI
	(cd MPI ; \
	$(MAKE) -f $(MPI_MAKEFILE) -j 1 "VPATH=$(VPATH)/MPI" \
	"MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
	"FPPFLAGS=$(FPPFLAGS)" \
	"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)


#----------------------------------------------------------------------------
#
clean:
	rm -f *.o *.*d  
	rm -f pi3 blacs_prb pblas_prb newcomm
	rm -rf MPI $(MPI_INTERFACE)

# Manual dependencies
m_mpi_utils.o: sys.o
sys.o: parallel.o
