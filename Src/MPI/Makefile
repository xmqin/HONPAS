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
.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .a 
#
default: module_built 
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90
#
INCFLAGS=-I$(MPI_INCLUDE)
#
Interfaces.f90 V_S.uses VS.uses: $(TEMPLATES)
	if [ -z "$(KINDS)" ] ; then  $(MAKE) kind_explorer ; fi
	@echo "The kind numbers for single and double precision reals follow"
	sh ${<D}/generate.sh "$(KINDS)"
#
kind_explorer: kind_explorer.o
	$(FC) -o kind_explorer $(LDFLAGS) kind_explorer.o
#
mpi.o: mpi__include.o Interfaces.o V_S.uses VS.uses
#
module_built: libmpi_f90.a
	@cp libmpi_f90.a ..
	@cp *.mod ..
	@touch module_built
#
libmpi_f90.a: mpi.o mpi__include.o Interfaces.o
	@$(AR) $(ARFLAGS_EXTRA) cru libmpi_f90.a mpi.o mpi__include.o Interfaces.o
	-$(RANLIB) libmpi_f90.a
#
clean:
	@rm -f Interfaces.f90 *.o *.mod
	@rm -f module_built *.uses kind_explorer




