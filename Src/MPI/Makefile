# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
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
INCFLAGS:= $(INCFLAGS) -I$(MPI_INCLUDE)
#
Interfaces.f90 V_S.uses VS.uses: $(TEMPLATES)
	if [ -z "$(KINDS)" ] ; then  $(MAKE) kind_explorer ; fi
	$(MAKE) int_explorer
	sh ${<D}/generate.sh "$(KINDS)"
#
kind_explorer: kind_explorer.o
	$(FC) $(FFLAGS) -o kind_explorer $(LDFLAGS) kind_explorer.o
#
int_explorer: int_explorer.o
	$(FC) $(FFLAGS) -o int_explorer $(LDFLAGS) int_explorer.o
#
mpi_interfaces.o: mpi__include.o Interfaces.o V_S.uses VS.uses
#
module_built: libmpi_f90.a
	@cp libmpi_f90.a ..
	@cp *.mod ..
	@touch module_built
#
libmpi_f90.a: mpi_interfaces.o mpi__include.o Interfaces.o mpi_siesta.o timer_mpi.o
	@$(AR) $(ARFLAGS_EXTRA) cru libmpi_f90.a mpi_interfaces.o \
               mpi_siesta.o mpi__include.o Interfaces.o timer_mpi.o
	-$(RANLIB) libmpi_f90.a
#
clean:
	@rm -f Interfaces.f90 *.o *.mod
	@rm -f module_built *.uses kind_explorer int_explorer
	@rm -f libmpi_f90.a

mpi_siesta.o: timer_mpi.o
Interfaces.o: timer_mpi.o


