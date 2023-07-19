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
SIESTA_ARCH=cscs-cray-mpi
#
# For Cray XT-3 at CSCS (Palu)
#
FC=ftn -target=catamount
FPP=linux-pgf90 -F
FC_ASIS=$(FC)
#
FFLAGS=-fastsse
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=            #  -DCDF
#
KINDS=4 8
MPI_LIBS=
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.
FPPFLAGS_MPI=-DMPI
#
LIBS= 
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
#
# Important 
# Compile atom.f and electrostatic.f without optimization.
# Make sure that the dependency is explicit, so that these
# lines work with VPATH
#
atom.o: atom.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
electrostatic.o: electrostatic.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
.F.o:
	$(FPP) $(FPPFLAGS) $<  ; mv $*.f aux_$*.f
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) -o $*.o aux_$*.f
	rm -f aux_$*.f
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FPP) $(FPPFLAGS) $<  ; mv $*.f aux_$*.f90
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) -o $*.o aux_$*.f90
	rm -f aux_$*.f90
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
