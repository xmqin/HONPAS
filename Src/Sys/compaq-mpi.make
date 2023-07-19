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
# Makefile include file for parallel Compaq machine at Grenoble

SIESTA_ARCH=compaq-mpi

FC=f90 

FFLAGS=  -O2
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -Rabc -ei

NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=

MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
FPPFLAGS_MPI=-DMPI

LIBS = -L/usr/local/lib/scalapack -lscalapack -lpblas -ltools -lblacsF77 \
        -lblacs -lblacsF77 \
        -lblacs -ldxml -lfmpi -lmpi -lelan

SYS=bsd
RANLIB=echo
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) 

# Actual compilation recipes for siesta code.

.F.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $(FPPFLAGS) $<
.f90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $<

