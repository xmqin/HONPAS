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
SIESTA_ARCH=nano-intel-serial
#
# Intel fortran compiler 9 for linux with mkl 8 optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
# Note: The -mp option is necessary to recover IEEE floating point precision.
#
FC=ifort
#
FFLAGS= -w  -O3 -mp
EXTRA_LIBS=-lpthread -lsvml
FFLAGS_DEBUG= -g 
FFLAGS_CHECKS= -g -O0 -debug full -traceback -C
LDFLAGS= -static
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#

LIBS=-L/opt/intel/mkl/10.0.3.020/lib/em64t -lmkl_lapack -lmkl_em64t \
      -lguide $(EXTRA_LIBS)
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#








