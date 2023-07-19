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
SIESTA_ARCH=intel9-mkl8
#
# Intel fortran compiler 9 for linux with mkl 8 optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
# Note: The -mp1 option is necessary to recover IEEE floating point precision.
#
FC=ifort
#
FFLAGS= -O2 -mp
EXTRA_LIBS=-lpthread -lsvml
FFLAGS_DEBUG= -g 
LDFLAGS= -static
COMP_LIBS=
RANLIB=echo
#
NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2-ifort
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
FPPFLAGS_CDF= -DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
GUIDE=/opt/intel/mkl/8.0.1/lib/32/libguide.a
LAPACK=/opt/intel/mkl/8.0.1/lib/32/libmkl_lapack.a
BLAS=/opt/intel/mkl/8.0.1/lib/32/libmkl_ia32.a
LIBS=$(LAPACK) $(BLAS)  $(GUIDE) $(EXTRA_LIBS) $(NETCDF_LIBS)
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








