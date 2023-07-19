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
.SUFFIXES: .f .F .o .a .f90 .F90
#
# For arina (Itanium cluster) at the SGI-UPV - new libraries
#
SIESTA_ARCH=ia64-linux-ifort-mlib

FPP=
FPP_OUTPUT= 
FC=mpif90    
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=-O0 -ftrapuv -CB
FPPFLAGS= -DWXML_INIT_FIX -DHAS_DLAMC3 -DALLOC_SAVE_BUG -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT
LDFLAGS= -openmp

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

SCALAPACK_LIBS=-L/opt/mlib/intel_9.0/hpmpi_2.2/lib/64 -lscalapack 
LAPACK_LIBS=-L/opt/mlib/intel_9.0/hpmpi_2.2/lib/64 -llapack -lveclib -lsolvers

COMP_LIBS=scalapack_extra.o dc_lapack.a 

NETCDF_LIBS=
NETCDF_INTERFACE=

LIBS=$(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(NETCDF_LIBS)

#SIESTA needs an F90 interface to MPI
#This will give you SIESTA's own implementation
#If your compiler vendor offers an alternative, you may change
#to it here.
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.

#Dependency rules are created by autoconf according to whether
#discrete preprocessing is necessary or not.
.F.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<
