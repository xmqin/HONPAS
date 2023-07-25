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

DUMMY_FOX= --enable-dummy
SIESTA_ARCH=x86_64

FPP=
FPP_OUTPUT= 
FC=mpifort
RANLIB=ranlib


SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=  -O2 -fPIC  -ftree-vectorize
# -traceback  -fp-model precise
FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 
LDFLAGS= 

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

BLAS_LIBS = -L/public/home/xmqin/MathLibs/lapack/3.11 -lblas
LAPACK_LIBS = -L/public/home/xmqin/MathLibs/lapack/3.11 -llapack
SCALAPACK_LIBS= -L/public/home/xmqin/MathLibs/lapack/3.11 -lscalapack

COMP_LIBS=dc_lapack.a

#Please install C++ GTO-ERIs library LIBINT (version = 1.1.x) and link it appropriately !
# In static link,  libderiv.a must be given first because libderiv.a depends on libint.a !
# You will get " undefined reference errors" if the order is reversed in static link !
#
# You can use linker flags "--start-group --end-group" to solve circular dependency problems.
#  -W1 --start-group $(LIBINT_LIBS) -Wl,--end-group -lstdc++
#
# -lsrdc++ for c++ library !
#
LIBINT_LIBS=./libint_install/lib/libderiv.a ./libint_install/lib/libint.a


NETCDF_LIBS=
NETCDF_INTERFACE=

LIBS= $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBINT_LIBS) -lpthread -lstdc++ $(NETCDF_LIBS)
#

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

