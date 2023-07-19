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
SIESTA_ARCH=x86_64-unknown-linux-gnu--unknown

FPP=
FPP_OUTPUT= 
FC=mpiifort
RANLIB=ranlib


SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS= -O2 -g  -assume byterecl -w -fPIC -fp-model source -heap-arrays
# -traceback  -fp-model precise
FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 
LDFLAGS= -static-intel

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

MKL_ROOT=/public/software/intel/2019/mkl/lib/intel64
#MKL_LIBS=-L/public/software/intel/2011/composer_xe_2011_sp1.8.273/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
MKL_LIBS=-L${MKL_ROOT} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#LAPACK_LIBS=-L${MKL_ROOT} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#LAPACK_LIBS=-L/public/software/intel/2011/composer_xe_2011_sp1.8.273/mkl/lib/intel64 -lmkl_lapack95_lp64
BLACS_LIBS=-L${MKL_ROOT} -lmkl_blacs_intelmpi_lp64
SCALAPACK_LIBS=-L{MKL_ROOT} -lmkl_scalapack_lp64

#FFTW_LIB = -L/public/home/xmqin/honpas_dgx/Library/fftw/3.3.8/intel/lib -lfftw3
#FFTW_INC = -I/public/home/xmqin/honpas_dgx/Library/fftw/3.3.8/intel/include


COMP_LIBS=dc_lapack.a

#Please install C++ GTO-ERIs library LIBINT (version = 1.1.4) and link it appropriately !
# In static link,  libderiv.a must be given first because libderiv.a depends on libint.a !
# You will get " undefined reference errors" if the order is reversed in static link !
#
# You can use linker flags "--start-group --end-group" to solve circular dependency problems.
#  -W1 --start-group $(LIBINT_LIBS) -Wl,--end-group -lstdc++
#
# -lsrdc++ for c++ library !
#

LIBINT_LIBS = /public/home/xmqin/MathLibs/libint/1.1.5/intel/lib/libderiv.a /public/home/xmqin/MathLibs/libint/1.1.5/intel/lib/libint.a
#ISDF_LIBS=/public/home/xmqin/honpas_dgx/Library/isdf_kmeans_weight/SRC/libisdf.a

#ISDF_LIBS=/raid/home/xmqin/honpas/honpas_kmeans_weight/Library/isdf_kmeans_weight/SRC/libisdf.a

#POISON_LIBS=/public/home/xmqin/honpas_dgx/Library/ISF/libpsolver_isf.a -lsvml


NETCDF_LIBS=
NETCDF_INTERFACE=

LIBS= $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(MKL_LIBS) $(LIBINT_LIBS) -lpthread -lstdc++ $(NETCDF_LIBS)
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

