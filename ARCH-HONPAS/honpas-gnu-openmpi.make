# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#-------------------------------------------------------------------
# arch.make file for gfortran compiler.
# To use this arch.make file you should rename it to
#   arch.make
# or make a sym-link.
# For an explanation of the flags see DOCUMENTED-TEMPLATE.make

.SUFFIXES:
.SUFFIXES: .f .F .o .c .a .f90 .F90

SIESTA_ARCH = x86_64

CC = gcc
FPP = $(FC) -E -P -x c
FC = mpif90
FC_SERIAL = gfortran

WITH_LIBINT = 1

FFLAGS = -O2 -fPIC  -ftree-vectorize

AR = ar
RANLIB = ranlib

SYS = nag

SP_KIND = 4
DP_KIND = 8
KINDS = $(SP_KIND) $(DP_KIND)


#MKLROOT= 
#

BLAS_LIBS = -L/public/home/xmqin/MathLibs/lapack/3.11 -lblas
LAPACK_LIBS = -L/public/home/xmqin/MathLibs/lapack/3.11 -llapack
SCALAPACK_LIBS= -L/public/home/xmqin/MathLibs/lapack/3.11 -lscalapack

LIBINT_LIBS=./libint_install/lib/libderiv.a ./libint_install/lib/libint.a

LDFLAGS =

COMP_LIBS = 
#libsiestaLAPACK.a libsiestaBLAS.a

FPPFLAGS = -DMPI $(DEFS_PREFIX)-DFC_HAVE_ABORT $(DEFS_PREFIX)-DHAVE_LIBINT

LIBS = $(COMP_LIBS) $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBINT_LIBS) -lpthread -lstdc++ $(NETCDF_LIBS)

MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.
# Dependency rules ---------

FFLAGS_DEBUG = -g  -O1 # -ggdb -fbacktrace -fcheck=bounds,do,mem,pointer,recursion -fimplicit-none -finit-real=snan -fno-automatic -fno-omit-frame-pointer -O0

# The atom.f code is very vulnerable. Particularly the Intel compiler
# will make an erroneous compilation of atom.f with high optimization
# levels.
atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $< 
state_analysis.o: state_analysis.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $< 

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

