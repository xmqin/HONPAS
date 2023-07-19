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
#-------------------------------------------------------------------
# DOCUMENTED arch.make
#
# The most useful makefile symbols are explained. Use this file as
# a guide when you are looking at the .make files in this directory,
# or after 'configure' has produced a draft arch.make for you.
#
# This block tells make to consider only these suffixes in its operation
# It is included in most makefiles in the source tree, but it does not
# hurt to have it here too.
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

#
# This string will be copied to the executable, so it pays to
# use something meaningful
SIESTA_ARCH=cscs-ibm-blanc-mpi

#
# In case your compiler does not understand the special meaning of 
# the .F and .F90 extensions ("files in need of preprocessing"), you
# will need to use an explicit preprocessing step.
FPP=/opt/ibmcmp/xlf/10.1/exe/cpp

#
# FC is typically the name of your fortran compiler. It is not always
# a good idea to add options here, except when they are essential for
# a proper operation. For example, here we request 64 bits and a special
# sub-compiler.
FC=mpfort -m64 -compiler xlf95_r
#
# The FC_SERIAL symbol is useful in at least two cases:
#   1. When the "MPI compiler environment" is so complex that it might
#      trick the configure scripts (for FoX at least).
#   2. When executables compiled with a (parallel) FC are flagged by 
#      the computer centers as "queuing-system-only". 
# Most utilities are thus compiled with FC_SERIAL, which in practice
# defaults to FC if it is not defined.
FC_SERIAL=xlf95_r -q64
#
# Here we should put mainly optimization flags
FFLAGS=-O3 -qstrict -qarch=pwr5 -qtune=pwr5 -qcache=auto -qessl
#
# Some systems do not have 'ranlib'. If so, use "echo" instead of "ranlib"
RANLIB=ranlib

# A compiler-specific file holding special versions of some routines
# For most f95 compilers, "nag" should work. (The name is historical)
SYS=xlf   # for IBM machines
#SYS=nag

# These symbols should not need to be specified. They will be detected
# automatically at the time of compiling the MPI interface. Set them
# only if the automatic detection fails and you are sure of their values.
SP_KIND=1
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

# Some compilers (notably IBM's) are not happy with the standard 
# syntax for definition of preprocessor symbols ( -DSOME_SYMBOL),
# and thy need a prefix (i.e. WF,-DSOME_SYMBOL). This is used
# in some utility makefiles.
DEFS_PREFIX=-WF,

#
# Used only at the linking stage. For example, you might neeed "-static"
LDFLAGS=

# Extra flags for library creation by the 'ar' command
ARFLAGS_EXTRA=
# Note that the 'ar' command can itself be specified by
# defining the AR variable. In most 'make' programs, AR is a
# built-in variable

#
# These symbols help to keep the building rules concise
# (they are generated automatically by the 'configure' script
# in some cases) 
FCFLAGS_fixed_f=-qfixed -qsuffix=cpp=f
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=-qfixed -qsuffix=cpp=F
FPPFLAGS_free_F90=

#
# This is the most installation-dependent part
# We can make things a bit easier by grouping symbols, and maybe
# using the -L flag to define search directories (see examples
# in this directory).
#
BLAS_LIBS=-lessl
LAPACK_LIBS=/apps/lapack/3.1.1/XL/lib/liblapack_ppc64.a -lessl
BLACS_LIBS=/apps/scalapack/64/lib/blacsCinit_MPI-ppc64-0.a /apps/scalapack/64/lib/blacsF77init_MPI-ppc64-0.a /apps/scalapack/64/lib/blacs_MPI-ppc64-0.a
SCALAPACK_LIBS=/apps/scalapack/64/lib/blacsCinit_MPI-ppc64-0.a /apps/scalapack/64/lib/blacs_MPI-ppc64-0.a -L/apps/scalapack/64/lib -lscalapack
#
# If you are using a "wrapper compiler" such as mpif90, MPI_LIBS can
# be left empty. If not, you might need something like -L somepath/ -l mpi ....
#
MPI_LIBS=
#
# Even if you have an optimized system library (such as ESSL), you might
# not have all of LAPACK. In particular, the divide_and_conquer routines
# might not be available, or they might be buggy. 
# In this case, you need to compile them from source 
# (COMP stands for "compiled")
# If you do not have any optimized linear algebra library, you can
# specify COMP_LIBS=linalg.a
#
##COMP_LIBS=dc_lapack.a 

#
# For netCDF support. Make sure you get a version compatible
# with the other options (for example, 32/64 bit). Don't forget
# to set -DCDF below.
#
NETCDF_ROOT=/apps/netcdf/64
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf

#
# This (as well as the -DMPI definition) is essential for MPI support
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.

# Preprocessor definitions or flags.
# Here we use FPPFLAGS (as 'configure' calls them), but historically
# it was very common to use DEFS. Try to use only FPPFLAGS from now on,
# converting any old arch.make files you might have lying around, and
# remember that you have to change the final building rules at the end
# to use only FPPFLAGS. DEFS is deprecated.

# CDF and MPI are self-explanatory
# Other definitions might be needed to work around some glitch in the compiler
#
# For old versions of gfortran, add -DGFORTRAN
#
FPPFLAGS= -WF,-DMPI -WF,-DFC_HAVE_ABORT -WF,-DCDF
#
# We put here all the neeeded libraries.
# Sometimes the BLAS are included in LAPACK (or it could be that everything
# is included in SCALAPACK...). You might need to experiment if you find 
# duplicate symbols. See examples in this directory.
#
LIBS=$(NETCDF_LIBS) $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) \
     $(MPI_LIBS) $(COMP_LIBS)

# Dependency rules ---------

# Some compilers are not able to compile certain files with full optimization,
# or they produce wrong results if they do. For example, the PGI compiler 
# has trouble with atom.f and electrostatic.f. In these cases, we need to
# insert extra lines. Use exactly the format shown, as it is general enough
# to work with VPATH.
#
FFLAGS_DEBUG=-g -O0   # your appropriate flags here...
#
atom.o: atom.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
electrostatic.o: electrostatic.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#

# Finally, the default building rules which will be used everywhere,
# unless overriden.
# These were created by a former run of 'configure'.
# See other examples in this directory. If you cut and paste, 
# MAKE SURE that there are TABS, not spaces, at the beginning.
#
# Important points to note:
#  - INCFLAGS must be present. It is used in several utility makefiles
#  - Either FPPFLAGS (preferred) or DEFS (deprecated) must be present
#    (see above) -- Note that the use of DEFS might break Util compilations.
#  - If your compiler does not recognize .F and .F90 extensions as in
#    need of preprocessing, you will need to use an intermediate
#    preprocessing step (see above about FPP). For example:
##
#.F90.o:
#        $(FPP) $(FPPFLAGS) $< > tmp_$*.f90
#        $(FC) -c $(FFLAGS) $(INCFLAGS) tmp_$*.f90
#        @mv tmp_$*.o $*.o
#        @rm -f tmp_$*.f90
#
#---------------- Example of actual rules  
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

