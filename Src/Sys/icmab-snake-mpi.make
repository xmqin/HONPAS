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
SIESTA_ARCH=IFC91_EM64T-ICT3.0-bin64-MPD
### Need . /apl/hyperbird/FC91_EM64T/bin/ifortvars.sh
### Need . /apl/INTEL/ICT3/mpi/3.0/bin64/mpivars.sh

#
# Intel fortran compiler 10.0 for linux with mkl 9.0 optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
# Note: The -mp1 option is necessary to recover IEEE floating point precision.
#
FC=mpiifort
FC_ASIS=$(FC)
#
FFLAGS= -w -O2 -mp1
FFLAGS_DEBUG= -g
LDFLAGS=
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/apl/INTEL/ICT3/mpi/3.0/include64
FPPFLAGS_MPI=-DMPI
#
MKLPATH=/apl/INTEL/ICT3/cmkl/9.0
LIBS=-L$(MKLPATH)/lib/em64t -lmkl_scalapack -lmkl_blacs_intelmpi20 -lmkl_lapack -lmkl_em64t -lguide -lpthread -lsvml
#LIBS=-L$(MKLPATH)/lib/32 -lmkl_scalapack -lmkl_blacs_intelmpi20 -lmkl_lapack -lmkl -lguide -lpthread -lsvml
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DWXML_INIT_FIX
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


