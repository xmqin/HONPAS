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
SIESTA_ARCH=nano-intel-mpi
#
# To run in parallel, make sure that
# /opt/intel/impi/3.1/bin64 is in your path, and do
#
#   mpirun -r ssh -np NPROCS --mpd=/opt/intel/impi/3.1/bin64/mpd siesta ....
#
#--------------------------------------------------------------------------
# Note: The -mpX option is necessary to recover IEEE floating point precision.
#
FC=/opt/intel/impi/3.1/bin64/mpiifort
#
#  You can play with other optimization options
#  I am not sure whether the compiler attempts to multithread the code
#
FFLAGS= -w  -O3 -mp
FFLAGS_CHECKS=-g -O0 -debug full -traceback -C
FFLAGS_DEBUG= -g 
LDFLAGS= -static    # do not remove just yet
EXTRA_LIBS=-lpthread -lsvml
COMP_LIBS=
RANLIB=echo
#
NETCDF_ROOT=/share/apps/netcdf-3.6.2-ifort
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
LIBS=-L/opt/intel/mkl/10.0.3.020/lib/em64t \
     -lmkl_scalapack -lmkl_blacs_intelmpi20_lp64 \
     -lmkl_lapack -lmkl_em64t -lguide $(EXTRA_LIBS) $(NETCDF_LIBS)
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








