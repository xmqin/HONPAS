SIESTA_ARCH=intel9-cmkl8-mpi
#
# arch.make created by Lucas Fernandez Seivane, quevedin@gmail.com
# You may need to change the name of the compiler, location of libraries...
# Modified by Alberto Garcia to suit cryst at the UPV.
#
LANG=
FC=mpiifort
FC_ASIS=$(FC)
#
FFLAGS=-O3 -xP -mp
FFLAGS_DEBUG= -g -O0
RANLIB=echo 
MPI_INCLUDE=/opt/intel/mpi/2.0/include
MPI_INTERFACE=libmpi_f90.a
FPPFLAGS_MPI=-DMPI
#
NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2-ifort
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF         # Note this !!
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#
LIBS=-L/opt/intel/cmkl/8.0/lib/32  -lmkl_scalapacktesting_intel80  \
      -lmkl_scalapack -lmkl_blacs_intelmpi20  \
      -lmkl_lapack -lmkl_ia32 -lguide -lpthread -lrt -lsvml $(NETCDF_LIBS)
SYS=nag
FPPFLAGS= $(FPPFLAGS_MPI) $(FPPFLAGS_CDF)
#
.F.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)   $<
#
