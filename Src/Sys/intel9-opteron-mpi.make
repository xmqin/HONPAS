SIESTA_ARCH=intel9-opteron-mpi
#
# arch.make based on one created by Lucas Fernandez Seivane, quevedin@gmail.com
# You may need to change the name of the compiler, location of libraries...
# Modified by Alberto Garcia to suit the opteron boxes at ICMAB
#
FC=mpiifort
FC_ASIS=$(FC)
#
# You can experiment with more aggressive optimizations (?)
#
FFLAGS=-O2 -mp
FFLAGS_DEBUG= -g -O0
RANLIB=echo 
MPI_INCLUDE=/apl/INTEL/ICT3/mpi/3.0/include64
MPI_INTERFACE=libmpi_f90.a
FPPFLAGS_MPI=-DMPI
#
# MPI_LIBS not needed if a "mpiXXX" form of the compiler is used
# (make sure to source the 
#  /somepath/to/ICT3/mpi/3.0/bin64/mpivars.sh file to enable it)
#
#MPI_LIBS=-L/apl/INTEL/ICT3/mpi/3.0/lib64 -lmpiif -lmpi  # -ldl
#
# Scalapack and Blacs were home-compiled
#
SCALAPACK=-L/apl/OPTERON/SCALAPACK/ -lscalapack
BLACS_ROOT=/apl/OPTERON/BLACS/LIB
BLACS=$(BLACS_ROOT)/blacsF77init.a $(BLACS_ROOT)/blacsCinit.a \
      $(BLACS_ROOT)/blacs.a
#
# Use AMD's library for lapack and blas references
# (make sure to put /somepath/to/ACML3.6.0/ifort64/lib 
#
LAPACK=-L/apl/OPTERON/ACML3.6.0/ifort64/lib -lacml
LIBS=$(SCALAPACK) $(BLACS)  $(LAPACK)  $(MPI_LIBS)

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
