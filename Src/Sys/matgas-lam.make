SIESTA_ARCH=matgas-lam
#
# arch.make based on one created by Lucas Fernandez Seivane, quevedin@gmail.com
# You may need to change the name of the compiler, location of libraries...
# Modified by Alberto Garcia to suit the opteron boxes at ICMAB
#
FC=ifort
FC_ASIS=$(FC)
#
# You can experiment with more aggressive optimizations (?)
#
FFLAGS=-O2 -mp
FFLAGS_DEBUG= -g -O0
RANLIB=echo 
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/apl/lam/include
MPI_LIBS=-L/apl/lam/lib -llamf77mpi -lmpi -llam -lguide -lpthread -lsvml
FPPFLAGS_MPI=-DMPI 
#
MKLPATH=/apl/INTEL/ICT3/cmkl/9.0/lib/32
GUIDE=$(MKLPATH)/libguide.a
BLACS=$(MKLPATH)/libmkl_blacs_intelmpi20.a
SCALAPACK=$(MKLPATH)/libmkl_scalapack.a
BLAS=$(MKLPATH)/libmkl_ia32.a
LIBS=-L$(MKLPATH) -lmkl_scalapack -lmkl_blacs_intelmpi -lmkl_lapack -lmkl_ia32 $(MPI_LIBS)

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
