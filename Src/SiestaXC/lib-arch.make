# 
# given in BSD_LICENSE
#
SIESTA_ARCH=macosx-openmpi
#
FC=/opt/openmpi.gcc.gfortran/bin/mpif90
FC_SERIAL=gfortran -m64
#
FFLAGS=-O0 -g -fbacktrace
FFLAGS_DEBUG= -g -Wall -Wextra
LDFLAGS=
RANLIB=echo
#
DEFS=-DGFORTRAN  -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DDEBUG_XC    # Note this !!
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
DEFS_MPI=-DMPI
#
DEFS:=$(DEFS_MPI) $(DEFS)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
