SIESTA_ARCH=XLF 64bits PARALLEL
#
FC= mpif90   #xlf90_r
#
FFLAGS_DEBUG= -g
FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q64
FFLAGS_parse=-qsuffix=f=f -qfree #-qfixed
LDFLAGS= -q64
COMP_LIBS=
RANLIB=echo
DEFS_PREFIX=-WF,
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=    libmpi_f90.a
MPI_INCLUDE=.
MPI_LIBS=         #-lblacsgm
FPPFLAGS_MPI=         -WF,-DMPI



#BLAS = -L/gpfs/apps/LAPACK/lib64 -llapack \
#	     -L/gpfs/apps/SCALAPACK/lib64 -lblas


BLAS= -L/gpfs/apps/SCALAPACK/lib64 -lscalapack \
        /gpfs/apps/SCALAPACK/lib64/blacsF77init_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib64/blacsCinit_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib64/blacs_MPI-PPC-0.a \
      -L/gpfs/apps/LAPACK/lib64 -llapack \
      -L/gpfs/apps/SCALAPACK/lib64 -lblas

##MPITRACER= -L/gpfs/apps/CEPBATOOLS/lib/64 -lmpitrace
#
##LIBS=  $(MPITRACER) $(BLAS)
LIBS=  $(BLAS)
SYS=xlf
FPPFLAGS= $(FPPFLAGS_MPI) $(FPPFLAGS_CDF)
FREE_F90=-qsuffix=f=f90
#
.F90.o:
	$(FC) -qsuffix=cpp=F90 -c $(INCFLAGS) $(FFLAGS) $(FPPFLAGS) $<
.f90.o:
	$(FC) -qsuffix=f=f90 -c $(INCFLAGS) $(FFLAGS)   $<
.F.o:
	$(FC) -qsuffix=cpp=F -c $(INCFLAGS) -qfixed $(FFLAGS) $(FPPFLAGS) $<
.f.o:
	$(FC) -qsuffix=f=f -qfixed -c $(INCFLAGS) $(FFLAGS)   $<
#
