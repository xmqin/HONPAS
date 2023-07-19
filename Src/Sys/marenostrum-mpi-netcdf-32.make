SIESTA_ARCH=mn-xlf-32bits-PARALLEL-netcdf
#
FC= mpif90   #xlf90_r
#
FFLAGS_DEBUG= -g
FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q32
FFLAGS_parse=-qsuffix=f=f -qfree #-qfixed
LDFLAGS= -q32
COMP_LIBS=
RANLIB=echo
DEFS_PREFIX=-WF,
#
NETCDF_ROOT=/gpfs/apps/NETCDF/3.6.2/32
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
#
NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
NETCDF_INTERFACE=
FPPFLAGS_CDF=-WF,-DCDF
#
MPI_INTERFACE=    libmpi_f90.a
MPI_INCLUDE=.
MPI_LIBS=         #-lblacsgm
FPPFLAGS_MPI=         -WF,-DMPI

BLAS= -L/gpfs/apps/SCALAPACK/lib32 -lscalapack \
        /gpfs/apps/SCALAPACK/lib32/blacsF77init_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib32/blacsCinit_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib32/blacs_MPI-PPC-0.a \
      -L/gpfs/apps/LAPACK/lib32 -llapack \
      -L/gpfs/apps/SCALAPACK/lib32 -lblas

##MPITRACER= -L/gpfs/apps/CEPBATOOLS/lib/32 -lmpitrace
#
##LIBS=  $(MPITRACER) $(BLAS) $(NETCDF_LIBS)
#
LIBS=  $(BLAS) $(NETCDF_LIBS)
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
