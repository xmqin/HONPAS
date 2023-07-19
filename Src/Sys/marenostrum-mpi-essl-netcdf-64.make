# --- Host ---
SIESTA_ARCH = marenostrum-mpi-essl-netcdf-64
SYS         = xlf

# --- Compiler ---
FC          = mpif90
FFLAGS      = -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto -q64
# FFLAGS_DEBUG = -g -C -qflttrap -qfloat=nans -qnosave -qundef

# --- BLAS & LAPACK ; BLACS & SCALAPACK ---
SCALAPACK64 = /gpfs/apps/SCALAPACK/lib64
SCALAPACK   = -lscalapack
BLACS       = $(SCALAPACK64)/blacsF77init_MPI-PPC-0.a \
              $(SCALAPACK64)/blacsCinit_MPI-PPC-0.a \
              $(SCALAPACK64)/blacs_MPI-PPC-0.a
LAPACK64    = /gpfs/apps/LAPACK/lib64
COMP_LIBS   = essl_lapack.a
BLAS        = -lessl
MATH_LIBS   = -L$(SCALAPACK64) -L$(LAPACK64) \
              $(SCALAPACK) $(BLACS) $(COMP_LIBS) $(BLAS)

# --- MPI ---
MPI_INTERFACE = libmpi_f90.a
MPI_INCLUDE = .
FPPFLAGS_MPI = -WF,-DMPI

# --- NETCDF ---
NETCDF_ROOT = /gpfs/apps/NETCDF/4.1.1/64
NETCDF_INCFLAGS = -I$(NETCDF_ROOT)/include
NETCDF_LIBS = -L$(NETCDF_ROOT)/lib -lnetcdf
FPPFLAGS_CDF = -WF,-DCDF

# --- FoX is off ---
DUMMY_FOX=--enable-dummy

# --- Preprocessor & linker ---
FPPFLAGS    = $(FPPFLAGS_MPI) $(FPPFLAGS_CDF)
LIBS        = $(MATH_LIBS) $(NETCDF_LIBS)
LDFLAGS     = -q64
RANLIB      = echo

# --- Compilation rules ---
.F90.o:
	$(FC) -qsuffix=cpp=F90       -c $(INCFLAGS) $(FFLAGS) $(FPPFLAGS) $<
.f90.o:
	$(FC) -qsuffix=f=f90         -c $(INCFLAGS) $(FFLAGS)             $<
.F.o:
	$(FC) -qsuffix=cpp=F -qfixed -c $(INCFLAGS) $(FFLAGS) $(FPPFLAGS) $<
.f.o:
	$(FC) -qsuffix=f=f   -qfixed -c $(INCFLAGS) $(FFLAGS)             $<

