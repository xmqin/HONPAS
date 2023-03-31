FC=mpif90
FC_SERIAL=gfortran

# Optimal flags
FFLAGS=-O3 -m64 -fPIC -fno-second-underscore -ftree-vectorize -fexpensive-optimizations \
	-finline-functions-called-once -funroll-loops -fvariable-expansion-in-unroller \
	-ftree-loop-optimize -frename-registers -fprefetch-loop-arrays -finline-small-functions \
	-fipa-pure-const -foptimize-sibling-calls -fipa-cp

# Default flags
FFLAGS=-O2 -m64 -fPIC

# This is for debugging purposes
#FFLAGS = -g -O0 -Warray-bounds -Wunused

PP = cpp -E -P -C -nostdinc

C_V = gnu-4.7.2
MPI_PATH=/opt/openmpi/1.8.3/$(C_V)
NCDF_PATH=/opt/netcdf/4.3.2/$(C_V)
Z_PATH=/opt/zlib/1.2.8/$(C_V)
HDF5_PATH=/opt/hdf5/1.8.12/$(C_V)
PNCDF_PATH=/opt/pnetcdf/1.5.0/$(C_V)

INC = -I$(Z_PATH)/include \
         -I$(HDF5_PATH)/include \
         -I$(PNCDF_PATH)/include \
         -I$(NCDF_PATH)/include \
         -I$(MPI_PATH)/include

LIB_PATH= -L$(Z_PATH)/lib \
	  -L$(HDF5_PATH)/lib \
	  -L$(PNCDF_PATH)/lib \
	  -L$(NCDF_PATH)/lib \
	  -Wl,-rpath=$(Z_PATH)/lib \
	  -Wl,-rpath=$(HDF5_PATH)/lib \
	  -Wl,-rpath=$(PNCDF_PATH)/lib \
	  -Wl,-rpath=$(NCDF_PATH)/lib \
	  -L$(MPI_PATH)/lib \
	  -Wl,-rpath=$(MPI_PATH)/lib


LIBS = $(LIB_PATH) -lnetcdff -lnetcdf -lpnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz

#FPPFLAGS += -DNCDF_4
FPPFLAGS += -DNCDF_PARALLEL -DNCDF_4
#FPPFLAGS += -DNCDF_PARALLEL

.F90.o:
	$(FC) -c $(INC) $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(INC) $(FFLAGS) $<

