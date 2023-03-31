FC=mpif90
FC_SERIAL=ifort

# Optimal flags
FFLAGS=-O3 -m64 -fPIC -xHost -fp-model strict

# Default flags
FFLAGS=-O2 -m64 -fPIC

# This is for debugging purposes
#FFLAGS = -g -O0 -check bounds traceback -warn unused -fp-model strict

PP = fpp -P

C_V = intel-13.1.1
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

AR=xiar

#FPPFLAGS += -DNCDF_4
FPPFLAGS += -DNCDF_PARALLEL -DNCDF_4
#FPPFLAGS += -DNCDF_PARALLEL

.F90.o:
	$(FC) -c $(INC) $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(INC) $(FFLAGS) $<

