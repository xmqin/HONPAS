#
# Architecture
#
FDF_ARCH=XLF 64bits PARALLEL

#
# Compiling, linking flags
#
FC= mpif90
FFLAGS= -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q32
LDFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q32

#
# Library utils
#
#AR=
RANLIB= echo

#
# MPI flags
#
MPIFLAGS= $(MPI_PATH) $(MPI_INCLUDE) $(MPI_LIBS)

MPI_PATH=
MPI_INCLUDE=
MPI_LIBS=

#
# Macro definitions
#
# DEBUG:    Full debugging information in library
# _MPI_:    Runs FDF under a MPI environment
# CLUSTER:  FDF on a Cluster (non-shared filesystem)
# BLOCKING: Blocking reading on a huge MPI execution
#
DEFS= $(DEFS_DEBUG) $(DEFS_MPI) $(DEFS_IO)

DEFS_DEBUG= #-WF,-DDEBUG
DEFS_MPI=   -WF,-D_MPI_
DEFS_IO=    #-WF,-DCLUSTER or #-WF,-DBLOCKING
