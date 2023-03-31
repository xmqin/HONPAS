# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
OBJDIR=Obj
#
default: swarm simplex
#
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
# Optimization flags are not really necessary, but do include
# the OpenMP compiler flag if your system supports it and it is worth it.
# Examples:
#
##FC=ifort -openmp    # -g -O1 -CB -ftrapuv
#
##FC=gfortran  -fopenmp 
#
#
SIMPLEX_OBJS=io.o parse.o minimizer.o vars_module.o simplex.o \
     sys.o precision.o amoeba.o
SWARM_OBJS=io.o parse.o minimizer.o vars_module.o swarm.o \
     sys.o precision.o 
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
simplex: $(SIMPLEX_OBJS)
	$(FC) -o $@ $(SIMPLEX_OBJS)
#
swarm: $(SWARM_OBJS)
	$(FC) -o $@ $(SWARM_OBJS)


clean:
	rm -f *.o *.mod swarm simplex

# DO NOT DELETE THIS LINE - used by make depend
amoeba.o: precision.o vars_module.o
minimizer.o: vars_module.o
simplex.o: amoeba.o minimizer.o precision.o vars_module.o
swarm.o: minimizer.o precision.o vars_module.o
vars_module.o: parse.o precision.o sys.o
m_amoeba.o: amoeba.o
