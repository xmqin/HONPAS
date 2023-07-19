# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
##VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
#
all: orbmol_proj
#
include ../../$(OBJDIR)/arch.make
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
##FFLAGS=$(FFLAGS_DEBUG)
#
#------------------
orbmol_proj: orbmol_proj.o
	$(FC) $(LDFLAGS) -o $@  orbmol_proj.o
#
#------------------
clean: 
	rm -f *.o orbmol_proj *.o *.*d
#










