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
# Makefile for Client-Server example
#
# JM Soler, A Garcia
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH:=$(shell pwd)/../../../../Src
#
#       driver: Mixed molecular dynamics
#       simple: Simple test
#       para  : Test driver for parallel calculation
#
default: simple driver para
#
ARCH_MAKE=../../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
FSIESTA_OBJ=../../../../Src/fsiesta.o      # Client helper routines
FSIESTA_MOD=../../../../Src/fsiesta.mod    # Client helper module
PXF_OBJ=../../../../Src/pxf.o              # Portable flush

#
driver.o: fmixmd.o sample.o
fmixmd.o: sample.o
#
OBJS=fsiesta.o sample.o fmixmd.o driver.o pxf.o
driver : $(OBJS)
	$(FC) -o driver  $(LDFLAGS) $(OBJS) 
#
simple: fsiesta.o simple.o  pxf.o
	$(FC) $(LDFLAGS) -o simple simple.o fsiesta.o pxf.o
#
para: fsiesta.o para.o pxf.o
	$(FC) $(LDFLAGS) -o para para.o fsiesta.o pxf.o
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f driver simple para *.o *.mod



