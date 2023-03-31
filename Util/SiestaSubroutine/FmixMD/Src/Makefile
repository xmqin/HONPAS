# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
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
#
driver.o: fmixmd.o sample.o
fmixmd.o: sample.o
#
OBJS=fsiesta_pipes.o sample.o fmixmd.o driver.o pxf.o
driver : $(OBJS)
	$(FC) -o driver  $(LDFLAGS) $(OBJS) 
#
simple: fsiesta_pipes.o simple.o  pxf.o
	$(FC) $(LDFLAGS) -o simple simple.o fsiesta_pipes.o pxf.o
#
para: fsiesta_pipes.o para.o pxf.o
	$(FC) $(LDFLAGS) -o para para.o fsiesta_pipes.o pxf.o
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f driver simple para *.o *.mod



