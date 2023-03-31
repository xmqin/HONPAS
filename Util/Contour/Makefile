# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for Contour utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj

all: grid1d grid2d

include ../../$(OBJDIR)/arch.make

# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
##FFLAGS=$(FFLAGS_DEBUG)

grid1d: grid1d.o
	$(FC) $(LDFLAGS) -o $@  grid1d.o
grid2d: grid2d.o
	$(FC) $(LDFLAGS) -o $@  grid2d.o

clean: 
	rm -f *.o grid1d grid2d
