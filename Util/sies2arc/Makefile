# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for sies2arc utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj

all: sies2arc

include ../../$(OBJDIR)/arch.make

# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
##FFLAGS=$(FFLAGS_DEBUG)

OBJS = sies2arc.o linepro.o wtof.o cell.o uncell.o
$(OBJS) : constants

sies2arc: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

clean: 
	rm -f *.o sies2arc
