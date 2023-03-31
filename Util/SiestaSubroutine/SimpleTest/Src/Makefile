# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Client-Server example
#
# JM Soler, A Garcia, E Anglada
#
WORK_DIR=$(shell pwd)
SIESTA_DIR=$(WORK_DIR)/../../../..
SRC_DIR=$(SIESTA_DIR)/Src
OBJDIR=Obj
OBJ_DIR=$(SIESTA_DIR)/$(OBJDIR)
#
.SUFFIXES:
.SUFFIXES: .c .f .F .o .a .f90 .F90
#
VPATH=$(SRC_DIR)
#
#  simple_pipes_serial     : Simple serial test with pipes version
#  simple_pipes_parallel   : Simple parallel test with pipes version
#  simple_sockets_serial   : Simple serial test with sockets version
#  simple_sockets_parallel : Simple parallel test with sockets version
#  simple_mpi_serial       : Simple serial test with MPI version
#  simple_mpi_parallel     : Simple parallel test with MPI version
#
# default: simple_pipes_serial simple_sockets_serial clean

all: default
default: simple_pipes_serial simple_pipes_parallel \
	simple_sockets_serial simple_sockets_parallel \
	simple_mpi_serial simple_mpi_parallel

ARCH_MAKE=$(OBJ_DIR)/arch.make
include $(ARCH_MAKE)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
# Routines used by pipes and sockets versions (VPATH knows where they are)
#
FSIESTA_PIPES= fsiesta_pipes.o pxf.o
FSIESTA_SOCKETS= fsockets.o sockets.o fsiesta_sockets.o 
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
WXML=libwxml.a
XC=libSiestaXC.a
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(LIBS)
#
libs_collected:
	(cd $(OBJ_DIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WORK_DIR) )
	touch libs_collected
#
# Pipes version
#
simple_pipes_serial: $(FSIESTA_PIPES) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_pipes_serial $(FSIESTA_PIPES) simple_serial.o 
#
simple_pipes_parallel: $(FSIESTA_PIPES) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_pipes_parallel \
	      $(FSIESTA_PIPES) simple_parallel.o
#
# Sockets version
#
simple_sockets_serial: $(FSIESTA_SOCKETS) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_sockets_serial \
	      $(FSIESTA_SOCKETS) simple_serial.o 
#
simple_sockets_parallel: $(FSIESTA_SOCKETS) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_sockets_parallel \
	      $(FSIESTA_SOCKETS) simple_parallel.o
#
# MPI version
#
phonons: libs_collected phonons.o 
	$(FC) $(LDFLAGS) -o phonons phonons.o $(ALL_LIBS)
#
simple_mpi_serial: libs_collected simple_serial.o 
	$(FC) $(LDFLAGS) -o simple_mpi_serial simple_serial.o $(ALL_LIBS)
#
simple_mpi_parallel: libs_collected simple_mpi_parallel.o 
	$(FC) $(LDFLAGS) -o simple_mpi_parallel \
	      simple_mpi_parallel.o $(ALL_LIBS)
#
intermediate: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make clean to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: clean
clean: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f simple_pipes_serial simple_pipes_parallel 
	rm -f simple_sockets_serial simple_sockets_parallel 
	rm -f simple_mpi_serial simple_mpi_parallel
	rm -f phonons
	rm -f *.a *.o *.mod libs_collected

