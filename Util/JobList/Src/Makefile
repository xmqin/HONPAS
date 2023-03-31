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
default: countJobs runJobs getResults horizontal
#
include ../../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
COUNTJOBS_OBJS=jobList.o countJobs.o
RUNJOBS_OBJS=jobList.o runJobs.o
GETRESULTS_OBJS=jobList.o getResults.o
HORIZONTAL_OBJS=jobList.o horizontal.o
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
countJobs: $(COUNTJOBS_OBJS)
	$(FC) -o $@ $(COUNTJOBS_OBJS)
#
runJobs: $(RUNJOBS_OBJS)
	$(FC) -o $@ $(RUNJOBS_OBJS)
#
getResults: $(GETRESULTS_OBJS)
	$(FC) -o $@ $(GETRESULTS_OBJS)
#
horizontal: $(HORIZONTAL_OBJS)
	$(FC) -o $@ $(HORIZONTAL_OBJS)
#
clean:
	rm -f *.o *.mod countJobs runJobs getResults horizontal

# DO NOT DELETE THIS LINE - used by make depend
countJobs.o: jobList.o
getResults.o: jobList.o
runJobs.o: jobList.o
