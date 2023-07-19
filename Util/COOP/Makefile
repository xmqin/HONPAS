#
# The VPATH directive below helps to re-use f2kcli.F90 from the top Src.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
#
all: mprop dm_creator fat
#
include ../../$(OBJDIR)/arch.make
#
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
INCFLAGS:=$(INCFLAGS) $(NETCDF_INCFLAGS)
#
MPROP_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys_local.o read_curves.o io_hs.o orbital_set.o main_vars.o
FAT_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys_local.o read_curves.o io_hs.o orbital_set.o main_vars.o
DM_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys_local.o iodm_netcdf.o io_hs.o  main_vars.o write_dm.o
#
default: mprop
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
mprop: $(MPROP_OBJS) mprop.o
	$(FC) -o $@ $(LDFLAGS) $(MPROP_OBJS) mprop.o
fat: $(FAT_OBJS) fat.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) fat.o
dm_creator: $(DM_OBJS) dm_creator.o
	$(FC) -o $@ $(LDFLAGS) $(DM_OBJS) dm_creator.o $(NETCDF_LIBS)
#
clean:
	rm -f *.o *.mod mprop  dm_creator

# DO NOT DELETE THIS LINE - used by make depend
dm_creator.o: io_hs.o iodm_netcdf.o main_vars.o
io.o: sys_local.o
io_hs.o: main_vars.o precision.o
iodm.o: precision.o sys_local.o
iodm_netcdf.o: precision.o sys_local.o
m_getopts.o: f2kcli.o
main_vars.o: m_getopts.o precision.o subs.o units.o
mprop.o: io_hs.o main_vars.o orbital_set.o read_curves.o
old_dm.o: io_hs.o iodm_netcdf.o main_vars.o orbital_set.o
orbital_set.o: main_vars.o subs.o
read_curves.o: orbital_set.o precision.o
subs.o: precision.o
units.o: precision.o
write_dm.o: precision.o sys_local.o
fat.o: mprop.o
