#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Use first sources here, then in the top level
VPATH:=$(VPATH_ROOT)/SiestaXC:$(VPATH_ROOT)
#
default: module
#
# Note: The copy operation below is only a possible way to
#       make the module and library files accessible.
#
module:  libSiestaXC.a
	cp *.mod *.a ..
#
SIESTA_SRCS= alloc.F90 bessph.f cellsubs.f  debugmpi.F \
             chkgmx.f fft1d.F interpolation.f90  \
             minvec.f m_io.f moreParallelSubs.F90 \
             class_OrbitalDistribution.F90 \
             m_timer.F90 timer_tree.F90 m_walltime.f90 parallel.F \
	     precision.F pxf.F90 radfft.f m_fft_gpfa.F sorting.f  timer.f90

SIESTA_OBJS:= $(SIESTA_SRCS:.f=.o)
SIESTA_OBJS:= $(SIESTA_OBJS:.F=.o)
SIESTA_OBJS:= $(SIESTA_OBJS:.f90=.o)
SIESTA_OBJS:= $(SIESTA_OBJS:.F90=.o)
#
SIESTA_SOURCES:=$(addprefix ../, $(SIESTA_SRCS))

LOCAL_SRCS= atomxc.F90 cellxc.F90 debugxc.f90 ggaxc.F ldaxc.f \
            siestaxc.F90 vdwxc.F90 xcmod.f90  fft3d.F90 fftr.F90 \
            mesh1d.f90 mesh3d.F90 am05.f90 vv_vdwxc.F90 array.F90 local_sys.F
ifdef WITH_LIBINT
LOCAL_SRCS += hsexc.f pbe0xc.f xc_xwpbe.F90
endif
LOCAL_OBJS:= $(LOCAL_SRCS:.f=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.f90=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F90=.o)
#
ALL_OBJS=$(SIESTA_OBJS) $(LOCAL_OBJS)
ALL_SOURCES=$(SIESTA_SOURCES) $(LOCAL_SRCS)
#
libSiestaXC.a: $(ALL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libSiestaXC.a $(ALL_OBJS) 
	@if [ ! -z "$(MPI_INTERFACE)" ] ; then \
            $(AR) $(ARFLAGS_EXTRA) cru libSiestaXC.a ../MPI/*.o ; fi
	-$(RANLIB) libSiestaXC.a
#
OTHER_FILES=lib-makefile lib-arch.make mpi-makefile lib-00-README test-makefile
#
.PHONY: dist
dist: siestaXC.tgz
siestaXC.tgz: $(ALL_SOURCES) $(OTHER_FILES) 
	rm -rf SiestaXC
	mkdir -p SiestaXC
	cp -fp $(ALL_SOURCES) SiestaXC
	cp -fp lib-makefile SiestaXC/makefile
	cp -fp lib-arch.make SiestaXC
	cp -frp ../MPI SiestaXC
	cp -fp mpi-makefile SiestaXC/MPI/Makefile
	cp -frp Testers SiestaXC
	cp -fp test-makefile SiestaXC/Testers/makefile
	cp BSD_LICENSE SiestaXC
	cp lib-00-README SiestaXC/00_README
	tar -cvzf siestaXC.tgz SiestaXC
	rm -rf SiestaXC
#
clean:
	rm -fr *.o *.*d *.a siestaXC.tgz
#
# DO NOT DELETE THIS LINE - used by make depend
alloc.o: debugmpi.o m_io.o parallel.o precision.o local_sys.o
am05.o: precision.o local_sys.o
atomxc.o: alloc.o ggaxc.o ldaxc.o mesh1d.o precision.o radfft.o
atomxc.o: local_sys.o vdwxc.o xcmod.o
bessph.o: precision.o local_sys.o
cellsubs.o: precision.o
cellxc.o: alloc.o cellsubs.o chkgmx.o debugxc.o fftr.o ggaxc.o ldaxc.o
cellxc.o: m_timer.o mesh3d.o moreParallelSubs.o parallel.o precision.o local_sys.o
cellxc.o: vdwxc.o xcmod.o
chkgmx.o: cellsubs.o minvec.o precision.o local_sys.o
debugxc.o: m_io.o moreParallelSubs.o parallel.o
debugmpi.o: parallel.o
fft1d.o: parallel.o precision.o local_sys.o
fft3d.o: alloc.o debugxc.o fft1d.o m_timer.o mesh3d.o precision.o local_sys.o
fftr.o: alloc.o fft3d.o mesh3d.o precision.o local_sys.o
ggaxc.o: ldaxc.o precision.o local_sys.o am05.o
interpolation.o: precision.o
ldaxc.o: precision.o local_sys.o
m_io.o: local_sys.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o local_sys.o
mesh1d.o: precision.o
mesh3d.o: alloc.o debugxc.o parallel.o precision.o sorting.o array.o local_sys.o
minvec.o: cellsubs.o precision.o sorting.o local_sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o local_sys.o
parallel.o: class_OrbitalDistribution.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
siestaxc.o: atomxc.o cellxc.o mesh3d.o precision.o xcmod.o
local_sys.o: parallel.o 
siesta_options.o: parallel.o precision.o units.o
timer.o: local_sys.o m_timer.o timer_tree.o 
timer_tree.o: m_walltime.o
vdwxc.o: alloc.o debugxc.o interpolation.o ggaxc.o ldaxc.o mesh1d.o
vdwxc.o: precision.o radfft.o vv_vdwxc.o local_sys.o
vv_vdwxc.o: alloc.o debugxc.o interpolation.o mesh1d.o
vv_vdwxc.o: precision.o radfft.o local_sys.o
xcmod.o: precision.o local_sys.o vdwxc.o
ifdef WITH_LIBINT
hsexc.o pbe0xc.o: xc_xwpbe.o
endif
