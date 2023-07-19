#
#  Makefile for fdf example and fdf module.
#  
.PHONY: check_flags
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
INCFLAGS:=-I$(VPATH) $(INCFLAGS)    # For VPATH operation
#
# Include copying operations in rule to make sure that
# they are always performed.
#
default: module
module:  libfdf.a
	cp libfdf.a ..
	@cp -p *.*d ..

check_flags:
	@echo "In fdf, INCFLAGS is: $(INCFLAGS)"
#
libfdf.a: check_flags fdf.o fdf_mod.o parse.o
	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a fdf.o fdf_mod.o parse.o
	-$(RANLIB) libfdf.a
#
fdf_mod.o: fdf_mod.f parse.o
fdf.o:		fdf.h
#
#---------------------------------------------------
OBJS=	sample.o fdf.o io_for_fdf_tests.o
#
sample: $(OBJS)
	$(FC) $(LDFLAGS) -o sample $(OBJS)
#
check:	sample sample.fdf
	sample < sample.fdf
#
io_sample: io_sample.o io_for_fdf_tests.o
	$(FC) $(LDFLAGS) -o io_sample io_sample.o io.o
#
io_check: io_sample
	io_sample
clean:
	rm -f *.o libfdf.a *.mod fdf_mod.o io_sample.o sample io_sample
	rm -f OUTPUT FDF.debug file2 file3 sample.out
	rm -f *.pc *.pcl
#
sample.o:	fdfdefs.h



