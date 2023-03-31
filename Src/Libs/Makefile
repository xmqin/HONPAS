# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.


default: libsiestaLAPACK.a libsiestaBLAS.a

ARCH_MAKE?=../arch.make
include $(ARCH_MAKE)

###########
# These are now the preferred complete BLAS and LAPACK
# sources.
# They are extracted using the linalg2file.py script.
#
# See create.sh for how the lapack.f and blas.f files are created.
#
# Note that the python script is extremely basic, but functions well for
# F77 code.

libsiestaBLAS.a: blas.o
	$(AR) $(ARFLAGS_EXTRA) cru $@ $<
	$(RANLIB) $@

libsiestaLAPACK.a: lapack.o
	$(AR) $(ARFLAGS_EXTRA) cru $@ $<
	$(RANLIB) $@

clean:
	rm -f lapack.o blas.o *.a
