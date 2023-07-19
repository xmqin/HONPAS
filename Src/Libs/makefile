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
#
#  Makefile for Libraries
#  Some systems only need the stuff in dc_lapack.f (Divide and Conquer)
#  TranSiesta needs zgesv_lapack.f
#
default: linalg.a
#
ARCH_MAKE?=../arch.make
include $(ARCH_MAKE)
#
ALL_OBJS= lapack.o dc_lapack.o zgesv_lapack.o zgeev_lapack.o blas.o machine.o
# Add machine.o here only if needed...
DC_OBJS=  dc_lapack.o   # machine.o
# Used in TranSIESTA
ZGESV_OBJS=  zgesv_lapack.o   # machine.o
# Used in TBTrans
ZGEEV_OBJS= zgeev_lapack.o
#
ESSL_OBJS= essl_lapack.o
#
default: linalg.a
#
libblas.a: blas.o
	$(AR) $(ARFLAGS_EXTRA) cru libblas.a blas.o
	$(RANLIB) libblas.a
	mv libblas.a ..
liblapack.a: lapack.o machine.o
	$(AR) $(ARFLAGS_EXTRA) cru liblapack.a lapack.o machine.o
	$(RANLIB) liblapack.a
	mv liblapack.a ..
linalg.a: $(ALL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru linalg.a $(ALL_OBJS)
	$(RANLIB) linalg.a
	mv linalg.a ..
dc_lapack.a: $(DC_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru dc_lapack.a $(DC_OBJS)
	$(RANLIB) dc_lapack.a
	mv dc_lapack.a ..
#
zgesv_lapack.a: $(ZGESV_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru zgesv_lapack.a $(ZGESV_OBJS)
	$(RANLIB) zgesv_lapack.a
	mv zgesv_lapack.a ..
#
zgeev_lapack.a: $(ZGEEV_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru zgeev_lapack.a $(ZGEEV_OBJS)
	$(RANLIB) zgeev_lapack.a
	mv zgeev_lapack.a ..
#
essl_lapack.a: $(ESSL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru essl_lapack.a $(ESSL_OBJS)
	$(RANLIB) essl_lapack.a
	mv essl_lapack.a ..
#
clean:
	rm -f $(ALL_OBJS) *.a




