#
#
FC= f77
CFLAGS=
FFLAGS=
LDFLAGS=
#
SRCS= bhspot.f potran.f symbols.f symbols.h set2.h plot.h makefile README \
      BHSINP
OBJS= bhspot.o potran.o symbols.o
#
bhs:	$(OBJS)
	$(FC) -o bhs $(LDFLAGS) $(OBJS)
#
symbols.o: set2.h
bhspot.o:  symbols.h
potran.o:  plot.h
#
dist:  
	tar cvf bhs.tar $(SRCS)
	gzip -vf bhs.tar
#
clean:
	/bin/rm -rf $(OBJS) bhs
#
