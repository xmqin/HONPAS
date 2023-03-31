#
#  Makefile for Sies2xcf set of tools
#
# Compiler and options:
   FC     = ifort
   FFLAGS = -mp1 -O3 -pc80 -prec_div -w
#  FC=f95
#  FFLAGS=-g -C=all -C=undefined
  LFLAGS = $(FFLAGS)
#
  EXECS = eig2bxsf xv2xsf md2axsf rho2xsf vib2xsf
  .f.o      : $(FC) -c $(FFLAGS) $<

all: $(EXECS)
   eig_obj = eig2bxsf.o inver3.o opnout.o 
   xv_obj  = xv2xsf.o  opnout.o
   md_obj  = md2axsf.o makebox.o  fillbox.o inver3.o  hit.o \
             opnout.o  test_ani.o test_md.o wraxsf1.o wraxsf2.o
   rho_obj = rho2xsf.o read_xv.o  makebox.o fillbox.o inver3.o  hit.o \
             intpl04.o opnout.o 
   vib_obj = vib2xsf.o read_xv.o  makebox.o fillbox.o inver3.o  hit.o \
	     displa.o  opnout.o   read_ev.o itochar.o w_arrow.o w_movie.o

 eig2bxsf : $(eig_obj)
	$(FC) $(LFLAGS) -o $@ $(eig_obj)
  xv2xsf   : $(xv_obj)
	$(FC) $(LFLAGS) -o $@ $(xv_obj)
  md2axsf  : $(md_obj)
	$(FC) $(LFLAGS) -o $@ $(md_obj)
  rho2xsf  : $(rho_obj)
	$(FC) $(LFLAGS) -o $@ $(rho_obj)
  vib2xsf  : $(vib_obj)
	$(FC) $(LFLAGS) -o $@ $(vib_obj)

clean:
	@echo " Cleaning up"
	rm -f $(EXECS) *.o core

