.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
# Choose compiler here.
# Optimization flags are not really necessary, but do include
# the OpenMP compiler flag if your system supports it and it is worth it.
#
##FC=ifort -openmp    # -g -O1 -CB -ftrapuv
#
##FC=gfortran  -fopenmp 
#
#
SIMPLEX_OBJS=io.o parse.o minimizer.o vars_module.o simplex.o \
     sys.o precision.o amoeba.o
SWARM_OBJS=io.o parse.o minimizer.o vars_module.o swarm.o \
     sys.o precision.o 
#
default: swarm simplex
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
simplex: $(SIMPLEX_OBJS)
	$(FC) -o $@ $(SIMPLEX_OBJS)
#
swarm: $(SWARM_OBJS)
	$(FC) -o $@ $(SWARM_OBJS)


clean:
	rm -f *.o *.mod swarm simplex

.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#


# DO NOT DELETE THIS LINE - used by make depend
amoeba.o: precision.o vars_module.o
minimizer.o: vars_module.o
simplex.o: amoeba.o minimizer.o precision.o vars_module.o
swarm.o: minimizer.o precision.o vars_module.o
v_test.o: vars_module.o
vars_module.o: parse.o precision.o sys.o
m_amoeba.o: amoeba.o
