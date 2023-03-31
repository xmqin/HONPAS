#
# Single-test makefile template for script usage
#
# We do not have popd and pushd in "sh" scripts, thus force the SHELL to be bash
SHELL=/bin/bash

MPI=mpirun -np 2
TS = ../../../../siesta


# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(TS))),"mpirun")
MPI=
endif
ifeq ($(strip $(firstword $(TS))),"mpiexec")
MPI=
endif

completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d work ] ; then rm -rf work ; fi; mkdir work
	@if [ -f script.sh ] ; then cp -f script.sh work ; fi
	@echo "    ==> Running script with TranSIESTA as $(TS)"
	@(cd work ; $(SHELL) script.sh "$(MPI) $(TS)" )
	@if [ -f completed ] ; then \
           echo "    ===> Script finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf work completed *.out 
