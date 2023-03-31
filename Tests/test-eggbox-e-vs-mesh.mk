#
# Single-test makefile template
#

GridMin=50
GridMax=500
GridIncrement=75

$(label)-completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for i in `cat $(name).pseudos` ; do \
		echo "    ==> Copying pseudopotential file for $$i..." ;\
		ln ../Pseudos/$$i.psf $(label)/$$i.psf ;\
	done 

	@(cd $(label) ; ../../Scripts/e-vs-mesh.sh ${SIESTA} $(name) $(GridMin) $(GridMax) $(GridIncrement) ) \
		&& touch $(label)-completed
	@if [ -f $(label)-completed ]; then \
	   echo "    ===> SIESTA finished";\
	   else \
	   echo " **** Test $(name) did not complete successfully";\
	fi

xmlcheck: pseudos $(label)-completed
#${SIESTA} 2>&1 > $(name).out < ../$(name).fdf
