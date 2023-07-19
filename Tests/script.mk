#
# Single-test makefile template for script usage
#
SIESTA=../../../siesta
#
completed:
	@echo ">>>> Running $(name) test..."
	@if [ -d work ] ; then rm -rf work ; fi; mkdir work
	@if [ -f script.sh ] ; then cp -f script.sh work ; fi
	@echo "    ==> Running script with SIESTA as ${SIESTA}"
	@(cd work ; sh script.sh "${SIESTA}")  && touch completed
	@if [ -f completed ] ; then \
           echo "    ===> Script finished successfully";\
         else \
           echo " **** Test $(name) did not complete successfully";\
         fi
#
xmlcheck: completed
	@echo "---- xmllint check $(name).xml ..."
	xmllint $(name).xml > /dev/null
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf work completed $(name).out $(name).xml
