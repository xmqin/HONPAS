#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                      --- Run all tests
#          make tests=fe_broyden     --- Run a single test
#          make batch                --- Dispatch to batch system
#                                        (site-dependent)
#                                        (needs extra configuration)
#
#          By default, the version of SIESTA used is the one in Src.
#          To use another one:
#
#          make SIESTA=/some/path/to/siesta-program ....
#
#          It is also possible to have separate working directories,
#          by using the a "label". For example:
#
#          make label=finer 
#     
#          The default label is "work".
#
# See also test.mk for details, and the scripts in Scripts for 
# examples of other options for batch processing.
#
#
REFERENCE=../../ref-version/Tests
#
label=work
#
ALL= h2o h2o_reparam h2o_basis h2o_dos h2o-pop      		\
       floating benzene bessel mgco3 si2x1h force_2 born born_spin	\
       var_cell constant_volume batio3 fe fe_broyden sih si64		\
       h2oZ sih_op_broyden h2o_op_broyden zmatrix md_pr md_npr		\
       md_anneal md_verlet md_nose si_bandpoints sih_fire		\
       graphite_c6 oxyn h2o_findp_bug h2o_radialgrid ch4	\
       si_coop born pb_filter_cutoff pb_filter_tol fire_benzene		\
       fe_clust_noncollinear h_chain n_chain si64_coop si_pdos_kgrid    \
       si2x1h-dipole sic-slab si-optical
#
#
#  Need to add the new tests and generate the reference outputs and xml files.
#
#
all:
	for dir in $(ALL); do\
	 $(MAKE) -C $$dir completed; \
	done

clean:
	@for i in $(ALL) ; do \
		(cd $$i ; make clean) ; \
	done
#


