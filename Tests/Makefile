#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                     --- Run all tests
#          make check               --- Run all tests and check with reference file
#          make check-only          --- Run all checks with reference file (run make before)
#          make show                --- Show the list of tests
#
#          make tests=fe_broyden { check | }     --- Run a single test
#          (Not implemented) make batch                --- Dispatch to batch system
#                                        (site-dependent)
#                                        (needs extra configuration)
#
#          By default, the version of SIESTA used is the one in the
#          parent directory (Obj or custom OBJDIR), and is executed
#          in parallel mode with 2 cores.
#
#          To change the MPI run command simply do:
#
#          make MPI="mpirun -np 4"
#
#          for 4 cores. Or
#
#          make MPI=
#
#          for a serial run.
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

# Default target
default: all

label = work

tests = 32_h2o ag anneal-cont ar2_vdw batio3 benzene bessel bessel-rich \
	born born_spin carbon_nanoscroll ch4 chargeconf-h2o \
	constant_volume dipole_correction fe fe_broyden \
	fe_clust_noncollinear fe_clust_noncollinear-gga fe_cohp \
	fen fe_noncol_kp fire_benzene floating \
	force_2 force_constants gate_G_charge gate_G_hartree ge111 \
	ge_fatbands_so	graphite_c6 graphite_c6_full graphite_vdw_df \
	h2_bessel h2o h2o_2 h2o_am05 h2o_bands h2o_bands_nc \
	h2o_bands_polarized h2o_basis h2o_coop h2o_dipole h2o_dipole2 \
	h2o_dos h2o_filteret_basis h2o_findp_bug h2o_netcdf h2o_op_broyden \
	h2o_orderN h2o-pop h2o_radialgrid h2o_reparam h2o-rhog-mixing h2oZ \
	h3po4_1 h3po4_2 h3po4_3 h_chain h_chain2 md_anneal md-anneal-cont \
	md-cont md-cont-no-restart-file md_nose md_npr md_pr md-single-run \
	md_verlet mgc-force mgco3 mix_broyden mix_linear mix_pulay mno_ldau \
	nanotube-c-5-0 n_chain OMM_h2o OMM_si oxyn partial pb_bulk pb_bulk_1000 \
	pb_filter_cutoff pb_filter_tol pd_ice ptcda ptcda-au scfmix \
	si001 si001-diags si001-ldos si111-spinpol si2x1h si2x1h-dipole \
	si2x1h-dipole-gcs si2x1h-quench si64 si64_coop si_bandpoints \
	si_coop	sic-slab si_fatbands sih sih_fire sih-mrrr sih_op_broyden \
	sinw sinw_2 si-optical si_pdos_gamma si_pdos_kgrid \
	var_cell wannier zmatrix

# Tests that won't work with more than a few processors
tests_mpi_max2 = h2_bessel h_chain h_chain2 \
	born born_spin \
	force_2 md_anneal n_chain

# Tests that have been disabled, currently
#  h2o_4 (ON, does not work in parallel)
#  ipi-test (requires ipi)
#  si64-pexsi
#  sih-pexsi
#  sih-pexsi-spin


# These tests are extremely time consuming and
# should only be runned sometimes
# Currently they may be runned individually.
#	SOC_FePt_xx SOC_FePt_xz SOC_FePt_zy SOC_FePt_zz
tests_soc  = SOC_Pt2_xx SOC_Pt2_xz SOC_Pt2_zy SOC_Pt2_zz
tests_soc += SOC_FePt_xx SOC_FePt_xz SOC_FePt_zy SOC_FePt_zz

# Tests only applicable for LUA
tests_lua  = lua_si111 lua_h2o

# Tests only applicable for FLOS
tests_flos  = flos_h2o_lbfgs flos_h2o_fire flos_h2o_meshcutoff flos_h2o_neb

# Accumulated variable for ALL tests
tests_all = $(tests) $(tests_soc) $(tests_lua) $(tests_flos)

all: $(tests)

# Define test specific targets:
tests: $(tests)
tests-soc: $(tests_soc)
tests-lua: $(tests_lua)
tests-flos: $(tests_flos)
tests-proc-limit: $(tests_mpi_max2)

.NOTPARALLEL: tests-ts
.PHONY: tests-ts
tests-ts:
	(cd TranSiesta-TBTrans ; MPI="$(MPI)" bash ./run_tests.sh)

# Ensure that tests may be called individually
# (and that they are run, PHONY is required)
# This ensures one can call:
#   make h2o
# and the h2o test will be runned.
.NOTPARALLEL: $(tests_all)
.PHONY: $(tests_all)
$(tests_all):
	$(MAKE) -C $@ completed

.PHONY: check
check:
	@for dir in $(tests) ; do\
	 $(MAKE) -C $$dir check; \
	done

check-only:
	@for dir in $(tests) ; do\
	 $(MAKE) -C $$dir check-only ; \
	done

.PHONY: show
show:
	@echo $(tests)

clean:
	@for i in $(tests_all) ; do \
		(cd $$i ; make clean) ; \
	done
