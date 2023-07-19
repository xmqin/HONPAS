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
# Makefile for siesta
#
# Remove implicit rules (including for Modula-2, which just confuses matters)
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH=.
#
default: what honpas #  fsiesta.o
#
configure: configure.ac aclocal.m4
	autoconf ${<D}/configure.ac
#
M4_FILES = m4/ac_config_aux_dirs.m4  m4/TW_CHECK_FC_FPP.m4      m4/TW_FC_ID_FLAGS.m4  \
 m4/ACX_MPI.m4             m4/TW_CHECK_FC_TR15580.m4  m4/TW_FC_ID.m4        \
 m4/fortran.m4             m4/TW_CHECK_FC_TR15581.m4  m4/TW_FIND_FC_BLAS.m4 \
 m4/TW_CHECK_BLACS.m4      m4/TW_CHECK_SCALAPACK.m4   m4/TW_FIND_LAPACK.m4  \
 m4/TW_CHECK_FC_90.m4      m4/TW_FC_CHECK_ABORT.m4    m4/TW_PATH_NETCDF.m4  \
 m4/TW_CHECK_FC_95.m4      m4/TW_FC_CHECK_DCFUNS.m4   m4/TW_TRY_DC_LAPACK.m4\
 m4/TW_CHECK_FC_FPP_90.m4  m4/TW_FC_CHECK_FLUSH.m4    m4/TW_FC_KINDS.m4
#
aclocal.m4: $(M4_FILES)
	cat ${<D}/$(M4_FILES)  > aclocal.m4
dep:
	@if [ $(VPATH) = "." ] ; then \
	sfmakedepend --depend=obj  --modext=o \
                              *.f *.f90 *.F *.F90 ; \
        else \
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 ; \
        fi
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
#
include arch.make
INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Uncomment the following line for debugging support
#
##FFLAGS=$(FFLAGS_DEBUG)
#
.PHONY: .honpas what version
#
what:
	@echo
	@if [ "$(VPATH)" = . ]; then \
              echo "** You can no longer build SIESTA in Src." ;\
              echo "** Go to the Obj directory and see the README file." ;\
              false ;\
         fi
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Src/Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o pxf.o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS =  automatic_cell.o atom_options.o \
        arw.o  atomlwf.o bands.o basis_enthalpy.o bessph.o bonds.o \
        born_charge.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o constr.o coxmol.o cross.o compute_norm.o \
	m_partial_charges.o \
	denmat.o denmatlomem.o detover.o dfscf.o diagon.o digcel.o fft3d.o \
	diagg.o diagk.o diagk_file.o diagkp.o diag2g.o diag2k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o doping.o dot.o  \
	m_efield.o egandd.o ener3.o ener3lomem.o extrapolon.o \
	fixed.o flib_spline.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o inver.o iodm_netcdf.o \
	iodmhs_netcdf.o iogrid_netcdf.o iolwf.o iozm.o ipack.o iopipes.o \
	iowfs_netcdf.o kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o metaforce.o \
        m_check_supercell.o m_filter.o mulliken.o minvec.o m_pulay.o naefs.o \
	new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o \
	ran3.o reclat.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o cgvc_zmatrix.o m_convergence.o \
        iocg.o ioeig.o iofa.o iokp.o iomd.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o kpoint_grid.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_init.o struct_init.o siesta_options.o siesta_geom.o \
        siesta_analysis.o siesta_end.o siesta_forces.o \
        siesta.o io.o sparse_matrices.o \
        coor.o atm_transfer.o \
        eggbox.o dsyevds.o zheevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o show_distribution.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o obc.o \
        precision.o sys.o m_cell.o recipes.o files.o \
        spatial.o parallel.o parallelsubs.o  parsing.o chemical.o xcmod.o atmparams.o \
        m_mpi_utils.o m_fdf_global.o m_history.o m_iorho.o  \
        atmfuncs.o listsc.o memoryinfo.o m_memory.o sorting.o \
        atomlist.o old_atmfuncs.o radial.o \
        alloc.o phonon.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o \
        molecularmechanics.o  zm_broyden_optim.o \
        cell_broyden_optim.o remove_intramol_pressure.o \
        m_ntm.o m_dipol.o  m_kinetic.o  m_rmaxh.o \
        m_forces.o m_stress.o m_eo.o m_spin.o m_gamma.o m_hsx.o m_fire.o \
        get_target_stress.o fire_optim.o zm_fire_optim.o \
        cell_fire_optim.o m_fire_para.o m_fire_mixing.o write_raw_efs.o \
        m_iodm.o m_sparse.o kpoint_pdos.o \
        domain_decom.o schecomm.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
        setup_ordern_indexes.o \
        atom.o atm_types.o dhscf.o initatom.o spher_harm.o xc.o \
        xc_xwpbe.o cell_pbc.o  broadcast_basis.o \
        nao2gto_transfer.o nao2gto_kinds.o nao2gto_types.o nao2gto_index.o \
        nao2gto_gamma_F.o nao2gto_expint.o nao2gto_erf_fn.o \
        nao2gto_mathconstants.o nao2gto_libint_wrapper_types.o \
        nao2gto_libint_wrapper.o nao2gto_update_hfx.o  nao2gto_update_hfx_force.o \
        nao2gto_full_DM.o nao2gto_hfx_potential.o nao2gto_hfx_gradient.o \
        nao2gto_prescreen.o nao2gto_eri.o nao2gto_gradient.o \
        powell.o nao2gto_util.o t_c_g0.o \
        nao2gto_contract_eri.o nao2gto_primitive_eri.o \
        nao2gto_nonlin.o
#--------------------------------------------------------------------------------------
# FoX whatnot
#
# FoX build targets:
FoX_configured=FoX/.config
FoX_built=FoX/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`FoX/FoX-config --fcflags`
FoX_LIBS=`FoX/FoX-config --libs --wcml`
#
# And add them to global compiler flags:
INCFLAGS:=$(INCFLAGS) $(FoX_FCFLAGS)
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
$(FoX_configured):
	(cd FoX; touch arch.make ; \
         CONFIGURE="$(VPATH)/FoX/configure"; \
         $$CONFIGURE VPATH="$(VPATH)/FoX" \
         FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS)" \
         --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd FoX; $(MAKE))
#
# Now ensure that building FoX is inserted into the dependency chain.
# Routines using FoX:
# version (compinfo.o) does not actually need FoX, but making it build after
# FoX is done avoids a warning as we use the FoX-config script everywhere.
#
version $(SYSOBJ) $(OBJS) \
basis_io.o m_efield.o initatom.o ioeig.o kpoint_grid.o meshsubs.o \
mulliken.o ordern.o outcell.o pdos.o post_scf_work.o radial.o   \
scfconvergence_test.o setup_hamiltonian.o siesta_analysis.o     \
siesta_cml.o siesta_cmlsubs.o siesta_end.o siesta_forces.o      \
siesta_init.o siesta_move.o siesta_options.o state_analysis.o   \
state_init.o struct_init.o write_subs.o xml.o: $(FoX_built)
#
#--------------------------------------------------------------
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) "VPATH=$(VPATH)/xmlparser" module)
#
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" module)
#
# Routines using fdf calls.
#
dhscf.o iolwf.o  grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_io.o basis_specs.o: $(FDF)
#
nao2gto_transfer.o nao2gto_hfx_potential.o : $(FDF)
#
ALL_OBJS=$(OBJS) $(SYSOBJ)
#HONPAS_OBJS=$(OBJS) $(SYOBJ)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) "VPATH=$(VPATH)/MPI")
#
# Libraries that might need to be compiled
#
libblas.a: Libs/blas.f
	@echo "==> Compiling libblas.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" libblas.a)
liblapack.a: Libs/lapack.f
	@echo "==> Compiling liblapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" liblapack.a)
linalg.a: Libs/blas.f Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" linalg.a)
dc_lapack.a:  Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" dc_lapack.a)
#
essl_lapack.a: Libs/essl_lapack.f
	@echo "==> Compiling essl_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" essl_lapack.a)
#

version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;s'SIESTA_VERSION'\
                                       $$(cat ${<D}/../version.info)'g" $<  > compinfo.F90
	@($(MAKE) "FPPFLAGS=$(FPPFLAGS)" "DEFS=$(DEFS)" compinfo.o)
	@mv compinfo.o version.o
	#@rm -f compinfo.F90
	@echo
#
check-siesta:
	@if [ -f .transiesta ]; then \
              echo "** You were building transiesta in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .siesta
#

honpas: check-siesta what version $(MPI_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) \
                $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o honpas \
               $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) $(MPI_INTERFACE)\
               $(COMP_LIBS) $(FoX_LIBS) $(LIBS)

siesta: check-siesta what version $(MPI_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) \
                $(COMP_LIBS) $(ALL_OBJS) 
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) $(MPI_INTERFACE)\
               $(COMP_LIBS) $(FoX_LIBS) $(LIBS) 
#
#------------------------------------------------------------------------
#----------------------------------------------------------------------------
# TranSIESTA...
#
TS_OBJS= m_ts_global_vars.o m_ts_options.o m_ts_aux_rout.o m_ts_kpoints.o \
         m_ts_iodm.o m_ts_io.o m_ts_scattering.o m_ts_contour.o \
         m_ts_electrode.o m_ts_in_siesta.o
#
transiesta: FPPFLAGS+=$(DEFS_PREFIX)-DTRANSIESTA 
transiesta: DEFS+=$(DEFS_PREFIX)-DTRANSIESTA
check-ts:
	@if [ -f .siesta ]; then \
              echo "** You were building plain siesta in this directory." ;\
              echo "** Please type 'make clean' first" ;\
              false ;\
	fi
	@touch .transiesta
#
# Library that might need to be compiled
#
zgesv_lapack.a: Libs/zgesv_lapack.f
	@echo "==> Compiling zgesv_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" zgesv_lapack.a)
#
transiesta: check-ts what version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF)  \
                    $(COMP_LIBS) $(ALL_OBJS) $(TS_OBJS)
	$(FC) -o transiesta \
               $(LDFLAGS) $(ALL_OBJS) $(TS_OBJS) $(FDF) \
               $(MPI_INTERFACE) $(NETCDF_INTERFACE) \
               $(COMP_LIBS) $(FoX_LIBS) $(LIBS)
#
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta honpas transiesta *.o  *.a *.pcl *.pc
	rm -f .siesta .honpas .transiesta
	rm -f *.mod
	rm -f *.bck *~
	rm -f aux_*.f aux_*.f90
	(cd fdf ; $(MAKE) clean)
	(cd FoX ; touch arch.make ; $(MAKE) distclean)
	(cd xmlparser ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd Libs ; $(MAKE) clean )
#
##pristine:  clean
##	cp -fp Include/constr.f .
#

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: parallel.o precision.o sys.o

arw.o: parallel.o sys.o

atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o

atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o

atom_options.o: sys.o

atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o

atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o

automatic_cell.o: atm_types.o atmfuncs.o parallel.o precision.o units.o

bands.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o m_fdf_global.o
bands.o: m_gamma.o parallel.o parallelsubs.o parsing.o precision.o
bands.o: siesta_geom.o sys.o writewave.o

basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o

basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o sys.o xml.o

basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o

basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o

bessph.o: parallel.o precision.o sys.o

bonds.o: alloc.o atmfuncs.o mneighb.o precision.o recipes.o sorting.o units.o

born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o

broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_memory.o m_mpi_utils.o
broyden_optim.o: parallel.o precision.o units.o


cart2frac.o: sys.o

cdiag.o: alloc.o m_spin.o memoryinfo.o parallel.o precision.o sys.o

cell_broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_memory.o
cell_broyden_optim.o: m_mpi_utils.o parallel.o precision.o sys.o units.o
cell_broyden_optim.o: zmatrix.o

cell_fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_memory.o parallel.o
cell_fire_optim.o: precision.o siesta_options.o sys.o zmatrix.o

cellxc.o: alloc.o meshsubs.o parallel.o parallelsubs.o precision.o sys.o
cellxc.o: xcmod.o

cgvc.o: conjgr_old.o m_fdf_global.o m_mpi_utils.o parallel.o precision.o
cgvc.o: units.o

cgvc_zmatrix.o: conjgr.o m_fdf_global.o m_mpi_utils.o parallel.o precision.o
cgvc_zmatrix.o: units.o zmatrix.o

cgwf.o: egandd.o onmod.o precision.o sys.o

chemical.o: precision.o sys.o

chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o

chkdim.o: sys.o

chkgmx.o: minvec.o parallel.o precision.o sys.o

coceri.o: files.o periodic_table.o precision.o

compute_dm.o: atomlist.o iodmhs_netcdf.o kpoint_grid.o m_energies.o m_eo.o
compute_dm.o: m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_ts_global_vars.o m_ts_in_siesta.o m_ts_options.o ordern.o
compute_dm.o: parallel.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o sys.o units.o

compute_norm.o: m_mpi_utils.o precision.o sparse_matrices.o

conjgr.o: precision.o

conjgr_old.o: precision.o


coor.o: alloc.o m_mpi_utils.o parallel.o precision.o siesta_geom.o sys.o
coor.o: zmatrix.o

coxmol.o: files.o periodic_table.o precision.o



denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o

denmatlomem.o: alloc.o globalise.o onmod.o precision.o

densematrix.o: precision.o

detover.o: alloc.o parallel.o parallelsubs.o precision.o

dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
dfscf.o: meshsubs.o meshsubs.o parallel.o precision.o sys.o

diag2g.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o

diag2k.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o

diagg.o: fermid.o parallel.o parallelsubs.o precision.o sys.o writewave.o

diagk.o: compute_norm.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagk.o: writewave.o

diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o: sys.o writewave.o

diagkp.o: fermid.o parallel.o parallelsubs.o precision.o sys.o

diagon.o: alloc.o densematrix.o m_memory.o parallel.o parallelsubs.o
diagon.o: precision.o siesta_options.o sys.o

diagpol.o: precision.o sys.o

diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o


dipole.o: parallel.o precision.o sys.o


dnaefs.o: atmfuncs.o mneighb.o precision.o

domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o

doping.o: alloc.o m_fdf_global.o meshsubs.o meshsubs.o parallel.o precision.o
doping.o: sys.o



dynamics.o: atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o precision.o
dynamics.o: sys.o units.o

egandd.o: alloc.o denmat.o ener3.o globalise.o gradient.o m_mpi_utils.o
egandd.o: on_subs.o onmod.o onmod.o precision.o sys.o

eggbox.o: parallel.o parsing.o precision.o

electrostatic.o: atm_types.o atmfuncs.o atmparams.o precision.o radfft.o
electrostatic.o: radial.o recipes.o sys.o

ener3.o: globalise.o m_mpi_utils.o onmod.o precision.o

ener3lomem.o: globalise.o m_mpi_utils.o onmod.o precision.o

extrapolon.o: parallel.o precision.o sys.o


fermid.o: parallel.o precision.o recipes.o sys.o

fft3d.o: alloc.o parallel.o parallelsubs.o precision.o sys.o


find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o

fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_mpi_utils.o parallel.o
fire_optim.o: precision.o siesta_options.o units.o

fixed.o: alloc.o parallel.o parsing.o precision.o

flib_spline.o: precision.o

forhar.o: alloc.o meshsubs.o precision.o



get_target_stress.o: m_mpi_utils.o parallel.o precision.o units.o

globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o

gradient.o: globalise.o m_mpi_utils.o onmod.o precision.o

gradientlomem.o: globalise.o onmod.o precision.o

grdsam.o: m_partial_charges.o parallel.o parsing.o precision.o siesta_options.o
grdsam.o: sys.o

hsparse.o: alloc.o atm_types.o atmfuncs.o listsc.o mneighb.o parallel.o
hsparse.o: parallelsubs.o precision.o siesta_options.o sorting.o sys.o



idiag.o: parallel.o sys.o

initparallel.o: alloc.o domain_decom.o m_fdf_global.o parallel.o parallelsubs.o
initparallel.o: precision.o siesta_options.o spatial.o sys.o


io.o: sys.o

iocg.o: files.o precision.o sys.o

iodm.o: files.o parallel.o parallelsubs.o precision.o sys.o

iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o

iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o

ioeig.o: files.o precision.o siesta_cml.o units.o

iofa.o: files.o precision.o

iogrid_netcdf.o: alloc.o parallel.o precision.o sys.o

iokp.o: files.o precision.o

iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o

iomd.o: files.o precision.o

iopipes.o: parallel.o precision.o sys.o

iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o

ioxv.o: files.o parallel.o precision.o

iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o

ipack.o: sys.o

kgrid.o: parallel.o precision.o units.o

kgridinit.o: minvec.o parallel.o precision.o units.o

kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o

kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o units.o

kpoint_pdos.o: find_kgrid.o minvec.o parallel.o precision.o siesta_options.o
kpoint_pdos.o: units.o

ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o parsing.o
ksv.o: precision.o sys.o

ksvinit.o: alloc.o parallel.o parsing.o precision.o



listsc.o: alloc.o

local_DOS.o: atomlist.o files.o kpoint_grid.o m_energies.o m_eo.o m_forces.o
local_DOS.o: m_gamma.o m_mpi_utils.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o sys.o

m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o

m_broyddj_nocomm.o: alloc.o m_fdf_global.o parallel.o precision.o sys.o

m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o

m_cell.o: precision.o siesta_geom.o units.o

m_check_supercell.o: minvec.o parallel.o precision.o


m_dipol.o: precision.o

m_efield.o: atmfuncs.o m_fdf_global.o m_mpi_utils.o parallel.o parallelsubs.o
m_efield.o: precision.o siesta_cml.o siesta_geom.o sys.o units.o

m_energies.o: precision.o

m_eo.o: precision.o

m_fdf_global.o: m_mpi_utils.o parallel.o precision.o


m_fire.o: parallel.o precision.o

m_fire_mixing.o: alloc.o m_fdf_global.o m_fire_para.o m_mpi_utils.o parallel.o
m_fire_mixing.o: precision.o

m_fire_para.o: m_mpi_utils.o parallel.o precision.o

m_forces.o: precision.o


m_getopts.o: f2kcli.o


m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o sys.o

m_iodm.o: alloc.o files.o m_mpi_utils.o parallel.o parallelsubs.o precision.o
m_iodm.o: sys.o

m_iorho.o: parallel.o parallelsubs.o precision.o sys.o

m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o

m_kinetic.o: precision.o

m_memory.o: memoryinfo.o parallel.o sys.o

m_mpi_utils.o: precision.o sys.o


m_options.o: precision.o

m_partial_charges.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o parallel.o
m_partial_charges.o: precision.o siesta_geom.o siesta_options.o sys.o

m_pulay.o: alloc.o atomlist.o m_spin.o parallel.o precision.o siesta_options.o
m_pulay.o: sparse_matrices.o sys.o

m_rmaxh.o: precision.o

m_sparse.o: alloc.o precision.o

m_spin.o: alloc.o m_fdf_global.o parallel.o precision.o


m_stress.o: precision.o

m_ts_aux_rout.o: precision.o sys.o

m_ts_contour.o: files.o m_ts_aux_rout.o m_ts_options.o precision.o sys.o

m_ts_electrode.o: files.o m_ts_aux_rout.o m_ts_io.o m_ts_kpoints.o
m_ts_electrode.o: m_ts_options.o precision.o

m_ts_global_vars.o: precision.o

m_ts_in_siesta.o: find_kgrid.o m_ts_contour.o m_ts_electrode.o m_ts_io.o
m_ts_in_siesta.o: m_ts_kpoints.o m_ts_options.o m_ts_scattering.o parallel.o
m_ts_in_siesta.o: parallelsubs.o precision.o sys.o

m_ts_io.o: files.o m_spin.o m_ts_kpoints.o parallel.o parallelsubs.o
m_ts_io.o: precision.o sys.o

m_ts_iodm.o: alloc.o files.o m_mpi_utils.o m_ts_global_vars.o parallel.o
m_ts_iodm.o: parallelsubs.o precision.o sys.o

m_ts_kpoints.o: minvec.o parallel.o precision.o siesta_cml.o siesta_options.o
m_ts_kpoints.o: sys.o units.o

m_ts_options.o: m_fdf_global.o m_ts_global_vars.o parallel.o precision.o
m_ts_options.o: siesta_options.o sys.o units.o

m_ts_scattering.o: m_ts_aux_rout.o m_ts_options.o parallelsubs.o precision.o
m_ts_scattering.o: sys.o

m_wallclock.o: m_walltime.o


madelung.o: parallel.o precision.o sys.o

matel.o: alloc.o atmfuncs.o precision.o radfft.o recipes.o spher_harm.o

md_out.o: files.o m_energies.o precision.o sys.o units.o

md_utils.o: precision.o

memory.o: memoryinfo.o parallel.o precision.o

memoryinfo.o: precision.o

meshmatrix.o: alloc.o atomlist.o meshsubs.o parallel.o parallelsubs.o
meshmatrix.o: precision.o

meshsubs.o: atm_types.o atmfuncs.o parallel.o parallelsubs.o precision.o
meshsubs.o: radial.o siesta_cml.o

metaforce.o: parallel.o parsing.o precision.o

minvec.o: precision.o sorting.o sys.o

mixer.o: atomlist.o iodm_netcdf.o m_broyden_mixing.o m_energies.o
mixer.o: m_fire_mixing.o m_iodm.o m_mpi_utils.o m_pulay.o m_spin.o m_steps.o
mixer.o: m_ts_global_vars.o m_ts_iodm.o m_ts_options.o parallel.o siesta_geom.o
mixer.o: siesta_options.o sparse_matrices.o


mneighb.o: alloc.o precision.o sys.o

molecularmechanics.o: alloc.o parallel.o parsing.o precision.o sys.o units.o

mulliken.o: atmfuncs.o parallel.o parallelsubs.o precision.o siesta_cml.o

naefs.o: atmfuncs.o mneighb.o precision.o


new_dm.o: alloc.o atomlist.o m_energies.o m_fdf_global.o m_iodm.o m_mpi_utils.o
new_dm.o: m_sparse.o m_spin.o m_steps.o m_ts_global_vars.o m_ts_iodm.o
new_dm.o: parallel.o parallelsubs.o parsing.o precision.o siesta_geom.o
new_dm.o: siesta_options.o sparse_matrices.o sys.o units.o

nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o

normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: sparse_matrices.o

obc.o: alloc.o precision.o

ofc.o: alloc.o files.o precision.o units.o

old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o

on_subs.o: alloc.o onmod.o onmod.o

onmod.o: precision.o

optical.o: alloc.o atomlist.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o parsing.o precision.o sys.o units.o

ordern.o: alloc.o cgwf.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
ordern.o: precision.o siesta_cml.o sys.o units.o

outcell.o: precision.o siesta_cml.o units.o

outcoor.o: alloc.o atmfuncs.o precision.o siesta_geom.o sys.o units.o

overfsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o

overlap.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o


parallelsubs.o: domain_decom.o parallel.o spatial.o sys.o



pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
pdos.o: parallelsubs.o precision.o siesta_geom.o units.o xml.o

pdosg.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o

pdosk.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o

pdoskp.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o

periodic_table.o: precision.o sys.o

phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o

phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o

phonon.o: files.o parallel.o precision.o sys.o

pixmol.o: files.o periodic_table.o precision.o

plcharge.o: files.o old_atmfuncs.o parallel.o precision.o

poison.o: alloc.o parallel.o precision.o sys.o

post_scf_work.o: atomlist.o m_energies.o m_spin.o mneighb.o setup_hamiltonian.o
post_scf_work.o: siesta_geom.o siesta_options.o sparse_matrices.o


printmatrix.o: alloc.o

projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o
projected_DOS.o: m_gamma.o m_mpi_utils.o m_spin.o parallel.o precision.o
projected_DOS.o: siesta_geom.o siesta_options.o sparse_matrices.o sys.o

propor.o: precision.o sys.o

proximity_check.o: m_ts_options.o mneighb.o parallel.o precision.o
proximity_check.o: siesta_geom.o siesta_options.o

pseudopotential.o: alloc.o atom_options.o flib_spline.o precision.o sys.o




radfft.o: alloc.o precision.o recipes.o

radial.o: precision.o xml.o


rdiag.o: alloc.o memoryinfo.o parallel.o precision.o sys.o

readsp.o: parallel.o precision.o sys.o

recipes.o: sys.o


redcel.o: precision.o sys.o

reinit.o: files.o parallel.o precision.o

remove_intramol_pressure.o: precision.o zmatrix.o

reoptical.o: parallel.o parsing.o precision.o sys.o

reord.o: alloc.o precision.o

rhooda.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o precision.o

rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofd.o: meshsubs.o meshsubs.o parallel.o precision.o sys.o

rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofdsp.o: meshsubs.o meshsubs.o precision.o sys.o

savepsi.o: alloc.o parallel.o parallelsubs.o precision.o


scfconvergence_test.o: m_convergence.o m_energies.o m_pulay.o
scfconvergence_test.o: m_ts_global_vars.o m_wallclock.o parallel.o siesta_cml.o
scfconvergence_test.o: siesta_options.o write_subs.o

schecomm.o: alloc.o

setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o

setspatial.o: alloc.o parallel.o precision.o spatial.o

setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o files.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_forces.o m_gamma.o m_hsx.o m_mpi_utils.o
setup_hamiltonian.o: m_ntm.o m_partial_charges.o m_spin.o m_steps.o m_stress.o
setup_hamiltonian.o: m_ts_global_vars.o m_ts_io.o m_ts_kpoints.o m_ts_options.o
setup_hamiltonian.o: metaforce.o molecularmechanics.o nlefsm.o normalize_dm.o
setup_hamiltonian.o: overfsm.o siesta_geom.o siesta_options.o sparse_matrices.o
setup_hamiltonian.o: sys.o

setup_kscell.o: minvec.o parallel.o precision.o units.o

setup_ordern_indexes.o: alloc.o domain_decom.o parallel.o spatial.o


shaper.o: alloc.o atmfuncs.o mneighb.o precision.o

show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o

siesta.o: m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o

siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dipol.o
siesta_analysis.o: m_energies.o m_eo.o m_fdf_global.o m_forces.o m_gamma.o
siesta_analysis.o: m_mpi_utils.o m_ntm.o m_partial_charges.o m_spin.o m_steps.o
siesta_analysis.o: parallel.o projected_DOS.o siesta_cml.o siesta_geom.o
siesta_analysis.o: siesta_options.o sparse_matrices.o write_subs.o writewave.o
siesta_analysis.o: zmatrix.o


siesta_cmlsubs.o: files.o parallel.o siesta_cml.o timestamp.o 

siesta_end.o: alloc.o m_wallclock.o parallel.o siesta_cml.o siesta_cmlsubs.o
siesta_end.o: timestamp.o

siesta_forces.o: compute_dm.o m_convergence.o m_steps.o m_ts_global_vars.o
siesta_forces.o: m_ts_options.o mixer.o parallel.o post_scf_work.o
siesta_forces.o: scfconvergence_test.o setup_hamiltonian.o siesta_cml.o
siesta_forces.o: siesta_options.o state_analysis.o state_init.o sys.o

siesta_geom.o: precision.o

siesta_init.o: alloc.o atomlist.o bands.o densematrix.o files.o ioxv.o
siesta_init.o: kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_fdf_global.o m_forces.o m_gamma.o
siesta_init.o: m_iostruct.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
siesta_init.o: m_ts_in_siesta.o m_ts_options.o m_wallclock.o metaforce.o
siesta_init.o: molecularmechanics.o parallel.o parallelsubs.o phonon.o
siesta_init.o: projected_DOS.o siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timestamp.o writewave.o
siesta_init.o: xcmod.o zmatrix.o

siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o iopipes.o ioxv.o
siesta_move.o: m_energies.o m_forces.o m_kinetic.o m_steps.o m_stress.o
siesta_move.o: parallel.o phonon.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_move.o: sys.o write_subs.o zm_broyden_optim.o zm_fire_optim.o zmatrix.o

siesta_options.o: files.o m_fdf_global.o m_mpi_utils.o memoryinfo.o parallel.o
siesta_options.o: precision.o siesta_cml.o sys.o units.o


sparse_matrices.o: precision.o

spatial.o: precision.o

state_analysis.o: atomlist.o born_charge.o m_energies.o m_forces.o m_iostruct.o
state_analysis.o: m_ntm.o m_spin.o m_stress.o m_wallclock.o parallel.o
state_analysis.o: remove_intramol_pressure.o siesta_cml.o siesta_geom.o
state_analysis.o: siesta_options.o sparse_matrices.o write_subs.o zmatrix.o

state_init.o: alloc.o atomlist.o domain_decom.o files.o hsparse.o iodm_netcdf.o
state_init.o: iodmhs_netcdf.o kpoint_grid.o m_check_supercell.o m_energies.o
state_init.o: m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o m_spin.o
state_init.o: m_steps.o m_ts_global_vars.o m_ts_io.o m_ts_options.o new_dm.o
state_init.o: overlap.o parallel.o proximity_check.o siesta_cml.o siesta_geom.o
state_init.o: siesta_options.o sparse_matrices.o sys.o write_subs.o zmatrix.o

struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_fdf_global.o m_iostruct.o m_mpi_utils.o parallel.o
struct_init.o: periodic_table.o siesta_cml.o siesta_geom.o siesta_options.o
struct_init.o: zmatrix.o

sys.o: parallel.o


timer.o: m_walltime.o parallel.o precision.o sys.o


transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o

typecell.o: precision.o

uncell.o: precision.o

units.o: precision.o


vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmat.o: meshsubs.o meshsubs.o parallel.o precision.o

vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmatsp.o: meshsubs.o meshsubs.o precision.o

vmb.o: parallel.o precision.o sys.o


write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o

write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o

write_subs.o: atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o m_ts_global_vars.o parallel.o precision.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o units.o zmatrix.o

writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o kpoint_grid.o
writewave.o: m_fdf_global.o parallel.o parallelsubs.o parsing.o precision.o
writewave.o: siesta_geom.o sys.o units.o

#xc.o: alloc.o precision.o sys.o xcmod.o

xcmod.o: parallel.o parsing.o precision.o sys.o


xml.o: precision.o


zm_broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_mpi_utils.o parallel.o
zm_broyden_optim.o: precision.o units.o zmatrix.o

zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o

zmatrix.o: alloc.o m_cell.o parallel.o parsing.o precision.o siesta_geom.o
zmatrix.o: sys.o units.o

band.o: bands.o
diagmemory.o: memoryinfo.o
listsc_module.o: listsc.o
m_born_charge.o: born_charge.o
m_broyden_optim.o: broyden_optim.o
m_cell_broyden_optim.o: cell_broyden_optim.o
m_cell_fire_optim.o: cell_fire_optim.o
m_cgwf.o: cgwf.o
m_chempot.o: chempot.o
m_compute_dm.o: compute_dm.o
m_conjgr.o: conjgr.o
m_conjgr_old.o: conjgr_old.o
m_denmat.o: denmat.o
m_doping.o: doping.o
m_dynamics.o: dynamics.o
m_eandg.o: egandd.o
m_ener3.o: ener3.o
m_fermid.o: fermid.o
m_find_kgrid.o: find_kgrid.o
m_fire_optim.o: fire_optim.o
m_forhar.o: forhar.o
m_gradient.o: gradient.o
m_hsparse.o: hsparse.o
m_intramol_pressure.o: remove_intramol_pressure.o
m_ioxv.o: ioxv.o
m_ksv.o: ksv.o
m_ksvinit.o: ksvinit.o
m_local_dos.o: local_DOS.o
m_minvec.o: minvec.o
m_mixer.o: mixer.o
m_new_dm.o: new_dm.o
m_nlefsm.o: nlefsm.o
m_norm.o: compute_norm.o
m_normalize_dm.o: normalize_dm.o
m_on_subs.o: on_subs.o
m_ordern.o: ordern.o
m_overfsm.o: overfsm.o
m_overlap.o: overlap.o
m_post_scf_work.o: post_scf_work.o
m_projected_dos.o: projected_DOS.o
m_proximity_check.o: proximity_check.o
m_radfft.o: radfft.o
m_recipes.o: recipes.o
m_scf_options.o: m_options.o
m_scfconvergence_test.o: scfconvergence_test.o
m_setup_hamiltonian.o: setup_hamiltonian.o
m_siesta_analysis.o: siesta_analysis.o
m_siesta_end.o: siesta_end.o
m_siesta_forces.o: siesta_forces.o
m_siesta_init.o: siesta_init.o
m_siesta_move.o: siesta_move.o
m_state_analysis.o: state_analysis.o
m_state_init.o: state_init.o
m_struct_init.o: struct_init.o
m_timestamp.o: timestamp.o
m_xyz_phiatm.o: matel.o
m_zm_broyden_optim.o: zm_broyden_optim.o
m_zm_fire_optim.o: zm_fire_optim.o
mesh.o: meshsubs.o
meshdscf.o: meshmatrix.o
meshphi.o: meshsubs.o
neighbour.o: mneighb.o
on_core.o: onmod.o
on_main.o: onmod.o
printmat.o: printmatrix.o
version_info.o: version.o
write_subs_energies.o: write_subs.o
write_subs_positions.o: write_subs.o
write_subs_pressure.o: write_subs.o

# HONPAS PART

atm_types.o: precision.o radial.o

atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o m_filter.o
atom.o: old_atmfuncs.o periodic_table.o precision.o pseudopotential.o recipes.o
atom.o: sys.o xcmod.o

broadcast_basis.o: atm_types.o parallel.o radial.o

dhscf.o: alloc.o atmfuncs.o doping.o files.o forhar.o iogrid_netcdf.o
dhscf.o: m_efield.o m_iorho.o m_partial_charges.o m_ts_global_vars.o
dhscf.o: m_ts_in_siesta.o m_ts_options.o meshsubs.o parallel.o parsing.o
dhscf.o: precision.o siesta_options.o sys.o units.o xcmod.o

initatom.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
initatom.o: electrostatic.o old_atmfuncs.o precision.o

spher_harm.o: alloc.o precision.o sys.o

xc.o: alloc.o precision.o sys.o xcmod.o xc_xwpbe.o
xc_xwpbe.o: nao2gto_erf_fn.o nao2gto_kinds.o
xc_xwpbe.o: nao2gto_mathconstants.o nao2gto_expint.o

cell_pbc.o: precision.o

t_c_g0.o: nao2gto_kinds.o

powell.o: nao2gto_kinds.o nao2gto_mathconstants.o

nao2gto_util.o: nao2gto_kinds.o nao2gto_mathconstants.o atm_types.o sys.o
nao2gto_transfer.o: atm_types.o atmfuncs.o sys.o nao2gto_mathconstants.o
nao2gto_transfer.o: parsing.o chemical.o nao2gto_util.o alloc.o
nao2gto_transfer.o: radial.o nao2gto_nonlin.o nao2gto_types.o

nao2gto_types.o: nao2gto_kinds.o m_fdf_global.o parallel.o xcmod.o

nao2gto_index.o: alloc.o
nao2gto_mathconstants.o: nao2gto_kinds.o
nao2gto_gamma_F.o: nao2gto_kinds.o nao2gto_mathconstants.o
nao2gto_erf_fn.o: nao2gto_kinds.o
nao2gto_expint.o: nao2gto_kinds.o

nao2gto_libint_wrapper.o: nao2gto_kinds.o  atm_types.o
nao2gto_libint_wrapper.o: nao2gto_libint_wrapper_types.o

nao2gto_primitive_eri.o: nao2gto_gamma_F.o nao2gto_libint_wrapper.o
nao2gto_primitive_eri.o: nao2gto_mathconstants.o nao2gto_kinds.o
nao2gto_primitive_eri.o: nao2gto_libint_wrapper_types.o
nao2gto_primitive_eri.o: nao2gto_types.o t_c_g0.o

nao2gto_contract_eri.o: nao2gto_kinds.o nao2gto_primitive_eri.o
nao2gto_contract_eri.o: nao2gto_libint_wrapper.o cell_pbc.o
nao2gto_contract_eri.o: nao2gto_libint_wrapper_types.o nao2gto_index.o

nao2gto_prescreen.o: atm_types.o atomlist.o sorting.o recipes.o \
                     atmfuncs.o nao2gto_kinds.o parallel.o alloc.o \
                     mneighb.o nao2gto_contract_eri.o listsc.o \
                     nao2gto_libint_wrapper_types.o nao2gto_index.o \
                     cell_pbc.o nao2gto_types.o powell.o \
                     nao2gto_primitive_eri.o nao2gto_full_DM.o

nao2gto_eri.o: nao2gto_kinds.o atmfuncs.o atm_types.o atomlist.o nao2gto_index.o
nao2gto_eri.o: alloc.o parallel.o parallelsubs.o parsing.o
nao2gto_eri.o: nao2gto_kinds.o nao2gto_gamma_F.o nao2gto_contract_eri.o
nao2gto_eri.o: nao2gto_mathconstants.o nao2gto_libint_wrapper_types.o
                         
nao2gto_hfx_potential.o: nao2gto_kinds.o atmfuncs.o parsing.o alloc.o
nao2gto_hfx_potential.o: atm_types.o atomlist.o nao2gto_index.o 
nao2gto_hfx_potential.o: parallel.o parallelsubs.o m_fdf_global.o
nao2gto_hfx_potential.o: nao2gto_kinds.o nao2gto_gamma_F.o 
nao2gto_hfx_potential.o: nao2gto_prescreen.o nao2gto_eri.o
nao2gto_hfx_potential.o: nao2gto_libint_wrapper_types.o
nao2gto_hfx_potential.o: nao2gto_libint_wrapper.o nao2gto_types.o
nao2gto_hfx_potential.o: atomlist.o parallelsubs.o

nao2gto_full_DM.o: precision.o atm_types.o atmfuncs.o atomlist.o listsc.o
nao2gto_update_hfx.o: precision.o atmfuncs.o parallel.o parallelsubs.o
nao2gto_update_hfx.o: alloc.o atm_types.o nao2gto_hfx_potential.o
nao2gto_update_hfx.o: nao2gto_full_DM.o

nao2gto_gradient.o: atmfuncs.o atm_types.o atomlist.o nao2gto_index.o
nao2gto_gradient.o: alloc.o parallel.o parallelsubs.o parsing.o
nao2gto_gradient.o: nao2gto_kinds.o nao2gto_gamma_F.o nao2gto_contract_eri.o
nao2gto_gradient.o: nao2gto_mathconstants.o nao2gto_libint_wrapper_types.o

nao2gto_hfx_gradient.o: atmfuncs.o parsing.o alloc.o
nao2gto_hfx_gradient.o: atm_types.o atomlist.o nao2gto_index.o
nao2gto_hfx_gradient.o: parallel.o parallelsubs.o m_fdf_global.o
nao2gto_hfx_gradient.o: nao2gto_kinds.o nao2gto_gamma_F.o
nao2gto_hfx_gradient.o: nao2gto_prescreen.o nao2gto_gradient.o
nao2gto_hfx_gradient.o: nao2gto_libint_wrapper_types.o
nao2gto_hfx_gradient.o: nao2gto_libint_wrapper.o nao2gto_types.o
nao2gto_hfx_gradient.o: atomlist.o nao2gto_full_DM.o

nao2gto_update_hfx_force.o: precision.o atmfuncs.o parallel.o parallelsubs.o
nao2gto_update_hfx_force.o: alloc.o atm_types.o nao2gto_hfx_gradient.o
nao2gto_update_hfx_force.o: nao2gto_full_DM.o

nao2gto_nonlin.o: precision.o nao2gto_types.o
