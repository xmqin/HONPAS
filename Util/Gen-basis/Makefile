# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Gen-basis and ioncat
#
# The idea is to use the code in the top Src directory as much as possible.
# This is achieved by the VPATH directive below.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH:=$(shell pwd)/../../Src

default: gen-basis ioncat

ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:= $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX) -UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

SYSOBJ=$(SYS).o

#
# Note that machine-specific files are now in top Src directory.
#
OBJS_GEN-BASIS=m_io.o io.o alloc.o basis_types.o precision.o  parallel.o \
          parsing.o basis_io.o  chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o interpolation.o \
	      bessph.o radfft.o sorting.o m_filter.o \
          basis_specs.o atom.o memoryinfo.o memory.o periodic_table.o\
          pseudopotential.o pxf.o dot.o  atom_options.o dftu_specs.o arw.o  \
          timer.o  xml.o m_walltime.o read_xc_info.o gen-basis.o $(SYSOBJ) \
          m_cite.o

OBJS_IONCAT=f2kcli.o m_getopts.o basis_types.o precision.o  parallel.o \
          parsing.o m_io.o io.o alloc.o atom_options.o basis_io.o atm_transfer.o atm_types.o \
          atmparams.o radial.o  memory.o chkdim.o \
	      pseudopotential.o chemical.o m_filter.o \
          basis_specs.o atom.o memoryinfo.o  m_memory.o periodic_table.o pxf.o \
          atmfuncs.o spher_harm.o interpolation.o old_atmfuncs.o arw.o \
          xml.o m_walltime.o ioncat.o dftu_specs.o m_cite.o $(SYSOBJ)


#
# Compile XC in a private directory, without MPI
#
XC=libSiestaXC.a
XC_MAKEFILE=$(VPATH)/SiestaXC/makefile
XC_INCFLAGS:=-I $(VPATH)/SiestaXC ##$(INCFLAGS)
$(XC): 
	(mkdir -p SiestaXC ; cd SiestaXC ; $(MAKE) -f $(XC_MAKEFILE) \
                          "FC=$(FC)" "VPATH_ROOT=$(VPATH)" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(XC_INCFLAGS)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "FFLAGS=$(FFLAGS)"  module )

#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)

#
# Simple hack to link to any compiled lapack/blas code
# If COMP_LIBS is empty, we do a no-op, else link the library file here.
#
link_linalg:
	for llibs in $(COMP_LIBS); do \
   [ -f ../../$(OBJDIR)/$$llibs ] && ln -sf ../../$(OBJDIR)/$$llibs . ; \
	done

#
# SiestaXC contains already the 'sys' module and the stand-alone 'die' routine.
#
$(OBJS_GEN-BASIS): INCFLAGS:= -I./SiestaXC $(INCFLAGS)
$(OBJS_IONCAT): INCFLAGS:= -I./SiestaXC $(INCFLAGS)
#
gen-basis: $(FDF) $(XC) $(OBJS_GEN-BASIS) link_linalg
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF) $(XC) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)

ioncat:  $(FDF) $(XC) $(OBJS_IONCAT)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(OBJS_IONCAT) $(FDF) $(XC) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod
	rm -rf fdf
	rm -rf SiestaXC
	rm -f _tmp_deps deps.list  protomake*

ALL_OBJS=$(OBJS_IONCAT) $(OBJS_GEN-BASIS)

dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		*.f *.f90 *.F *.F90
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
	@sed "/Makefile/! s/sys.o//g"  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
#
#       The last sed command above will remove all mentions to 'sys.o' from this Makefile,
#       (except in lines containing 'Makefile'). The 'sys' module is now imported from
#       the SiestaXC subdirectory (see below for redefinition of INCFLAGS) and libSiestaXC.a
#       is then linked in explicitly for 'ioncat'.
#       This seems to be necessary to avoid duplication of symbols if another 'local_sys'
#       file is used. (Recall that the 'local_sys' trick was used to lighten the dependency
#       cascade from the top-level 'sys', but for Gen-basis both the local program AND SiestaXC
#       were using it.)
#       Empty dependency lines below correspond to those files which use 'sys'.
#

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: m_io.o parallel.o precision.o 
amn.o: alloc.o atmfuncs.o atomlist.o m_orderbands.o m_trialorbitalclass.o
amn.o: matel_registry.o mneighb.o new_matel.o parallel.o parallelsubs.o
amn.o: precision.o siesta2wannier90.o siesta_geom.o 
arw.o: alloc.o parallel.o precision.o 
atm_transfer.o: atm_types.o atmparams.o dftu_specs.o old_atmfuncs.o
atm_transfer.o: periodic_table.o radial.o 
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o 
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o interpolation.o
atom.o: m_filter.o old_atmfuncs.o periodic_table.o precision.o
atom.o: pseudopotential.o 
atom_graph.o: alloc.o atm_types.o atmfuncs.o class_OrbitalDistribution.o
atom_graph.o: class_SpData1D.o class_SpData2D.o class_SpData2D.o
atom_graph.o: class_Sparsity.o dftu_specs.o intrinsic_missing.o mneighb.o
atom_graph.o: precision.o radial.o sorting.o 
atom_options.o: 
atomlist.o: alloc.o atm_types.o atmfuncs.o dftu_specs.o parallel.o precision.o
atomlist.o: radial.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o 
automatic_cell.o: atm_types.o atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o atmfuncs.o atomlist.o densematrix.o diag.o diag_option.o
bands.o: files.o get_kpoints_scale.o m_spin.o parallel.o parallelsubs.o
bands.o: precision.o siesta_geom.o  writewave.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o  xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o 
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o 
bessph.o: precision.o 
bloch_unfold.o: units.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: precision.o siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o chemical.o dftu_specs.o m_mpi_utils.o parallel.o
broadcast_basis.o: radial.o
broadcast_projections.o: m_trialorbitalclass.o parallel.o siesta2wannier90.o
broyden_optim.o: m_broyddj_nocomm.o m_memory.o parallel.o precision.o
broyden_optim.o: siesta_options.o  units.o
bsc_cellxc.o: alloc.o mesh.o moremeshsubs.o parallel.o parallelsubs.o
bsc_cellxc.o: precision.o 
cart2frac.o: precision.o 
cell_broyden_optim.o: m_broyddj_nocomm.o parallel.o precision.o  units.o
cell_broyden_optim.o: zmatrix.o
cell_fire_optim.o: alloc.o m_fire.o m_memory.o parallel.o precision.o
cell_fire_optim.o: siesta_options.o  zmatrix.o
cellsubs.o: precision.o
cgvc.o: alloc.o conjgr_old.o m_mpi_utils.o parallel.o precision.o  units.o
cgvc_zmatrix.o: alloc.o conjgr.o m_mpi_utils.o parallel.o precision.o 
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: egandd.o onmod.o parallel.o precision.o 
chemical.o: parallel.o precision.o 
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o 
chkdim.o: 
chkgmx.o: cellsubs.o minvec.o precision.o 
class_Fstack_Data1D.o: class_Data1D.o class_Data1D.o
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Fstack_Pair_SpData1D.o: class_Pair_SpData1D.o class_Pair_SpData1D.o
class_Geometry.o: alloc.o
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_Pair_SpData1D.o: class_SpData1D.o class_SpData1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o
class_Sparsity.o: alloc.o
coceri.o: files.o periodic_table.o precision.o units.o
compute_dm.o: atomlist.o class_SpData1D.o compute_ebs_shift.o diagon.o
compute_dm.o: iodmhs_netcdf.o kpoint_grid.o m_dminim.o m_energies.o m_eo.o
compute_dm.o: m_hsx.o m_pexsi_driver.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_transiesta.o m_ts_global_vars.o m_zminim.o normalize_dm.o
compute_dm.o: ordern.o parallel.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o  units.o
compute_ebs_shift.o: m_mpi_utils.o parallel.o precision.o
compute_energies.o: atomlist.o class_SpData1D.o class_SpData2D.o dhscf.o
compute_energies.o: files.o m_dipol.o m_energies.o m_mpi_utils.o m_ntm.o
compute_energies.o: m_rhog.o m_spin.o m_ts_global_vars.o precision.o
compute_energies.o: siesta_geom.o siesta_options.o sparse_matrices.o
compute_energies.o: ts_energies.o
compute_max_diff.o: m_mpi_utils.o precision.o
compute_norm.o: m_mpi_utils.o m_spin.o precision.o sparse_matrices.o
compute_pw_matrix.o: alloc.o m_planewavematrix.o m_planewavematrixvar.o
compute_pw_matrix.o: parallel.o precision.o siesta2wannier90.o
compute_pw_matrix.o: sparse_matrices.o
compute_rhog.o: alloc.o atomlist.o dhscf.o files.o m_dipol.o m_energies.o
compute_rhog.o: m_ntm.o m_spin.o m_stress.o siesta_geom.o siesta_options.o
compute_rhog.o: sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o  units.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o units.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
create_Sparsity_Union.o: class_OrbitalDistribution.o class_Sparsity.o
create_Sparsity_Union.o: m_region.o parallel.o precision.o
debugmpi.o: parallel.o
delk.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o
delk.o: m_planewavematrixvar.o mesh.o meshdscf.o meshphi.o parallel.o
delk.o: parallelsubs.o precision.o siesta_geom.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o linpack.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
dfscf.o: meshphi.o parallel.o parallelsubs.o precision.o 
dftu.o: alloc.o atm_types.o atmfuncs.o compute_max_diff.o dftu_specs.o
dftu.o: m_mpi_utils.o mneighb.o new_matel.o parallel.o parallelsubs.o
dftu.o: precision.o radial.o
dftu_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
dftu_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o parallel.o
dftu_specs.o: precision.o pseudopotential.o radial.o  units.o
dhscf.o: alloc.o atmfuncs.o delk.o dfscf.o dipole.o doping_uniform.o files.o
dhscf.o: forhar.o iogrid_netcdf.o m_charge_add.o m_efield.o m_hartree_add.o
dhscf.o: m_iorho.o m_mesh_node.o m_ncdf_io.o m_ncdf_siesta.o
dhscf.o: m_partial_charges.o m_rhog.o m_spin.o m_ts_electype.o
dhscf.o: m_ts_global_vars.o m_ts_hartree.o m_ts_options.o m_ts_voltage.o mesh.o
dhscf.o: meshdscf.o meshsubs.o moremeshsubs.o parallel.o parsing.o precision.o
dhscf.o: rhofft.o rhoofd.o rhoofdsp.o siesta_options.o  units.o vmat.o
dhscf.o: vmatsp.o
diag.o: alloc.o diag_option.o parallel.o precision.o 
diag2g.o: fermid.o intrinsic_missing.o parallel.o parallelsubs.o precision.o
diag2g.o:  writewave.o
diag2k.o: fermid.o parallel.o parallelsubs.o precision.o  writewave.o
diag2kp.o: alloc.o fermid.o intrinsic_missing.o m_spin.o parallel.o
diag2kp.o: parallelsubs.o precision.o 
diag2kspiral.o: fermid.o parallel.o parallelsubs.o precision.o 
diag3g.o: fermid.o intrinsic_missing.o parallel.o parallelsubs.o precision.o
diag3g.o:  writewave.o
diag3k.o: fermid.o parallel.o parallelsubs.o precision.o  writewave.o
diag3kp.o: alloc.o fermid.o intrinsic_missing.o m_spin.o parallel.o
diag3kp.o: parallelsubs.o precision.o 
diag_option.o: parallel.o precision.o
diagg.o: alloc.o fermid.o intrinsic_missing.o m_spin.o parallel.o
diagg.o: parallelsubs.o precision.o  writewave.o
diagk.o: compute_norm.o fermid.o m_spin.o parallel.o parallelsubs.o precision.o
diagk.o:  writewave.o
diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o:  writewave.o
diagkp.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o 
diagon.o: alloc.o densematrix.o diag_option.o m_memory.o m_spin.o parallel.o
diagon.o: parallelsubs.o precision.o siesta_options.o 
diagonalizeHk.o: alloc.o atomlist.o densematrix.o m_orderbands.o m_spin.o
diagonalizeHk.o: parallel.o parallelsubs.o precision.o siesta2wannier90.o
diagonalizeHk.o: sparse_matrices.o units.o
diagpol.o: parallel.o precision.o 
die.o: parallel.o siesta_cml.o
die.o: parallel.o
dipole.o: atmfuncs.o intrinsic_missing.o m_char.o m_cite.o mesh.o parallel.o
dipole.o: precision.o siesta_cml.o  units.o
dm_charge.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
dm_charge.o: class_Sparsity.o m_mpi_utils.o m_spin.o precision.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o 
doping_uniform.o: alloc.o m_ntm.o mesh.o parallel.o precision.o 
dynamics.o: alloc.o atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o
dynamics.o: precision.o  units.o
egandd.o: alloc.o denmat.o ener3.o globalise.o gradient.o on_subs.o onmod.o
egandd.o: onmod.o precision.o 
eggbox.o: parallel.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o bessph.o interpolation.o
electrostatic.o: precision.o radfft.o radial.o 
ener3.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
extrapolateSpData2D.o: class_Data2D.o class_OrbitalDistribution.o
extrapolateSpData2D.o: class_SpData2D.o class_Sparsity.o restructSpData2D.o
extrapolon.o: parallel.o precision.o 
fdf_extra.o: m_region.o
fermid.o: errorf.o parallel.o precision.o 
fft.o: alloc.o fft1d.o m_timer.o mesh.o parallel.o parallelsubs.o precision.o
fft.o: 
fft1d.o: parallel.o precision.o 
final_H_f_stress.o: alloc.o atomlist.o class_SpData1D.o class_SpData2D.o
final_H_f_stress.o: compute_max_diff.o dftu.o dftu_specs.o dnaefs.o files.o
final_H_f_stress.o: grdsam.o kinefsm.o m_dipol.o m_energies.o m_forces.o
final_H_f_stress.o: m_hsx.o m_mpi_utils.o m_ncdf_siesta.o m_ntm.o m_spin.o
final_H_f_stress.o: m_stress.o m_ts_io.o m_ts_kpoints.o m_ts_options.o
final_H_f_stress.o: metaforce.o molecularmechanics.o naefs.o nlefsm.o overfsm.o
final_H_f_stress.o: parallel.o siesta_geom.o siesta_options.o sparse_matrices.o
final_H_f_stress.o: spinorbit.o  units.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fire_optim.o: alloc.o m_fire.o m_mpi_utils.o parallel.o precision.o
fire_optim.o: siesta_options.o units.o
flook_siesta.o: m_os.o parallel.o siesta_dicts.o siesta_geom.o siesta_options.o
flook_siesta.o: 
forhar.o: alloc.o mesh.o moremeshsubs.o parallel.o precision.o
fsiesta_mpi.o: precision.o siesta_analysis.o siesta_end.o siesta_forces.o
fsiesta_mpi.o: siesta_init.o siesta_master.o siesta_move.o 
fsiesta_sockets.o: fsockets.o
geom_helper.o: intrinsic_missing.o
get_kpoints_scale.o: parallel.o siesta_geom.o units.o
get_target_stress.o: parallel.o precision.o  units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o dhscf.o files.o m_mpi_utils.o m_partial_charges.o m_spin.o
grdsam.o: parallel.o precision.o siesta_geom.o siesta_options.o  units.o
hsparse.o: alloc.o atm_types.o atmfuncs.o atomlist.o dftu_specs.o listsc.o
hsparse.o: mneighb.o parallel.o parallelsubs.o precision.o radial.o sorting.o
hsparse.o: sparse_matrices.o
idiag.o: parallel.o 
init_output.o: files.o
initatom.o: atmparams.o atom.o atom_options.o basis_io.o basis_specs.o
initatom.o: basis_types.o chemical.o dftu_specs.o electrostatic.o m_spin.o
initatom.o: old_atmfuncs.o precision.o pseudopotential.o
initparallel.o: alloc.o class_OrbitalDistribution.o domain_decom.o parallel.o
initparallel.o: parallelsubs.o precision.o siesta_options.o sparse_matrices.o
initparallel.o: spatial.o 
io.o: m_io.o
iocg.o: files.o precision.o 
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o 
iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o 
iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o 
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o units.o
iogrid_netcdf.o: alloc.o parallel.o precision.o 
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o 
iomd.o: files.o precision.o units.o
iopipes.o: parallel.o precision.o 
iosockets.o: cellsubs.o fsockets.o m_mpi_utils.o parallel.o precision.o 
iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o 
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: 
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o  units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o new_matel.o parallel.o parallelsubs.o
kinefsm.o: precision.o
kpoint_convert.o: precision.o  units.o
kpoint_grid.o: find_kgrid.o m_spin.o minvec.o parallel.o precision.o
kpoint_grid.o: siesta_cml.o siesta_options.o  units.o
kpoint_pdos.o: find_kgrid.o m_spin.o minvec.o parallel.o precision.o
kpoint_pdos.o: siesta_options.o  units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o precision.o 
ksvinit.o: alloc.o parallel.o precision.o 
listsc.o: alloc.o
local_DOS.o: atomlist.o dhscf.o diagon.o files.o kpoint_grid.o m_energies.o
local_DOS.o: m_eo.o m_forces.o m_ntm.o m_partial_charges.o m_spin.o parallel.o
local_DOS.o: siesta_geom.o siesta_options.o sparse_matrices.o 
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o 
m_broyddj_nocomm.o: alloc.o parallel.o precision.o 
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_charge_add.o: intrinsic_missing.o m_cite.o m_geom_objects.o m_mesh_node.o
m_charge_add.o: parallel.o precision.o units.o
m_check_walltime.o: m_walltime.o
m_digest_nnkp.o: alloc.o atomlist.o files.o m_mpi_utils.o m_noccbands.o
m_digest_nnkp.o: m_spin.o m_trialorbitalclass.o parallel.o precision.o
m_digest_nnkp.o: siesta_geom.o siesta_options.o sparse_matrices.o  units.o
m_diis.o: class_Data1D.o class_Fstack_Pair_Data1D.o class_Pair_Data1D.o m_svd.o
m_diis.o: parallel.o precision.o
m_dipol.o: precision.o
m_dminim.o: densematrix.o files.o parallel.o parallelsubs.o precision.o
m_dminim.o: siesta_options.o 
m_dscfcomm.o: alloc.o parallel.o parallelsubs.o precision.o schecomm.o
m_efield.o: atmfuncs.o intrinsic_missing.o mesh.o parallel.o precision.o
m_efield.o: siesta_cml.o siesta_geom.o  units.o
m_energies.o: m_ts_global_vars.o precision.o
m_eo.o: precision.o
m_exp_coord.o: files.o m_os.o parallel.o precision.o units.o
m_filter.o: bessph.o precision.o radfft.o 
m_fire.o: parallel.o precision.o
m_fire_mixing.o: alloc.o m_fire_para.o m_mpi_utils.o parallel.o precision.o
m_fire_para.o: m_mpi_utils.o parallel.o precision.o
m_fixed.o: fdf_extra.o intrinsic_missing.o m_region.o parallel.o precision.o
m_forces.o: precision.o
m_geom_box.o: intrinsic_missing.o m_geom_aux.o
m_geom_coord.o: intrinsic_missing.o m_geom_aux.o parallel.o
m_geom_objects.o: m_geom_aux.o m_geom_box.o m_geom_coord.o m_geom_plane.o
m_geom_objects.o: m_geom_square.o
m_geom_plane.o: intrinsic_missing.o m_geom_aux.o
m_geom_square.o: intrinsic_missing.o m_geom_aux.o m_geom_plane.o
m_getopts.o: f2kcli.o
m_handle_sparse.o: class_Data2D.o class_OrbitalDistribution.o class_SpData1D.o
m_handle_sparse.o: class_SpData2D.o class_Sparsity.o geom_helper.o m_iodm.o
m_handle_sparse.o: m_os.o m_region.o m_ts_io.o parallel.o precision.o units.o
m_hartree_add.o: intrinsic_missing.o m_cite.o m_geom_objects.o m_mesh_node.o
m_hartree_add.o: parallel.o precision.o units.o
m_hs_matrix.o: alloc.o cellsubs.o geom_helper.o precision.o 
m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o 
m_integrate.o: precision.o
m_io.o: 
m_io_s.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_io_s.o: class_Sparsity.o parallel.o precision.o
m_io_yaml.o: m_energies.o precision.o version.o
m_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_iodm.o: m_io_s.o m_os.o parallel.o
m_iodm_old.o: alloc.o files.o m_matio.o m_mpi_utils.o parallel.o parallelsubs.o
m_iodm_old.o: precision.o sparse_matrices.o 
m_iorho.o: alloc.o parallel.o parallelsubs.o precision.o 
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o  units.o
m_kinetic.o: precision.o
m_mat_invert.o: intrinsic_missing.o m_pivot_array.o precision.o
m_matio.o: alloc.o parallel.o
m_memory.o: memoryinfo.o parallel.o 
m_mesh_node.o: intrinsic_missing.o parallel.o precision.o units.o
m_mixing.o: class_Data1D.o class_Fstack_Data1D.o parallel.o precision.o
m_mixing_scf.o: class_Fstack_Data1D.o m_mixing.o parallel.o precision.o
m_monitor.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
m_monitor.o: parallel.o precision.o
m_mpi_utils.o: precision.o 
m_ncdf_io.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ncdf_io.o: class_Sparsity.o m_io_s.o parallel.o precision.o
m_ncdf_siesta.o: atm_types.o atmparams.o atomlist.o class_Sparsity.o files.o
m_ncdf_siesta.o: kpoint_grid.o m_energies.o m_forces.o m_ncdf_io.o m_spin.o
m_ncdf_siesta.o: m_stress.o m_ts_electype.o m_ts_kpoints.o m_ts_options.o
m_ncdf_siesta.o: parallel.o precision.o radial.o siesta_geom.o siesta_options.o
m_ncdf_siesta.o: sparse_matrices.o timestamp.o
m_new_dm.o: alloc.o atomlist.o class_Data2D.o class_Fstack_Data1D.o
m_new_dm.o: class_Fstack_Pair_Geometry_SpData2D.o class_Geometry.o
m_new_dm.o: class_OrbitalDistribution.o class_Pair_Geometry_SpData2D.o
m_new_dm.o: class_SpData2D.o class_Sparsity.o files.o m_energies.o
m_new_dm.o: m_handle_sparse.o m_iodm.o m_mixing.o m_mixing_scf.o m_spin.o
m_new_dm.o: m_spin.o m_steps.o m_svd.o m_ts_electype.o m_ts_global_vars.o
m_new_dm.o: m_ts_iodm.o m_ts_method.o m_ts_options.o parallel.o parsing.o
m_new_dm.o: precision.o restructSpData2D.o siesta_geom.o siesta_options.o
m_new_dm.o: sparse_matrices.o  units.o
m_noccbands.o: alloc.o atmfuncs.o m_spin.o parallel.o precision.o siesta_geom.o
m_noccbands.o: 
m_options.o: precision.o
m_orderbands.o: alloc.o parallel.o parallelsubs.o siesta2wannier90.o 
m_overkkneig.o: alloc.o atomlist.o m_orderbands.o parallel.o parallelsubs.o
m_overkkneig.o: precision.o siesta_geom.o sparse_matrices.o
m_partial_charges.o: atmfuncs.o atomlist.o mesh.o meshphi.o parallel.o
m_partial_charges.o: precision.o siesta_geom.o siesta_options.o 
m_pexsi.o: precision.o
m_pexsi_dos.o: alloc.o class_Distribution.o m_mpi_utils.o m_redist_spmatrix.o
m_pexsi_dos.o: parallel.o precision.o units.o
m_pexsi_driver.o: alloc.o class_Distribution.o extrae_module.o m_convergence.o
m_pexsi_driver.o: m_mpi_utils.o m_pexsi.o m_redist_spmatrix.o parallel.o
m_pexsi_driver.o: precision.o siesta_options.o units.o
m_pexsi_local_dos.o: alloc.o atomlist.o class_Distribution.o dhscf.o files.o
m_pexsi_local_dos.o: m_energies.o m_forces.o m_mpi_utils.o m_ntm.o
m_pexsi_local_dos.o: m_redist_spmatrix.o m_spin.o parallel.o precision.o
m_pexsi_local_dos.o: siesta_geom.o siesta_options.o sparse_matrices.o 
m_pexsi_local_dos.o: units.o
m_pivot.o: class_Sparsity.o m_pivot_methods.o m_region.o
m_pivot_methods.o: class_Sparsity.o m_region.o precision.o
m_planewavematrix.o: atomlist.o dhscf.o m_planewavematrixvar.o precision.o
m_planewavematrix.o: siesta_geom.o sparse_matrices.o
m_planewavematrixvar.o: precision.o
m_pulay.o: alloc.o atomlist.o m_mpi_utils.o m_spin.o m_svd.o parallel.o
m_pulay.o: precision.o siesta_options.o sparse_matrices.o 
m_redist_spmatrix.o: alloc.o class_Distribution.o
m_region.o: class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
m_region.o: intrinsic_missing.o
m_rhog.o: alloc.o atomlist.o class_Data1D.o class_Fstack_Pair_Data1D.o
m_rhog.o: class_Pair_Data1D.o m_diis.o m_mpi_utils.o m_spin.o parallel.o
m_rhog.o: precision.o siesta_options.o sorting.o
m_rmaxh.o: precision.o
m_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData2D.o
m_sparse.o: class_Sparsity.o geom_helper.o intrinsic_missing.o parallel.o
m_sparse.o: precision.o
m_sparsity_handling.o: class_Data1D.o class_Data2D.o
m_sparsity_handling.o: class_OrbitalDistribution.o class_SpData1D.o
m_sparsity_handling.o: class_SpData2D.o class_Sparsity.o geom_helper.o
m_sparsity_handling.o: intrinsic_missing.o m_interpolate.o m_region.o
m_sparsity_handling.o: precision.o
m_spin.o: alloc.o get_kpoints_scale.o m_cite.o parallel.o precision.o 
m_stress.o: precision.o
m_supercell.o: atom_graph.o class_OrbitalDistribution.o class_SpData2D.o
m_supercell.o: parallel.o parallelsubs.o precision.o
m_svd.o: precision.o
m_target_stress.o: parallel.o precision.o  units.o
m_test_io.o: class_SpData2D.o class_Sparsity.o m_iodm.o m_matio.o parallel.o
m_test_io.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o 
m_transiesta.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_transiesta.o: class_SpData2D.o class_Sparsity.o files.o m_energies.o
m_transiesta.o: m_interpolate.o m_os.o m_ts_contour_eq.o m_ts_contour_neq.o
m_transiesta.o: m_ts_electype.o m_ts_fullg.o m_ts_fullk.o m_ts_gf.o
m_transiesta.o: m_ts_kpoints.o m_ts_method.o m_ts_mumpsg.o m_ts_mumpsk.o
m_transiesta.o: m_ts_options.o m_ts_sparse.o m_ts_tri_common.o m_ts_tri_init.o
m_transiesta.o: m_ts_trig.o m_ts_trik.o parallel.o precision.o ts_charge.o
m_transiesta.o: ts_dq.o units.o
m_trialorbitalclass.o: precision.o units.o
m_trimat_invert.o: class_TriMat.o intrinsic_missing.o m_pivot_array.o
m_trimat_invert.o: precision.o
m_ts_aux.o: precision.o
m_ts_cctype.o: m_gauss_fermi_inf.o m_ts_io_ctype.o precision.o
m_ts_chem_pot.o: m_ts_cctype.o m_ts_io_ctype.o parallel.o precision.o units.o
m_ts_contour.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_contour_eq.o
m_ts_contour.o: m_ts_contour_neq.o m_ts_electype.o m_ts_io_contour.o
m_ts_contour.o: precision.o
m_ts_contour_eq.o: m_gauss_fermi_17.o m_gauss_fermi_18.o m_gauss_fermi_19.o
m_ts_contour_eq.o: m_gauss_fermi_20.o m_gauss_fermi_22.o m_gauss_fermi_24.o
m_ts_contour_eq.o: m_gauss_fermi_26.o m_gauss_fermi_28.o m_gauss_fermi_30.o
m_ts_contour_eq.o: m_gauss_fermi_inf.o m_gauss_quad.o m_integrate.o m_io.o
m_ts_contour_eq.o: m_ts_aux.o m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_ts_contour_eq.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_ts_contour_eq.o: units.o
m_ts_contour_neq.o: m_gauss_quad.o m_integrate.o m_io.o m_ts_aux.o
m_ts_contour_neq.o: m_ts_cctype.o m_ts_chem_pot.o m_ts_electype.o
m_ts_contour_neq.o: m_ts_io_contour.o m_ts_io_ctype.o parallel.o precision.o
m_ts_contour_neq.o: units.o
m_ts_debug.o: class_Sparsity.o class_TriMat.o geom_helper.o parallel.o
m_ts_debug.o: precision.o
m_ts_dm_update.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
m_ts_dm_update.o: class_SpData2D.o class_Sparsity.o geom_helper.o
m_ts_dm_update.o: intrinsic_missing.o precision.o
m_ts_elec_se.o: intrinsic_missing.o m_ts_cctype.o m_ts_electype.o precision.o
m_ts_elec_se.o: units.o
m_ts_electrode.o: alloc.o class_SpData1D.o class_SpData2D.o class_Sparsity.o
m_ts_electrode.o: geom_helper.o intrinsic_missing.o m_iterator.o m_mat_invert.o
m_ts_electrode.o: m_pivot_array.o m_ts_elec_se.o m_ts_electype.o parallel.o
m_ts_electrode.o: precision.o  units.o
m_ts_electype.o: bloch_unfold.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_electype.o: class_SpData2D.o class_Sparsity.o create_Sparsity_SC.o
m_ts_electype.o: geom_helper.o intrinsic_missing.o m_char.o m_geom_box.o
m_ts_electype.o: m_handle_sparse.o m_iodm.o m_os.o m_region.o m_ts_chem_pot.o
m_ts_electype.o: m_ts_io.o m_ts_io_ctype.o m_ts_iodm.o parallel.o precision.o
m_ts_electype.o: units.o
m_ts_full_scat.o: intrinsic_missing.o m_ts_cctype.o m_ts_electype.o
m_ts_full_scat.o: m_ts_method.o precision.o
m_ts_fullg.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_fullg.o: class_SpData2D.o class_Sparsity.o intrinsic_missing.o
m_ts_fullg.o: m_iterator.o m_ts_cctype.o m_ts_contour_eq.o m_ts_contour_neq.o
m_ts_fullg.o: m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o m_ts_electype.o
m_ts_fullg.o: m_ts_full_scat.o m_ts_gf.o m_ts_method.o m_ts_options.o
m_ts_fullg.o: m_ts_sparse.o m_ts_sparse_helper.o m_ts_weight.o parallel.o
m_ts_fullg.o: precision.o ts_dq.o units.o
m_ts_fullk.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_fullk.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
m_ts_fullk.o: intrinsic_missing.o m_iterator.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_fullk.o: m_ts_contour_neq.o m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o
m_ts_fullk.o: m_ts_electype.o m_ts_full_scat.o m_ts_gf.o m_ts_kpoints.o
m_ts_fullk.o: m_ts_method.o m_ts_options.o m_ts_sparse.o m_ts_sparse_helper.o
m_ts_fullk.o: m_ts_weight.o parallel.o precision.o ts_dq.o units.o
m_ts_gf.o: m_os.o m_ts_cctype.o m_ts_contour_eq.o m_ts_contour_neq.o
m_ts_gf.o: m_ts_electrode.o m_ts_electype.o parallel.o precision.o 
m_ts_gf.o: units.o
m_ts_global_vars.o: densematrix.o parallel.o
m_ts_hartree.o: m_mesh_node.o m_ts_tdir.o parallel.o precision.o units.o
m_ts_io.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_io.o: class_SpData2D.o class_Sparsity.o geom_helper.o m_io_s.o m_ncdf_io.o
m_ts_io.o: m_os.o m_sparse.o parallel.o precision.o 
m_ts_io_contour.o: precision.o units.o
m_ts_io_ctype.o: m_io.o parallel.o precision.o units.o
m_ts_iodm.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_ts_iodm.o: m_io_s.o m_os.o parallel.o precision.o
m_ts_kpoints.o: files.o find_kgrid.o kpoint_grid.o m_ts_global_vars.o
m_ts_kpoints.o: m_ts_tdir.o minvec.o parallel.o precision.o siesta_cml.o
m_ts_kpoints.o: siesta_options.o 
m_ts_method.o: alloc.o fdf_extra.o geom_helper.o m_region.o m_ts_electype.o
m_ts_mumps_init.o: class_OrbitalDistribution.o class_Sparsity.o
m_ts_mumps_init.o: create_Sparsity_Union.o m_ts_elec_se.o m_ts_electype.o
m_ts_mumps_init.o: m_ts_method.o m_ts_sparse.o precision.o
m_ts_mumps_scat.o: m_ts_electype.o precision.o
m_ts_mumpsg.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_mumpsg.o: class_SpData2D.o class_Sparsity.o intrinsic_missing.o
m_ts_mumpsg.o: m_iterator.o m_ts_cctype.o m_ts_contour_eq.o m_ts_contour_neq.o
m_ts_mumpsg.o: m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o m_ts_electype.o
m_ts_mumpsg.o: m_ts_gf.o m_ts_method.o m_ts_mumps_init.o m_ts_mumps_scat.o
m_ts_mumpsg.o: m_ts_options.o m_ts_sparse.o m_ts_sparse_helper.o m_ts_weight.o
m_ts_mumpsg.o: parallel.o precision.o ts_dq.o units.o
m_ts_mumpsk.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_mumpsk.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
m_ts_mumpsk.o: intrinsic_missing.o m_iterator.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_mumpsk.o: m_ts_contour_neq.o m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o
m_ts_mumpsk.o: m_ts_electype.o m_ts_gf.o m_ts_kpoints.o m_ts_method.o
m_ts_mumpsk.o: m_ts_mumps_init.o m_ts_mumps_scat.o m_ts_options.o m_ts_sparse.o
m_ts_mumpsk.o: m_ts_sparse_helper.o m_ts_weight.o parallel.o precision.o
m_ts_mumpsk.o: ts_dq.o units.o
m_ts_options.o: atomlist.o files.o intrinsic_missing.o m_cite.o m_mixing.o
m_ts_options.o: m_mixing_scf.o m_os.o m_ts_chem_pot.o m_ts_contour.o
m_ts_options.o: m_ts_contour_eq.o m_ts_contour_neq.o m_ts_electype.o
m_ts_options.o: m_ts_global_vars.o m_ts_hartree.o m_ts_method.o
m_ts_options.o: m_ts_mumps_init.o m_ts_tdir.o m_ts_weight.o parallel.o
m_ts_options.o: precision.o siesta_options.o ts_dq.o units.o
m_ts_pivot.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_pivot.o: create_Sparsity_Union.o geom_helper.o intrinsic_missing.o
m_ts_pivot.o: m_char.o m_pivot.o m_region.o m_sparsity_handling.o m_ts_debug.o
m_ts_pivot.o: m_ts_electype.o parallel.o precision.o
m_ts_rgn2trimat.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_rgn2trimat.o: create_Sparsity_Union.o geom_helper.o m_region.o
m_ts_rgn2trimat.o: m_ts_electype.o m_ts_method.o m_ts_tri_common.o parallel.o
m_ts_rgn2trimat.o: precision.o
m_ts_sparse.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse.o: class_Sparsity.o create_Sparsity_SC.o geom_helper.o
m_ts_sparse.o: intrinsic_missing.o m_region.o m_sparsity_handling.o
m_ts_sparse.o: m_ts_debug.o m_ts_electype.o m_ts_method.o parallel.o
m_ts_sparse.o: precision.o
m_ts_sparse_helper.o: class_OrbitalDistribution.o class_SpData1D.o
m_ts_sparse_helper.o: class_SpData1D.o class_SpData2D.o class_SpData2D.o
m_ts_sparse_helper.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
m_ts_sparse_helper.o: m_region.o m_ts_electype.o m_ts_method.o precision.o
m_ts_tri_common.o: intrinsic_missing.o m_region.o m_ts_electype.o precision.o
m_ts_tri_init.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o
m_ts_tri_init.o: create_Sparsity_SC.o create_Sparsity_Union.o fdf_extra.o
m_ts_tri_init.o: m_pivot.o m_pivot_methods.o m_region.o m_sparsity_handling.o
m_ts_tri_init.o: m_ts_debug.o m_ts_electype.o m_ts_method.o m_ts_pivot.o
m_ts_tri_init.o: m_ts_rgn2trimat.o m_ts_sparse.o m_ts_tri_common.o parallel.o
m_ts_tri_init.o: precision.o
m_ts_tri_scat.o: alloc.o class_TriMat.o m_region.o m_trimat_invert.o
m_ts_tri_scat.o: m_ts_electype.o m_ts_method.o m_ts_trimat_invert.o precision.o
m_ts_trig.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_trig.o: class_SpData2D.o class_Sparsity.o class_TriMat.o
m_ts_trig.o: intrinsic_missing.o m_iterator.o m_mat_invert.o m_region.o
m_ts_trig.o: m_trimat_invert.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_trig.o: m_ts_contour_neq.o m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o
m_ts_trig.o: m_ts_electype.o m_ts_gf.o m_ts_method.o m_ts_options.o
m_ts_trig.o: m_ts_sparse.o m_ts_sparse_helper.o m_ts_tri_common.o
m_ts_trig.o: m_ts_tri_init.o m_ts_tri_scat.o m_ts_trimat_invert.o m_ts_weight.o
m_ts_trig.o: parallel.o precision.o ts_dq.o units.o
m_ts_trik.o: alloc.o class_OrbitalDistribution.o class_SpData1D.o
m_ts_trik.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o class_TriMat.o
m_ts_trik.o: intrinsic_missing.o m_iterator.o m_mat_invert.o m_region.o
m_ts_trik.o: m_trimat_invert.o m_ts_cctype.o m_ts_contour_eq.o
m_ts_trik.o: m_ts_contour_neq.o m_ts_debug.o m_ts_dm_update.o m_ts_elec_se.o
m_ts_trik.o: m_ts_electype.o m_ts_gf.o m_ts_kpoints.o m_ts_method.o
m_ts_trik.o: m_ts_options.o m_ts_sparse.o m_ts_sparse_helper.o
m_ts_trik.o: m_ts_tri_common.o m_ts_tri_init.o m_ts_tri_scat.o
m_ts_trik.o: m_ts_trimat_invert.o m_ts_weight.o parallel.o precision.o ts_dq.o
m_ts_trik.o: units.o
m_ts_trimat_invert.o: class_TriMat.o m_pivot_array.o m_region.o
m_ts_trimat_invert.o: m_trimat_invert.o m_ts_electype.o m_ts_method.o
m_ts_trimat_invert.o: precision.o
m_ts_voltage.o: iogrid_netcdf.o m_geom_box.o m_mesh_node.o m_ncdf_io.o
m_ts_voltage.o: m_ts_electype.o m_ts_options.o m_ts_tdir.o parallel.o
m_ts_voltage.o: precision.o units.o
m_ts_weight.o: class_OrbitalDistribution.o class_SpData2D.o class_Sparsity.o
m_ts_weight.o: geom_helper.o intrinsic_missing.o m_ts_chem_pot.o
m_ts_weight.o: m_ts_contour_neq.o m_ts_electype.o parallel.o precision.o
m_ts_weight.o: units.o
m_wallclock.o: m_walltime.o
m_writedelk.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o
m_writedelk.o: parallelsubs.o precision.o siesta_geom.o 
m_zminim.o: files.o parallel.o parallelsubs.o precision.o siesta_options.o
m_zminim.o: 
madelung.o: parallel.o precision.o 
matel_registry.o: m_trialorbitalclass.o precision.o radial.o spher_harm.o 
md_out.o: files.o m_energies.o precision.o  units.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
memory_all.o: rusage.o
memory_snapshot.o: rusage.o
memoryinfo.o: precision.o
mesh.o: precision.o
meshcomm.o: mesh.o precision.o
meshdscf.o: alloc.o atomlist.o m_dscfcomm.o meshphi.o parallel.o parallelsubs.o
meshdscf.o: precision.o
meshphi.o: alloc.o precision.o
meshsubs.o: alloc.o atm_types.o atmfuncs.o cellsubs.o chkgmx.o fft1d.o mesh.o
meshsubs.o: meshphi.o moremeshsubs.o parallel.o parallelsubs.o precision.o
meshsubs.o: radial.o siesta_cml.o 
metaforce.o: alloc.o parallel.o precision.o 
minvec.o: cellsubs.o precision.o sorting.o 
mixer.o: atomlist.o m_mixing.o m_mixing_scf.o m_spin.o parallel.o precision.o
mixer.o: siesta_geom.o siesta_options.o sparse_matrices.o
mmn.o: alloc.o atomlist.o m_digest_nnkp.o m_noccbands.o m_planewavematrixvar.o
mmn.o: parallel.o precision.o siesta2wannier90.o siesta_geom.o
mmn.o: sparse_matrices.o
mneighb.o: alloc.o precision.o 
molecularmechanics.o: alloc.o parallel.o precision.o  units.o
moments.o: atmfuncs.o parallel.o parallelsubs.o precision.o spinorbit.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o 
moremeshsubs.o: alloc.o mesh.o parallel.o precision.o schecomm.o 
mulliken.o: alloc.o atmfuncs.o m_spin.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o new_matel.o precision.o
new_matel.o: alloc.o errorf.o interpolation.o matel_registry.o parallel.o
new_matel.o: precision.o radfft.o siesta_options.o spher_harm.o 
nlefsm.o: alloc.o atm_types.o atmfuncs.o atomlist.o chemical.o mneighb.o
nlefsm.o: new_matel.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: siesta_options.o sparse_matrices.o 
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o 
on_subs.o: alloc.o onmod.o onmod.o
onmod.o: precision.o
optical.o: alloc.o atomlist.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o precision.o  units.o
ordern.o: alloc.o cgwf.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
ordern.o: precision.o siesta_cml.o  units.o
outcell.o: precision.o siesta_cml.o units.o
outcoor.o: alloc.o atmfuncs.o precision.o siesta_geom.o  units.o
overfsm.o: alloc.o atmfuncs.o m_spin.o mneighb.o new_matel.o parallel.o
overfsm.o: parallelsubs.o precision.o
overlap.o: alloc.o atmfuncs.o atomlist.o m_iodm_old.o m_matio.o mneighb.o
overlap.o: new_matel.o parallel.o parallelsubs.o precision.o
parallelsubs.o: domain_decom.o parallel.o spatial.o 
pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o diag_option.o files.o
pdos.o: m_energies.o m_spin.o parallel.o parallelsubs.o precision.o
pdos.o: siesta_geom.o units.o xml.o
pdos2g.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdos2k.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdos3g.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdos3k.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdosg.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdosk.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
pdoskp.o: alloc.o parallel.o parallelsubs.o precision.o  units.o
periodic_table.o: precision.o 
phirphi.o: alloc.o atmfuncs.o mneighb.o new_matel.o parallel.o parallelsubs.o
phirphi.o: precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o new_matel.o
phirphi_opt.o: parallel.o parallelsubs.o precision.o spher_harm.o 
pixmol.o: files.o periodic_table.o precision.o units.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o cellsubs.o chkgmx.o fft.o parallel.o precision.o 
post_scf_work.o: atomlist.o class_Fstack_Pair_Geometry_SpData2D.o
post_scf_work.o: class_Geometry.o class_Pair_Geometry_SpData2D.o compute_dm.o
post_scf_work.o: diagon.o final_H_f_stress.o kpoint_grid.o m_dminim.o
post_scf_work.o: m_energies.o m_eo.o m_spin.o m_steps.o m_zminim.o mneighb.o
post_scf_work.o: parallel.o siesta_geom.o siesta_options.o sparse_matrices.o
print_spin.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
print_spin.o: siesta_cml.o sparse_matrices.o
printmatrix.o: alloc.o
projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o m_spin.o
projected_DOS.o: parallel.o precision.o siesta_geom.o siesta_options.o
projected_DOS.o: sparse_matrices.o  units.o
propor.o: precision.o 
proximity_check.o: chemical.o m_ts_global_vars.o mneighb.o parallel.o
proximity_check.o: precision.o siesta_geom.o siesta_options.o units.o
pseudopotential.o: alloc.o atom_options.o interpolation.o precision.o 
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
read_options.o: files.o m_charge_add.o m_cite.o m_hartree_add.o m_mixing_scf.o
read_options.o: m_spin.o m_target_stress.o parallel.o precision.o siesta_cml.o
read_options.o: siesta_options.o  units.o
read_xc_info.o: parallel.o precision.o 
redcel.o: parallel.o precision.o 
register_rfs.o: atm_types.o atmfuncs.o matel_registry.o parallel.o precision.o
register_rfs.o: radial.o
reinit.o: files.o parallel.o siesta_master.o 
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: m_spin.o parallel.o precision.o 
reord.o: alloc.o precision.o
reordpsi.o: alloc.o m_orderbands.o parallel.o parallelsubs.o precision.o
reordpsi.o: siesta2wannier90.o
restructSpData2D.o: class_Data2D.o class_OrbitalDistribution.o class_SpData2D.o
restructSpData2D.o: class_Sparsity.o
rhofft.o: cellsubs.o fft.o precision.o 
rhooda.o: atmfuncs.o atomlist.o mesh.o meshphi.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o m_spin.o mesh.o
rhoofd.o: meshdscf.o meshphi.o parallel.o parallelsubs.o precision.o 
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o
rhoofdsp.o: meshdscf.o meshphi.o parallel.o precision.o 
save_density_matrix.o: atomlist.o class_SpData2D.o files.o iodm_netcdf.o
save_density_matrix.o: m_energies.o m_iodm.o m_matio.o m_ncdf_siesta.o m_spin.o
save_density_matrix.o: m_steps.o m_ts_global_vars.o m_ts_iodm.o m_ts_options.o
save_density_matrix.o: precision.o siesta_geom.o siesta_options.o
save_density_matrix.o: sparse_matrices.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: atomlist.o dftu.o dftu_specs.o m_convergence.o
scfconvergence_test.o: m_energies.o m_spin.o m_ts_global_vars.o m_ts_options.o
scfconvergence_test.o: m_wallclock.o parallel.o precision.o siesta_cml.o
scfconvergence_test.o: siesta_geom.o siesta_options.o sparse_matrices.o units.o
scfconvergence_test.o: write_subs.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o 
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o class_SpData1D.o class_SpData2D.o
setup_H0.o: dhscf.o dnaefs.o kinefsm.o m_energies.o m_mpi_utils.o m_ntm.o
setup_H0.o: m_spin.o metaforce.o molecularmechanics.o naefs.o nlefsm.o
setup_H0.o: siesta_geom.o siesta_options.o sparse_matrices.o spinorbit.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o class_SpData1D.o
setup_hamiltonian.o: class_SpData2D.o dftu.o dftu_specs.o dhscf.o files.o
setup_hamiltonian.o: m_dipol.o m_energies.o m_hsx.o m_mpi_utils.o m_ntm.o
setup_hamiltonian.o: m_partial_charges.o m_rhog.o m_spin.o m_steps.o m_stress.o
setup_hamiltonian.o: metaforce.o molecularmechanics.o parallel.o siesta_geom.o
setup_hamiltonian.o: siesta_options.o sparse_matrices.o 
setup_ordern_indexes.o: alloc.o domain_decom.o parallel.o spatial.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o 
siesta.o: alloc.o extrae_module.o m_io_yaml.o m_mpi_utils.o m_steps.o
siesta.o: parallel.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta2wannier90.o: files.o m_digest_nnkp.o m_spin.o m_trialorbitalclass.o
siesta2wannier90.o: parallel.o precision.o siesta_options.o 
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o dhscf.o files.o
siesta_analysis.o: flook_siesta.o kpoint_grid.o ksv.o ksvinit.o local_DOS.o
siesta_analysis.o: m_dipol.o m_energies.o m_eo.o m_forces.o m_iodm_old.o
siesta_analysis.o: m_mpi_utils.o m_ntm.o m_partial_charges.o m_pexsi_dos.o
siesta_analysis.o: m_pexsi_local_dos.o m_spin.o m_steps.o m_ts_global_vars.o
siesta_analysis.o: parallel.o projected_DOS.o siesta2wannier90.o siesta_cml.o
siesta_analysis.o: siesta_geom.o siesta_options.o sparse_matrices.o units.o
siesta_analysis.o: write_subs.o writewave.o zmatrix.o
siesta_cmlsubs.o: files.o m_uuid.o parallel.o siesta_cml.o timestamp.o
siesta_cmlsubs.o: 
siesta_dicts.o: atomlist.o class_SpData1D.o files.o kpoint_grid.o m_energies.o
siesta_dicts.o: m_forces.o m_mixing_scf.o m_steps.o m_stress.o precision.o
siesta_dicts.o: siesta_geom.o siesta_options.o sparse_matrices.o
siesta_end.o: alloc.o atomlist.o bands.o densematrix.o diag.o
siesta_end.o: extrae_eventllist.o flook_siesta.o listsc.o m_cite.o m_dscfcomm.o
siesta_end.o: m_fixed.o m_io.o m_mixing_scf.o m_rhog.o m_wallclock.o meshdscf.o
siesta_end.o: meshphi.o moremeshsubs.o new_matel.o parallel.o precision.o
siesta_end.o: siesta_cml.o siesta_cmlsubs.o siesta_master.o siesta_options.o
siesta_end.o: sparse_matrices.o timestamp.o
siesta_forces.o: atomlist.o class_Fstack_Data1D.o class_SpData2D.o compute_dm.o
siesta_forces.o: compute_energies.o compute_max_diff.o densematrix.o
siesta_forces.o: dm_charge.o files.o final_H_f_stress.o flook_siesta.o
siesta_forces.o: kpoint_grid.o m_check_walltime.o m_convergence.o m_energies.o
siesta_forces.o: m_forces.o m_iodm_old.o m_mixing.o m_mixing_scf.o
siesta_forces.o: m_mpi_utils.o m_ncdf_siesta.o m_pexsi.o m_pexsi_driver.o
siesta_forces.o: m_rhog.o m_spin.o m_steps.o m_stress.o m_transiesta.o
siesta_forces.o: m_ts_electype.o m_ts_global_vars.o m_ts_method.o
siesta_forces.o: m_ts_options.o mixer.o parallel.o post_scf_work.o precision.o
siesta_forces.o: save_density_matrix.o scfconvergence_test.o setup_H0.o
siesta_forces.o: setup_hamiltonian.o siesta_cml.o siesta_dicts.o siesta_geom.o
siesta_forces.o: siesta_master.o siesta_options.o sparse_matrices.o
siesta_forces.o: state_analysis.o state_init.o  timer.o ts_dq.o units.o
siesta_forces.o: write_subs.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o class_Fstack_Pair_Geometry_SpData2D.o
siesta_init.o: diag_option.o files.o flook_siesta.o ioxv.o kpoint_grid.o
siesta_init.o: kpoint_pdos.o ksvinit.o m_check_walltime.o m_cite.o m_energies.o
siesta_init.o: m_eo.o m_fixed.o m_forces.o m_iostruct.o m_mpi_utils.o
siesta_init.o: m_new_dm.o m_rmaxh.o m_spin.o m_steps.o m_supercell.o m_timer.o
siesta_init.o: m_wallclock.o metaforce.o molecularmechanics.o object_debug.o
siesta_init.o: parallel.o parallelsubs.o projected_DOS.o siesta_cmlsubs.o
siesta_init.o: siesta_dicts.o siesta_geom.o siesta_options.o sparse_matrices.o
siesta_init.o: struct_init.o  timer.o timestamp.o ts_init.o units.o
siesta_init.o: writewave.o zmatrix.o
siesta_master.o: iopipes.o iosockets.o precision.o 
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o flook_siesta.o ioxv.o
siesta_move.o: m_check_walltime.o m_energies.o m_exp_coord.o m_forces.o
siesta_move.o: m_kinetic.o m_mpi_utils.o m_steps.o m_stress.o m_target_stress.o
siesta_move.o: parallel.o siesta_cml.o siesta_dicts.o siesta_geom.o
siesta_move.o: siesta_master.o siesta_options.o  units.o write_subs.o
siesta_move.o: zm_broyden_optim.o zm_fire_optim.o zmatrix.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_Sparsity.o precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o 
spinorbit.o: atm_types.o atmfuncs.o atmparams.o basis_types.o m_mpi_utils.o
spinorbit.o: parallel.o parallelsubs.o precision.o pseudopotential.o 
state_analysis.o: atomlist.o born_charge.o flook_siesta.o m_energies.o
state_analysis.o: m_fixed.o m_forces.o m_ntm.o m_spin.o m_stress.o
state_analysis.o: m_ts_global_vars.o m_ts_options.o m_wallclock.o parallel.o
state_analysis.o: remove_intramol_pressure.o siesta_cml.o siesta_geom.o
state_analysis.o: siesta_options.o sparse_matrices.o ts_charge.o units.o
state_analysis.o: write_subs.o zmatrix.o
state_init.o: alloc.o atomlist.o class_Data2D.o
state_init.o: class_Fstack_Pair_Geometry_SpData2D.o
state_init.o: class_Pair_Geometry_SpData2D.o class_SpData1D.o class_SpData2D.o
state_init.o: class_Sparsity.o create_Sparsity_SC.o dftu_specs.o domain_decom.o
state_init.o: files.o hsparse.o iodm_netcdf.o iodmhs_netcdf.o ioxv.o
state_init.o: kpoint_grid.o m_energies.o m_eo.o m_handle_sparse.o m_mixing.o
state_init.o: m_mixing_scf.o m_mpi_utils.o m_ncdf_siesta.o m_new_dm.o m_os.o
state_init.o: m_pivot_methods.o m_rmaxh.o m_sparse.o m_sparsity_handling.o
state_init.o: m_spin.o m_steps.o m_supercell.o m_test_io.o m_ts_electype.o
state_init.o: m_ts_global_vars.o m_ts_io.o m_ts_kpoints.o m_ts_options.o
state_init.o: m_ts_sparse.o m_ts_tri_init.o normalize_dm.o overlap.o parallel.o
state_init.o: proximity_check.o restructSpData2D.o siesta_cml.o siesta_dicts.o
state_init.o: siesta_geom.o siesta_options.o sparse_matrices.o  ts_dq.o
state_init.o: units.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o ioxv.o m_exp_coord.o
struct_init.o: m_iostruct.o m_mpi_utils.o m_steps.o parallel.o periodic_table.o
struct_init.o: siesta_cml.o siesta_geom.o siesta_master.o siesta_options.o
struct_init.o: units.o zmatrix.o
: parallel.o siesta_cml.o
timer.o: extrae_eventllist.o extrae_module.o m_timer.o parallel.o 
timer.o: timer_tree.o
timer_tree.o: m_walltime.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o 
ts_charge.o: class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
ts_charge.o: m_ts_electype.o m_ts_method.o parallel.o precision.o
ts_dq.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
ts_dq.o: m_interpolate.o m_ts_cctype.o m_ts_chem_pot.o m_ts_contour_eq.o
ts_dq.o: m_ts_electype.o m_ts_method.o parallel.o precision.o sorting.o
ts_dq.o: ts_charge.o units.o
ts_energies.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
ts_energies.o: class_Sparsity.o geom_helper.o m_energies.o m_spin.o
ts_energies.o: m_ts_electype.o m_ts_method.o m_ts_options.o parallel.o
ts_energies.o: precision.o sparse_matrices.o
ts_init.o: m_fixed.o m_os.o m_ts_cctype.o m_ts_electrode.o m_ts_electype.o
ts_init.o: m_ts_gf.o m_ts_global_vars.o m_ts_kpoints.o m_ts_method.o
ts_init.o: m_ts_options.o parallel.o siesta_options.o
ts_show_regions.o: m_region.o m_ts_electype.o m_ts_method.o parallel.o
ts_show_regions.o: precision.o units.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vacuum_level.o: moreParallelSubs.o precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o m_spin.o mesh.o
vmat.o: meshdscf.o meshphi.o parallel.o parallelsubs.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmatsp.o: meshphi.o parallel.o precision.o
vmb.o: m_fixed.o parallel.o precision.o 
write_inp_wannier.o: alloc.o atmfuncs.o atomlist.o m_ntm.o m_orderbands.o
write_inp_wannier.o: mneighb.o parallel.o parallelsubs.o precision.o
write_inp_wannier.o: siesta2wannier90.o siesta_geom.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o units.o
write_orb_indx.o: atmfuncs.o cellsubs.o files.o precision.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: atomlist.o iofa.o m_energies.o m_forces.o m_iostruct.o m_region.o
write_subs.o: m_spin.o m_steps.o m_stress.o m_ts_electype.o m_ts_global_vars.o
write_subs.o: m_ts_method.o m_ts_options.o parallel.o precision.o siesta_cml.o
write_subs.o: siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o diag.o diag_option.o
writewave.o: get_kpoints_scale.o kpoint_grid.o m_spin.o m_spin.o parallel.o
writewave.o: parallelsubs.o precision.o siesta_geom.o  units.o
xml.o: precision.o
zm_broyden_optim.o: m_broyddj_nocomm.o parallel.o precision.o  units.o
zm_broyden_optim.o: zmatrix.o
zm_fire_optim.o: m_fire.o parallel.o precision.o siesta_options.o units.o
zm_fire_optim.o: zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o  units.o
gen-basis.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
gen-basis.o: chemical.o old_atmfuncs.o parallel.o precision.o 
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
band.o: bands.o
bloch_unfold_m.o: bloch_unfold.o
class_cdata1d.o: class_Data1D.o
class_cdata2d.o: class_Data2D.o
class_cspdata1d.o: class_SpData1D.o
class_cspdata2d.o: class_SpData2D.o
class_ctrimat.o: class_TriMat.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_dtrimat.o: class_TriMat.o
class_fstack_ddata1d.o: class_Fstack_Data1D.o
class_fstack_gdata1d.o: class_Fstack_Data1D.o
class_fstack_pair_ddata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_dspdata1d.o: class_Fstack_Pair_SpData1D.o
class_fstack_pair_gdata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_fstack_pair_gspdata1d.o: class_Fstack_Pair_SpData1D.o
class_fstack_pair_sdata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_sspdata1d.o: class_Fstack_Pair_SpData1D.o
class_fstack_sdata1d.o: class_Fstack_Data1D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_itrimat.o: class_TriMat.o
class_ldata1d.o: class_Data1D.o
class_ldata2d.o: class_Data2D.o
class_lspdata1d.o: class_SpData1D.o
class_lspdata2d.o: class_SpData2D.o
class_ltrimat.o: class_TriMat.o
class_pair_ddata1d.o: class_Pair_Data1D.o
class_pair_dspdata1d.o: class_Pair_SpData1D.o
class_pair_gdata1d.o: class_Pair_Data1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_pair_gspdata1d.o: class_Pair_SpData1D.o
class_pair_sdata1d.o: class_Pair_Data1D.o
class_pair_sspdata1d.o: class_Pair_SpData1D.o
class_sdata1d.o: class_Data1D.o
class_sdata2d.o: class_Data2D.o
class_sspdata1d.o: class_SpData1D.o
class_sspdata2d.o: class_SpData2D.o
class_strimat.o: class_TriMat.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
class_ztrimat.o: class_TriMat.o
diagmemory.o: memoryinfo.o
dipole_m.o: dipole.o
f90sockets.o: fsockets.o
fsiesta.o: fsiesta_sockets.o
iofa_m.o: iofa.o
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_born_charge.o: born_charge.o
m_broyden_optim.o: broyden_optim.o
m_cell_broyden_optim.o: cell_broyden_optim.o
m_cell_fire_optim.o: cell_fire_optim.o
m_cgwf.o: cgwf.o
m_chempot.o: chempot.o
m_chkgmx.o: chkgmx.o
m_compute_dm.o: compute_dm.o
m_compute_ebs_shift.o: compute_ebs_shift.o
m_compute_energies.o: compute_energies.o
m_compute_max_diff.o: compute_max_diff.o
m_compute_rhog.o: compute_rhog.o
m_conjgr.o: conjgr.o
m_conjgr_old.o: conjgr_old.o
m_delk.o: delk.o
m_denmat.o: denmat.o
m_dfscf.o: dfscf.o
m_dftu.o: dftu.o
m_dhscf.o: dhscf.o
m_diag.o: diag.o
m_diag_option.o: diag_option.o
m_diagon.o: diagon.o
m_dm_charge.o: dm_charge.o
m_dnaefs.o: dnaefs.o
m_doping_uniform.o: doping_uniform.o
m_dynamics.o: dynamics.o
m_eandg.o: egandd.o
m_ener3.o: ener3.o
m_errorf.o: errorf.o
m_fermid.o: fermid.o
m_fft.o: fft.o
m_final_h_f_stress.o: final_H_f_stress.o
m_find_kgrid.o: find_kgrid.o
m_fire_optim.o: fire_optim.o
m_forhar.o: forhar.o
m_get_kpoints_scale.o: get_kpoints_scale.o
m_gradient.o: gradient.o
m_grdsam.o: grdsam.o
m_hsparse.o: hsparse.o
m_intramol_pressure.o: remove_intramol_pressure.o
m_ioxv.o: ioxv.o
m_kinefsm.o: kinefsm.o
m_ksv.o: ksv.o
m_ksvinit.o: ksvinit.o
m_linpack.o: linpack.o
m_local_dos.o: local_DOS.o
m_matel_registry.o: matel_registry.o
m_minvec.o: minvec.o
m_mixer.o: mixer.o
m_naefs.o: naefs.o
m_new_matel.o: new_matel.o
m_nlefsm.o: nlefsm.o
m_norm.o: compute_norm.o
m_normalize_dm.o: normalize_dm.o
m_object_debug.o: object_debug.o
m_on_subs.o: on_subs.o
m_ordern.o: ordern.o
m_overfsm.o: overfsm.o
m_overlap.o: overlap.o
m_pexsi_solver.o: m_pexsi_driver.o
m_post_scf_work.o: post_scf_work.o
m_projected_dos.o: projected_DOS.o
m_proximity_check.o: proximity_check.o
m_radfft.o: radfft.o
m_restruct_spdata2d.o: restructSpData2D.o
m_rhofft.o: rhofft.o
m_rhoofd.o: rhoofd.o
m_rhoofdsp.o: rhoofdsp.o
m_rusage.o: rusage.o
m_save_density_matrix.o: save_density_matrix.o
m_scf_options.o: m_options.o
m_scfconvergence_test.o: scfconvergence_test.o
m_setup_h0.o: setup_H0.o
m_setup_hamiltonian.o: setup_hamiltonian.o
m_siesta2wannier90.o: siesta2wannier90.o
m_siesta_analysis.o: siesta_analysis.o
m_siesta_end.o: siesta_end.o
m_siesta_forces.o: siesta_forces.o
m_siesta_init.o: siesta_init.o
m_siesta_move.o: siesta_move.o
m_state_analysis.o: state_analysis.o
m_state_init.o: state_init.o
m_struct_init.o: struct_init.o
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
m_ts_init.o: ts_init.o
m_vmat.o: vmat.o
m_vmatsp.o: vmatsp.o
m_zm_broyden_optim.o: zm_broyden_optim.o
m_zm_fire_optim.o: zm_fire_optim.o
mtprng.o: m_uuid.o
neighbour.o: mneighb.o
on_core.o: onmod.o
on_main.o: onmod.o
printmat.o: printmatrix.o
t_spin.o: m_spin.o
timer_options.o: timer.o
trialorbitalclass.o: m_trialorbitalclass.o
ts_charge_m.o: ts_charge.o
ts_dq_m.o: ts_dq.o
ts_energies_m.o: ts_energies.o
version_info.o: version.o
write_subs_energies.o: write_subs.o
write_subs_positions.o: write_subs.o
write_subs_pressure.o: write_subs.o
