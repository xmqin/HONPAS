! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE siesta_options

#ifdef SIESTA__FLOOK
  use flook, only : luaState
#endif

  implicit none
  
  integer, parameter, private :: dp = selected_real_kind(10,100)

  PUBLIC
  save

  ! Compatibility options
  ! -- pre 4.0 DM and H flow logic
  logical :: compat_pre_v4_DM_H      ! General switch
  logical :: mix_after_convergence ! Mix DM or H even after convergence
  logical :: recompute_H_after_scf ! Update H while computing forces

  ! -- pre 4.0 coordinate output logic -- to be implemented
  logical :: compat_pre_v4_dynamics      ! General switch

  
  logical :: mix_scf_first ! Mix first SCF step?
  logical :: mix_scf_first_force ! Mix first SCF step? and force it!
  logical :: mix_charge    ! New: mix fourier components of rho
  logical :: mixH          ! Mix H instead of DM
  logical :: h_setup_only  ! H Setup only
  logical :: chebef        ! Compute the chemical potential in ordern?
  logical :: dumpcharge    ! Write electron density?
  logical :: fire_mix      ! SCF mixing with FIRE method
  logical :: fixspin       ! Keep the total spin fixed?
  logical :: init_anti_ferro ! Antiferro spin ordering in initdm?
  logical :: initdmaux     ! Re-initialize DM when auxiliary supercell changes?        
  logical :: allow_dm_reuse! Allow re-use of the previous geometry DM ? (with possible extrapolation)
  logical :: allow_dm_extrapolation ! Allow the extrapolation of previous geometries' DM ?
  logical :: change_kgrid_in_md ! Allow k-point grid to change in MD calculations
  logical :: use_aux_cell  ! Force the use of the auxiliary cell
  logical :: negl          ! Neglect hamiltonian matrix elements without overlap?
  logical :: noeta         ! Use computed chemical potential instead of eta in ordern?
  integer :: diag_wfs_cache! WFS cache used in diagonalization routine (0=none, 1=cdf)
  logical :: outlng        ! Long output in the output file?
  logical :: pulfile       ! Use file to store Pulay info in pulayx? (Obsolete)
  logical :: RelaxCellOnly ! Relax only lattice vectors, not atomic coordinates
  logical :: RemoveIntraMolecularPressure   ! Remove molecular virial contribution to p
  logical :: savehs        ! Write file with Hamiltonian electrostatic potential?
  logical :: savevh        ! Write file with Hartree electrostatic potential?
  logical :: savevna       ! Write file with neutral-atom potential?
  logical :: savevt        ! Write file with total effective potential?
  logical :: savedrho      ! Write file with diff. between SCF and atomic density?
  logical :: saverho       ! Write file with electron density?
  logical :: saverhoxc     ! Write file with electron density including nonlinear core correction?
  logical :: savepsch      ! Write file with ionic (local pseudopotential) charge?
  logical :: savetoch      ! Write file with total charge?
  logical :: savebader     ! Write file with charge for Bader analysis?
  logical :: usesaveddata  ! Default for UseSavedData flag
  logical :: usesavecg     ! Use continuation file for CG geometry relaxation?
  logical :: usesavelwf    ! Use continuation file for Wannier functions?
  logical :: usesavedm     ! Use cont. file for density matrix?
  logical :: usesavedmloc  ! Temporary to keep usesavedm value
  logical :: usesavexv     ! Use cont. file for atomic positions and velocities?
  logical :: usesavezm     ! Use cont. file for Z-matrix?
  logical :: writeig       ! Write eigenvalues?
  logical :: writbk        ! Write k vectors of bands?
  logical :: writmd        
  logical :: writpx        ! Write atomic coordinates at every geometry step?
  logical :: writb         ! Write band eigenvalues?
  logical :: writec        ! Write atomic coordinates at every geometry step?
  logical :: write_coop    ! Write information for COOP/COHP analysis ?
  logical :: save_ORB_INDX ! Write orbital information to ORB_INDX file ?

  ! Create graphviz information to visualize connectivity graph
  integer :: write_GRAPHVIZ
!----------------------------------------------------
! Wannier90 interface
!
  logical :: w90_processing   ! Will we call the interface with Wannier90
  logical :: w90_write_mmn    ! Write the Mmn matrix for the interface with Wannier
  logical :: w90_write_amn    ! Write the Amn matrix for the interface with Wannier
  logical :: w90_write_eig    ! Write the eigenvalues or the interface with Wannier
  logical :: w90_write_unk    ! Write the unks for the interface with Wannier
  logical :: hasnobup         ! Is the number of bands with spin up for 
                              !   wannierization defined?
  logical :: hasnobdown       ! Is the number of bands with spin down for 
                              !   wannierization defined?
  logical :: hasnob           ! Is the number of bands for wannierization defined?
                              !   (for non spin-polarized calculations).
  integer :: nobup            ! Number of bands with spin up for wannierization
  integer :: nobdown          ! Number of bands with spin down for wannierization
  integer :: nob              ! Number of bands for wannierization
                              !   (for non spin-polarized calculations).

!----------------------------------------------------
  logical :: writef        ! Write atomic forces at every geometry step?
  logical :: writek        ! Write the k vectors of the BZ integration mesh?
  logical :: writic        ! Write the initial atomic ccordinates?
  logical :: varcel        ! Change unit cell during relaxation or dynamics?
  logical :: do_pdos       ! Compute the projected density of states?
  logical :: write_tshs_history ! Write the MD track of Hamiltonian and overlap matrices in transiesta format
  logical :: write_hs_history ! Write the MD track of Hamiltonian and overlap matrices
  logical :: writedm       ! Write file with density matrix?
  logical :: write_dm_at_end_of_cycle ! Write DM at end of SCF cycle? (converged or not)
  logical :: writeH        ! Write file with Hamiltonian? (in "DM" format)
  logical :: write_H_at_end_of_cycle ! Write H at end of SCF cycle? 
  logical :: writedm_cdf   ! Write file with density matrix in netCDF form?
#ifdef NCDF_4
  logical :: write_cdf     ! Write file with all information attached
  integer :: cdf_comp_lvl  ! The compression level of the Netcdf-4 file
  logical :: cdf_w_parallel  ! Allows writing NetCDF files in parallel
  logical :: cdf_r_parallel  ! Allows reading NetCDF files in parallel, parallel read does not impose the same requirements as w_parallel
#endif
  logical :: writedm_cdf_history   ! Write file with SCF history of DM in netCDF form?
  logical :: writedmhs_cdf ! Write file with DM_in, H, DM_out, and S in netCDF form?
  logical :: writedmhs_cdf_history   ! Write file with SCF history in netCDF form?
  logical :: read_charge_cdf   ! Read charge density from file in netCDF form?
  logical :: read_deformation_charge_cdf   ! Read deformation charge density from file in netCDF form?
!
  logical :: save_initial_charge_density ! Just save the initial charge density used
  logical :: analyze_charge_density_only ! Exit dhscf after processing charge

  logical :: atmonly       ! Set up pseudoatom information only?
  logical :: harrisfun     ! Use Harris functional?
  logical :: muldeb        ! Write Mulliken populations at every SCF step?
  logical :: spndeb        ! Write spin-polarization information at every SCF step?
  logical :: orbmoms       ! Write orbital moments?

  ! Convergence options
  logical :: converge_FreeE   ! free Energy conv. to finish SCF iteration?
  real(dp):: tolerance_FreeE  ! Free-energy tolerance
  logical :: converge_Eharr   ! to finish SCF iteration?
  real(dp):: tolerance_Eharr  ! Harris tolerance
  logical :: converge_EDM     ! to finish SCF iteration?
  real(dp):: tolerance_EDM    ! Tolerance in change of EDM elements to finish SCF iteration
  logical :: converge_DM      ! to finish SCF iteration?
  real(dp):: dDtol            ! Tolerance in change of DM elements to finish SCF iteration
  logical :: converge_H       ! to finish SCF iteration?
  real(dp):: dHtol            ! Tolerance in change of H elements to finish SCF iteration
  
  logical :: broyden_optim ! Use Broyden method to optimize geometry?
  logical :: fire_optim    ! Use FIRE method to optimize geometry?
  logical :: struct_only   ! Output initial structure only?
  logical :: use_struct_file ! Read structural information from a special file?
  logical :: bornz          ! Calculate Born polarization charges?
  logical :: SCFMustConverge ! Do we have to converge for each SCF calculation?
  logical :: GeometryMustConverge ! Do we *have to* converge the relaxation?
  logical :: want_domain_decomposition ! Use domain decomposition for orbitals 
  logical :: want_spatial_decomposition ! Use spatial decomposition for orbitals
  logical :: hirshpop        ! Perform Hirshfeld population analysis?
  logical :: voropop         ! Perform Voronoi population analysis?
  logical :: partial_charges_at_every_geometry
  logical :: partial_charges_at_every_scf_step
  logical :: monitor_forces_in_scf ! Compute forces and stresses at every step

  logical :: minim_calc_eigenvalues ! Use diagonalization at the end of each MD step to find eigenvalues for OMM

  integer :: ia1           ! Atom index
  integer :: ia2           ! Atom index
  integer :: ianneal       ! Annealing option read in redata and passed to anneal
  integer :: idyn          ! Geommetry relaxation/dynamics option
  integer :: ifinal        ! Last geommetry iteration step for some types of dynamics
  integer :: ioptlwf       ! Order-N functional option read in redata used in ordern
  integer :: iquench       ! Quenching option, read in redata, used in dynamics routines
  integer :: isolve        ! Option to find density matrix: 0=>diag, 1=>order-N
  integer :: istart        ! First geommetry iteration step for certain types of dynamics
  integer :: DM_history_depth ! Number of previous density matrices used in extrapolation and reuse
  integer :: maxsav        ! Number of previous density matrices used in Pulay mixing
  integer :: broyden_maxit ! Max. iterations in Broyden geometry relaxation
  integer :: mullipop      ! Option for Mulliken population level of detail
  integer :: ncgmax        ! Max. number of conjugate gradient steps in order-N minim.
  integer :: nkick         ! Period between 'kick' steps in SCF iteration
  integer :: nmove         ! Number of geometry iterations
  integer :: nscf          ! Number of SCF iteration steps
  integer :: min_nscf      ! Minimum number of SCF iteration steps
  integer :: pmax          
  integer :: neigwanted    ! Wanted number of eigenstates (per k point)
  integer :: level         ! Option for allocation report level of detail
  integer :: call_diagon_default    ! Default number of SCF steps for which to use diagonalization before OMM
  integer :: call_diagon_first_step ! Number of SCF steps for which to use diagonalization before OMM (first MD step)

  real(dp) :: beta          ! Inverse temperature for Chebishev expansion.
  real(dp) :: bulkm         ! Bulk modulus
  real(dp) :: charnet       ! Net electric charge
  real(dp) :: rijmin        ! Min. permited interatomic distance without warning
  real(dp) :: dm_normalization_tol    ! Threshold for DM normalization mismatch error
  logical  :: normalize_dm_during_scf ! Whether we normalize the DM 
  real(dp) :: dt            ! Time step in dynamics
  real(dp) :: dx            ! Atomic displacement used to calculate Hessian matrix
  real(dp) :: dxmax         ! Max. atomic displacement allowed during geom. relaxation
  real(dp) :: eta(2)        ! Chemical-potential param. Read by redata used in ordern
  real(dp) :: etol          ! Relative tol. in CG minim, read by redata, used in ordern
  real(dp) :: ftol          ! Force tolerance to stop geometry relaxation
  real(dp) :: g2cut         ! Required planewave cutoff of real-space integration mesh
  real(dp) :: mn            ! Mass of Nose thermostat
  real(dp) :: mpr           ! Mass of Parrinello-Rahman variables
  real(dp) :: occtol        ! Occupancy threshold to build DM
  real(dp) :: rcoor         ! Cutoff radius of Localized Wave Functions in ordern
  real(dp) :: rcoorcp       ! Cutoff radius to find Fermi level by projection in ordern
  real(dp) :: rmax_bonds    ! Cutoff length for bond definition
  real(dp) :: strtol        ! Stress tolerance in relaxing the unit cell
  real(dp) :: taurelax      ! Relaxation time to reach desired T and P in anneal
  real(dp) :: temp          
  real(dp) :: tempinit      ! Initial ionic temperature read in redata
  real(dp) :: threshold     ! Min. size of arrays printed by alloc_report
  real(dp) :: tp            ! Target pressure. Read in redata. Used in dynamics routines
  real(dp) :: total_spin    ! Total spin used in spin-polarized calculations
  real(dp) :: tt            ! Target temperature. Read in redata. Used in dynamics rout.
  real(dp) :: wmix          ! Mixing weight for DM in SCF iteration
  real(dp) :: wmixkick      ! Mixing weight for DM in special 'kick' SCF steps

  ! Matrix element compatibility variable
  integer :: matel_NRTAB = 1024
  
  character(len=164) :: sname   ! System name, used to initialise read

  integer,  parameter :: SOLVE_DIAGON = 0
  integer,  parameter :: SOLVE_ORDERN = 1
  integer,  parameter :: SOLVE_TRANSI = 2
  integer,  parameter :: SOLVE_MINIM  = 3
  integer,  parameter :: SOLVE_PEXSI  = 4
  integer,  parameter :: MATRIX_WRITE = 5
  
#ifdef SIESTA__FLOOK
  ! LUA-handle
  type(luaState) :: LUA
#endif

#ifdef HAVE_LIBINT
  ! Hartree-Fock / hybrid XC calculations
  logical :: hfx_wanted = .false.
#endif

END MODULE siesta_options
