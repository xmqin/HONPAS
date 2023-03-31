! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE Kpoint_pdos
!
! Contains data structures and routines to deal with the kpoint-grid
! for the projected density of states calculation. 
! Modelled on the equivalent module for the SCF calculation.
!
  USE precision, only : dp

  implicit none

  private
  
  logical, public, save   :: pdos_kgrid_first_time = .true.
  logical, public, save   :: gamma_pdos
  integer, public, save   :: maxk_pdos         ! 
  integer, public, save   :: nkpnt_pdos        ! Total number of k-points
  real(dp)                :: eff_kgrid_cutoff  ! Effective kgrid_cutoff

  real(dp), pointer, public, save :: kweight_pdos(:) 
  real(dp), pointer, public, save :: kpoints_pdos(:,:)

  integer,  dimension(3,3), save  :: kscell = 0
  real(dp), dimension(3), save    :: kdispl = 0.0_dp

  logical, save     :: user_requested_mp = .false.
  logical, save     :: user_requested_cutoff = .false.

  logical, save     :: time_reversal_symmetry = .true.
  logical, save     :: firm_displ = .false.

  public :: setup_kpoint_pdos

  CONTAINS

  subroutine setup_Kpoint_pdos( ucell, different_pdos_grid )
  USE parallel, only  : Node
  USE fdf, only       : fdf_defined, fdf_get
  USE m_find_kgrid, only : find_kgrid
  use m_spin, only: TrSym

  implicit none
  real(dp), intent(in)  :: ucell(3,3)
  logical,  intent(out) :: different_pdos_grid

  logical :: spiral
  if (pdos_kgrid_first_time) then
    nullify(kweight_pdos,kpoints_pdos)
    spiral = fdf_defined('SpinSpiral')
      ! Allow the user to control the use of time-reversal-symmetry
      ! By default, it is on, except for "spin-spiral" calculations
      ! and/or non-collinear calculations
    time_reversal_symmetry = fdf_get(             &
           "TimeReversalSymmetryForKpoints",   &
           (.not. spiral) .and. TrSym)

    call setup_pdos_kscell(ucell, firm_displ)

    pdos_kgrid_first_time = .false.

  else
    if ( user_requested_mp    ) then
      ! no need to set up the kscell again
    else
      ! This was wrong in the old code
      call setup_pdos_kscell(ucell, firm_displ)
    endif
  endif

  if (user_requested_mp.or.user_requested_cutoff) then
    different_pdos_grid = .true.
  else
    different_pdos_grid = .false.
  endif

! If the grid hasn't been explicit specified then just set dummy values 
  if (different_pdos_grid) then
    call find_kgrid(ucell,kscell,kdispl,firm_displ, &
                    time_reversal_symmetry,         &
                    nkpnt_pdos,kpoints_pdos,kweight_pdos,eff_kgrid_cutoff)

    maxk_pdos = nkpnt_pdos
    gamma_pdos = (nkpnt_pdos == 1 .and. dot_product(kpoints_pdos(:,1),kpoints_pdos(:,1)) < 1.0e-20_dp)

    if (Node .eq. 0) call siesta_write_k_points_pdos()
  else
    nkpnt_pdos = 0
    maxk_pdos = nkpnt_pdos
    gamma_pdos = .true.
  endif

  end subroutine setup_Kpoint_pdos

!--------------------------------------------------------------------
  subroutine setup_pdos_kscell( cell, firm_displ )

! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! ***************** OUTPUT *********************************************
! logical firm_displ   : User-specified displacements (firm)?

!   The relevant fdf labels are PDOS.kgrid_cutoff and PDOS.kgrid_Monkhorst_Pack.
!   If both are present, PDOS.kgrid_Monkhorst_Pack has priority. If neither is
!   present, the cutoff default is zero, producing only the gamma point.
!   Examples of fdf data specifications:
!     PDOS.kgrid_cutoff  50. Bohr
!     %block PDOS.kgrid_Monkhorst_Pack  # Defines kscell and kdispl
!     4  0  0   0.50                    # (kscell(i,1),i=1,3), kdispl(1)
!     0  4  0   0.50                    # (kscell(i,2),i=1,3), kdispl(2)
!     0  0  4   0.50                    # (kscell(i,3),i=1,3), kdispl(3)
!     %endblock PDOS.kgrid_Monkhorst_Pack
! **********************************************************************

!  Modules

    use precision,  only : dp
    use m_minvec,   only : minvec
    use sys,        only : die
    use fdf

    implicit          none

! Passed variables
    real(dp), intent(in)   :: cell(3,3)
    logical, intent(out)   :: firm_displ

! Internal variables
    integer           i, j,  factor(3,3), expansion_factor

    real(dp)          scmin(3,3),  vmod, cutoff
    real(dp)          ctransf(3,3)
    logical           mp_input

    real(dp), parameter :: defcut = 0.0_dp
    integer, dimension(3,3), parameter :: unit_matrix =  &
                         reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      mp_input = fdf_block('PDOS.kgrid_Monkhorst_Pack',bfdf)
      if ( mp_input ) then
         user_requested_mp = .true.
         do i= 1, 3
            if (.not. fdf_bline(bfdf,pline))            &
              call die('setup_pdos_kscell: ERROR in ' // &
                       'PDOS.kgrid_Monkhorst_Pack block')
            kscell(1,i) = fdf_bintegers(pline,1)
            kscell(2,i) = fdf_bintegers(pline,2)
            kscell(3,i) = fdf_bintegers(pline,3)
            if ( fdf_bnvalues(pline) > 3 ) then
               kdispl(i) = mod(fdf_bvalues(pline,4), 1._dp)
            else
               kdispl(i) = 0._dp
            end if
          enddo
          call fdf_bclose(bfdf)
         firm_displ = .true.

      else

         cutoff = fdf_physical('PDOS.kgrid_cutoff',defcut,'Bohr')
         if (cutoff /= defcut) then
         !!  write(6,"(a,f10.5)") "PDOS Kgrid cutoff input: ", cutoff
            user_requested_cutoff = .true.
         endif

         kdispl(1:3) = 0.0_dp  ! Might be changed later
         firm_displ = .false.  ! In future we might add new options
                               ! for user-specified displacements
         
         ! Find equivalent rounded unit-cell
         call minvec( cell, scmin, ctransf )

         expansion_factor = 1
         do j = 1,3
            factor(j,1:3) = 0
            vmod = sqrt(dot_product(scmin(1:3,j),scmin(1:3,j)))
            factor(j,j) = int(2.0_dp*cutoff/vmod) + 1
            expansion_factor = expansion_factor * factor(j,j)
         enddo
         ! Generate actual supercell skeleton
         kscell = matmul(ctransf, factor)
         ! Avoid confusing permutations
         if (expansion_factor == 1) then
            kscell = unit_matrix
         endif
      endif

  end subroutine setup_pdos_kscell

  subroutine siesta_write_k_points_pdos()

    USE siesta_options, only: writek
    USE units, only: Ang

    implicit none

    integer  :: ik, ix, i

    if ( writek ) then
      write(6,'(/,a)') 'siesta: k-point coordinates (Bohr**-1) and weights for PDOS:'
      write(6,'(a,i4,3f12.6,3x,f12.6)')                          &
              ('siesta: ', ik, (kpoints_pdos(ix,ik),ix=1,3), kweight_pdos(ik), ik=1,nkpnt_pdos)
    endif

    write(6,'(/a,i6)')  'siesta: PDOS.k-grid: Number of k-points =', nkpnt_pdos
    write(6,'(a,f10.3,a)')  'siesta: PDOS.k-grid: Cutoff (effective) =', eff_kgrid_cutoff/Ang, ' Ang'
    write(6,'(a)') 'siesta: PDOS.k-grid: Supercell and displacements'
    write(6,'(a,3i4,3x,f8.3)') 'siesta: PDOS.k-grid: ', (kscell(i,1),i=1,3), kdispl(1)
    write(6,'(a,3i4,3x,f8.3)') 'siesta: PDOS.k-grid: ', (kscell(i,2),i=1,3), kdispl(2)
    write(6,'(a,3i4,3x,f8.3)') 'siesta: PDOS.k-grid: ', (kscell(i,3),i=1,3), kdispl(3)

  end subroutine siesta_write_k_points_pdos

END MODULE Kpoint_pdos
