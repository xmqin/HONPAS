! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE Kpoint_grid
!
! Contains data structures and routines to deal with the kpoint-grid
! for the self-consistent calculation
! Other uses (bands, optical, polarization) have their own structures.
!
  USE precision, only : dp

  implicit none

  public :: setup_kpoint_grid, scf_kgrid_first_time, gamma_scf, &
            nkpnt, kweight, kpoint, kscell, kdispl
  public :: eff_kgrid_cutoff

  private
  
  logical                  :: scf_kgrid_first_time = .true.
  logical                  :: gamma_scf
  integer                  :: nkpnt             ! Total number of k-points
  real(dp)                 :: eff_kgrid_cutoff  ! Effective kgrid_cutoff

  real(dp),        pointer :: kweight(:) => null()
  real(dp),        pointer :: kpoint(:,:) => null()

  integer,  dimension(3,3) :: kscell = 0
  real(dp), dimension(3)   :: kdispl = 0.0_dp
  logical                  :: user_requested_mp = .false.
  logical                  :: user_requested_cutoff = .false.

  logical, save            :: time_reversal_symmetry = .true.
  logical                  :: firm_displ = .false.

  CONTAINS

  subroutine setup_Kpoint_grid( ucell )
  USE parallel, only  : Node
  USE fdf, only       : fdf_defined, fdf_get
  USE m_find_kgrid, only : find_kgrid
  use m_spin, only: TrSym

    implicit none
    real(dp) :: ucell(3,3)
    logical  :: spiral

    if (scf_kgrid_first_time) then
       spiral = fdf_defined('SpinSpiral')
          ! Allow the user to control the use of time-reversal-symmetry
          ! By default, it is on, except for "spin-spiral" calculations
          ! and/or non-collinear calculations
       time_reversal_symmetry = fdf_get(             &
            "TimeReversalSymmetryForKpoints",   &
            (.not. spiral) .and. TrSym)
       call setup_scf_kscell(ucell, firm_displ)

       scf_kgrid_first_time = .false.

    else
       if ( user_requested_mp    ) then
          ! no need to set up the kscell again
       else
          ! This was wrong in the old code
          call setup_scf_kscell(ucell, firm_displ)
       endif
    endif
    
    call find_kgrid(ucell,kscell,kdispl,firm_displ,     &
                    time_reversal_symmetry,             &
                    nkpnt,kpoint,kweight, eff_kgrid_cutoff)

    gamma_scf =  (nkpnt == 1 .and.  &
                  dot_product(kpoint(:,1),kpoint(:,1)) < 1.0e-20_dp)

    if (Node .eq. 0) call siesta_write_k_points()

  end subroutine setup_Kpoint_grid

!--------------------------------------------------------------------
  subroutine setup_scf_kscell( cell, firm_displ )

! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! ***************** OUTPUT *********************************************
! logical firm_displ   : User-specified displacements (firm)?

!   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
!   If both are present, kgrid_Monkhorst_Pack has priority. If none is
!   present, the cutoff default is zero, producing only the gamma point.
!   Examples of fdf data specifications:
!     kgrid_cutoff  50. Bohr
!     %block kgrid_Monkhorst_Pack  # Defines kscell and kdispl
!     4  0  0   0.50               # (kscell(i,1),i=1,3), kdispl(1)
!     0  4  0   0.50               # (kscell(i,2),i=1,3), kdispl(2)
!     0  0  4   0.50               # (kscell(i,3),i=1,3), kdispl(3)
!     %endblock kgrid_Monkhorst_Pack
! **********************************************************************

!  Modules

      use precision,  only : dp
      use parallel,   only : Node
      use m_minvec,   only : minvec
      use sys,        only : die
      use fdf

      implicit          none

! Passed variables
      real(dp), intent(in)   :: cell(3,3)
      logical, intent(out)   :: firm_displ

! Internal variables
      integer           i, j, factor(3,3), expansion_factor
      real(dp)          scmin(3,3),  vmod, cutoff
      real(dp)          ctransf(3,3)
      logical           mp_input

      real(dp), parameter :: defcut = 0.0_dp
      integer, dimension(3,3), parameter :: unit_matrix =  &
                         reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline


      mp_input = fdf_block('kgrid_Monkhorst_Pack',bfdf)
      if ( mp_input ) then
         user_requested_mp = .true.
         do i= 1, 3
            if (.not. fdf_bline(bfdf,pline))            &
              call die('setup_scf_kscell: ERROR in ' // &
                       'kgrid_Monkhorst_Pack block')
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

         cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
         if (cutoff /= defcut) then
         !!  write(6,"(a,f10.5)") "Kgrid cutoff input: ", cutoff
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

    end subroutine setup_scf_kscell

    subroutine siesta_write_k_points()
      USE siesta_options, only: writek
      USE units, only: Ang
      USE siesta_cml

      implicit none

      integer  :: ik, ix, i
      external :: iokp

      if ( writek ) then
         write(6,'(/,a)') 'siesta: k-point coordinates (Bohr**-1) and weights:'
         write(6,'(a,i4,3f12.6,3x,f12.6)')                          &
              ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik), &
              ik=1,nkpnt)
      endif
      ! Always write KP file
      call iokp( nkpnt, kpoint, kweight )
      write(6,'(/a,i6)')  'siesta: k-grid: Number of k-points =', nkpnt
      write(6,'(a,f10.3,a)')  'siesta: k-grid: Cutoff (effective) =',  &
           eff_kgrid_cutoff/Ang, ' Ang'
      write(6,'(a)') 'siesta: k-grid: Supercell and displacements'
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',        &
           (kscell(i,1),i=1,3), kdispl(1)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',        &
           (kscell(i,2),i=1,3), kdispl(2)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',        &
           (kscell(i,3),i=1,3), kdispl(3)
      if (cml_p) then
         call cmlStartPropertyList(xf=mainXML, title="k-points", &
                                      dictRef="siesta:kpoints")
         call cmlAddProperty(xf=mainXML, value=nkpnt, dictref='siesta:nkpnt', &
                             units="cmlUnits:countable")
         do ik=1, nkpnt
           call cmlAddKPoint(xf=mainXML, coords=kpoint(:,ik), weight=kweight(ik))
         enddo
         call cmlAddProperty(xf=mainXML, value=eff_kgrid_cutoff/Ang,     & 
              dictref='siesta:kcutof', units='siestaUnits:angstrom')
         call cmlEndPropertyList(mainXML)
         call cmlAddProperty(xf=mainXML, value=kscell, dictref='siesta:kscell', &
                             units="siestaUnits:Ang")
         call cmlAddProperty(xf=mainXML, value=kdispl, dictref='siesta:kdispl', &
                             units="siestaUnits:Ang")
      endif
    end subroutine siesta_write_k_points

END MODULE Kpoint_grid
