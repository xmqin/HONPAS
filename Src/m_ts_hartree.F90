!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_ts_hartree

! Module for fixing the Hartree potential so that the potential fluctuations
! does not go wild.
! This is necessary to get a stable SCF solution
!
! Created and copyrighted by: Nick Papior Andersen, 2014
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
  
  use precision, only : dp

  implicit none
  
  private
  save

  ! The idea is to have sub routines in this module to do
  ! various Hartree potential fixes
  public :: read_ts_hartree_options
  public :: ts_hartree_fix

  !< Define which planes to use for the average
  logical, public :: TS_HA_PLANES(2, 3) = .false.

  ! The fraction of the actual fix
  real(dp), public :: TS_HA_frac = 1._dp

  ! The Hartree offset potential as determined from the electrode calculation
  real(dp), public :: TS_HA_offset = 0._dp

contains

  subroutine read_ts_hartree_options()

    use fdf, only: fdf_get, leqi
    use parallel, only: IONode
    use m_ts_tdir, only : ts_tidx

    character(len=64) :: def, c
    integer :: i, pn, N

    ! Predefine the default location (should not be overwritten
    ! unless you are a developer).
    select case ( ts_tidx )
    case ( 1 )
      def = 'A'
    case ( 2 )
      def = 'B'
    case ( 3 )
      def = 'C'
    case default
      def = 'USER-REQUIRED'
    end select
    
    c = fdf_get('TS.Hartree.Fix', def)

    if ( ts_tidx > 0 .and. .not. leqi(c, def) ) then
      if ( IONode ) then
        write(*,'(a)') 'WARNING ******************************'
        write(*,'(a)') 'ts: Requesting a different hartree plane than the default!'
        write(*,'(a)') 'ts: This is now YOUR responsibility!'
        write(*,'(a)') 'WARNING ******************************'
      end if
    end if
    
    if ( leqi(c, 'USER-REQUIRED') ) then
      call die('TS.Hartree.Fix: invalid option for fixing the Hartree potential &
          &you have to specify at least one plane "[-+][ABC]"!')
    end if

    ! Now parse it
    N = len_trim(c)
    i = 1
    do while ( i <= N )
      if ( c(i:i) == '-' ) then
        pn = 1
      else if ( c(i:i) == '+' ) then
        pn = 2
      else
        pn = 0
      end if

      ! Skip to read plane
      if ( pn > 0 ) then
        i = i + 1
        if ( i > N ) call die('TS.Hartree.Fix, error in option, could not parse it!')
      else
        pn = 1
      end if
      
      select case ( c(i:i) )
      case ( 'A' )
        TS_HA_PLANES(pn, 1) = .true.
      case ( 'B' )
        TS_HA_PLANES(pn, 2) = .true.
      case ( 'C' )
        TS_HA_PLANES(pn, 3) = .true.
      case default
        call die('TS.Hartree.Fix: Could not parse plane, must be one of [ABC]!')
      end select
      i = i + 1
    end do
      
    TS_HA_frac = fdf_get('TS.Hartree.Fix.Frac',1._dp)
    TS_HA_offset = fdf_get('TS.Hartree.Offset',0._dp, 'Ry')

  end subroutine read_ts_hartree_options

  ! Fix the potential
  subroutine ts_hartree_fix(nmesh, nmeshl, Vscf)

    use parallel, only: IONode
    use precision, only : grid_p
    use units, only: eV
#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum
    use mpi_siesta, only : MPI_Comm_World, MPI_integer
    use mpi_siesta, only : MPI_double_precision
#endif
    use m_mesh_node, only: mesh_correct_idx
    use m_mesh_node, only : offset_i, offset_r, dMesh, dL

    integer, intent(in) :: nmesh(3), nmeshl(3)
    real(grid_p), intent(inout) :: Vscf(nmeshl(1),nmeshl(2),nmeshl(3))

    ! Internal variables
    integer :: nlp
    integer :: ip, i
#ifdef MPI
    integer :: MPIerror
#endif
    real(dp) :: Vav, Vtot

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE TS_VH_fix' )
#endif

    ! Quick return if no planes are fixed
    if ( .not. any(TS_HA_PLANES) ) return

    ! Initialize summation
    Vtot = 0._dp

    ! Initialize counters
    nlp = 0

    do ip = 1, 3
      
      ! Lower plane
      if ( TS_HA_PLANES(1,ip) ) then
        i = 1 - offset_i(ip)
        if ( 0 < i .and. i <= nmeshl(ip) ) then
          select case ( ip )
          case ( 1 )
            Vtot = Vtot + sum(Vscf(i, :, :))
            nlp = nlp + nmeshl(2) * nmeshl(3)
          case ( 2 )
            Vtot = Vtot + sum(Vscf(:, i, :))
            nlp = nlp + nmeshl(1) * nmeshl(3)
          case ( 3 )
            Vtot = Vtot + sum(Vscf(:, :, i))
            nlp = nlp + nmeshl(1) * nmeshl(2)
          end select
        end if
      end if

      ! Upper plane
      if ( TS_HA_PLANES(2,ip) ) then
        i = nmesh(ip) - offset_i(ip)
        if ( 0 < i .and. i <= nmeshl(ip) ) then
          select case ( ip )
          case ( 1 )
            Vtot = Vtot + sum(Vscf(i, :, :))
            nlp = nlp + nmeshl(2) * nmeshl(3)
          case ( 2 )
            Vtot = Vtot + sum(Vscf(:, i, :))
            nlp = nlp + nmeshl(1) * nmeshl(3)
          case ( 3 )
            Vtot = Vtot + sum(Vscf(:, :, i))
            nlp = nlp + nmeshl(1) * nmeshl(2)
          end select
        end if
      end if

    end do
    
    ! Scale the contribution
    Vtot = Vtot * TS_HA_frac

#ifdef MPI
    call MPI_AllReduce(Vtot,Vav,1,MPI_double_precision,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    Vtot = Vav
    call MPI_AllReduce(nlp,i,1,MPI_integer,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    nlp = i
#endif

    if ( nlp == 0 ) then
       call die('The partitioning of the basal plane went wrong. &
            &No points are encapsulated.')
    end if

    Vav = Vtot / nlp - TS_HA_offset

    if ( IONode ) then
       write(*,'(a,e12.5,a)')'ts-Vha: ',Vav / eV,' eV'
    end if
    
    ! Align potential
    Vscf(:,:,:) = Vscf(:,:,:) - Vav

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS TS_VH_fix' )
#endif

  end subroutine ts_hartree_fix

end module m_ts_hartree
