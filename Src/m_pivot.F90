!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This module can process a sparsity pattern and/or as subset
! of the sparsity pattern and create a pivot table which can
! be determined by an algorithm.

module m_pivot

  use m_region

  implicit none

  public 

  integer, parameter :: PVT_CUTHILL_MCKEE     = 1 
  integer, parameter :: PVT_REV_CUTHILL_MCKEE = 2
  integer, parameter :: PVT_GPS = 3
  integer, parameter :: PVT_REV_GPS = 4
  integer, parameter :: PVT_GGPS = 5
  integer, parameter :: PVT_REV_GGPS = 6
  integer, parameter :: PVT_PCG = 7
  integer, parameter :: PVT_REV_PCG = 8
#ifdef SIESTA__METIS
  integer, parameter :: PVT_METIS_NODEND = 9
#endif
  integer, parameter :: PVT_CONNECT = 10
  integer, parameter :: PVT_REV_CONNECT = 11
#ifdef SIESTA__METIS
  integer, parameter :: PVT_METIS_PARTGRAPHKWAY = 12
  integer, parameter :: PVT_METIS_PARTGRAPHRECURSIVE = 13
#endif

contains

  ! Main routine to create pivot table for the sparsity pattern
  ! The sparse pattern MUST not be distributed.
  subroutine sp_pvt(n,sp,pvt,method,sub,start,priority, only_sub)
    use class_Sparsity
    use m_pivot_methods
    
    ! The size of the full sparsity pattern
    integer, intent(in) :: n
    type(Sparsity), intent(inout) :: sp
    ! The pivot table returned
    type(tRgn), intent(inout) :: pvt
    integer, intent(in) :: method
    ! The subset-region where the pivot table should be considered
    type(tRgn), intent(in), optional :: sub
    ! Certain algorithms (Cuthill-Mckee can choose rather
    ! arbitrary starting points, hence we can force it to start
    ! amongst any of these)
    type(tRgn), intent(in), optional :: start
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    logical, intent(in), optional :: only_sub

    type(tRgn) :: lsub
    integer :: n_nzs
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    ! the sparse pattern is an intrinsically good way
    ! to deal with pivoting routines
    ! It intrinsically holds the degree of each graph point
    ! and swapping/queing becomes easy

    ! Check that we have a UC sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nnzs = n_nzs)
    if ( maxval(l_col) > n ) then
       call die('sp_pvt: Sparsity pattern is not an UC &
            &sparse pattern, several matrix blocks are appended.')
    end if

    if ( present(sub) ) then
       call rgn_copy(sub,lsub)
       call rgn_sort(lsub)
    else
       ! This will automatically be sorted
       call rgn_range(lsub,1,n)
    end if

    ! Call the appropriate routine
    if      ( method == PVT_CUTHILL_MCKEE     ) then
       call Cuthill_Mckee(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt, &
            start = start , priority = priority, only_sub=only_sub )
       pvt%name = 'Cuthill-Mckee'
    else if ( method == PVT_REV_CUTHILL_MCKEE ) then
       call rev_Cuthill_Mckee(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt, &
            start = start , priority = priority, only_sub=only_sub )
       pvt%name = 'rev-Cuthill-Mckee'
    else if ( method == PVT_GPS               ) then
       call GPS(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority )
       pvt%name = 'Gibbs-Poole-Stockmeyer'
    else if ( method == PVT_REV_GPS           ) then
       call rev_GPS(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority )
       pvt%name = 'rev-Gibbs-Poole-Stockmeyer'
    else if ( method == PVT_CONNECT           ) then
       call connectivity_graph(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt,start, &
            priority = priority, only_sub=only_sub )
       pvt%name = 'Connect-Graph ('//trim(start%name)//')'
    else if ( method == PVT_REV_CONNECT       ) then
       call rev_connectivity_graph(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt,start, &
            priority = priority, only_sub=only_sub )
       pvt%name = 'rev-Connect-Graph ('//trim(start%name)//')'
    else if ( method == PVT_PCG               ) then
       call PCG(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority, only_sub=only_sub )
       pvt%name = 'Peripheral-Connect-Graph'
    else if ( method == PVT_REV_PCG           ) then
       call rev_PCG(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority, only_sub=only_sub )
       pvt%name = 'rev-Peripheral-Connect-Graph'
    else if ( method == PVT_GGPS              ) then
       call GGPS(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority )
       pvt%name = 'General-Gibbs-Poole-Stockmeyer'
    else if ( method == PVT_REV_GGPS          ) then
       call rev_GGPS(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt , priority = priority )
       pvt%name = 'rev-General-Gibbs-Poole-Stockmeyer'
#ifdef SIESTA__METIS
    else if ( method == PVT_METIS_NODEND ) then
       call metis_NodeND_pvt(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt, priority = priority)
       pvt%name = 'metis-NodeND'
    else if ( method == PVT_METIS_PARTGRAPHKWAY ) then
       call metis_PartGraphKway_pvt(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt, priority = priority)
       pvt%name = 'metis-PartGraphKway'
    else if ( method == PVT_METIS_PARTGRAPHRECURSIVE ) then
       call metis_PartGraphRecursive_pvt(n,n_nzs,ncol,l_ptr,l_col,lsub,pvt, priority = priority)
       pvt%name = 'metis-PartGraphRecursive'
#endif
    else
       call die('m_pivot: Programming error, unknown method')
    end if

    call rgn_delete(lsub)
    
  end subroutine sp_pvt
    
end module m_pivot
