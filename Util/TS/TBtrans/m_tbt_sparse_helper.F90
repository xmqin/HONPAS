! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

! Creating a Hamiltonian/Overlap matrix at a specified k-point
! in a global UC sparsity pattern is enabled using these routines:

!  - create_region_HS
!    1. accepts a distributed matrix
!    2. requires the output matrix to be globalized in the sense
!       that the sparsity pattern is a UC sparsity pattern
!       NO dublicate entries. 
!       REQUIREMENT: each row MUST be sorted in column index
!    3. Creates the transposed matrix in the output matrix
!       A simple argument is that any traversal in the sparsity
!       pattern will be faster if it is following fortran order.
!       This reduces cache-misses
!    4. The matrix is forced Hermitian
!    5. The k-points are based on the m_tbt_k 'kRegion' data-type.

module m_tbt_sparse_helper

  use precision, only : dp
  use m_ts_method, only : orb_type, TYP_BUFFER, TYP_DEVICE
#ifdef MPI
  use m_ts_sparse_helper, only: AllReduce_SpData
#endif
  use m_ts_sparse_helper, only: symmetrize_HS

  implicit none

  private

  interface create_region_HS
     module procedure create_HS_kpt
  end interface 
  public :: create_region_HS

contains

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
  subroutine create_HS_kpt(dit,sp, &
       Ef, cell, n_k, r_k, na_u, lasto, &
       N_Elec, Elecs, no_u, n_s, &
       n_nzs, H, S, sc_off, SpArrH, SpArrS, &
       nwork, work)

    use parallel, only : Node
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use m_region

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use m_ts_electype

    use m_tbt_kregions, only : kRegion, kregion_k

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef, cell(3,3)
    ! Number of different k-regions
    integer, intent(in) :: n_k
    ! Different k-regions (0 is Gamma)
    type(kRegion), intent(in) :: r_k(0:n_k)
    ! Number of atoms, and the last orbital per atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: no_u, n_s
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(n_nzs),S(n_nzs)
    ! The supercell offsets
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    ! The arrays we will save in...
    type(zSpData1D), intent(inout) :: SpArrH, SpArrS
    ! we pass a work array
    integer, intent(in) :: nwork
    ! work-array
    complex(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: k_ncol(:), k_ptr(:), k_col(:)
    real(dp) :: bk(3), k(3), rcell(3,3)
    complex(dp), pointer :: zH(:), zS(:)
    complex(dp) :: ph(0:n_s-1)
    type(tRgn) :: ro
    type(Sparsity), pointer :: sp_k
    integer :: no_l, lio, io, ind, jo, ind_k, kn, i, il
     
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &tbtrans went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_k => spar(SpArrH)
    call attach(sp_k, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! obtain the value arrays...
    zH => val(SpArrH)
    zS => val(SpArrS)

    call reclat(cell,rcell,1)

    zH(:) = dcmplx(0._dp,0._dp)
    zS(:) = dcmplx(0._dp,0._dp)

!$OMP parallel default(shared), &
!$OMP&private(il,i,io,lio,kn,ind,jo,ind_k,bk,k,ro)

! No data race condition as each processor takes a separate row
    do il = 0 , n_k

       ! Update k-point
       if ( il == 0 ) then
          ! Gamma-region
          bk(:) = 0._dp
          k(:)  = 0._dp
       else
          call kregion_k(il,bk)
          call kpoint_convert(rcell,bk,k,-2)
       end if
       
!$OMP do 
       do i = 0 , n_s - 1
          ph(i) = cdexp(dcmplx(0._dp, &
               k(1) * sc_off(1,i) + &
               k(2) * sc_off(2,i) + &
               k(3) * sc_off(3,i)))
       end do
!$OMP end do nowait

!$OMP single
       ! Convert to orbital space
       call rgn_Atom2Orb(r_k(il)%atm,na_u,lasto,ro)
!$OMP end single

!$OMP do
    do i = 1 , ro%n
       io = ro%r(i)
       ! obtain the global index of the orbital.
       if ( orb_type(io) /= TYP_BUFFER ) then

       lio = index_global_to_local(dit,io,Node)
       if ( lio > 0 ) then

       kn = k_ncol(io)
       ! if there is no contribution in this row
       if ( kn /= 0 ) then

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          if ( orb_type(jo) == TYP_BUFFER ) cycle

          ! Do a check whether we have connections
          ! across the junction...
          ! This is the same as removing all electrode connections
          if ( count(OrbInElec(Elecs,io) .neqv. OrbInElec(Elecs,jo)) > 1 ) cycle
           
          ! find the equivalent position in the sparsity pattern
          ! of the full unit cell
          ind_k = k_ptr(io)

          ! Notice that SFIND REQUIRES that the sparsity pattern
          ! is SORTED!
          ! Thus it will only work for UC sparsity patterns.
          ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+kn),jo)
          ! if ( ind_k <= k_ptr(io) ) &
          ! call die('Could not find k-point index')
          if ( ind_k <= k_ptr(io) ) cycle

          jo = (l_col(ind)-1) / no_u
          
          zH(ind_k) = zH(ind_k) + H(ind) * ph(jo)
          zS(ind_k) = zS(ind_k) + S(ind) * ph(jo)

       end do

       end if

       end if

       end if

    end do
!$OMP end do nowait

!$OMP single
    call rgn_delete(ro)
!$OMP end single

    end do
!$OMP end parallel
     
#ifdef MPI
    if ( dist_nodes(dit) > 1 ) then
       ! Note that zH => val(SpArrH)
       ! Note that zS => val(SpArrS)
       call AllReduce_SpData(SpArrH,nwork,work)
       call AllReduce_SpData(SpArrS,nwork,work)
    end if
#endif

    ! We symmetrize AND shift
    call symmetrize_HS(N_Elec,Elecs,Ef,SpArrH,SpArrS)
     
    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_kpt

end module m_tbt_sparse_helper

