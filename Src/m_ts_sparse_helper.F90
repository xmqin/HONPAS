!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

! Creating a Hamiltonian/Overlap matrix at a specified k-point
! in a global UC sparsity pattern is enabled using these routines:

!  - create_HS
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
!  - create_U
!    1. accepts a globalized sparsity pattern of a matrix
!    2. Creates a matrix in Upper Triangular form
!       which directly can be inserted in LAPACK [dz]spgvd/[dz]spev
!       routines.
!    3. Symmetrization is not enforced by this routine
!  - create_Full
!    1. accepts a globalized sparsity pattern of a matrix
!    2. creates the full matrix (with all the zeroes) which can
!       be directly used for BLAS/LAPACK matrix operations
!    3. Symmetrization is not enforced by this routine

!  - AllReduce_SpData
!    1. Accepts a matrix which will be reduced in the
!       MPI_Comm_World communicator.
!    2. It takes a work-array as an additional argument
!       (however we should consider the MPI_INPLACE)

module m_ts_sparse_helper

  use precision, only : dp
  use m_ts_method, only : orb_type, TYP_BUFFER, TYP_DEVICE

  implicit none

  private

#ifdef MPI
  interface AllReduce_SpData
     module procedure AllReduce_dSpData1D
     module procedure AllReduce_zSpData1D
     module procedure AllReduce_dSpData2D
     module procedure AllReduce_zSpData2D
  end interface 
  public :: AllReduce_SpData
#endif

  interface create_HS
     module procedure create_HS_Gamma
     module procedure create_HS_kpt
  end interface 
  public :: create_HS

  interface symmetrize_HS
     module procedure symmetrize_HS_Gamma
     module procedure symmetrize_HS_kpt
  end interface 
  public :: symmetrize_HS

  interface create_U
     module procedure create_Gamma_U
     module procedure create_kpt_U
  end interface 
  public :: create_U

  interface create_Full
     module procedure create_Gamma_Full
     module procedure create_kpt_Full
  end interface 
  public :: create_Full

contains

   ! Helper routine to create and distribute the sparse 
   ! k-point Hamiltonian.
  subroutine create_HS_kpt(dit,sp, &
       Ef, &
       N_Elec, Elecs, no_u, n_s, &
       n_nzs, H, S, sc_off, SpArrH, SpArrS, k, &
       nwork, work)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use m_ts_electype

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
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
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
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
    complex(dp), pointer :: zH(:), zS(:)
    complex(dp) :: ph(0:n_s-1)
    type(Sparsity), pointer :: sp_k
    integer :: no_l, lio, io, io_T, ind, jo, jo_T, ind_k, kn
     
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_k => spar(SpArrH)
    call attach(sp_k, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! obtain the value arrays...
    zH => val(SpArrH)
    zS => val(SpArrS)

    ! Pre-calculate phases
    do jo = 0 , n_s - 1
       ph(jo) = cdexp(dcmplx(0._dp, &
            k(1) * sc_off(1,jo) + &
            k(2) * sc_off(2,jo) + &
            k(3) * sc_off(3,jo)))
    end do

    zH(:) = dcmplx(0._dp,0._dp)
    zS(:) = dcmplx(0._dp,0._dp)

! No data race condition as each processor takes a separate row
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,io_T,kn,ind,jo,jo_T,ind_k)
    do lio = 1 , no_l

       ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio)
      io_T = orb_type(io)
      kn = k_ncol(io)
      ! if there is no contribution in this row
      if ( kn /= 0 ) then
        
#ifndef TS_NOCHECKS
       ! The io orbitals are in the range [1;no_u_TS]
       ! This should be redundant as it is catched by kn==0
       if ( io_T == TYP_BUFFER ) then
         call die('Error in code, &
             &please contact Nick Papior Andersen nickpapior@gmail.com')
       end if
#endif

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

         ! as the local sparsity pattern is a super-cell pattern,
         ! we need to check the unit-cell orbital
         ! The unit-cell column index
         jo = UCORB(l_col(ind),no_u)
         jo_T = orb_type(jo)

         ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
         if ( jo_T == TYP_BUFFER ) cycle

         ! Do a check whether we have connections
         ! across the junction...
         ! This is the same as removing all electrode connections
         if ( io_T > 0 .and. jo_T > 0 .and. io_T /= jo_T ) cycle
           
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

    end do
!$OMP end parallel do
     
#ifdef MPI
    if ( dist_nodes(dit) > 1 ) then
       ! Note that zH => val(SpArrH)
       ! Note that zS => val(SpArrS)
       call AllReduce_SpData(SpArrH,nwork,work)
       call AllReduce_SpData(SpArrS,nwork,work)
    end if
#endif

    ! We symmetrize AND shift
    call symmetrize_HS_kpt(N_Elec,Elecs,Ef,SpArrH,SpArrS)
     
    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_kpt

  subroutine symmetrize_HS_kpt(N_Elec, Elecs, Ef,SpArrH, SpArrS)
    use class_Sparsity
    use class_zSpData1D

    use m_ts_electype

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in) :: Ef
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(zSpData1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    type(Sparsity), pointer :: s
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), pointer :: zH(:), zS(:)
    integer :: iEl, jEl, nr, io, ind, jo, rin, rind
    real(dp) :: E_Ef(0:N_Elec)

    ! create the overlap electrode fermi-level
    E_Ef(:) = Ef
    do iEl = 1 , N_Elec
      ! Note that for bulk V_frac_CT will be set to 0.
      E_Ef(iEl) = Ef - Elecs(iEl)%V_frac_CT * Elecs(iEl)%mu%mu
    end do

    s  => spar(SpArrH)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=nr)
    zH => val(SpArrH)
    zS => val(SpArrS)

    ! This loop is across the local rows...
! No data race condition as each processor takes a separate row
!$OMP parallel do default(shared), &
!$OMP&private(io,iEl,jEl,ind,jo,rin,rind)
    do io = 1 , nr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(io) /= 0 ) then

       iEl = orb_type(io)
       if ( iEl /= TYP_BUFFER ) then ! this means a buffer row
       jEl = 0

       ! Now we loop across the update region
       ! This one must *per definition* have less elements.
       ! Hence, we can exploit this, and find equivalent
       ! super-cell orbitals.
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          jo = l_col(ind)

          ! As we symmetrize we do not need
          ! to cycle all points through two times...
          if ( jo < io ) cycle

          ! an electrode will not connect with
          ! other electrodes, hence, we only need to check
          ! whether the electrode connects with itself
          if ( iEl > 0 ) then
             jEl = orb_type(jo)
             if ( jEl /= iEl ) then
                if ( jEl /= TYP_DEVICE ) &
                     call die('Programming error, &
                     &contact Nick Papior Andersen, nickpapior@gmail.com')
             else if ( Elecs(jEl)%bulk ) then
                ! we must refrain from shifting the fermi-level
                ! it has already been done... This should never occur though
                jEl = 0
             end if
          end if

          ! We will find the Hermitian part:
          ! The fact that we have a SYMMETRIC
          ! update region makes this *tricky* part easy...
          rin  = l_ptr(jo)
          ! TODO, this REQUIRES that l_col(:) is sorted
          rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
#ifndef TS_NOCHECKS
          ! We do a check, just to be sure...
          if ( rind <= rin ) &
               call die('ERROR symmetrization orbital does not &
               &exist.')
#endif

          ! Symmetrize (notice that we transpose here!)
          ! See prep_GF
          zS(rind) = 0.5_dp * ( zS(ind) + dconjg(zS(rind)) )
          zS(ind)  = dconjg(zS(rind))

          zH(rind) = 0.5_dp * ( zH(ind) + dconjg(zH(rind)) ) &
               - E_Ef(jEl) * zS(rind)
          zH(ind)  = dconjg(zH(rind))

          if ( ind == rind ) then
             ! This is the diagonal matrix elements
             zS(ind) = dreal(zS(ind))
             zH(ind) = dreal(zH(ind))
          end if
                      
       end do

       end if
       end if

    end do
!$OMP end parallel do

  end subroutine symmetrize_HS_kpt


  ! Helper routine to create and distribute the sparse 
  ! k-point Hamiltonian.
  subroutine create_HS_Gamma(dit,sp, &
       Ef, &
       N_Elec, Elecs, no_u, &
       n_nzs, H, S, SpArrH, SpArrS, &
       nwork, work)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    use m_ts_electype

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Fermi-level
    real(dp), intent(in) :: Ef
    ! The electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: no_u
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    ! The hamiltonian and overlap sparse matrices 
    real(dp), intent(in) :: H(n_nzs),S(n_nzs)
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(dSpData1D), intent(inout) :: SpArrH, SpArrS
    ! we pass a work array
    integer, intent(in) :: nwork
    ! work-array
    real(dp), intent(in out) :: work(nwork)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer  :: k_ncol(:), k_ptr(:), k_col(:)
    real(dp), pointer :: dH(:), dS(:)
    type(Sparsity), pointer :: sp_G
    integer :: no_l, lio, io, ind, jo, ind_k, io_T, jo_T
    
    ! obtain the local number of rows and the global...
    no_l = nrows(sp)
    if ( no_u /= nrows_g(sp) ) then
       call die('Creating the k-&point matrix in &
            &transiesta went wrong. Please TODO...')
    end if

    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! obtain the full sparsity unit-cell
    sp_G => spar(SpArrH)
    call attach(sp_G, n_col=k_ncol,list_ptr=k_ptr,list_col=k_col)

    ! obtain the value arrays...
    dH => val(SpArrH)
    dS => val(SpArrS)

    dH(:) = 0._dp
    dS(:) = 0._dp

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,io_T,ind,jo,jo_T,ind_k)
    do lio = 1 , no_l

      ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio)
      io_T = orb_type(io)
      ! if there is no contribution in this row
      if ( k_ncol(io) /= 0 ) then
        
#ifndef TS_NOCHECKS
       ! The io orbitals are in the range [1;no_u]
       if ( io_T == TYP_BUFFER ) then
         call die('Error in programming: Contact &
             &Nick Papior Andersen: nickpapior@gmail.com')
       end if
#endif

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

         ! as the local sparsity pattern is a super-cell pattern,
         ! we need to check the unit-cell orbital
         ! The unit-cell column index
         jo = UCORB(l_col(ind),no_u)
         jo_T = orb_type(jo)

         ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
         if ( jo_T == TYP_BUFFER ) cycle

         ! Do a check whether we have connections
         ! across the junction...
         ! This is the same as removing all electrode connections
         if ( io_T > 0 .and. jo_T > 0 .and. io_T /= jo_T ) cycle
         
         ! find the equivalent position in the sparsity pattern
         ! of the full unit cell
         ind_k = k_ptr(io)
         ind_k = ind_k + SFIND(k_col(ind_k+1:ind_k+k_ncol(io)),jo)
         ! if ( ind_k <= k_ptr(io) ) &
         ! call die('Could not find k-point index')
         if ( ind_k <= k_ptr(io) ) cycle
         
         dH(ind_k) = dH(ind_k) + H(ind)
         dS(ind_k) = dS(ind_k) + S(ind)
         
       end do

      end if

    end do
!$OMP end parallel do
     
#ifdef MPI
    ! Note that dH => val(SpArrH)
    ! Note that dS => val(SpArrS)
    call AllReduce_SpData(SpArrH,nwork,work)
    call AllReduce_SpData(SpArrS,nwork,work)
#endif

    ! We need to do symmetrization AFTER reduction as we need the full
    ! Hamiltonian before we can do anything
    call symmetrize_HS_Gamma(N_Elec,Elecs,Ef,SpArrH,SpArrS)

    ! It could be argued that MPI reduction provides
    ! numeric fluctuations.
    ! However, the block-cyclic distribution ensures that
    ! there are no two elements accessed by two or more processors.
    ! This makes all non-local elements ZERO, and there should not
    ! be flucuations on adding ZEROS as they are *only* dependent
    ! on order of summation.

  end subroutine create_HS_Gamma

  subroutine symmetrize_HS_Gamma(N_elec, Elecs, Ef, SpArrH, SpArrS)
    use class_Sparsity
    use class_dSpData1D

    use m_ts_electype

    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in) :: Ef
    ! The arrays we will save in... these are the entire TS-region sparsity
    type(dSpData1D), intent(inout) :: SpArrH, SpArrS

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    type(Sparsity), pointer :: s
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: dH(:), dS(:)
    real(dp) :: E_Ef(0:N_Elec)
    integer :: iEl, jEl, nr, io, ind, jo, rin, rind
   
    ! create the overlap electrode fermi-level
    E_Ef(:) = Ef
    do iEl = 1 , N_Elec
      E_Ef(iEl) = Ef - Elecs(iEl)%V_frac_CT * Elecs(iEl)%mu%mu
    end do
    
    s  => spar(SpArrH)
    call attach(s, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=nr)
    dH => val(SpArrH)
    dS => val(SpArrS)

    ! This loop is across the local rows...
!$OMP parallel do default(shared), &
!$OMP&private(io,iEl,jEl,ind,jo,rin,rind)
    do io = 1 , nr

       ! Quickly go past the empty regions... (we have nothing to update)
       if ( l_ncol(io) /= 0 ) then

       iEl = orb_type(io)
       if ( iEl /= TYP_BUFFER ) then ! this means a buffer row
       jEl = 0

       ! Now we loop across the update region
       ! This one must *per definition* have less elements.
       ! Hence, we can exploit this, and find equivalent
       ! super-cell orbitals.
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          jo = l_col(ind)

          ! As we symmetrize we do not need
          ! to cycle all points through two times...
          if ( jo < io ) cycle

          ! an electrode will not connect with
          ! other electrodes, hence, we only need to check
          ! whether the electrode connects with itself
          if ( iEl > 0 ) then
             jEl = orb_type(jo)
             if ( jEl /= iEl ) then
#ifndef TS_NOCHECKS
                if ( jEl /= TYP_DEVICE ) &
                     call die('Programming error, &
                     &contact Nick Papior Andersen, nickpapior@gmail.com')
#endif
             else if ( Elecs(jEl)%bulk ) then
                ! we must refrain from shifting the fermi-level
                ! it has already been done... This should never occur though
                jEl = 0
             end if
          end if

          ! We will find the Hermitian part:
          ! The fact that we have a SYMMETRIC
          ! update region makes this *tricky* part easy...
          rin  = l_ptr(jo)
          ! TODO, this REQUIRES that l_col(:) is sorted
          rind = rin + SFIND(l_col(rin+1:rin+l_ncol(jo)),io)
#ifndef TS_NOCHECKS
          ! We do a check, just to be sure...
          if ( rind <= rin ) &
               call die('ERROR symmetrization orbital does not &
               &exist.')
#endif

          ! Symmetrize (for Gamma, transposed is the same!)
          dS(ind) = 0.5_dp * ( dS(ind) + dS(rind) )
          dH(ind) = 0.5_dp * ( dH(ind) + dH(rind) ) &
               - E_Ef(jEl) * dS(ind)

          ! we have a real Matrix
          dH(rind) = dH(ind)
          dS(rind) = dS(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
    
  end subroutine symmetrize_HS_Gamma


  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_Gamma_U(dit,sp, &
       no, r, &
       n_nzs, A, A_UT)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: A(n_nzs)
    ! The UT format matrix
    real(dp), intent(out) :: A_UT(no*(no+1)/2)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, j, i, idx
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    A_UT(:) = 0._dp

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,j,idx)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          idx = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          j = rgn_pivot(r,idx)
          if ( j <= 0 ) cycle

          ! Calculate position
          ! For Gamma, we do not need the complex conjugate...
          if ( i > j ) then
             idx = j + (i -1)*i /2
!$OMP atomic
             A_UT(idx) = A_UT(idx) + 0.5_dp * A(ind)
          else if ( i < j ) then
             idx =  i + (j-1)*j/2
!$OMP atomic
             A_UT(idx) = A_UT(idx) + 0.5_dp * A(ind)
          else
             idx =  i + (j-1)*j/2
!$OMP atomic
             A_UT(idx) = A_UT(idx) +          A(ind)
          end if

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_Gamma_U

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_Gamma_Full(dit,sp, &
       no, r, &
       n_nzs, A, A_full)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs
    real(dp), intent(in) :: A(n_nzs)
    ! The matrix
    real(dp), intent(out) :: A_full(no,no)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, jo, i
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    A_full(:,:) = 0._dp

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,jo)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          jo = rgn_pivot(r,jo)
          if ( jo <= 0  ) cycle

          ! Calculate position
          ! For Gamma, we do not need the complex conjugate...
          A_full(jo,io) = A_full(jo,io) + A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_Gamma_Full

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_kpt_U(dit,sp, &
       no, r, &
       n_nzs, n_s, A, sc_off, A_UT, k)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs, n_s
    real(dp), intent(in) :: A(n_nzs), sc_off(3,0:n_s-1)
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
    ! The UT format matrix
    complex(dp), intent(out) :: A_UT(no*(no+1)/2)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, j, i, idx, is
    complex(dp) :: ph
    real(dp) :: w
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    A_UT(:) = dcmplx(0._dp,0._dp)

    w = log(0.5)

!$OMP parallel do default(shared), firstprivate(w), &
!$OMP&private(i,io,lio,ind,j,is,ph,idx)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          is = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          j = rgn_pivot(r,is)
          if ( j <= 0  ) cycle

          is = (l_col(ind)-1)/no_u

          ! Calculate position
          if ( i > j ) then
             ph = cdexp(dcmplx(w, - &
                  k(1) * sc_off(1,is) - &
                  k(2) * sc_off(2,is) - &
                  k(3) * sc_off(3,is)))
             idx = j + (i - 1)* i/2
          else if ( i < j ) then
             ph = cdexp(dcmplx(w, &
                  k(1) * sc_off(1,is) + &
                  k(2) * sc_off(2,is) + &
                  k(3) * sc_off(3,is)))
             idx = i  + (j-1)*j/2
          else
             ! diagonal elements are not "double" counted
             ph = cdexp(dcmplx(0._dp, &
                  k(1) * sc_off(1,is) + &
                  k(2) * sc_off(2,is) + &
                  k(3) * sc_off(3,is)))
             idx = i  + (j-1)*j/2
          end if

!$OMP atomic
          A_UT(idx) = A_UT(idx) + ph * A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_kpt_U

  ! Helper routine to create and distribute an upper
  ! tri-angular matrix of the hamiltonian in a specified
  ! region.
  subroutine create_kpt_full(dit,sp, &
       no, r, &
       n_nzs, n_s, A, sc_off, A_full, k)

    use class_OrbitalDistribution
    use class_Sparsity

    use m_region

    use geom_helper,       only : UCORB

! *********************
! * INPUT variables   *
! *********************
    ! the distribution that the H and S live under
    type(OrbitalDistribution), intent(inout) :: dit
    ! The (local) sparsity pattern that H, S lives by
    type(Sparsity), intent(inout) :: sp
    ! Number of orbitals that form this array
    integer, intent(in) :: no
    ! The region that we wish to create the UT matrix in
    type(tRgn), intent(in) :: r
    ! The number of elements in the sparse arrays
    integer, intent(in) :: n_nzs, n_s
    real(dp), intent(in) :: A(n_nzs), sc_off(3,0:n_s-1)
    ! The k-point we will create
    real(dp), intent(in) :: k(3)
    ! The UT format matrix
    complex(dp), intent(out) :: A_full(no,no)

! *********************
! * LOCAL variables   *
! *********************
    ! Create loop-variables for doing stuff
    integer, pointer  :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind, jo, i, is
    complex(dp) :: ph(0:n_s-1)
    
    ! Create all the local sparsity super-cell
    call attach(sp, n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    do is = 0 , n_s - 1
       ph(is) = cdexp(dcmplx(0._dp, &
            k(1) * sc_off(1,is) + &
            k(2) * sc_off(2,is) + &
            k(3) * sc_off(3,is)))
    end do

    A_full(:,:) = dcmplx(0._dp,0._dp)

    ! Loop over region orbitals
!$OMP parallel do default(shared), private(i,io,lio,ind,jo,is)
    do i = 1 , no
       
       ! Global orbital
       io = r%r(i)
       ! obtain the global index of the orbital.
       lio = index_global_to_local(dit,io)
       if ( lio > 0 ) then
       
       ! if there is no contribution in this row
       if ( l_ncol(lio) > 0 ) then

       ! Loop number of entries in the row... (in the index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)

          ! If the orbital is not in the region, we skip it
          jo = rgn_pivot(r,jo)
          if ( jo <= 0  ) cycle

          is = (l_col(ind)-1)/no_u

          A_full(io,jo) = A_full(io,jo) + ph(is) * A(ind)

       end do

       end if
       end if

    end do
!$OMP end parallel do
     
  end subroutine create_kpt_full


! ************************************************
! * Routines for handling the sparsity pattern.  *
! * We supply routines for initialization and    *
! * broadcasting values.                         *
! ************************************************

#ifdef MPI

  ! **** Double precision complex ****
  subroutine AllReduce_z1D(nnzs,arr,nwork,work)
    use mpi_siesta
    integer, intent(in) :: nnzs
    complex(dp), intent(inout) :: arr(nnzs)
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    integer :: MPIerror, i
    i = 0
    do while ( i + nwork <= nnzs )
      call zcopy(nwork,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nwork, &
          MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
      i = i + nwork
    end do
    if ( i < nnzs ) then
      call zcopy(nnzs-i,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nnzs-i, &
          MPI_Double_Complex, MPI_Sum, MPI_Comm_World, MPIerror)
    end if
  end subroutine AllReduce_z1D

  subroutine AllReduce_zSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_zSpData1D
    type(zSpData1D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    complex(dp), pointer :: arr(:)
    integer :: n_nzs
    n_nzs = nnzs(sp_arr)
    arr => val(sp_arr)
    call AllReduce_z1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_zSpData1D

  subroutine AllReduce_zSpData2D(sp_arr,nwork,work,dim2_count)
    use mpi_siesta
    use class_zSpData2D
    type(zSpData2D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    integer, intent(in), optional :: dim2_count
    complex(dp), pointer :: arr(:,:)
    integer :: n_nzs, d2
    arr => val(sp_arr)
    d2 = size(arr,dim=2)
    if ( present(dim2_count) ) d2 = dim2_count
    n_nzs = nnzs(sp_arr) * d2
    call AllReduce_z1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_zSpData2D



  ! **** Double precision ****
  subroutine AllReduce_d1D(nnzs,arr,nwork,work)
    use mpi_siesta
    integer, intent(in) :: nnzs
    real(dp), intent(inout) :: arr(nnzs)
    integer, intent(in) :: nwork
    real(dp), intent(inout) :: work(nwork)
    integer :: MPIerror, i
    i = 0
    do while ( i + nwork <= nnzs )
      call dcopy(nwork,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nwork, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
      i = i + nwork
    end do
    if ( i < nnzs ) then
      call dcopy(nnzs-i,arr(i+1),1,work(1),1)
      call MPI_AllReduce(work(1),arr(i+1),nnzs-i, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
    end if
  end subroutine AllReduce_d1D
  
  subroutine AllReduce_dSpData1D(sp_arr,nwork,work)
    use mpi_siesta
    use class_dSpData1D
    type(dSpData1D), intent(inout) :: sp_arr
    integer, intent(in)     :: nwork
    real(dp), intent(inout) :: work(nwork)
    real(dp), pointer :: arr(:)
    integer :: n_nzs
    n_nzs = nnzs(sp_arr)
    arr => val(sp_arr)
    call AllReduce_d1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_dSpData1D

  subroutine AllReduce_dSpData2D(sp_arr,nwork,work,dim2_count)
    use mpi_siesta
    use class_dSpData2D
    type(dSpData2D), intent(inout) :: sp_arr
    integer, intent(in) :: nwork
    real(dp), intent(inout) :: work(nwork)
    integer, intent(in), optional :: dim2_count
    real(dp), pointer :: arr(:,:)
    integer :: n_nzs, d2
    arr => val(sp_arr)
    d2 = size(arr,dim=2)
    if ( present(dim2_count) ) d2 = dim2_count
    n_nzs = nnzs(sp_arr) * d2
    call AllReduce_d1D(n_nzs,arr,nwork,work)
  end subroutine AllReduce_dSpData2D

#endif

end module m_ts_sparse_helper

#ifdef MPI
subroutine my_full_G_reduce(sp_arr,nwork,work,dim2_count)
  use precision, only : dp
  use class_dSpData2D
  use m_ts_sparse_helper, only : AllReduce_SpData
  type(dSpData2D), intent(inout) :: sp_arr
  integer, intent(in)     :: nwork
  real(dp), intent(inout) :: work(nwork)
  integer, intent(in) :: dim2_count
  call AllReduce_SpData(sp_arr,nwork,work,dim2_count)
end subroutine my_full_G_reduce
#endif

