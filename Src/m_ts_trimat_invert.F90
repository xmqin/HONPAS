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

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.
! It is customized to deal with calculating the inverted matrix
! only in a certain region.

module m_ts_trimat_invert

  ! We use the general inversion module to obtain the
  ! same routines etc.

  use precision, only : dp
  use class_zTrimat
  use m_trimat_invert, only : calc_Mnn_inv
  use m_trimat_invert, only : calc_Xn_div_Cn_p1, calc_Yn_div_Bn_m1
  use m_trimat_invert, only : Xn_div_Cn_p1, Yn_div_Bn_m1

  use m_pivot_array, only : Npiv, ipiv, init_pivot
  
  implicit none

  private

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

  public :: invert_BiasTriMat_prep
#ifdef TBTRANS
  public :: BiasTriMat_prep
#endif
  public :: invert_BiasTriMat_col
  public :: invert_BiasTriMat_rgn
  public :: TriMat_Bias_idxs

contains

  subroutine invert_BiasTriMat_prep(M,Minv, part_cols, all_nn)

    type(zTriMat), intent(inout) :: M, Minv
    ! This is a 3x... list of
    !  1 == part
    !  2 == smallest column in the part
    !  3 == largest column in the part
    integer, intent(in), optional :: part_cols(:,:)
    logical, intent(in), optional :: all_nn

    complex(dp), pointer :: Mpinv(:)

    integer :: n, np, sNm1, sNp1
    integer :: sCol, eCol, i
    logical :: piv_initialized

    np = parts(M)
#ifndef TS_NOCHECKS
    if ( np /= parts(Minv) ) then
       call die('Could not calculate the inverse on non equal sized &
            &matrices')
    end if
    if ( np == 1 ) then
       call die('This matrix is not tri-diagonal')
    end if
#endif
    piv_initialized = .true.
#ifndef TS_NOCHECKS
    do n = 1 , np
       if ( Npiv < nrows_g(M,n) ) piv_initialized = .false.
    end do
    if ( .not. piv_initialized ) then
       call die('Pivoting array for inverting matrix not set.')
    end if
#endif

    call timer('V_TM_Pinv',1)

    ! Calculate all Xn/Cn+1
    do n = np - 1 , 1 , -1 
       Mpinv => val(Minv,n+1,n+1)
       sNp1  =  nrows_g(M,n+1)
       call calc_Xn_div_Cn_p1(M,Minv, n, Mpinv, sNp1**2 )
    end do
    ! Calculate all Yn/Bn-1
    do n = 2 , np
       Mpinv => val(Minv,n-1,n-1)
       sNm1  =  nrows_g(M,n-1)
       call calc_Yn_div_Bn_m1(M,Minv, n, Mpinv, sNm1**2 )
    end do

    ! At this point the status is:
    !   - M contains the original matrix
    !   - Minv contains all Xn/Cn+1 and Yn/Bn-1 in all parts

    piv_initialized = .false.
    if ( present(all_nn) ) piv_initialized = all_nn
    if ( .not. present(part_cols) ) piv_initialized = .true.

    if ( piv_initialized ) then
       
       ! We need to calculate all Mnn
       do n = 1 , np
          call calc_Mnn_inv(M,Minv,n)
       end do

    else if ( present(part_cols) ) then

       ! This loops an array of parts to invert
       do i = 1 , size(part_cols,2)

          ! Note that these *must* be separate parts
          n = part_cols(1,i)
          sCol = part_cols(2,i)
          eCol = part_cols(3,i)

          call calc_Mnn_inv_cols(M,Minv,n,sCol,eCol)
          
       end do

    else
       call die('inv-BTM-prep: Wrong options.')
    end if

    ! At this point we have calculated the 
    !  Mnn matrices for the overlapping regions for the
    !  electrodes.
    ! This is now saved in Minv.
    ! Minv contains all information to calculate
    ! all Gf columns

    call timer('V_TM_Pinv',2)

  end subroutine invert_BiasTriMat_prep

#ifdef TBTRANS
  subroutine BiasTrimat_prep(M,N_Elec,Elecs,has_El,part_cols)

    use m_ts_electype

    type(zTriMat), intent(inout) :: M
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    logical, intent(in) :: has_El(N_Elec)
    integer, intent(inout), allocatable :: part_cols(:,:)

    integer :: iEl, off, sCol, eCol
    integer :: n, np, sN, Ne, sNm1, sNp1
    ! We need to allow the electrodes to be split
    ! across two blocks.
    integer :: tmp(3,N_Elec*2)

    np = parts(M)

    Ne = 0
    off = 0
    do n = 1 , np
       ! Get part information
       sN   = nrows_g(M,n)
       
       sCol = huge(1)
       eCol = -1
       
       do iEl = 1 , N_Elec
          
          if ( .not. has_El(iEl) ) cycle
          
          ! We know that any diagonal region with the
          ! electrodes will only occupy at most 2 parts
          ! Hence, it makes no sense to loop them individually
          ! Also inDpvt is a sorted array with device indices
          sNm1 = Elecs(iEl)%inDpvt%r(1)
          sNp1 = Elecs(iEl)%inDpvt%r(Elecs(iEl)%inDpvt%n)
          if ( which_part(M,sNm1) <= n .and. &
               n <= which_part(M,sNp1) ) then
             
             ! get number of columns that belongs to
             ! the electrode in the 'n' diagonal part
             ! this means we only calculate "what is needed"
             sCol = min(sCol, max(sNm1 - off ,  1) )
             eCol = max(eCol, min(sNp1 - off , sN) )
             
          end if
          
       end do
       
       if ( eCol /= -1 ) then
          Ne = Ne + 1
          tmp(1,Ne) = n
          tmp(2,Ne) = sCol
          tmp(3,Ne) = eCol
       end if
       
       off = off + sN
       
    end do

    allocate(part_cols(3,Ne))
    part_cols(:,:) = tmp(:,1:Ne)
    
  end subroutine BiasTrimat_prep
#endif

  subroutine invert_BiasTriMat_col(M,Minv,El,calc_parts)

    use m_ts_electype
    use m_ts_method, only : orb_offset

    type(zTriMat), intent(inout) :: M, Minv
    type(Elec), intent(in) :: El
    logical, intent(in) :: calc_parts(:)

    complex(dp), pointer :: Mpinv(:), Mp(:)
    complex(dp), pointer :: Xn(:), Yn(:)
    complex(dp), pointer :: z(:)

    integer :: nr, np, no
    integer :: idx_o
    integer :: sPart, ePart, lsPart, lePart
    integer :: sColF, eColF, sIdxF, eIdxF
    integer :: sColT, eColT, sIdxT, eIdxT
    integer :: sN, sNc, sNm1, sNp1, n, s
    integer :: off
    logical :: piv_initialized

    ! In this routine M should have been processed through invert_PrepTriMat
    ! So M contains all *needed* inv(Mnn) and all Xn/Cn+1 and Yn/Bn-1.
    ! So we will save the result in Minv
#ifndef TS_NOCHECKS
    if ( parts(M) /= parts(Minv) ) then
       call die('Could not calculate the inverse on non equal sized &
            &matrices')
    end if
    if ( parts(M) == 1 ) then
       call die('This matrix is not tri-diagonal')
    end if
    piv_initialized = .true.
    do n = 1 , parts(M) 
       if ( Npiv < nrows_g(M,n) ) piv_initialized = .false.
    end do
    if ( .not. piv_initialized ) then
       call die('Pivoting array for inverting matrix not set.')
    end if

    if ( parts(M) /= size(calc_parts) ) then
       call die('Error in code, calc_parts, not consistent')
    end if
#endif

    call timer('V_TM_inv',1)

    nr = nrows_g(M)
    np = parts(M)

    idx_o = El%idx_o - orb_offset(El%idx_o)
    no = TotUsedOrbs(El)

    sPart = which_part(M,idx_o)
    ePart = which_part(M,idx_o+no-1)
#ifndef TS_NOCHECKS
    if ( sPart < 1 ) call die('Error in the Bias inversion, sPart')
    if ( ePart - sPart + 1 > 2 ) call die('Error in trimat partition')
    if ( ePart > parts(M) ) call die('Error in the Bias inversion, ePart')
#endif

    ! Point to the matrices
    z => val(Minv,all=.true.)

    ! First we need to copy over the Mnn with the electrode part!

    ! with this offset we can calculate the column offset for
    ! the current part
    off = 0
    do n = 1 , sPart - 1
       off = off + nrows_g(M,n)
    end do
    idx_o = idx_o - off
    if ( idx_o <= 0 ) call die('Error in electrode setup')
    do n = sPart , ePart

       ! current count of orbitals in the tri-diagonal segment
       if ( n > 1  ) sNm1 = nrows_g(M,n-1)
                     sN   = nrows_g(M,n  )
       if ( n < np ) sNp1 = nrows_g(M,n+1)
       
       ! placement of the already inverted matrix
       Mp => val(M,n,n)

       ! get number of columns that belongs to
       ! the electrode in the 'n' diagonal part
       sColF = max(idx_o          ,  1)
       eColF = min(idx_o + no - 1 , sN)
#ifndef TS_NOCHECKS
       if ( eColF < sColF ) &
            call die('Here: Something went wrong')
#endif
       sIdxF = (sColF-1) * sN + 1
       eIdxF =  eColF    * sN

       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Minv,no,n,sIdxT,eIdxT)
       Mpinv => z(sIdxT:eIdxT)

       ! get the placement in the inversed column
       if ( 1 <= idx_o ) then
          ! we are taking the first part of the inversed matrix
          sColT = 1
          eColT = min(eColF-sColF+1,no)
       else
          sColT = -idx_o + 2 ! we have to pass zero
          eColT = no
       end if
       sIdxT = (sColT-1) * sN + 1
       eIdxT =  eColT    * sN

#ifndef TS_NOCHECKS
       if ( eIdxT - sIdxT /= eIdxF - sIdxF ) & 
            call die('Error in determining column')
#endif
       call zcopy(eIdxT-sIdxT+1,Mp(sIdxF),1,Mpinv(sIdxT),1)

       if ( sPart == ePart ) cycle ! we have everything! :)

       ! We need to calculate the remaining 
       ! inverted matrix (they currently reside 
       ! with the neighbouring cells)

       ! first calculate the missing number of columns
       sNc = no - (eColF - sColF + 1)

       if ( n == sPart ) then
          ! we miss the right part of Mnn
          ! we thus need the
          ! Mmn = -Ym+1/Bm * Mm+1n

          Yn => Yn_div_Bn_m1(M,n+1)

          ! placement of the inverted matrix
          Mp => val(M,n+1,n+1)

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sN,sNc,sNp1, &
               zm1, Yn, sN, Mp(1),sNp1,z0, Mpinv(eIdxT+1),sN)

       else if ( n == ePart ) then
          ! we miss the left part of Mnn
          ! we thus need the
          ! Mmn = -Xm-1/Cm * Mm-1n

          Xn => Xn_div_Cn_p1(M,n-1)

          ! placement of the inverted matrix
          Mp => val(M,n-1,n-1)
          s = sNm1 ** 2

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sN,sNc,sNm1, &
               zm1, Xn, sN, Mp(s-sNm1*sNc+1),sNm1,z0, Mpinv(1),sN)

       end if

       ! update offset on rows
       off = off + sN
       idx_o = idx_o - sN

    end do

    ! We now have inv(Mnn) in the correct place.
    ! figure out what we should calculate
    ! We can not save a lot of computations here
    ! as we can not be sure of the full extend of
    ! the "end"-blocks. Only if two consecutive 
    ! blocks are fully encompassed in one electrode
    ! can we reduce the computational work.
    do n = 1 , np
       if ( calc_parts(n) ) then
          lsPart = max(1,n-1)
          exit
       end if
    end do
    do n = np , 1 , -1
       if ( calc_parts(n) ) then
          lePart = min(n+1,np)
          exit
       end if
    end do

    ! We now calculate:
    !  Mmn = -Ym+1/Bm * Mm+1n, for m<n
    do n = sPart - 1 , lsPart , - 1
       
       sN   = nrows_g(M,n)
       sNp1 = nrows_g(M,n+1)

       ! get Ym+1/Bm
       Yn => Yn_div_Bn_m1(M,n+1)

       ! Get Mm+1n
       call TriMat_Bias_idxs(Minv,no,n+1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,no,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)

#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,no,sNp1, &
            zm1, Yn, sN, Mp(1),sNp1,z0, Mpinv(1),sN)
       
    end do

    ! We now calculate:
    !  Mmn = -Xm-1/Cm * Mm-1n, for m>n
    do n = ePart + 1 , lePart

       sNm1 = nrows_g(M,n-1)
       sN   = nrows_g(M,n)

       ! get Xm-1/Cm
       Xn => Xn_div_Cn_p1(M,n-1)
       
       ! Get Mm-1n
       call TriMat_Bias_idxs(Minv,no,n-1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,no,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)
       
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,no,sNm1, &
            zm1, Xn, sN, Mp(1),sNm1,z0, Mpinv(1),sN)
       
    end do

    ! At this point the total 
    ! inverted column is placed at the end of
    ! the tri-mat inversion.

    call timer('V_TM_inv',2)

  end subroutine invert_BiasTriMat_col

  subroutine invert_BiasTriMat_rgn(M,Minv,r,pvt,r_col,only_diag)

    use m_region

    type(zTriMat), intent(inout) :: M, Minv
    ! The pivoting table for the tri-diagonal matrix
    ! pvt is the pivoting back to original index (to remove
    ! rgn_pivot calls, i.e. pvt%r(r%r) == (/1,...,no_u/)
    type(tRgn), intent(in) :: r, pvt
    ! The siesta indices for the column we wish to calculate
    type(tRgn), intent(in) :: r_col
    ! Whether we should only calculate the diagonal Green function
    logical, intent(in), optional :: only_diag

    complex(dp), pointer :: Mpinv(:), Mp(:)
    complex(dp), pointer :: XY(:)
    complex(dp), pointer :: z(:)

    ! Used for BLAS calls (local variables)
    complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
    complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

    integer :: nr, np
    integer :: sPart, ePart
    integer :: i, i_Elec, idx_Elec
    integer :: sIdxF, eIdxF, sIdxT, eIdxT
    integer :: sN, n, in, s, sNo, nb
    integer, pointer :: crows(:)

    ! In this routine M should have been processed through invert_PrepTriMat
    ! So M contains all *needed* inv(Mnn) and all Xn/Cn+1 and Yn/Bn-1.
    ! So we will save the result in Minv
#ifndef TS_NOCHECKS
    if ( parts(M) /= parts(Minv) ) then
       call die('Could not calculate the inverse on non equal sized &
            &matrices')
    end if
    if ( parts(M) == 1 ) then
       call die('This matrix is not tri-diagonal')
    end if
#endif
    call timer('V_TM_inv',1)

    nr = nrows_g(M)
    np = parts(M)
    crows => cum_rows(M)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    sPart = huge(1)
    ePart = 0
    do n = 1 , r_col%n
       s = pvt%r(r_col%r(n))
       sPart = min(sPart,which_part(M,s))
       ePart = max(ePart,which_part(M,s))
    end do
#ifndef TS_NOCHECKS
    if ( sPart < 1 ) call die('Error in the Bias inversion, sPart')
    if ( ePart - sPart + 1 > 2 ) call die('Error in trimat partition')
    if ( ePart > np ) call die('Error in the Bias inversion, ePart')
#endif

    ! Point to the matrices
    z => val(Minv,all=.true.)

    ! CHECK
    ! This requires that the o_inD is sorted
    ! according to the device region.
    ! Check m_tbt_regions to assert this!

    i_Elec = 1
    do while ( i_Elec <= r_col%n ) 

       idx_Elec = pvt%r(r_col%r(i_Elec))

       ! We start by copying over the Mnn in blocks

       ! We start by creating a region of consecutive
       ! memory.
       n = which_part(M,idx_Elec)
       sN = nrows_g(M,n)
       nb = 1
       do while ( i_Elec + nb <= r_col%n )
          i = pvt%r(r_col%r(i_Elec+nb))
          ! In case it is not consecutive
          if ( i - idx_Elec /= nb ) exit
          ! In case the block changes, then
          ! we cut the block size here.
          if ( n /= which_part(M,i) ) exit
          nb = nb + 1
       end do
       
       ! Copy over this portion of the Mnn array
       
       ! Figure out which part we have Mnn in
       i = pvt%r(r_col%r(i_Elec))
       n = which_part(M,i)

       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxT,eIdxT)
       Mpinv => z(sIdxT:eIdxT)

       ! Correct the copied elements
       ! Figure out the placement in the copied to array
       ! First we calculate the starting index of the block
       sIdxT = ( i_Elec -    1  ) * sN + 1
       eIdxT = ( i_Elec + nb - 1) * sN

       ! *** Now we have the matrix that we can save the 
       !     block Mnn in

       ! We now need to figure out the placement of the 
       ! Mnn part that we have already calculated.
       Mp => val(M,n,n)
       i = pvt%r(r_col%r(i_Elec))
       sIdxF = (i-(crows(n)-sN)-1) * sN + 1
       i = pvt%r(r_col%r(i_Elec+nb-1))
       eIdxF = (i-(crows(n)-sN)) * sN

       ! Check that we have something correct...
!print *,trim(El%name),sIdxT,eIdxT,sIdxF,eIdxF
!print *,trim(El%name),eIdxT-sIdxT,eIdxF-sIdxF

#ifndef TS_NOCHECKS
       if ( eIdxT - sIdxT /= eIdxF - sIdxF ) & 
            call die('Error in determining column')
#endif

       ! Prepare for next segment
       Mp => Mp(sIdxF:)

       ! Copy over diagonal block
       call zcopy(eIdxT-sIdxT+1,Mp(1),1,Mpinv(sIdxT),1)

       ! Calculate the off-diagonal Green function in the regions
       ! of interest
       do in = n - 1 , sPart , - 1

!print *,'calculating Mn-1n',in,n, i_Elec
          ! Number of orbitals in the other segment
          sNo = nrows_g(M,in)
          ! grab the Yn matrix to perform the 
          ! Gf off-diagonal calculation.
          XY => Yn_div_Bn_m1(M,in+1)

          ! placement of the inverted matrix
          call TriMat_Bias_idxs(Minv,r_col%n,in,sIdxT,eIdxT)
          Mpinv => z(sIdxT:eIdxT)
          sIdxT = ( i_Elec - 1 ) * sNo + 1

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,nb,sN, &
               zm1, XY, sNo, Mp(1),sN,z0, Mpinv(sIdxT),sNo)

          Mp => Mpinv(sIdxT:)
          sN = sNo
          
       end do
          
       ! Reset arrays to just before off-diagonal
       sN = nrows_g(M,n)
       Mp => val(M,n,n)
       Mp => Mp(sIdxF:)

       ! Calculate the off-diagonal Green function in the regions
       ! of interest
       do in = n + 1 , ePart

!print *,'calculating Mn+1n',in,n, i_Elec
          ! Number of orbitals in the other segment
          sNo = nrows_g(M,in)
          ! grab the Xn matrix to perform the 
          ! Gf off-diagonal calculation.
          XY => Xn_div_Cn_p1(M,in-1)

          ! placement of the inverted matrix
          call TriMat_Bias_idxs(Minv,r_col%n,in,sIdxT,eIdxT)
          Mpinv => z(sIdxT:eIdxT)
          sIdxT = ( i_Elec - 1 ) * sNo + 1

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,nb,sN, &
               zm1, XY, sNo, Mp(1),sN,z0, Mpinv(sIdxT),sNo)

          Mp => Mpinv(sIdxT:)
          sN = sNo
          
       end do
       
       ! Update current segment of the electrode copied entries.
       i_Elec = i_Elec + nb

    end do

    if ( present(only_diag) ) then
       if ( only_diag ) then
          call timer('V_TM_inv',2)
          return
       end if
    end if
    
    ! Now we need to calculate the remaining column

    ! We now calculate:
    !  Mmn = -Ym+1/Bm * Mm+1n, for m<n
    do n = sPart - 1 , 1 , - 1
       
       sN  = nrows_g(M,n)
       sNo = nrows_g(M,n+1)

       ! get Ym+1/Bm
       XY => Yn_div_Bn_m1(M,n+1)

       ! Get Mm+1n
       call TriMat_Bias_idxs(Minv,r_col%n,n+1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)

#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,r_col%n,sNo, &
            zm1, XY, sN, Mp(1),sNo,z0, Mpinv(1),sN)
       
    end do

    ! We now calculate:
    !  Mmn = -Xm-1/Cm * Mm-1n, for m>n
    do n = ePart + 1 , np

       sNo = nrows_g(M,n-1)
       sN  = nrows_g(M,n)

       ! get Xm-1/Cm
       XY => Xn_div_Cn_p1(M,n-1)
       
       ! Get Mm-1n
       call TriMat_Bias_idxs(Minv,r_col%n,n-1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)
       
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,r_col%n,sNo, &
            zm1, XY, sN, Mp(1),sNo,z0, Mpinv(1),sN)
       
    end do

    ! At this point the total 
    ! inverted column is placed at the end of
    ! the tri-mat inversion.

    call timer('V_TM_inv',2)
    
  end subroutine invert_BiasTriMat_rgn


  ! We will partition the system by:
  ! 1. nrows_g(tri,1) x no
  ! 2. nrows_g(tri,2) x no
  ! ...
  subroutine TriMat_Bias_idxs(M,no,p,sIdx,eIdx)
    type(zTriMat), intent(in) :: M
    ! no is the number of orbitals we wish to take out
    ! p is the part that we wish to point to
    integer, intent(in) :: no, p
    integer, intent(out) :: sIdx, eIdx
    integer, pointer :: crows(:)

    integer :: cum

    crows => cum_rows(M)
    
    cum = nrows_g(M,p)
    eIdx = no * cum - 1
    cum = nrows_g(M) - crows(p) + cum
    
    ! This is the number of elements already occupied
    sIdx = elements(M, all=.true.) - no * cum + 1
    eIdx = sIdx + eIdx

  end subroutine TriMat_Bias_idxs

  
  subroutine calc_Mnn_inv_cols(M,Minv,n,sCol,eCol)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, sCol, eCol
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Xn(:), Yn(:), Cn(:), Bn(:)
    integer :: sNm1, sN, sNp1, ierr, i

    if ( 1 < n )        sNm1 = nrows_g(M,n-1)
                        sN   = nrows_g(M,n)
    if ( n < parts(M) ) sNp1 = nrows_g(M,n+1)

    if ( sCol == 1 .and. eCol == sN ) then
       call calc_Mnn_inv(M,Minv,n)
       return
    end if
    
    ! Retrieve Ann
    Mp => val(M,n,n)
    if ( n == 1 ) then
       ! First we calculate M11^-1
       ! Retrieve the X1/C2 array
       Xn => Xn_div_Cn_p1(Minv,n)
       ! The C2 array
       Cn => val(M,n,n+1)
       ! Calculate: A1 - X1
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)

    else if ( n == parts(M) ) then

       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Yn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)

    else
       ! Retrieve the Xn/Cn+1 array
       Xn => Xn_div_Cn_p1(Minv,n)
       ! The Cn+1 array
       Cn => val(M,n,n+1)
       ! Calculate: An - Xn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)
       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Xn - Yn
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)
  
    end if

    ! Retrieve the position in the inverted matrix
    Mpinv => val(Minv,n,n)
    Mpinv((sCol-1)*sN+1:eCol*sN) = dcmplx(0._dp,0._dp)
    do i = sCol - 1 , eCol - 1
       Mpinv(i * sN + i + 1) = dcmplx(1._dp,0._dp)
    end do

    i = eCol - sCol + 1
    call zgesv(sN,i,Mp,sN,ipiv,Mpinv((sCol-1)*sN+1),sN,ierr)
    if ( ierr /= 0 ) then
       call die('Error in inverting the partial bias block')
    end if

  end subroutine calc_Mnn_inv_cols

end module m_ts_trimat_invert
    
