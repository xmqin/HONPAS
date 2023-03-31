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

module m_ts_tri_scat

  use precision, only : dp
  use m_ts_method

  implicit none

  private

  public :: GF_Gamma_GF
  public :: dir_GF_Gamma_GF
  public :: has_full_part
  public :: insert_Self_Energies

  complex(dp), parameter :: z0 = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1 = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1= dcmplx(-1._dp, 0._dp)
  complex(dp), parameter :: zi = dcmplx( 0._dp, 1._dp)

contains

  ! The problem of this routine is that we wish not to
  ! overwrite the old half-inverted matrix, that would mean
  ! that we need to do the full calculation for each electrode
  ! (which can be quite time-consuming!)
  subroutine GF_Gamma_GF(Gf_tri, El, no, calc_parts, nwork, work)

    use alloc, only : re_alloc, de_alloc

    use class_zTriMat
    use m_ts_trimat_invert, only : TriMat_Bias_idxs
    use m_ts_electype

    implicit none

! *********************
! * INPUT variables   *
! *********************
    ! The Green function column
    type(zTriMat), intent(inout) :: Gf_tri
    type(Elec), intent(in) :: El ! contains: (Sigma - Sigma^dagger) ^T
    integer, intent(in) :: no ! The dimension of (Sigma - Sigma^dagger) ^T
    logical, intent(in) :: calc_parts(:)

! *********************
! * OUTPUT variables  *
! *********************
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)


    ! local variables
    complex(dp), pointer :: fGf(:), Gf(:), GGG(:)
    integer :: nr, np
    integer :: sIdx, eIdx
    integer :: cp, n
    integer :: sN, sNc

    integer :: lsPart, lePart
    integer :: BsPart, BePart

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

#ifndef TBTRANS
    call timer("GFGGF",1)
#else
#ifdef TBTRANS_TIMING
    call timer("GFGGF",1)
#endif
#endif

    ! tri-diagonal parts information
    nr = nrows_g(Gf_tri)
    np = parts(Gf_tri)

    ! Which parts are needed
    ! In this case we need to calculate till the end
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

    ! Capture the full elements
    fGf => val(Gf_tri,all=.true.)

    do n = lsPart , lePart

       if ( .not. calc_parts(n) ) cycle

       ! Calculate the \Gamma Gf^\dagger n,1
       sN = nrows_g(Gf_tri,n)
       if ( nwork < sN * no ) then
          print *,nwork,sN*no
          call die('Work size not big enough')
       end if

       ! correct to the quantities that is available
       BsPart = max(n-1,lsPart)
       BePart = min(n+1,lePart)

       call TriMat_Bias_idxs(Gf_tri,no,n,sIdx,eIdx)
       ! obtain the Gf in the respective column
       Gf => fGf(sIdx:eIdx)
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'T','C',no,sN,no, z1, El%Gamma, no, &
            Gf, sN, z0, work, no)
       
       ! Now we are ready to perform the multiplication
       ! for the requested region

#ifdef TRANSIESTA_DEBUG
       write(*,'(a,2(tr1,i0),a,2(tr1,i0))')'GfGGf at:',BsPart,n,' --',BePart,n
#endif
       
       ! this will populate in ascending column major order
       do cp = BsPart , BePart

          ! skip unneeded elements
          if ( .not. calc_parts(cp) ) cycle

          sNc = nrows_g(Gf_tri,cp)

          ! Retrieve Gf block
          call TriMat_Bias_idxs(Gf_tri,no,cp,sIdx,eIdx)
          Gf => fGf(sIdx:eIdx)

          ! retrieve the GGG block
          GGG => val(Gf_tri,cp,n)

          ! We need only do the product in the closest
          ! regions (we don't have information anywhere else)
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N', sNc, sN, no, z1, &
               Gf, sNc, work, no, z0, GGG, sNc)
          
       end do
       
    end do
       
#ifndef TBTRANS
    call timer("GFGGF",2)
#else
#ifdef TBTRANS_TIMING
    call timer("GFGGF",2)
#endif
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF


#ifdef TBTRANS
  subroutine dir_GF_Gamma_GF(Gf_tri, A_tri, r, pvt, El, calc_parts, TrGfG)
#else
  subroutine dir_GF_Gamma_GF(Gf_tri, A_tri, r, pvt, El, calc_parts)
#endif

    use alloc, only : re_alloc, de_alloc

    use class_zTriMat
    use m_ts_electype

    use m_trimat_invert, only: Xn_div_Cn_p1, Yn_div_Bn_m1

    implicit none

! *********************
! * INPUT variables   *
! *********************
    ! The BTD matrix with the diagonal blocks needed
    ! and _all_ X/B and Y/C
    type(zTriMat), intent(inout) :: Gf_tri
    ! This is on entry an empty matrix and will
    ! on exit it contains the spectral function
    ! for the electrode
    type(zTriMat), intent(inout) :: A_tri

    type(tRgn), intent(in) :: r ! the pivoting array for the sparse pattern
    type(tRgn), intent(in) :: pvt ! the pivoting of r back to the sparse pattern
    type(Elec), intent(in) :: El ! contains: (Sigma - Sigma^dagger) ^T
    logical, intent(in) :: calc_parts(:)
#ifdef TBTRANS
! **********************
! * OUTPUT variables   *
! **********************
    real(dp), intent(out), optional :: TrGfG
#endif

    
    ! local variables
    complex(dp), pointer :: fGf(:), fA(:), Gf(:), A(:)

    complex(dp), pointer :: ztmp(:)
    logical :: tmp_allocated
    integer :: nrtmp
    complex(dp), pointer :: zwork(:) => null()
    logical :: work_allocated 
    integer, pointer :: crows(:)

    integer :: i_Elec, idx_Elec, no, nb
#ifdef TBTRANS
    integer :: i_diag
#endif
    integer :: np
    integer :: in, n, i, off
    integer :: sPart, ePart, lPart
    integer :: idx, sIdx, eIdx, sIdxF, eIdxF
    integer :: sN, sNc, sNo
#ifndef TS_NOCHECKS
    integer :: oi
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE dGFGammaGF' )
#endif

#ifndef TBTRANS
    call timer("dGFGGF",1)
#else
#ifdef TBTRANS_TIMING
    call timer("dGFGGF",1)
#endif
#endif

    ! get electrode quantities
    no = El%o_inD%n
    
    ! intersecting parts
    sPart = huge(1)
    ePart = 0
#ifndef TS_NOCHECKS
    oi = 0
#endif
    do n = 1 , no
      i = pvt%r(El%o_inD%r(n))
#ifndef TS_NOCHECKS
      if ( oi > i ) then
        call rgn_print(r)
        call rgn_print(El%o_inD)
        print *,n,oi,i
        call die('BTD error in sorting of o_inD for electrode.')
      end if
      oi = i
#endif
      if ( i < sPart ) sPart = i
      if ( ePart < i ) ePart = i
    end do
    sPart = which_part(A_tri, sPart)
    ePart = which_part(A_tri, ePart)
    if ( ePart - sPart > 1 ) then
       print *,sPart, ePart
       call die('Gamma filling more than 2 blocks, how is that &
            &possible.')
    end if

    ! tri-diagonal parts information
    np = parts(Gf_tri)
    crows => cum_rows(A_tri)
    ! Capture the full elements
    fGf => val(Gf_tri, all=.true.)
    fA => val(A_tri, all=.true.)

    ! In the direct spectral function calculation
    ! we need to start from the diagonal part of
    ! the triple-product

    ! So we start by calculating the triple-product
    ! of the electrode diagonal elements
    sN = nrows_g(A_tri, sPart)
    ! get first column/row in this part
    i = crows(sPart) - sN + 1
    ! get starting index
    sIdx = index(A_tri,i,i) - 1
    
    ! do the same for the last part with
    sNc = nrows_g(A_tri, ePart)
    i = crows(ePart)
    ! get ending index
    eIdx = index(A_tri,i,i) + 1

    ! Figure out if the work array can be associated
    ! or needs to be allocated

    ! first calculate required work-array
    if ( sPart == ePart ) then
       ! only one block is needed
       sNc = sN
    else
       ! both blocks are needed
       sNc = sN + sNc
    end if
    ! get Gf.Gamma column size (maximum size)
    idx = no * sNc
    
    ! idx now contains the maximum number of
    ! elements needed to retain the part-column
    ! of the scattering matrix
    work_allocated = .false.
    ! get tmp work size
    sIdxF = sIdx ! not-used array in A matrix [start]
    eIdxF = size(fA) - eIdx + 1 ! not-used array in A matrix [end]
    if ( idx <= sIdxF ) then
       ! the work-array can be at the start of
       ! the A_tri array
       zwork => fA(1:idx)
       ! Associated the temporary work array
       if ( sIdxF - idx > eIdxF ) then
          ! the temporary array is largest in the
          ! beginning
          ztmp => fA(idx+1:sIdx)
       else
          ! anything
          ztmp => fA(eIdx:)
       end if
    else if ( idx <= eIdxF ) then
       zwork => fA(eIdx:eIdx+idx-1)
       ! Associated the temporary work array
       if ( sIdxF > eIdxF - idx ) then
          ! the temporary array is largest in the
          ! beginning
          ztmp => fA(1:sIdx)
       else
          ztmp => fA(eIdx+idx:)
       end if
    else
       work_allocated = .true.
       allocate(zwork(idx))
       ! Associated the temporary work array
       if ( sIdxF > eIdxF ) then
          ! the temporary array is largest in the
          ! beginning
          ztmp => fA(1:sIdx)
       else
          ztmp => fA(eIdx:)
       end if
    end if

    ! calculate number of rows we can simultaneously
    ! cobe with.
    ! Calculate a size according to 4 MB
    ! This corresponds roughly to a 500 X 500 matrix
    i = 4 * 1024._dp**2 / 16._dp / no
    i = max(2,i)
    nrtmp = size(ztmp) / no
    if ( nrtmp < i ) then
       nrtmp = i
       nullify(ztmp)
       allocate(ztmp(no*nrtmp))
       tmp_allocated = .true.
    else
       tmp_allocated = .false.
    end if


    ! Now calculate the column where we need G
    i_Elec = 1
    idx = 1
    n = sPart
    do while ( i_Elec <= no )

       ! get orbital index for the current column/row
       idx_Elec = pvt%r(El%o_inD%r(i_Elec))
       
       ! We start by copying over the Gfnn in blocks

       ! ... create a region of consecutive
       ! memory.
       if ( idx_Elec <= crows(sPart) ) then
         n = sPart
       else
         ! If we skip to a new block we must add sN
         ! to account for the new diagonal position
         idx = idx + nrows_g(A_tri, sPart) ! add the offset for the first
         n = ePart
       end if
       sN = nrows_g(A_tri, n)
       eIdx = crows(n)
       sIdx = eIdx - sN + 1

       ! we only require up to 'no' columns
       nb = 1
       do while ( i_Elec + nb <= no )
          i = pvt%r(El%o_inD%r(i_Elec+nb))
          ! In case it is not consecutive
          if ( i - idx_Elec /= nb ) exit
          ! In case the block changes, then
          ! we cut the block size here.
          if ( i < sIdx .or. eIdx < i ) exit
          nb = nb + 1
       end do
       
       ! Copy over this portion of the Gfnn array

       ! We now need to figure out the placement of the 
       ! Gf part that we have already calculated.
       Gf => val(Gf_tri,n,n)
       ! first column in the current diagonal block
       i = crows(n) - sN + 1
       sIdxF = (idx_Elec-i) * sN + 1
       eIdxF = (idx_Elec+nb-i) * sN

       ! Copy over diagonal block to its correct placement
       sIdx = idx
       eIdx = sN - 1
       do i = 1 , nb
          zwork(sIdx:sIdx+eIdx) = Gf(sIdxF:sIdxF+eIdx)
          sIdx = sIdx + sNc
          sIdxF = sIdxF + sN
       end do

       ! now we have the diagonal part of the Green function
       ! in the correct place

       ! Calculate the off-diagonal Green function in the regions
       ! -> above
       i = sN
       sIdx = idx ! next Gf
       eIdx = idx ! current Gf
       do in = n - 1 , sPart , - 1

          ! Number of orbitals in the other segment
          sNo = nrows_g(Gf_tri, in)
          sIdx = eIdx - sNo

          ! grab the Yn matrix to perform the 
          ! Gf off-diagonal calculation.
          Gf => Yn_div_Bn_m1(Gf_tri,in+1)

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,nb,i, zm1, Gf, sNo, &
               zwork(eIdx), sNc, z0, zwork(sIdx), sNc)

          ! correct next multiplication step
          eIdx = sIdx
          i = sNo
          
       end do
          
       ! Calculate the off-diagonal Green function in the regions
       ! -> below
       i = sN
       sIdx = idx ! current Gf
       eIdx = idx ! next Gf
       do in = n + 1 , ePart

          ! Number of orbitals in the other segment
          sNo = nrows_g(Gf_tri, in)
          eIdx = sIdx + i

          ! grab the Xn matrix to perform the 
          ! Gf off-diagonal calculation.
          Gf => Xn_div_Cn_p1(Gf_tri,in-1)

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,nb,i, zm1, Gf, sNo, &
               zwork(sIdx), sNc, z0, zwork(eIdx), sNc)

          sIdx = eIdx
          i = sNo
          
       end do
       
       ! Update current segment of the electrode copied entries.
       i_Elec = i_Elec + nb
       idx = idx + nb * sNc

       ! Revert offset in case the Gamma is split multiple times in the blocks
       if ( n /= sPart ) idx = idx - nrows_g(A_tri, sPart)

    end do

    ! now we have calculated the column for the blocks that
    ! is encompassed by the scattering matrix
#ifdef TBTRANS
    if ( present(TrGfG) ) TrGfG = 0._dp ! initialize
    ! initialize diagonal counter.
    ! Note that we know the columns are coinciding with
    ! the electrode scattering matrix.
    i_diag = 0
#endif
    
    ! Calculate the "diagonal" triple-matrix-product
    ! starting index:
    off = 0 ! offset from block start
    n = sPart
    sN = nrows_g(Gf_tri, sPart)
    idx_Elec = crows(sPart) - sN + 1 ! current row in BTD
    i_Elec = 1 ! loop row
    do while ( i_Elec <= sNc )

       ! Get current number of orbitals in this block
       sN = nrows_g(A_tri, n)

       ! the maximum size is nrtmp
       ! the minimum size is remaining elements in current block
       nb = min(nrtmp,sN - (i_Elec - off) + 1)

       ! do Gf.Gamma
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','T',nb,no,no, z1, zwork(i_Elec), sNc, &
            El%Gamma, no, z0, ztmp(1), nb)

#ifdef TBTRANS
       if ( present(TrGfG) ) then
          ! The idx_Elec counter tracks the row
          ! in the BTD matrix.
          ! But we need to figure out the diagonal
          ! element.
          ! Here we skip to the correct column, of
          ! the first index in the currently calculated
          ! Gf.G MM
          in = i_diag * nb + 1
          do i = 0 , nb - 1
             ! This should be fast because inDpvt is sorted
             if ( in_rgn(El%inDpvt,idx_Elec+i) ) then
                ! we have hit a diagonal element
                i_diag = i_diag + 1
                ! add the contribution
                TrGfG = TrGfG - aimag(ztmp(in))
                ! step to next column
                in = in + nb
             end if
             ! always step to the next row
             in = in + 1
          end do
       end if
#endif

       eIdxF = 1
       do in = sPart , ePart

          ! get number of columns of this part
          sNo = nrows_g(A_tri, in)
          fA => val(A_tri,n,in)
          
          ! do Gf.Gamma.Gf^dagger (with correct offset of fA)
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','C',nb,sNo,no, z1, ztmp(1), nb, &
               zwork(eIdxF), sNc, z0, fA(i_Elec-off), sN)

          eIdxF = eIdxF + sNo

       end do

       i_Elec = i_Elec + nb
       idx_Elec = idx_Elec + nb

       ! We are stepping to the next block
       if ( idx_Elec > crows(n) ) then
         off = off + sN
         n = n + 1
       end if

    end do

#ifdef TBTRANS
    if ( present(TrGfG) ) then
       ! Now we have:
       !   T = Tr[G \Gamma]
       ! Multiply by two to get the Tr[G^\dagger Gamma]
       ! contribution as well.
       ! See e.g. the discussion for GF_Gamma in
       ! m_tbt_tri_scat.F90
       TrGfG = TrGfG * 2._dp
    end if
#endif

    ! Deallocate ztmp if allocated
    if ( tmp_allocated ) then
       deallocate(ztmp)
    end if


    ! Now we have calculate the triple matrix product for
    ! the diagonal blocks where the scattering matrix
    ! lives.
    ! Now we simply need to propagate the result
    ! up and down the block matrix


    ! calculate above matrix
    lPart = max(1,sPart - 1)
    do n = 1 , sPart - 1
       if ( calc_parts(n) ) then
          lPart = max(1,n-1)
          exit
       end if
    end do

    ! Do calculation of spectral function
    ! in the above parts
    do n = sPart , lPart, -1
       
       ! we are now operating on block row n
       sN = nrows_g(A_tri,n)

       if ( n < sPart ) then

          ! this is the preceeding block
          sNo = nrows_g(A_tri,n+1)

          ! calculate n,n+1
          A => val(A_tri,n,n+1)
          ! get diagonal A
          Gf => val(A_tri,n+1,n+1)
          ! Get Y
          ztmp => Yn_div_Bn_m1(Gf_tri,n+1)
          
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sN,sNo,sNo, zm1, ztmp, sN, &
               Gf, sNo, z0, A, sN)

          ! skip calculating the diagonal block
          if ( n == lPart .and. .not. calc_parts(n) ) exit
          
          ! calculate n,n
          Gf => A(:)
          A => val(A_tri,n,n)
          
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','C',sN,sN,sNo, zm1, Gf, sN, &
               ztmp, sN, z0, A, sN)

       end if
       
       if ( n <= 1 ) exit
       
       sNo = nrows_g(A_tri,n-1)
       ! calculate n,n-1
       A => val(A_tri,n,n-1)
       ! get diagonal A
       Gf => val(A_tri,n,n)
       ! Get Y
       ztmp => Yn_div_Bn_m1(Gf_tri,n)
       
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','C',sN,sNo,sN, zm1, Gf, sN, &
            ztmp, sNo, z0, A, sN)
       
    end do


    ! calculate below matrx
    lPart = min(np,ePart + 1)
    do n = np , ePart + 1, -1
       if ( calc_parts(n) ) then
          lPart = min(np,n+1)
          exit
       end if
    end do

    ! Do calculation of spectral function
    ! in the above parts
    do n = ePart , lPart

       ! we are now operating on block row n
       sN  = nrows_g(A_tri,n)

       if ( ePart < n ) then

          ! this is the preceeding block
          sNo = nrows_g(A_tri,n-1)

          ! calculate n,n-1
          A => val(A_tri,n,n-1)
          ! get diagonal A
          Gf => val(A_tri,n-1,n-1)
          ! Get X
          ztmp => Xn_div_Cn_p1(Gf_tri,n-1)
          
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sN,sNo,sNo, zm1, ztmp, sN, &
               Gf, sNo, z0, A, sN)

          ! skip calculating the diagonal block
          if ( n == lPart .and. .not. calc_parts(n) ) exit

          ! and calculate n,n
          Gf => A(:)
          A => val(A_tri,n,n)
          
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','C',sN,sN,sNo, zm1, Gf, sN, &
               ztmp, sN, z0, A, sN)

       end if
       
       if ( n >= np ) exit
       
       sNo = nrows_g(A_tri,n+1)
       ! calculate n,n+1
       A => val(A_tri,n,n+1)
       ! get diagonal A
       Gf => val(A_tri,n,n)
       ! Get X
       ztmp => Xn_div_Cn_p1(Gf_tri,n)
       
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','C',sN,sNo,sN, zm1, Gf, sN, &
            ztmp, sNo, z0, A, sN)

    end do

    if ( work_allocated ) then
       deallocate(zwork)
    end if

#ifndef TBTRANS
    call timer("dGFGGF",2)
#else
#ifdef TBTRANS_TIMING
    call timer("dGFGGF",2)
#endif
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS dGFGammaGF' )
#endif

  end subroutine dir_GF_Gamma_GF

  function has_full_part(N_tri_part,tri_parts, &
       part,io1,io2) result(has)
    integer, intent(in) :: N_tri_part, tri_parts(N_tri_part), part, io1, io2
    logical :: has
    integer :: i, io

    io = 1
    do i = 1 , part - 1
       io = io + tri_parts(i)
    end do
    
    has = io1 <= io .and. &
         io + tri_parts(i) - 1 <= io2

  end function has_full_part

  ! Generic routine for inserting the self-energies in the 
  ! tri-diagonal matrices
  subroutine insert_Self_Energies(Gfinv_tri, Gfinv, pvt, El)
    use m_region
    use m_ts_electype
    use class_zTriMat
    type(zTriMat), intent(inout) :: GFinv_tri
    complex(dp), intent(inout) :: Gfinv(:)
    type(tRgn), intent(in) :: pvt
    type(Elec), intent(in) :: El

    integer :: no, off, i, j, ii, idx
    
    no = TotUsedOrbs(El)
    off = El%idx_o - 1

    if ( El%Bulk ) then
!$OMP do private(j,i,ii,idx)
       do j = off + 1 , off + no
          ii = (j-off-1) * no
          do i = 1 , no
             idx = index(GFinv_tri,pvt%r(i+off),pvt%r(j))
             Gfinv(idx) = El%Sigma(ii+i)
          end do
       end do
!$OMP end do
    else
!$OMP do private(j,i,ii,idx)
       do j = off + 1 , off + no
          ii = (j-off-1) * no
          do i = 1 , no
             idx = index(GFinv_tri,pvt%r(i+off),pvt%r(j))
             Gfinv(idx) = Gfinv(idx) - El%Sigma(ii+i)
          end do
       end do
!$OMP end do
    end if
    
  end subroutine insert_Self_Energies

end module m_ts_tri_scat
