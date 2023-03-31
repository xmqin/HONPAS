! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_tbt_tri_scat

  use precision, only : dp
  use units, only : Pi
  use m_region
  
  use class_zTriMat

  use m_ts_tri_scat, only : GF_Gamma_GF, dir_GF_Gamma_GF
  use m_ts_tri_common, only : GFGGF_needed_worksize

  use m_ts_electype

  implicit none

  private

  public :: A_DOS   ! Spectral function density of states
  public :: GF_DOS  ! Green's function density of states
  public :: A_Gamma ! Calculate the transmission from spectral function . Gamma
  public :: A_Gamma_Block ! Calculate the transmission from spectral function . Gamma (in block form)
  public :: TT_eigen ! Eigenvalue calculation of the transmission eigenvalues
  public :: GF_Gamma ! Calculate the transmission from Green function . Gamma (same-lead contribution)
  public :: GF_T, GF_T_solve
  
  public :: insert_Self_Energy
  public :: insert_Self_Energy_Dev

  ! From ts_tri_scat
  public :: GF_Gamma_GF
  public :: dir_GF_Gamma_GF
  public :: GFGGF_needed_worksize
#ifdef NCDF_4
  public :: GF_COP, A_COP
  public :: GF_COHP_add_dH, A_COHP_add_dH
  public :: orb_current
  public :: orb_current_add_dH
  public :: GF_DM, A_DM
#endif

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)
  complex(dp), parameter :: zi  = dcmplx( 0._dp, 1._dp)

contains

  ! Calculate the DOS from a non-fully calculated Green function.
  ! We assume that the diagonal Green function matrices are already calculated
  ! and the remaining \tilde X and \tilde Y matrices are present.
  !    all GF_nn are in Gf_tri
  !    all \tilde Yn and \tilde Xn are in Gf_tri
  ! After this routine, all off-diagonal Gf blocks are in work_tri (correctly
  ! positioned).
  ! I.e. the full Gf in the blocks can be extracted from Gf_tri and work_tri.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine GF_DOS(r,Gfd_tri,Gfo_tri,S_1D,pvt,DOS)
    use class_Sparsity
    use class_zSpData1D
    
    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(zSpData1D), intent(inout) :: S_1D ! (transposed S(k))
    type(tRgn), intent(in) :: pvt
    real(dp), intent(out) :: DOS(r%n)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: np, n, no_o, no_i
    integer :: br, io, ind, bc
    real(dp) :: lDOS

#ifdef TBTRANS_TIMING
    call timer('Gf-DOS',1)
#endif

    np = parts(Gfd_tri)
    
    ! First calculate all off-diagonal green-function elements
    no_o = nrows_g(Gfd_tri,1)
    no_i = nrows_g(Gfd_tri,2)
    call calc(2,1)
    do n = 2, np - 1
      no_o = nrows_g(Gfd_tri,n)
      no_i = nrows_g(Gfd_tri,n + 1)
      call calc(n+1,n)
      no_i = nrows_g(Gfd_tri,n - 1)
      call calc(n-1,n)
    end do
    no_o = nrows_g(Gfd_tri,np)
    no_i = nrows_g(Gfd_tri,np-1)
    call calc(np-1,np)

    ! At this point we have calculated all Green function matrices
    ! All diagonal elements are in Gfd_tri,
    ! all off-diagonal elements are in Gfo_tri

    ! The DOS per orbital is calculated like this (.=matrix multiplication):
    !   DOS(io) = - Im[ (Gf-Gf^\dagger) . S ](io,io) / Pi
    !           = - \sum_jo Im[ {Gf(io, jo)-Gf^\dagger(io,jo)} * S(jo, io)] / Pi
    !
    ! The fact that we need Gf - Gf^\dagger can be checked
    ! by a simple tight-binding calculation with large overlap matrices
    ! and *very* small eta values (and using sum_elec ADOS == DOS).
    ! In this case the k-resolved DOS is only correct if one uses
    ! the above equation. It should however be noted that the full
    ! DOS is independent on Gf or Gf - Gf^\dagger choice!

    sp => spar(S_1D)
    S => val(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

!$OMP parallel do default(shared), private(br,io,lDOS,ind,bc,GfGfd)
    do br = 1, r%n
      io = r%r(br)

      ! Loop columns in S(k)^T (actually the rows)
      lDOS = 0._dp
      do ind = l_ptr(io) + 1, l_ptr(io) + ncol(io)
        bc = pvt%r(l_col(ind))
        if ( bc > 0 ) then
          call calc_GfGfd(br, bc, GfGfd)
          lDOS = lDOS + dimag( GfGfd * S(ind) )
        end if
      end do
       
      DOS(br) = - lDOS / (2._dp * Pi)
       
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('Gf-DOS',2)
#endif

  contains

    subroutine calc(m,n)
      integer, intent(in) :: m,n
      complex(dp), pointer :: Gf(:), Mnn(:), XY(:)

      XY  => val(Gfd_tri,m,n)
      Mnn => val(Gfd_tri,n,n)
      Gf  => val(Gfo_tri,m,n)
      
      ! We need to calculate the 
      ! Mnm1n/Mnp1n Green's function
#ifdef USE_GEMM3M
      call zgemm3m( &
#else
      call zgemm( &
#endif
           'N','N',no_i,no_o,no_o, &
           zm1, XY,no_i, Mnn,no_o,z0, Gf,no_i)
      
    end subroutine calc

    subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G
      integer :: p_r, i_r, p_c, i_c, i

      call part_index(Gfo_tri, br, p_r, i_r)
      call part_index(Gfo_tri, bc, p_c, i_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfd(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        G = G - conjg(Gfd(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      else
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfo(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        i = index_block(Gfo_tri, p_c, p_r)
        G = G - conjg(Gfo(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      end if

    end subroutine calc_GfGfd

  end subroutine GF_DOS


  ! Calculate the DOS from a fully calculated spectral function.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine A_DOS(r,A_tri,S_1D,pvt,DOS)
    use class_Sparsity
    use class_zSpData1D

    type(tRgn), intent(in) :: r ! BTD matrix elements
    type(zTriMat), intent(inout) :: A_tri
    type(zSpData1D), intent(inout) :: S_1D ! (transposed S(k))
    type(tRgn), intent(in) :: pvt ! from sparse matrix to BTD
    real(dp), intent(out) :: DOS(r%n)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), A(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: io, ind, idx, br, bc
    real(dp) :: lDOS

#ifdef TBTRANS_TIMING
    call timer('A-DOS',1)
#endif

    ! Get data arrays
    A => val(A_tri)

    sp => spar(S_1D)
    S => val(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    ! The DOS per orbital is calculated like this:
    !   ADOS(io) = Re[ A . S ](io,io) / 2Pi
    !            = Re[A(io, jo) * S(jo, io)] / 2Pi

!$OMP parallel do default(shared), private(br,io,lDOS,ind,bc,idx)
    do br = 1, r%n
      io = r%r(br)

      ! Loop columns in S(k)^T (actually the rows)
      lDOS = 0._dp
      do ind = l_ptr(io) + 1, l_ptr(io) + ncol(io)
        bc = pvt%r(l_col(ind))
        if ( bc > 0 ) then
          idx = index(A_tri, br, bc)
          lDOS = lDOS + dreal( A(idx) * S(ind) )
        end if
      end do
      
      DOS(br) = lDOS / (2._dp * Pi)
      
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('A-DOS',2)
#endif

  end subroutine A_DOS
  
#ifdef NCDF_4

  ! Calculate the COOP contribution from a fully calculated Green function.
  ! We assume that the Green function distribution like this:
  !    all GF_nn are in Gfd_tri (diagonal)
  !    all GF_mn (m/=n) are in Gfo_tri (off-diagonal)
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine GF_COP(r,Gfd_tri,Gfo_tri,pvt,sp,M,sc_off,k,ph,COP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: M(:) ! S for COOP, H for COHP
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COP ! COOP or COHP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:), c_col(:)
    integer :: no_u, br, io, ind, iind, bc
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('Gf-COP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COOP/COHP by looping the sparse matrix
    ! The following discussion in concerning COOP, but
    ! there is no ambiguity in the two methods.
    
    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = - Im{ [Gf - Gf^\dagger](io,jo) * S(jo,io) * e^(ik.R) } / 2Pi
    ! Here we want:
    !   DOS(io) = \sum_jo COOP(io,jo)
    ! since we know that COOP(io,jo) is the io -> jo DOS.
    ! As COOP is interesting in the supercell picture we have
    ! to calculate it with the daggered component (Gf - Gf^\dagger) (also why we need /2).
    ! Note that this is not necessary if S is S(k). I.e. it is because
    ! we want the cross-cell COOP curves as well.

    ! Create the phases
    ! Since we have to do Gf.S we simply
    ! create S(-k) (which is S^T)
    ! and thus get the correct values.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, - &
          k(1) * sc_off(1,io) - &
          k(2) * sc_off(2,io) - &
          k(3) * sc_off(3,io))) / (2._dp * Pi)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COP)
    C => val(COP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

    C(:) = 0._dp

!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss,GfGfd)
    do br = 1, r%n
      io = r%r(br)
      
      ! Get lookup columns for the COOP
      call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))
      
      ! Loop on overlap entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)
        
        ! Check if the orbital exists in the region
        iind = cptr(io) + ssearch_find(ss, l_col(ind))
        
        ! if zero the element does not exist
        ! This is the case on the elements connecting out
        ! of the device region
        if ( iind <= cptr(io) ) cycle
        
        ! COOP(iind) = - Im[ (G(io,jo) - G^\dagger(io,jo)) * S(jo,io) ] / 2Pi
        bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
        call calc_GfGfd(br, bc, GfGfd)

        C(iind) = -aimag( GfGfd * M(ind) * ph( (l_col(ind)-1)/no_u ))

      end do
          
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('Gf-COP',2)
#endif

  contains

    subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G
      integer :: p_r, i_r, p_c, i_c, i
      
      call part_index(Gfo_tri, br, p_r, i_r)
      call part_index(Gfo_tri, bc, p_c, i_c)

      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfd(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        G = G - conjg(Gfd(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      else
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfo(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        i = index_block(Gfo_tri, p_c, p_r)
        G = G - conjg(Gfo(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      end if

    end subroutine calc_GfGfd

  end subroutine GF_COP

  subroutine Gf_COHP_add_dH(dH_1D,sc_off,k,ph,Gfd_tri,Gfo_tri,r,COHP,pvt)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData1D), intent(in) :: dH_1D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: COHP
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:)
    type(Sparsity), pointer :: c_sp
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd
    real(dp), pointer :: C(:)
    integer :: no_u, br, io, jo, ind, iind

#ifdef TBTRANS_TIMING
    call timer('COHP-Gf-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_1D)
    dH => val(dH_1D)
    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    ! Create the phases
    ! We are using the explicit H(j, i) and thus the phases are consistent with +
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, + &
          k(1) * sc_off(1,io) + &
          k(2) * sc_off(2,io) + &
          k(3) * sc_off(3,io))) / (2._dp * Pi)
    end do

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

!$OMP parallel do default(shared), private(br,io,iind,jo,ind,col,GfGfd)
    do br = 1, r%n
      io = r%r(br)
      
      ! Loop on the COHP indices
      do iind = cptr(io) + 1, cptr(io) + cncol(io)

        ! Here we will calculate the COHP contribution from dH
        !  COHP(iind) = -Im{ [Gf(io, jo) - Gf^\dagger(io,jo)] * dH(jo, io) } / 2pi

        ! Since we are looping the dH indices we have to 

        ! Get column Gf orbital
        jo = ucorb(ccol(iind), no_u)

        ! Check if the jo,io orbital exists in dH
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))

          ! Note that we here find the dH(jo,io) value (in the supercell picture)
          ind = l_ptr(jo) + SFIND(col, TO(ccol(iind)) + io)

          if ( ind > l_ptr(jo) ) then

            call calc_GfGfd(br, pvt%r(jo), GfGfd)
            ! COHP(iind) += - Im[ (G(io,jo) - G^\dagger(io,jo)) * dH(jo,io)] / 2Pi
            C(iind) = C(iind) &
                - aimag( GfGfd * dH(ind) * ph( (l_col(ind)-1)/no_u ))

          end if

        end if
        
      end do
      
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('COHP-Gf-dH',2)
#endif

  contains

    subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G
      integer :: p_r, i_r, p_c, i_c, i

      call part_index(Gfo_tri, br, p_r, i_r)
      call part_index(Gfo_tri, bc, p_c, i_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfd(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        G = G - conjg(Gfd(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      else
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfo(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        i = index_block(Gfo_tri, p_c, p_r)
        G = G - conjg(Gfo(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      end if

    end subroutine calc_GfGfd

    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i

      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)
        
        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('Gf_COHP_add_dH: could not find transpose supercell index')

    end function TO
    
  end subroutine Gf_COHP_add_dH

  ! Calculate the COOP contribution from a fully calculated spectral function.
  ! This routine utilizes the sparse matrix as a loop, instead of looping
  ! all BTD matrix elements.
  ! This turns out to be much faster for (at least tight-binding calculations).
  subroutine A_COP(r,A_tri,pvt,sp,M,sc_off,k,ph,COP)
    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: A_tri
    type(tRgn), intent(in) :: pvt
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: M(:) ! S for COOP, H for COHP
    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(dSpData1D), intent(inout) :: COP ! COOP or COHP

    type(Sparsity), pointer :: c_sp
    real(dp), pointer :: C(:)
    complex(dp), pointer :: A(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer :: no_u, br, io, ind, iind, bc
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('A-COP',1)
#endif

#ifdef TBT_PHONON
    call die('Currently not implemented for PHtrans')
#endif

    ! Extract COOP/COHP by looping the sparse matrix
    ! The following disôcussion in concerning COOP, but
    ! there is no ambiguity in the two methods.

    ! The COOP calculation can be written as
    !
    !   COOP(io,jo) = Re{ A(io,jo) * S(jo,io) * e^(ik.R) } / 2Pi
    ! Here we want:
    !   ADOS(io) = \sum_jo COOP(io,jo)
    ! since we know that COOP(io,jo) is the io -> jo ADOS.
    ! Note that this is not necessary if S is S(k). I.e. it is because
    ! we want the cross-cell COOP curves as well.

    ! Create the phases
    ! Since we have to do A.S we simply
    ! create the S(-k) (which is S^T)
    ! and thus get the correct values.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, - &
          k(1) * sc_off(1,io) - &
          k(2) * sc_off(2,io) - &
          k(3) * sc_off(3,io))) / (2._dp * Pi)
    end do

    call attach(sp,nrows_g=no_u, n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    c_sp => spar(COP)
    C => val(COP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)

    A => val(A_tri)

    C(:) = 0._dp

!$OMP parallel do default(shared), private(br,io,ind,iind,bc,ss)
    do br = 1, r%n
      io = r%r(br)

      ! Get lookup columns for the COOP
      call ssearch_init(ss, ccol(cptr(io)+1:cptr(io)+cncol(io)))

      ! Loop on overlap entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

        ! Check if the orbital exists in the region
        iind = cptr(io) + ssearch_find(ss, l_col(ind))

        ! if zero the element does not exist
        ! This is the case on the elements connecting out
        ! of the device region
        if ( iind <= cptr(io) ) cycle

        ! COOP(iind) = Re[ A(io,jo) * S(jo,io) ] / (2 pi)
        bc = pvt%r(ucorb(l_col(ind),no_u)) ! pivoted orbital index in tri-diagonal matrix
        bc = index(A_tri,br,bc)

        C(iind) = real(A(bc) * M(ind) * ph( (l_col(ind)-1)/no_u ), dp)

      end do
          
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('A-COP',2)
#endif

  end subroutine A_COP

  subroutine A_COHP_add_dH(dH_1D,sc_off,k,ph,A_tri,r,COHP,pvt)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData1D), intent(in) :: dH_1D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: COHP
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:)
    type(Sparsity), pointer :: c_sp
    integer, pointer :: cncol(:), cptr(:), ccol(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp), pointer :: A(:)
    real(dp), pointer :: C(:)
    integer :: no_u, iu, io, i, ind, iind, jo, iA

#ifdef TBTRANS_TIMING
    call timer('COHP-A-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_1D)
    dH => val(dH_1D)
    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    c_sp => spar(COHP)
    C => val(COHP)
    call attach(c_sp, n_col=cncol, list_ptr=cptr, list_col=ccol)
    
    ! Create the phases
    ! We are using the explicit H(j, i) and thus the phases are consistent with +
    do i = 1 , size(sc_off, dim=2)
      ph(i-1) = cdexp(dcmplx(0._dp, + &
          k(1) * sc_off(1,i) + &
          k(2) * sc_off(2,i) + &
          k(3) * sc_off(3,i))) / (2._dp * Pi)
    end do

    A => val(A_tri)
    
!$OMP parallel do default(shared), private(iu,io,iind,jo,ind,col,iA)
    do iu = 1, r%n
      io = r%r(iu)

      ! Loop on the COHP indices
      do iind = cptr(io) + 1, cptr(io) + cncol(io)

        ! Here we will calculate the COHP contribution from dH
        !  COHP(iind) == A(io, jo) * dH(jo, io) / 2pi

        ! Get column A orbital
        jo = ucorb(ccol(iind), no_u)

        ! Check if the jo,io orbital exists in dH
        if ( l_ncol(jo) > 0 ) then
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))

          ! Note that we here find the dH(jo,io) value (in the supercell picture)
          ind = l_ptr(jo) + SFIND(col, TO(ccol(iind)) + io)

          if ( ind > l_ptr(jo) ) then

            iA = index(A_tri,iu,pvt%r(jo)) ! A_ij

            ! COHP                    Aij  * Hji
            C(iind) = C(iind) + real(A(iA) * dH(ind) * ph( (l_col(ind)-1)/no_u ), dp)

          end if

        end if

      end do
       
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('COHP-A-dH',2)
#endif
    
  contains
    
    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i
      
      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)

        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('A_COHP_add_dH: could not find transpose supercell index')

    end function TO
    
  end subroutine A_COHP_add_dH


#ifdef NOT_WORKING
  ! A simple routine to calculate the DOS
  ! from a partially calculated GF
  ! When entering this routine Gf_tri
  ! should contain:
  ! all GF_nn
  ! all Yn/Bn-1 and all Xn/Cn+1
  ! This lets us calculate all entries
  subroutine GF_DOS_proj(r,Gf_tri,S_1D,N_mol,mols,DOS,bGfk,nwork,work)
    use class_Sparsity
    use class_zSpData1D

    use m_tbt_proj

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gf_tri
    type(zSpData1D), intent(inout) :: S_1D
    integer, intent(in) :: N_mol
    type(tProjMol), intent(in) :: mols(N_mol)
    real(dp), intent(out) :: DOS(r%n)
    complex(dp), intent(out) :: bGfk(:)
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), Gf(:), Mnn(:), XY(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: off1, off2, n, in
    integer :: jo, ii, i, j, no_o, no_i, ind, np, iD

    ! For looping the molecule projections
    integer :: Ns, Nl, Nm_dos
    integer :: im, ip, idx, no, i_o
    integer :: step_o

    S  => val(S_1D)
    sp => spar(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)
    
    ! Initialize DOS to 0
    DOS(:) = 0._dp

    off2 = 0
    np = parts(Gf_tri)

    ! Find maximum size of molecule orbitals
    no = maxval(mols(:)%orb%n)
    ! Calculate the size of the calculated bra at each index
    no_i = 0
    Nm_dos = 0
    do im = 1 , N_mol
      if ( .not. mols(im)%DOS ) cycle
      Nm_dos = Nm_dos + 1
      no_i = no_i + size(mols(im)%proj) * no
    end do
    ! Find maximum work size needed to retain the Gf
    no_o = 0
    do n = 1 , np - 1
      no_o = max(no_o,nrows_g(Gf_tri,n)*nrows_g(Gf_tri,n+1))
    end do
    ! Get the starting position of the projection matrices
    idx = no_o + 1
    ! Calculate maximum number of state orbitals we can
    ! accomodate simultaneously
    max_p = (nwork - no_o) / no_i
    if ( max_p < 1 ) then
      call die('Work size for projection of Gf not sufficient. &
          Try and use fewer projections, or simply do not calculate &
          the DOS projection.')
    end if

    do n = 1 , np

      no_o = nrows_g(Gf_tri,n)
      
      ! Calculate the step size for the projection
      ! on this column
      step_o = min(max_p,no_o)
      
      ! Loop over smaller group of columns in this block-column
      do i_o = 1 , no_o, step_o
        
        im = 0
        do i = 1 , N_mol
          if ( .not. mols(i)%DOS ) cycle
          no = mols(i)%orb%n
          Ns = size(mols(i)%proj)
          ! step calculated DOS for molecule
          im = im + 1
!$OMP parallel do default(shared), private(j,ip,ii)
          do j = 1 , step_o
            ! Calculate the projection matrix on these column
            ! indices
            do ip = 1 , Ns
              ! We have all molecules
              ii = idx + (((im-1)*step_o+j-1)*Ns+ip-1) * no + 1
              call proj_state_bra(mols(i),mols(i)%proj(ip), &
                  i_o+j, zwork(ii:ii+no-1) )
            end do
          end do
!$OMP end parallel do
        end do
        
        do in = max(1,n-1) , min(n+1,np)
          
          no_i = nrows_g(Gf_tri,in)
          
          if ( in < n ) then
            off1 = off2 - no_i
          else if ( n < in ) then
            off1 = off2 + no_o
          else
            off1 = off2
          end if
          
          if ( in == n ) then
            ! Retrieve the central part of the
            ! matrix
            Gf => val(Gf_tri,n,n)
            ! re-point
            Gf => Gf((i_o-1)*no_o+1:)
            
          else
            
            XY  => val(Gf_tri,in,n)
            Mnn => val(Gf_tri,n,n)
            ! re-point
            Mnn => Mnn((i_o-1)*no_o+1:)
            
            Gf  => work(1:no_o*no_i)
            
            ! We need to calculate the 
            ! Mnm1n/Mnp1n Green's function
#ifdef USE_GEMM3M
            call zgemm3m( &
#else
            call zgemm( &
#endif
                  'N','N',no_i,step_o,no_o, &
                  zm1, XY,no_i, Mnn,no_o,z0, Gf,no_i)

          end if

!$OMP parallel do default(shared), private(j,ii,jo,ind,i,ip,im,iD,lcol)
          do j = 1 , step_o
            ii = (j-1) * no_i
            iD = off2 + j
            jo = r%r(iD)
            lcol => l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo))
            ! get the equivalent one in the
            ! overlap matrix
            ! REMEMBER, S is transposed!
            ! Hence we do not need conjg :)
            do i = 1 , no_i
              ind = SFIND(lcol,r%r(off1+i))
              if ( ind == 0 ) cycle
              ind = l_ptr(jo) + ind
              DOS(iD) = DOS(iD) - dimag( Gf(ii+i) * S(ind) )
            end do
          end do
!$OMP end parallel do

        end do
        
        ! Update the offset
        off2 = off2 + no_o
        
      end do
      
    end do

    call dscal(r%n, 1._dp / Pi, DOS, 1)

  contains
    
    subroutine calc_state_Gf(N_mol,mols,Gf,step_o,zw,bGfk)
      
      iG = 0
      im = 0
      do i = 1 , N_mol
        if ( .not. mols(i)%DOS ) cycle
        no = mols(i)%orb%n
        Ns = size(mols(i)%proj)
        ! step calculated DOS for molecule
        im = im + 1
        do ip = 1 , Ns
          ! We have all molecules
          ii = idx + (((im-1)*step_o+j-1)*Ns+ip-1) * no + 1
          iG = iG + 1
          bGfk(iG) = bGfk(iG) + zw(
        end do
      end do
    end subroutine calc_state_Gf

  end subroutine GF_DOS_PROJ

#endif
#endif

  ! The simplest routine to do the transport calculation
  ! It takes the spectral function and multiplies it with
  ! the scattering matrix of the down-projected self-energy
  ! and calculates the transmission.
  ! We do this by taking advantage of the transposed scattering
  ! matrix: \Gamma
  subroutine A_Gamma(A_tri,El,T)

    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(Elec), intent(in) :: El ! contains Gamma == (Sigma - Sigma^\dagger)^T
    real(dp), intent(out) :: T

    ! Here we need a double loop
    integer :: no
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer :: o
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)

    ! External BLAS routine
    complex(dp), external :: zdotu
    
#ifdef TBTRANS_TIMING
    call timer('A-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region
    T = 0._dp

    ! Loop columns
    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(A_tri,El,i_Elec,in,ii)
      isN = nrows_g(A_tri,in)
      
      ! Get starting placement of column in the current block
      ! of the spectral function (zero based)
      if ( in == 1 ) then
        A_i = El%inDpvt%r(i_Elec) - 1
      else
        A_i = El%inDpvt%r(i_Elec) - crows(in-1) - 1
      end if

      if ( ii == no ) then

        ! The easy calculation, note that ii == no, only
        ! if the entire electrode sits in one block
        A => val(A_tri,in,in)
        do o = 0 , no - 1
          T = T + zdotu(no,A((A_i+o)*isN+A_i+1),1,El%Gamma(o*no+1),1)
        end do

        ! Quick break of loop
        exit

      end if
      
      ! Loop rows
      j_Elec = 1
      do while ( j_Elec <= no ) 

        ! We start by creating a region of consecutive memory.
        call consecutive_index(A_tri,El,j_Elec,jn,jj)
        jsN = nrows_g(A_tri,jn)

        ! Get the block with the spectral function
        A => val(A_tri,jn,in)

        if ( jn == 1 ) then
          A_j = El%inDpvt%r(j_Elec)
        else
          A_j = El%inDpvt%r(j_Elec) - crows(jn-1)
        end if

        do o = 0 , ii - 1
          T = T + zdotu(jj,A((A_i+o)*jsN+A_j),1, &
              El%Gamma((i_Elec-1+o)*no+j_Elec),1)
        end do

        j_Elec = j_Elec + jj

      end do

      i_Elec = i_Elec + ii

    end do

#ifdef TBTRANS_TIMING
    call timer('A-Gamma',2)
#endif
    
  end subroutine A_Gamma

  
  ! On entry A_tri is the spectral function
  ! on return the first El%o_inD%n x El%o_inD%n will be the
  ! G.Gamma.Gf.El%Gamma matrix
  ! This will enable eigenvalue calculators and possibly
  ! speed up the calculation of the transmission.
  subroutine A_Gamma_Block(A_tri,El,T,nwork,work)

    use intrinsic_missing, only : transpose, trace
    
    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(Elec), intent(inout) :: El
    real(dp), intent(out) :: T
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)

    ! Here we need a double loop
    integer :: no
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)
    complex(dp) :: z
    
#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n
    if ( no ** 2 > nwork ) then
       call die('A_Gamma_Block: Insufficient work-size')
    end if

    ! "sadly" Gamma is saved in transposed form, hence
    ! we transpose, and return it to original form, when returning
    call transpose(no,El%Gamma)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Loop columns
    i_Elec = 1
    ! The first column calculation initializes the result
    z = z0
    do while ( i_Elec <= no ) 
      
      ! We start by creating a region of consecutive memory.
      call consecutive_index(A_tri,El,i_Elec,in,ii)
      isN = nrows_g(A_tri,in)

      ! Get starting placement of column in the current block
      ! of the spectral function (zero based)
      if ( in == 1 ) then
        A_i = El%inDpvt%r(i_Elec) - 1
      else
        A_i = El%inDpvt%r(i_Elec) - crows(in-1) - 1
      end if

      if ( ii == no ) then
        ! The easy calculation, note that ii == no, only
        ! if the entire electrode sits in one block
        A => val(A_tri,in,in)

#ifdef USE_GEMM3M
        call zgemm3m( &
#else
        call zgemm( &
#endif
            'N','N',no,no,no, z1, A(A_i*(isN+1)+1), isN, &
            El%Gamma(1), no, z0, work(1), no)

        ! Quick break of loop
        exit

      end if

      ! Loop rows
      j_Elec = 1
      do while ( j_Elec <= no ) 

        ! We start by creating a region of consecutive memory.
        call consecutive_index(A_tri,El,j_Elec,jn,jj)
        jsN = nrows_g(A_tri,jn)

        ! Get the block with the spectral function
        A => val(A_tri,jn,in)

        if ( jn == 1 ) then
          A_j = El%inDpvt%r(j_Elec)
        else
          A_j = El%inDpvt%r(j_Elec) - crows(jn-1)
        end if

#ifdef USE_GEMM3M
        call zgemm3m( &
#else
        call zgemm( &
#endif
            'N','N',jj,no,ii, z1, A(A_i*jsN + A_j), jsN, &
            El%Gamma(i_Elec), no, z, work(j_Elec), no)

        j_Elec = j_Elec + jj

      end do
       
      i_Elec = i_Elec + ii
      ! Now we have already filled the first entries, sum...
      z = z1
      
    end do
    
    ! Calculate transmission
    T = dreal(trace(no,work))
    
    ! Now we have the square matrix product
    !   tt = G \Gamma_1 G^\dagger \Gamma_El
    
    call transpose(no,El%Gamma)
    
#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',2)
#endif
    
  end subroutine A_Gamma_Block


  subroutine TT_eigen(n,tt,nwork,work,eig)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: tt(n*n)
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    complex(dp), intent(inout) :: eig(n)

    real(dp) :: rwork(n*2)
    complex(dp) :: z
    integer :: i, j

#ifdef TBTRANS_TIMING
    call timer('TT-eig',1)
#endif

    ! To remove any singular values we add a 1e-3 to the diagonal
    do i = 1 , n
      tt((i-1)*n+i) = tt((i-1)*n+i) + 1.e-3_dp
    end do

    call zgeev('N','N',n,tt,n,eig,work(1),1,work(1),1, &
        work,nwork,rwork,i)
    if ( i /= 0 ) then
      print *,i
      call die('TT_eigen: Could not calculate eigenvalues.')
    end if

    ! Sort the eigenvalues, and simultaneously shift them back
    eig(1) = eig(1) - 1.e-3_dp
    do i = 2 , n
      eig(i) = eig(i) - 1.e-3_dp
      do j = 1 , i - 1
        if ( dreal(eig(j)) < dreal(eig(i)) ) then
          z = eig(j)
          eig(j) = eig(i)
          eig(i) = z
        end if
      end do
    end do

#ifdef TBTRANS_TIMING
    call timer('TT-eig',2)
#endif
    
  end subroutine TT_eigen
  
  subroutine GF_Gamma(Gfcol,El,T)

    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(Elec), intent(inout) :: El
    real(dp), intent(out) :: T

    complex(dp), pointer :: Gf(:)
    complex(dp), pointer :: z(:)

    integer :: no, np
    integer :: i, ii, i_Elec
    integer, pointer :: crows(:)

    integer :: sN, n, nb
    ! BLAS routines
    complex(dp), external :: zdotu, zdotc

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',1)
#endif

    no = El%inDpvt%n
    np = parts(Gfcol)
    crows => cum_rows(Gfcol)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    T = 0._dp

    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no,n,i,ii)

      i = i + El%inDpvt%r(i_Elec) - (crows(n)-sN) - 1
      Gf => z(i:ii)

#ifdef TBT_T_G_GAMMA_OLD
       
      ! Number of columns that we want to do product of
      ii = 1
      do i = 1 , no
        T = T - aimag( zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) ) ! G \Gamma
        ii = ii + sN
      end do
      ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
      !    Tr[(G \Gamma)^\dagger]
      ! Hence the below calculation shouldn't be necessary
      ii = (i_Elec - 1) * no + 1
      do i = 1 , nb
        T = T + aimag( zdotc(no,Gf(i),sN,El%Gamma(ii),1) )! G^\dagger \Gamma
        ii = ii + no
      end do

#else
      
      ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
      !    Tr[(G \Gamma)^\dagger]
      ! Hence we only calculate one of the contributions and double
      ! it after
      ! Indeed we actually need to calculate:
      !    i Tr[G \Gamma] and since Gamma is not having the i factor
      ! we may take the negative real part.
      ii = 1
      do i = 1 , no
        T = T - aimag( zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) )! G \Gamma
        ii = ii + sN
      end do
       
#endif

      i_Elec = i_Elec + nb

    end do

    ! Now we have:
    !   T = Tr[G \Gamma - G^\dagger \Gamma]
#ifndef TBT_T_G_GAMMA_OLD
    T = T * 2._dp
#endif

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',2)
#endif
    
  end subroutine GF_Gamma

  subroutine GF_T(Gfcol,El,T_Gf,T_self,nzwork,zwork)

    use intrinsic_missing, only: TRACE
    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(Elec), intent(inout) :: El
    real(dp), intent(out) :: T_Gf, T_self
    integer, intent(in) :: nzwork
    complex(dp), intent(out) :: zwork(nzwork)

    complex(dp), pointer :: z(:)

    integer :: no, np
    integer :: i, ii, i_Elec
    integer, pointer :: crows(:)
    integer :: sN, n
    integer :: nb
    ! BLAS routines
    complex(dp), external :: zdotu

#ifdef TBTRANS_TIMING
    call timer('Gf-T',1)
#endif

    no = El%inDpvt%n
    np = parts(Gfcol)
    crows => cum_rows(Gfcol)

#ifndef TS_NOCHECKS
    if ( no**2 > nzwork ) call die('GF_T: no**2 < nzwork')

    ! First we check that we can use the first elements
    ! of Gfcol as temporary storage
    call TriMat_Bias_idxs(Gfcol,no,1,i,ii)
    if ( i < no ** 2 ) then
      write(*,'(a)') 'Remove TBT.T.Gf from your fdf file. &
          &It is not possible in your current setup.'
      call die('GF_T: Size of temporary array not possible.')
    end if

#endif

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    i_Elec = 1
    do while ( i_Elec <= no ) 

      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no,n,i,ii)

      i = i + El%inDpvt%r(i_Elec) - (crows(n)-sN) - 1

      ! Calculate the G \Gamma
#ifdef USE_GEMM3M
      call zgemm3m( &
#else
      call zgemm( &
#endif
           'N','T',nb,no,no, z1, z(i),sN, &
           El%Gamma(1), no, z0, zwork(i_Elec), no)
       
      i_Elec = i_Elec + nb

    end do

    ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
    !    Tr[(G \Gamma)^\dagger]
    ! Now we have:
    !    = G \Gamma
    T_Gf = - aimag( TRACE(no,zwork) ) * 2._dp


    ! Now we need to correct for the current electrode
    
    ! Now we can calculate the spectral function for this
    ! electrode where it lives
    i_Elec = 1
    do while ( i_Elec <= no ) 
       
      ! We start by creating a region of consecutive memory.
      call consecutive_index(Gfcol,El,i_Elec,n,nb)
      sN = nrows_g(Gfcol,n)

      ! get placement of the diagonal block in the column
      call TriMat_Bias_idxs(Gfcol,no,n,i,ii)

      i = i + El%inDpvt%r(i_Elec) - (crows(n)-sN) - 1

      ii = ( i_Elec-1 ) * no + 1
      ! Calculate the G \Gamma G^\dagger
#ifdef USE_GEMM3M
      call zgemm3m( &
#else
      call zgemm( &
#endif
           'N','C',no,nb,no, z1, zwork(1),no, &
           z(i), sN, z0, z(ii), no)
       
      i_Elec = i_Elec + nb

    end do

    ! The remaining calculation is very easy as we simply
    ! need to do the sum of the trace
    T_self = - aimag( zdotu(no*no,z(1),1,El%Gamma(1),1) )

#ifdef TBTRANS_TIMING
    call timer('Gf-T',2)
#endif
    
  end subroutine GF_T

  subroutine GF_T_solve(N_Elec,T,has_all)
    integer, intent(in) :: N_Elec
    real(dp), intent(inout) :: T(N_Elec+1,N_Elec)
    ! Whether or not we have calculated all transmissions
    logical, intent(in) :: has_all

    real(dp) :: TT(3)

    select case ( N_Elec ) 
    case ( 1 )

      ! For one electrodes, we simply return immediately
      return

    case ( 2 )

      ! The simple case is when we have 2 electrodes
      
      T(2,1) = T(N_Elec+1,1) - T(1,1)
      if ( has_all ) then
        T(1,2) = T(N_Elec+1,2) - T(2,2)
      else
        T(1,2) = T(2,1)
      end if

      return

    case ( 3 )
      
      if ( .not. has_all ) then
        call die('GF_T_solve: Can not separate transmissions (need all bulk).')
      end if

    case default

      call die('Calculating transmission from underdetermined &
          &system is not allowed. Remove TBT.T.Gf.')
      
    end select

    ! RHS
    TT(1) = T(N_Elec+1,1) - T(1,1)
    TT(2) = T(N_Elec+1,2) - T(2,2)
    TT(3) = T(N_Elec+1,3) - T(3,3)

    ! Calculate them exactly
    ! We do not need LAPACK here as this can
    ! only be definitely solved for N_Elec == 3
    T(2,1) = (TT(1) - TT(3) + TT(2)) * 0.5_dp
    T(1,2) = T(2,1)
    T(3,1) = TT(1) - T(2,1)
    T(1,3) = T(3,1)
    T(3,2) = TT(3) - T(3,1)
    T(2,3) = T(3,2)

  end subroutine GF_T_solve

  subroutine consecutive_index(Tri,El,current,p,n)
    type(zTriMat), intent(inout) :: Tri
    type(Elec), intent(in) :: El
    integer, intent(in) :: current
    integer, intent(out) :: p, n

    ! Local variables
    integer :: idx_Elec, i, sIdx, eIdx
    integer, pointer :: crows(:)
    
    idx_Elec = El%inDpvt%r(current)
    p = which_part(Tri,idx_Elec)
    crows => cum_rows(Tri)
    eIdx = crows(p)
    sIdx = eIdx - nrows_g(Tri,p) + 1

    n = 1
    do while ( current + n <= El%inDpvt%n )
      i = El%inDpvt%r(current+n)
      ! In case it is not consecutive
      if ( i - idx_Elec /= n ) exit
      ! In case the block changes, then
      ! we cut the block size here.
      if ( i < sIdx .or. eIdx < i ) exit
      n = n + 1
    end do

  end subroutine consecutive_index


#ifdef NCDF_4
  subroutine orb_current(sp,H,S,sc_off,k,ph,cE,A_tri,r,orb_J,pvt)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use geom_helper,       only : UCORB
    use sorted_search_m, only: ssearch_t, ssearch_init, ssearch_find

    use m_ts_cctype, only: ts_c_idx

    type(Sparsity), intent(inout) :: sp
    ! We require that the input Hamiltonian is Hermitian
    real(dp), intent(in) :: H(:), S(:), sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(ts_c_idx) :: cE
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of orb_J
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: orb_J
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: i_sp
    integer, pointer :: i_ncol(:), i_ptr(:), i_col(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    complex(dp), pointer :: A(:)
    complex(dp) :: Hi
    real(dp), pointer :: J(:)
    real(dp) :: E
    integer :: no_u, iu, io, ind, iind, ju, jo
    type(ssearch_t) :: ss

#ifdef TBTRANS_TIMING
    call timer('orb-current',1)
#endif

    ! Retrieve energy
    E = real(cE%e,dp)

    call attach(sp, nrows_g=no_u, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    i_sp => spar(orb_J)
    J    => val (orb_J)
    call attach(i_sp, n_col=i_ncol, list_ptr=i_ptr, list_col=i_col)

    ! Create the phases
    ! We are using the symmetric H(j, i) = H(i, j) relation.
    ! So since we are taking the complex part on the first entry we retrieve the H(j,i) (in k-space)
    ! component.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, + &
          k(1) * sc_off(1,io) + &
          k(2) * sc_off(2,io) + &
          k(3) * sc_off(3,io)))
    end do

    A => val(A_tri)

    ! we need this in case the device region gets enlarged due to dH
    J(:) = 0._dp

!$OMP parallel do default(shared), private(iu,io,ju,jo,iind,ind,Hi,ss)
    do iu = 1, r%n
      io = r%r(iu)

#ifndef TS_NOCHECKS
      if ( i_ncol(io) == 0 ) call die('orb_current: J has zero columns &
          &for at least one row')
#endif

      ! Get lookup columns for the orbital current
      call ssearch_init(ss, i_col(i_ptr(io)+1:i_ptr(io)+i_ncol(io)))

      ! Loop on Hamiltonian entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

        ! Check if the orbital exists in the region
        iind = i_ptr(io) + ssearch_find(ss, l_col(ind))
        ! if zero the element does not exist
        ! This is the case on the elements connecting out
        ! of the device region
        if ( iind <= i_ptr(io) ) cycle

        ! H_ind == H_ij

        ! We may take the conjugate later as E is a real quantity
        Hi = (H(ind) - E * S(ind)) * ph( (l_col(ind)-1)/no_u )

        ! J(iind) = J(io,jo)
        jo = ucorb(l_col(ind),no_u)

        ! Get spectral function indices
        ju = pvt%r(jo) ! pivoted orbital index in tri-diagonal matrix
        jo = index(A_tri,iu,ju) ! A_ij
        ju = index(A_tri,ju,iu) ! A_ji

        ! We skip the pre-factors as the units are "never" used

        ! Jij                Hji    * Aij    Hij * Aji
        J(iind) = aimag( dconjg(Hi) * A(jo) - Hi * A(ju) )

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('orb-current',2)
#endif

  end subroutine orb_current
  
  subroutine orb_current_add_dH(dH_1D,sc_off,k,ph,A_tri,r,orb_J,pvt)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND
    use geom_helper,       only : UCORB

    type(zSpData1D), intent(in) :: dH_1D
    real(dp), intent(in) :: sc_off(:,:), k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of orb_J
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: orb_J
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: dH(:)
    type(Sparsity), pointer :: i_sp
    integer, pointer :: i_ncol(:), i_ptr(:), i_col(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:), col(:)

    complex(dp) :: p
    complex(dp), pointer :: A(:)
    real(dp), pointer :: J(:)
    integer :: no_u, iu, io, ind, iind, ju, jo, jj

#ifdef TBTRANS_TIMING
    call timer('orb-current-dH',1)
#endif

    ! Retrieve dH
    sp => spar(dH_1D)
    dH => val (dH_1D)
    call attach(sp, nrows_g=no_u, &
        n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    i_sp => spar(orb_J)
    J    => val (orb_J)
    call attach(i_sp, n_col=i_ncol, list_ptr=i_ptr, list_col=i_col)

    ! Create the phases
    ! We are using the explicit H(j, i) and thus the phases are consistent with +
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, + &
          k(1) * sc_off(1,io) + &
          k(2) * sc_off(2,io) + &
          k(3) * sc_off(3,io)))
    end do

    A => val(A_tri)

!$OMP parallel do default(shared), &
!$OMP&private(iu,io,iind,jo,ju,ind,col,jj,p)
    do iu = 1, r%n
      io = r%r(iu)

      ! Loop on the orbital current indices
      do iind = i_ptr(io) + 1, i_ptr(io) + i_ncol(io)

        ! Here we will calculate the orbital current from dH
        ! onto orbital:
        !  J(iind) == J(io, jo)

        ! Get jo orbital
        jo = ucorb(i_col(iind), no_u)
        ju = pvt%r(jo) ! pivoted orbital index in tri-diagonal matrix

        
        ! Check if the jo, io orbital exists in dH
        if ( l_ncol(jo) < 1 ) then
          ind = -1
        else
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+l_ncol(jo))
          ! Get transpose element
          jj = TO(i_col(iind)) + io
          ind = l_ptr(jo) + SFIND(col, jj)
        end if

        if ( ind > l_ptr(jo) ) then

          ! Add orbital current from ji
          p = ph( (l_col(ind)-1)/no_u )

          ! Check for the Hamiltonian element H_ji
          jj = index(A_tri,iu,ju) ! A_ij

          ! Jij                      Aij   * Hji
          J(iind) = J(iind) + aimag( A(jj) * dH(ind) * p )

        end if

        ! Check if the io, jo orbital exists in dH
        if ( l_ncol(io) < 1 ) then
          ind = -1
        else
          col => l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io))
          ind = l_ptr(io) + SFIND(col, i_col(iind))
        end if

        if ( ind > l_ptr(io) ) then

          ! Add orbital current from ij
          p = ph( (l_col(ind)-1)/no_u )

          ! Check for the Hamiltonian element H_ij
          jj = index(A_tri,ju,iu) ! A_ji

          ! Jij -=                   Aji   * Hij
          J(iind) = J(iind) - aimag( A(jj) * dH(ind) * p )

        end if

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('orb-current-dH',2)
#endif

  contains

    function TO(io) result(jo)
      integer, intent(in) :: io
      integer :: jo, isc, i

      ! Get the current supercell index
      isc = (io-1)/no_u + 1
      
      do i = 1, size(sc_off, dim=2)

        ! We have to check for the opposite super-cell to get the
        ! transpose element.
        ! 0.001 Bohr seems like a more than accurate difference for
        ! unit-cells.
        if ( all( abs(sc_off(:,i) + sc_off(:, isc)) < 0.001_dp) ) then
          jo = (i - 1) * no_u
          return
        end if

      end do

      jo = 0
      call die('orb_current_add_dH: could not find transpose supercell index')
      
    end function TO
    
  end subroutine orb_current_add_dH


  subroutine GF_DM(sc_off,k,ph,Gfd_tri,Gfo_tri,r,pvt,spDM)

    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB

    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: Gfd_tri, Gfo_tri
    ! The region that specifies the size of spDM
    type(tRgn), intent(in) :: r
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt
    type(dSpData1D), intent(inout) :: spDM

    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    type(Sparsity), pointer :: sp
    real(dp), pointer :: DM(:)
    complex(dp), pointer :: Gfd(:), Gfo(:)
    complex(dp) :: GfGfd

    integer :: no_u, iu, io, ind, ju

#ifdef TBTRANS_TIMING
    call timer('Gf-DM',1)
#endif

    sp => spar(spDM)
    DM => val(spDM)
    call attach(sp, nrows_g=no_u, n_col=ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    ! Since we have to do Gf.exp(ikR) we simply
    ! create exp(-ikR) for the supercell connections.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, - &
          k(1) * sc_off(1,io) - &
          k(2) * sc_off(2,io) - &
          k(3) * sc_off(3,io))) / (2._dp * Pi)
    end do

    Gfd => val(Gfd_tri)
    Gfo => val(Gfo_tri)

    ! we need this in case the device region gets enlarged due to dH
    DM(:) = 0._dp

!$OMP parallel do default(shared), private(iu,io,ind,ju,GfGfd)
    do iu = 1, r%n
      io = r%r(iu)

#ifndef TS_NOCHECKS
      if ( ncol(io) == 0 ) call die('Gf_DM: DM has zero columns &
          &for at least one row')
#endif

      ! Loop on DM entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

        ju = pvt%r(ucorb(l_col(ind), no_u))
        call calc_GfGfd(iu, ju, GfGfd)
        DM(ind) = - aimag( GfGfd * ph((l_col(ind) - 1) / no_u) )

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('Gf-DM',2)
#endif
    
  contains
    
    subroutine calc_GfGfd(br, bc, G)
      integer, intent(in) :: br, bc
      complex(dp), intent(inout) :: G
      integer :: p_r, i_r, p_c, i_c, i

      call part_index(Gfo_tri, br, p_r, i_r)
      call part_index(Gfo_tri, bc, p_c, i_c)
      
      if ( p_r == p_c ) then
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfd(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        G = G - conjg(Gfd(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      else
        i = index_block(Gfo_tri, p_r, p_c)
        G = Gfo(i + i_r + (i_c-1) * Gfo_tri%data%tri_nrows(p_r))
        i = index_block(Gfo_tri, p_c, p_r)
        G = G - conjg(Gfo(i + i_c + (i_r-1) * Gfo_tri%data%tri_nrows(p_c)))
      end if

    end subroutine calc_GfGfd
    
  end subroutine GF_DM

  subroutine A_DM(sc_off,k,ph,A_tri,r,pvt,spDM)

    use class_Sparsity
    use class_dSpData1D
    use geom_helper,       only : UCORB

    real(dp), intent(in) :: sc_off(:,:)
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: ph(0:)
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of spDM
    type(tRgn), intent(in) :: r
    ! The pivoting region that transfers r%r(iu) to io
    type(tRgn), intent(in) :: pvt
    type(dSpData1D), intent(inout) :: spDM

    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: A(:)
    real(dp), pointer :: DM(:)
    integer :: no_u, iu, io, ind, ju

#ifdef TBTRANS_TIMING
    call timer('A-DM',1)
#endif

    sp => spar(spDM)
    DM => val(spDM)
    call attach(sp, nrows_g=no_u, n_col=ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    ! Since we have to do Gf.exp(ikR) we simply
    ! create exp(-ikR) for the supercell connections.
    do io = 1 , size(sc_off, dim=2)
      ph(io-1) = cdexp(dcmplx(0._dp, - &
          k(1) * sc_off(1,io) - &
          k(2) * sc_off(2,io) - &
          k(3) * sc_off(3,io))) / (2._dp * Pi)
    end do

    A => val(A_tri)

    ! we need this in case the device region gets enlarged due to dH
    DM(:) = 0._dp

!$OMP parallel do default(shared), private(iu,io,ind,ju)
    do iu = 1, r%n
      io = r%r(iu)

#ifndef TS_NOCHECKS
      if ( ncol(io) == 0 ) call die('A_DM: DM has zero columns &
          &for at least one row')
#endif

      ! Loop on DM entries here...
      do ind = l_ptr(io) + 1 , l_ptr(io) + ncol(io)

        ju = pvt%r(ucorb(l_col(ind), no_u))
        ju = index(A_tri, iu, ju)
        DM(ind) = real(A(ju) * ph((l_col(ind) - 1) / no_u), dp)

      end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('A-DM',2)
#endif

  end subroutine A_DM

#endif

  subroutine insert_Self_energy(n1,n2,M,r,El,off1,off2)

    ! The sizes of the matrix
    integer, intent(in) :: n1, n2
    complex(dp), intent(inout) :: M(n1,n2)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! Electrodes...
    type(Elec), intent(inout) :: El
    ! The offsets of the matrix
    integer, intent(in) :: off1, off2

    ! local variables
    integer :: j, je, i, ie, no, idx

    idx = El%idx_o - 1
    no = TotUsedOrbs(El)

    ! We are dealing with the intrinsic electrode
    ! self energy
    ! Here we have two options,
    ! Bulk) We are dealing with a bulk electrode
    ! not bulk) A non-bulk electrode

    if ( El%Bulk ) then
!$OMP do private(j,je,i,ie)
      do j = 1 , n2
        je = r%r(off2+j) - idx
        if ( 1 <= je .and. je <= no ) then
          je = (je - 1) * no
          do i = 1 , n1
            ie = r%r(off1+i) - idx
            if ( ie < 1 ) cycle
            if ( no < ie ) cycle
             
            M(i,j) = El%Sigma(je + ie)
            
          end do
        end if
      end do
!$OMP end do
    else
!$OMP do private(j,je,i,ie)
      do j = 1 , n2
        je = r%r(off2+j) - idx
        if ( 1 <= je .and. je <= no ) then
          je = (je - 1) * no
          do i = 1 , n1
            ie = r%r(off1+i) - idx
            if ( ie < 1 ) cycle
            if ( no < ie ) cycle

            M(i,j) = M(i,j) - El%Sigma(je + ie)

          end do
        end if
      end do
!$OMP end do
    end if

  end subroutine insert_Self_energy


  subroutine insert_Self_energy_Dev(Gfinv_tri,Gfinv,r,El)

    type(zTriMat), intent(inout) :: GFinv_tri
    complex(dp), intent(inout) :: Gfinv(:)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    type(Elec), intent(in) :: El

    ! local variables
    integer :: j, je, i, ii, idx, no

    no = El%o_inD%n

    ! A down-folded self-energy, this
    ! is always considered to be "non-bulk" as
    ! we have it downfolded.

!$OMP do private(j,ii,je,i,idx)
    do j = 1 , no
      ii = (j-1)*no
      ! grab the index in the full tri-diagonal matrix
      je = El%inDpvt%r(j)
      do i = 1 , no
        
        idx = index(GFinv_tri,El%inDpvt%r(i),je)
        
        Gfinv(idx) = Gfinv(idx) - El%Sigma(ii+i)
        
      end do
    end do
!$OMP end do

  end subroutine insert_Self_energy_Dev

end module m_tbt_tri_scat
