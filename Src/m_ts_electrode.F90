! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_ts_electrode
!
! Routines that are used for Electrodes GFs calculations
! Heavily updated by Nick Papior Andersen, 2012
!

  use precision, only : dp

  implicit none

  public :: create_Green
  public :: init_Electrode_HS
  public :: calc_next_GS_Elec

  private

  ! BLAS parameters
  complex(dp), parameter :: z_1  = dcmplx(1._dp,0._dp)
  complex(dp), parameter :: z_m1 = dcmplx(-1._dp,0._dp)
  complex(dp), parameter :: z_0  = dcmplx(0._dp,0._dp)

  interface set_HS_transfer
     module procedure set_HS_Transfer_1d
     module procedure set_HS_Transfer_2d
  end interface set_HS_transfer

contains


  ! Calculates the surface Green function for the electrodes
  ! Handles both the left and right one
  ! this is the Sancho, Sancho and Rubio algorithm
  subroutine SSR_sGreen_DOS(no,ZE,H00,S00,H01,S01,accu, GS, &
       DOS, T, &
       nwork, zwork, &
       iterations, final_invert)
       
! ***************** INPUT **********************************************
! integer     no      : Number of orbitals in the electrode
! complex(dp) ZE      : The energy of the Green function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell (discarding T-direction)
! complex(dp) S00     : Overlap matrix within the first unit cell (discarding T-direction)
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell (in T-direction)
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell (in T-direction)
! real(dp) accu       : Define the accuracy needed for convergence.
! ***************** OUTPUT *********************************************
! complex(dp) GS      : Surface Green function of the electrode
! real(dp) DOS        : DOS of bulk electrode (additive)
! real(dp) T          : transmission of bulk electrode (additive)
! **********************************************************************
    use m_pivot_array, only : ipiv
    use m_mat_invert
    use precision, only: dp
    use units, only : Pi
    use intrinsic_missing, only: transpose

! ***********************
! * INPUT variables     *
! ***********************
    integer,     intent(in) :: no
    complex(dp), intent(in) :: ZE 
    complex(dp), intent(in) :: H00(no*no),S00(no*no)
    complex(dp), intent(in) :: H01(no*no),S01(no*no)
    real(dp), intent(in) :: accu

    integer, intent(in) :: nwork

    logical, intent(in), optional :: final_invert

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(inout), target :: GS(no*no)
    complex(dp), intent(inout), target :: zwork(nwork)
    real(dp), intent(inout) :: DOS(no)
    real(dp), intent(inout) :: T

    integer, intent(inout), optional :: iterations

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: nom1, no2, nosq
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2
    logical :: as_first

    real(dp) :: ro
    complex(dp) :: zij, zji

    complex(dp), dimension(:), pointer :: rh,rh1,w,alpha,beta,gb
    complex(dp), dimension(:), pointer :: gsL,gsR
    
    complex(dp), external :: zdotu, zdotc

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE SSR_sGreen_DOS' )
#endif

    ! Initialize counter
    if ( present(iterations) ) iterations = 0

!    call timer('ts_GS',1)

    nom1 = no - 1
    no2  = no * 2
    nosq = no * no

    if ( nwork < 9 * nosq ) call die('SSR_sGreen_DOS: &
         &Not enough work space')
    i = 0
    rh  => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    rh1 => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    alpha => zwork(i+1:i+nosq) 
    i = i + nosq
    beta => zwork(i+1:i+nosq) 
    i = i + nosq
    w => zwork(i+1:i+nosq)
    i = i + nosq
    GB => zwork(i+1:i+nosq) 
    i = i + nosq

    gsL => zwork(i+1:i+nosq) 
    gsR => GS

!$OMP parallel default(shared), private(i,j,ic,ic2)

! gb    =   Z*S00-H00
!$OMP do
    do i = 1 , nosq
       GB(i)    = ZE * S00(i) - H00(i)
    end do
!$OMP end do nowait
! alpha = -(Z*S01-H01)
!$OMP do
    do i = 1 , nosq
       alpha(i) = H01(i) - ZE * S01(i)
    end do
!$OMP end do nowait
    ! zero arrays
    ! We do not start with H00 as we then
    ! will still need to re-adjust when
    ! calculating the scattering matrices.
!$OMP do
    do i = 1 , nosq
       gsL(i) = z_0
    end do
!$OMP end do nowait
!$OMP do
    do i = 1 , nosq
       gsR(i) = z_0
    end do
!$OMP end do nowait

! beta = -(Z*S10-H10)
!$OMP do
    do j = 1 , no
       ic = no * (j-1)
       do i = 1 , no
          ic2 = no*(i-1) + j
          beta(ic+i) = dconjg(H01(ic2)) - ZE * dconjg(S01(ic2))
       end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! Initialize loop
    ro = accu + 1._dp
    as_first = .false.
    do while ( ro > accu ) 

       ! Increment iterations
       if ( present(iterations) ) &
            iterations = iterations + 1

! rh = -(Z*S01-H01) ,j<no
! rh = -(Z*S10-H10) ,j>no
!$OMP parallel default(shared), private(i)
       
!$OMP do
       do i = 1, nosq
          rh(i)      = alpha(i)
       end do
!$OMP end do nowait
!$OMP do
       do i = 1, nosq
          rh(nosq+i) = beta(i)
       end do
!$OMP end do nowait
!$OMP do
       do i = 1, nosq
          ! w = Z*S00-H00
          w(i) = GB(i)
       end do
!$OMP end do nowait
       
!$OMP end parallel

! rh =  rh1^(-1)*rh
! rh =  t0
       call zgesv(no, no2, w, no, ipiv, rh, no, ierr)
       if ( ierr /= 0 ) then
          write(*,*) 'ERROR: SSR_sGreen_DOS 1 MATRIX INVERSION FAILED'
          write(*,*) 'ERROR: LAPACK INFO = ',ierr
       end if

       ! switch pointers instead of copying elements
       call switch_alpha_beta_rh1(as_first)

! alpha = -(Z*S01-H01)*t0
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_1,rh1(1),no,rh(1),no,z_0,alpha,no)
! beta  = -(Z*S10-H10)*t0 ??
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_1,rh1(nosq+1),no,rh(nosq+1),no,z_0,beta,no)

! ba    = (Z*S10-H10)*t0b
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_0,w,no)
!$OMP parallel default(shared), private(i)

!$OMP do
       do i = 1 , nosq
          GB(i)  = GB(i) + w(i)
       end do
!$OMP end do nowait
!$OMP do 
       do i = 1 , nosq
          gsL(i) = gsL(i) + w(i)
       end do
!$OMP end do nowait
       
!$OMP end parallel

! ab    = (Z*S01-H01)*t0
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_m1,rh1(1),no,rh(nosq+1),no,z_0,w,no)

       ro = -1._dp
!$OMP parallel default(shared), private(i)
!$OMP do
       do i = 1 , nosq
          GB(i)  = GB(i) + w(i)
       end do
!$OMP end do nowait
!$OMP do
       do i = 1 , nosq
          gsR(i) = gsR(i) + w(i)
       end do
!$OMP end do nowait
!$OMP do reduction(max:ro)
       do i = 1 , nosq
          ! update the criteria
          ro = max(ro,abs(w(i)))
       end do
!$OMP end do nowait
!$OMP end parallel
       
    end do

    ! *** Initiate DOS and bulk transmission calculation...
    
    ! Invert to obtain the bulk Green function
    call mat_invert(GB,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: SSR_sGreen_DOS GB MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    ! Calculate scattering matrices for left-right self-energy
    ! and correct the self-energy with the bulk Hamiltonian
    ! to get the correct self-energy
!$OMP parallel do default(shared), private(i,j,ic,ic2,zij,zji)
    do j = 1 , no
       do i = 1 , j - 1
          ic  = (j-1)*no+i
          ic2 = (i-1)*no+j
          
          ! Calculate bulk Hamiltonian and overlap
          zij = ZE*S00(ic ) - H00(ic )
          zji = ZE*S00(ic2) - H00(ic2)
          
          ! left scattering states
          alpha(ic)  = gsL(ic ) - dconjg(gsL(ic2))
          alpha(ic2) = gsL(ic2) - dconjg(gsL(ic ))
          
          ! Correct for bulk self-energy + bulk Hamiltonian
          gsL(ic ) = zij + gsL(ic )
          gsL(ic2) = zji + gsL(ic2)

          ! right scattering states (transposed)
          beta(ic2) = gsR(ic ) - dconjg(gsR(ic2))
          beta(ic ) = gsR(ic2) - dconjg(gsR(ic ))

          ! Correct for bulk self-energy + bulk Hamiltonian
          gsR(ic ) = zij + gsR(ic )
          gsR(ic2) = zji + gsR(ic2)
          
       end do
       ic = (j-1)*no+j
       
       ! Add bulk Hamiltonian and overlap
       zij = ZE*S00(ic) - H00(ic)
       
       ! left scattering state
       alpha(ic) = gsL(ic) - dconjg(gsL(ic))
       
       ! Correct for bulk self-energy + bulk Hamiltonian
       gsL(ic) = zij + gsL(ic)
       
       ! right scattering state (transposed)
       beta(ic) = gsR(ic) - dconjg(gsR(ic))
       
       ! Correct for bulk self-energy + bulk Hamiltonian
       gsR(ic) = zij + gsR(ic)
       
    end do
!$OMP end parallel do

    
    ! gsR and GS are the same array
    ! Hence it we should not return the
    ! inverted matrix, copy gsR to other
    ! work-array
    if ( present(final_invert) ) then
       if ( .not. final_invert ) then

          ! we return a non-inverted matrix
          ! hence prohibit the inversion of the matrix
          ! by moving data to another work-array
!$OMP parallel do default(shared), private(i)
          do i = 1 , nosq
             rh1(i) = gsR(i)
          end do
!$OMP end parallel do
          gsR => rh1(1:nosq)

       end if
    end if
    

    ! Invert to get the Surface Green function (Left)
    call mat_invert(gsL,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: SSR_sGreen_DOS GSL MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    ! Invert to get the Surface Green function (Right)
    call mat_invert(gsR,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: SSR_sGreen_DOS GSR MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if
    

    ! Calculate bulk transmision
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',no,no,no,z_1,GB   ,no,alpha,no,z_0,w    ,no)
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','C',no,no,no,z_1,w   ,no,GB    ,no,z_0,alpha,no)

    ! Calculate bulk-transmission (same as matrix product + trace)
    T = T - real(zdotu(nosq,alpha,1,beta,1),dp)


    ! We now calculate the density of states...
!$OMP parallel default(shared), private(i) 
!$OMP do
    do i = 1 , nosq
       alpha(i) = H01(i) -        ZE  * S01(i)
    end do
!$OMP end do nowait
!$OMP do
    do i = 1 , nosq
       ! notice, we utilize the relation (H10-z*S10) = (H01-conjg(z)*S01)^H
       beta(i)  = H01(i) - dconjg(ZE) * S01(i)
    end do
!$OMP end do nowait
!$OMP end parallel

    ! Transpose GB to make the latter dot product easier
    call transpose(no,GB)
    i = 1
    j = nosq + 1 
    ! DOS = diag{ G_b * S00 + 
    !             G_l * (H01 - E * S01 ) * G_b * S10 +
    !             G_r * (H10 - E * S10 ) * G_b * S01   }
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',no,no,no,z_1,gsL  ,no,alpha,no,z_0,w    ,no)
    ! Note that GB is transposed, (so transpose back here)
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','T',no,no,no,z_1,w    ,no,GB   ,no,z_0,rh(i),no)
    ! Note that GB is transposed, (so transpose back here)
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'C','T',no,no,no,z_1,beta ,no,GB   ,no,z_0,rh(j),no)
    ! Transpose the matrix product to align the following
    ! dot product (reduces the need to calculate the total
    ! matrix before we take the trace)
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'T','T',no,no,no,z_1,rh(j),no,gsR  ,no,z_0,w    ,no)

    ! We resort to use the zdotu/zdotc to not calculate
    ! unnessary elements (for large systems this should
    ! speed it up!) For small systems this has little
    ! overhead
    i = 1
    do j = 1 , no
       
       ! Calculate for the bulk
       ro =    aimag(zdotu(no,GB(i),1,S00(i),1))
       ! For the left self-energy
       ro = ro+aimag(zdotc(no,S01(i),1,rh(i),1))
       ! For the right self-energy
       ro = ro+aimag(zdotu(no,w(i),1,S01(i),1))

       ! Calculate the total DOS
       DOS(j) = DOS(j) - ro / Pi
       
       i = i + no
       
    end do

!    call timer('ts_GS',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS SSR_sGreen_DOS' )
#endif

  contains

    ! We supply a routine to switch the pointer position of alpha,beta / rh1
    subroutine switch_alpha_beta_rh1(as_first)
      logical, intent(inout) :: as_first
      integer :: i 
      ! start
      i = 2 * nosq

      if ( as_first ) then
         rh1 => zwork(i+1:i+2*nosq) 
         i = i + 2*nosq
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
      else
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
         i = i + nosq
         rh1 => zwork(i+1:i+2*nosq) 
      end if
      as_first = .not. as_first

    end subroutine switch_alpha_beta_rh1

  end subroutine SSR_sGreen_DOS

  ! Calculates the surface Green function for the electrodes
  ! Handles both the left and right one
  ! this is the Sancho, Sancho and Rubio algorithm
  subroutine SSR_sGreen_NoDOS(no,ZE,H00,S00,H01,S01,accu, GS, &
       nwork, zwork, &
       iterations, final_invert)
       
! ***************** INPUT **********************************************
! integer     no      : Number of orbitals in the electrode
! complex(dp) ZE      : The energy of the Green function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell (discarding T-direction)
! complex(dp) S00     : Overlap matrix within the first unit cell (discarding T-direction)
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell (in T-direction)
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell (in T-direction)
! ***************** OUTPUT *********************************************
! complex(dp) GS      : Surface Green function of the electrode
! **********************************************************************
    use m_pivot_array, only : ipiv
    use m_mat_invert
    use precision, only: dp

! ***********************
! * INPUT variables     *
! ***********************
    integer,     intent(in) :: no
    complex(dp), intent(in) :: ZE 
    complex(dp), intent(in) :: H00(no*no),S00(no*no)
    complex(dp), intent(in) :: H01(no*no),S01(no*no)
    real(dp), intent(in) :: accu

    integer,     intent(in) :: nwork

    logical, intent(in), optional :: final_invert

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(inout), target :: GS(no*no)
    complex(dp), intent(inout), target :: zwork(nwork)

    integer, intent(inout), optional :: iterations

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: nom1, no2, nosq
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2
    logical :: as_first

    real(dp) :: ro

    complex(dp), dimension(:), pointer :: rh,rh1,w,alpha,beta,GB

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE SSR_sGreen_NoDOS' )
#endif

    ! Initialize counter
    if ( present(iterations) ) iterations = 0

!    call timer('ts_GS',1)

    nom1 = no - 1
    no2  = 2 * no
    nosq = no * no

    if ( nwork < 8 * nosq ) call die('SSR_sGreen_NoDOS: &
         &Not enough work space')
    i = 0
    rh  => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    rh1 => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    alpha => zwork(i+1:i+nosq) 
    i = i + nosq
    beta => zwork(i+1:i+nosq) 
    i = i + nosq
    w => zwork(i+1:i+nosq)
    i = i + nosq
    GB => zwork(i+1:i+nosq) 

!$OMP parallel default(shared), private(i,j,ic,ic2)

! gb    =   Z*S00-H00
!$OMP do 
    do i = 1 , nosq
       GB(i)    = ZE * S00(i) - H00(i)
    end do
!$OMP end do nowait
! gs  = gb
!$OMP do 
    do i = 1 , nosq
       GS(i)    = GB(i)
    end do
!$OMP end do nowait
! alpha = -(Z*S01-H01)
!$OMP do 
    do i = 1 , nosq
       alpha(i) = H01(i) - ZE * S01(i)
    end do
!$OMP end do nowait
! beta = -(Z*S10-H10)
!$OMP do
    do j = 1 , no
       ic = no * (j-1) + 1
       do i = 0 , nom1
          ic2 = j + no*i
          beta(ic+i) = dconjg(H01(ic2)) - ZE * dconjg(S01(ic2))
       end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! Initialize loop
    ro = accu + 1._dp
    as_first = .false.
    do while ( ro > accu ) 

       ! Increment iterations
       if ( present(iterations) ) &
            iterations = iterations + 1


!$OMP parallel default(shared), private(i)
! rh = -(Z*S01-H01) ,j<no
!$OMP do
       do i = 1, nosq
          rh(i)      = alpha(i)
       end do
!$OMP end do nowait
! rh = -(Z*S10-H10) ,j>no
!$OMP do
       do i = 1, nosq
          rh(nosq+i) = beta(i)
       end do
!$OMP end do nowait
!$OMP do
       do i = 1, nosq
          ! w = Z*S00-H00
          w(i)       = GB(i)
       end do
!$OMP end do nowait
!$OMP end parallel

! rh =  rh1^(-1)*rh
! rh =  t0
       call zgesv(no, no2, w, no, ipiv, rh, no, ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'ERROR: SSR_sGreen_NoDOS 1 MATRIX INVERSION FAILED'
          write(*,*) 'ERROR: LAPACK INFO = ',ierr
       end if

       ! switch pointers instead of copying elements
       call switch_alpha_beta_rh1(as_first)

! alpha = -(Z*S01-H01)*t0
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_1,rh1(1),no,rh(1),no,z_0,alpha,no)
! beta  = -(Z*S10-H10)*t0 ??
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_1,rh1(nosq+1),no,rh(nosq+1),no,z_0,beta,no)

! gb = gb + [ba    = (Z*S10-H10)*t0b]
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_1,GB,no)

! ab    = (Z*S01-H01)*t0
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',no,no,no,z_m1,rh1(1),no,rh(nosq+1),no,z_0,w,no)

       ro = -1._dp
!$OMP parallel default(shared), private(i)
!$OMP do
       do i = 1 , nosq
          GB(i) = GB(i) + w(i)
       end do
!$OMP end do nowait
!$OMP do
       do i = 1 , nosq
          GS(i) = GS(i) + w(i)
       end do
!$OMP end do nowait
!$OMP do reduction(max:ro)
       do i = 1 , nosq
          ! update the criteria
          ro = max(ro,abs(w(i)))
       end do
!$OMP end do nowait
!$OMP end parallel

    end do

    if ( present(final_invert) ) then
       if ( final_invert ) then
          ! Invert to get the Surface Green function
          call mat_invert(GS,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
       end if
    end if

    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: SSR_sGreen_NoDOS GS MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

!    call timer('ts_GS',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS SSR_sGreen_NoDOS' )
#endif

  contains

    ! We supply a routine to switch the pointer position of alpha,beta / rh1
    subroutine switch_alpha_beta_rh1(as_first)
      logical, intent(inout) :: as_first
      integer :: i 
      ! start
      i = 2 * nosq

      if ( as_first ) then
         rh1 => zwork(i+1:i+2*nosq) 
         i = i + 2*nosq
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
      else
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
         i = i + nosq
         rh1 => zwork(i+1:i+2*nosq) 
      end if
      as_first = .not. as_first

    end subroutine switch_alpha_beta_rh1

  end subroutine SSR_sGreen_NoDOS

  ! Print information regarding the Green function file
  ! In particular the size of the corresponding Green function
  ! file.
  subroutine print_Elec_Green(El, NE, nkpt)
    
    use precision,  only : dp
    use parallel, only: IONode
    use units,      only : eV

    use m_ts_electype

    ! Input variables
    type(Elec), intent(in) :: El
    ! Number of energy-points
    integer :: NE
    ! Number of k-points
    integer :: nkpt

    ! Local variables
    integer :: i, j
    real(dp) :: b
    ! Whether the electrode is pre-expanded...
    integer :: nq
    logical :: pre_expand

    if ( .not. IONode ) return

    nq = El%Bloch%size()
    pre_expand = El%pre_expand > 0 .and. nq > 1

    write(*,'(/,2a)') 'Calculating all surface Green functions for: ',trim(name(El))
    write(*,'(a,f14.5,1x,a)') &
         ' Fermi level shift in electrode (chemical potential) : ',El%mu%mu/eV,' eV'

    ! Show the number of used atoms and orbitals
    write(*,'(a,i6,'' / '',i6)') ' Atoms available    / used atoms   : ', &
         El%na_u,El%na_used
    write(*,'(a,i6,'' / '',i6)') ' Orbitals available / used orbitals: ', &
         El%no_u,El%no_used

    write(*,'(1x,a,i0)') 'Total self-energy calculations: ',nq*NE*nkpt

    ! We show them in units of reciprocal lattice vectors
    do i = 1 , 3
      if ( El%Bloch%B(i) > 1 ) then
        write(*,'(3(a,i0))') ' Bloch expansion k-points in A_',i, &
            ' direction [b_',i,']: ',El%Bloch%B(i)
      end if
    end do

    write(*,'(2a)') ' Saving surface Green functions in: ',trim(El%GFfile)

    if ( pre_expand .and. El%pre_expand == 1 ) then
       ! only the SSG is expanded
       b = El%nspin * nkpt * ( 2 + NE * nq )
    else
       b = El%nspin * nkpt * ( 2 + NE ) * nq
    end if

    ! Correct estimated file-size for fully expanded
    if ( pre_expand ) b = b * nq

    ! to complex double precision
    b = b * El%no_used ** 2 * 16._dp

    ! to MB
    b = b / 1024._dp ** 2

    if ( b > 2001._dp ) then
       b = b / 1024._dp
       write(*,'(a,f10.3,a)') ' Estimated file size: ',b,' GB'
    else
       write(*,'(a,f10.3,a)') ' Estimated file size: ',b,' MB'
    end if

  end subroutine print_Elec_Green


! ##################################################################
! ## Driver subroutine for calculating the (ideal)                ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                 Updated by : Nick Papior Andersen            ##
! ## It has now been parallelized to speed up electrode           ##
! ## surface Green function generation.                           ##
! ## It generates the surface Green function by handling          ##
! ## repetition as well.                                          ##
! ##################################################################
  subroutine create_Green(El, &
       ucell,nkpnt,kpoint,kweight, &
       NEn,ce, &
       DOS,T)

    use precision,  only : dp
    use parallel  , only : Node, Nodes, IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Sum, MPI_Max, MPI_integer
    use mpi_siesta, only : MPI_Wait,MPI_Status_Size
    use mpi_siesta, only : MPI_double_complex
    use mpi_siesta, only : MPI_double_precision
#endif
    use m_ts_electype
    use m_mat_invert

    use m_ts_elec_se, only : update_UC_expansion_A

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use m_iterator

! ***********************
! * INPUT variables     *
! ***********************
    type(Elec), intent(inout) :: El  ! The electrode 
    integer,  intent(in)      :: nkpnt ! Number of k-points
    real(dp), intent(in)      :: kpoint(3,nkpnt) ! k-points
    real(dp), intent(in)      :: kweight(nkpnt) ! weights of kpoints
    real(dp), dimension(3,3)  :: ucell ! The unit cell of the CONTACT
    integer, intent(in)       :: NEn ! Number of energy points
    complex(dp), intent(in)   :: ce(NEn) ! the energy points

! ***********************
! * OUTPUT variables    *
! ***********************
    real(dp), intent(inout), optional :: DOS(El%no_u,NEn,El%nspin)
    real(dp), intent(inout), optional :: T(NEn,El%nspin)

! ***********************
! * LOCAL variables     *
! ***********************
    ! Array for holding converted k-points
    real(dp) :: bkpt(3), kpt(3), kq(3), wq, rcell(3,3)
    real(dp), allocatable :: lDOS(:)
    
    ! Dimensions
    integer :: nq, nspin, n_s
    integer :: nuo_E, nS, nuou_E, nuS, no_X, n_X

    ! Electrode transfer and hamiltonian matrix
    complex(dp), pointer :: H00(:) => null()
    complex(dp), pointer :: S00(:) => null()
    complex(dp), pointer :: H01(:) => null()
    complex(dp), pointer :: S01(:) => null()
    complex(dp), pointer :: zwork(:) => null()
    complex(dp), pointer :: zHS(:) => null()
    real(dp), allocatable :: sc_off(:,:)

    ! Expanded arrays
    complex(dp), pointer :: X(:) => null()

    ! Green function variables
    complex(dp), pointer :: GS(:)
    complex(dp), pointer :: Hq(:), Sq(:), Gq(:)
    complex(dp) :: ZEnergy

    ! In order to print information about the recursize algorithm
    integer, allocatable :: iters(:,:,:,:)
    real(dp) :: i_mean, i_std

    integer :: uGF
    ! Big loop counters
    type(itt2) :: it2
    integer, pointer :: ispin, ikpt
    integer :: iEn, iqpt
    ! Counters
    integer :: i, j, io, jo, off

    logical :: CalcDOS, CalcT, pre_expand
    logical :: is_left, Gq_allocated, reduce_size

#ifdef MPI
    integer :: MPIerror, curNode
    integer :: req, status(MPI_Status_Size)
    integer, allocatable :: reqs(:)
#endif
    
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE create_Green' )
#endif

    call timer('TS_SE',1)

    CalcDOS = present(DOS)
    CalcT = present(T)

    ! Check input for what to do
    if( El%inf_dir == INF_NEGATIVE ) then
       is_left = .true.
    else if( El%inf_dir == INF_POSITIVE ) then
       is_left = .false.
    else
       call die("init electrode has received wrong job ID [L,R].")
    endif

    ! Initialize TSGF-file
    call init_TSGF()

    ! capture information from the electrode
    nspin  = El%nspin
    nuo_E  = El%no_u
    nS     = nuo_E ** 2
    nuou_E = El%no_used
    nuS    = nuou_E ** 2
    ! create expansion k-points (weight of q-points)
    nq     = El%Bloch%size()
    wq     = 1._dp / real(nq,dp)
    ! We also need to invert to get the contribution in the
    reduce_size = nuo_E /= nuou_E
    no_X = nuou_E * nq
    n_X  = no_X ** 2
    pre_expand = El%pre_expand > 0 .and. nq > 1

    ! Calculate offsets
    n_s = size(El%isc_off,dim=2)
    allocate(sc_off(3,n_s))
    sc_off = matmul(El%cell,El%isc_off)

    ! Print information on file-size and electrode type.
    call print_Elec_Green(El, NEn, nkpnt)
    
    ! Initialize Green function and Hamiltonian arrays
    nullify(GS)
    if ( nS /= nuS ) then
       allocate(GS(nS))
       call memory('A','Z',nS,'create_green')
    !else
    !  the regions are of same size, so we can just point
    !  to the correct memory segment
    end if

    ! Allocate work array
    i = max(nS*9,nuS*nq*2)
    if ( pre_expand ) then
       i = max(i,n_X)
    end if
    allocate(zwork(i))
    call memory('A','Z',i,'create_green')

    ! Point the hamiltonian and the overlap to the work array
    ! The work-array is only used for calculation the surface
    ! Green function and
    Hq => zwork(1:nuS*nq)
    Sq => zwork(nuS*nq+1:nuS*nq*2)
    if ( size(zwork) >= nS * 9 + nuS*nq ) then
       Gq => zwork(nS*9+1:nS*9+nuS*nq)
       Gq_allocated = .false.
    else
       nullify(Gq)
       allocate(Gq(nuS*nq))
       call memory('A','Z',size(Gq),'create_green')
       Gq_allocated = .true.
    end if

    if ( pre_expand ) then
       ! We allocate space for pre-expansion of the arrays
       allocate(X(n_X))
    end if

    ! all the Hamiltonian and overlaps
    allocate(zHS(nS * nq * 4))
    call memory('A','Z',nS * nq * 4,'create_green')

    ! Prepare for the inversion
    i = max(no_X,nuo_E)
    call init_mat_inversion(i)

    ! Reset bulk DOS
    if ( CalcDOS ) then
       allocate(lDOS(nuo_E))
       DOS(:,:,:) = 0._dp
    end if
    if ( CalcT ) then
       T = 0._dp
    end if

!******************************************************************
!           Start Green function calculation
!******************************************************************
    
#ifdef MPI
    if ( IONode ) then
       allocate(reqs(Nodes-1))
       call memory('A','I',Nodes-1,'create_green')
       ! Create request handles for communication
       ! This is a rather new feature which enhances communication times.
       ! However, this is perhaps overkill as we never have VERY many 
       ! contour points. Say NEn > 1000
       ! Look in the loop for MPI_Start(...) for where this is used
       do i = 1 , Nodes - 1
          if ( pre_expand ) then
             call MPI_Recv_Init(X(1),n_X,MPI_double_complex, &
                  i,i,MPI_Comm_World,reqs(i),MPIerror)
          else
             call MPI_Recv_Init(Gq(1),nuS*nq,MPI_double_complex, &
                  i,i,MPI_Comm_World,reqs(i),MPIerror)
          end if
       end do
    else
       ! Create request handles for communication
       if ( pre_expand ) then
          call MPI_Send_Init(X(1),n_X,MPI_double_complex, &
               0,Node,MPI_Comm_World,req,MPIerror)
       else
          call MPI_Send_Init(Gq(1),nuS*nq,MPI_double_complex, &
               0,Node,MPI_Comm_World,req,MPIerror)
       end if
    end if
#endif

    ! prepare the iteration counter
    allocate(iters(nq,NEn,nkpnt,2))
    if ( IONode ) then
       ! TODO when adding new surface-Green functions schemes, please update here
       write(*,'(1x,a)') 'Lopez Sancho, Lopez Sancho & Rubio recursive &
            &surface self-energy calculation...'
    end if

    ! start up the iterators
    call itt_init  (it2,end1=nspin,end2=nkpnt)
    call itt_attach(it2,cur1=ispin,cur2=ikpt)

    call reclat(El%cell,rcell,1)

    ! do spin and k-point loop in one go...
    do while ( .not. itt_step(it2) )
       
       if ( itt_stepped(it2,1) ) then
          ! Number of iterations
          iters(:,:,:,:) = 0
       end if
       
       ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
       call Elec_kpt(El,ucell,kpoint(:,ikpt),bkpt, opt = 2)
       ! We need to save the k-point for the "expanded" super-cell
       El%bkpt_cur = bkpt
       
       ! loop over the repeated cell...
       HSq_loop: do iqpt = 1 , nq
             
          ! point to the correct segment of memory
          H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
          S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
          H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
          S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)

          ! init qpoint in reciprocal lattice vectors
          kpt = bkpt(:) + q_exp(El,iqpt)
          ! Convert to 1/Bohr
          call kpoint_convert(rcell,kpt,kq,-2)

          ! Setup the transfer matrix and the intra cell at the k-point and q-point
          ! Calculate transfer matrices @Ef (including the chemical potential)
          call set_HS_Transfer(ispin, El, n_s,sc_off, kq, &
               nuo_E, H00,S00,H01,S01)

          i = (iqpt-1)*nuS
          if ( reduce_size ) then
             if( is_left ) then
                ! Left, we use the last orbitals
                off = nuo_E - nuou_E + 1
                do jo = off - 1 , nuo_E - 1
                   do io = off , nuo_E
                      i = i + 1
                      Hq(i) = H00(nuo_E*jo+io)
                      Sq(i) = S00(nuo_E*jo+io)
                   end do
                end do
             else
                ! Right, the first orbitals
                do jo = 0 , nuou_E - 1
                   do io = 1 , nuou_E
                      i = i + 1
                      Hq(i) = H00(nuo_E*jo+io)
                      Sq(i) = S00(nuo_E*jo+io)
                   end do   ! io
                end do      ! jo
             end if
          end if
          
       end do HSq_loop

       ! Save Hamiltonian and overlap
       call store_HS()
       
       Econtour_loop: do iEn = 1, NEn

#ifdef MPI
          ! Every node takes one energy point
          ! This asserts that IONode = Node == 0 will have iEn == 1
          ! Important !
          curNode = MOD(iEn-1,Nodes)
          E_Nodes: if ( curNode == Node ) then
#endif
             ! as we already have shifted H,S to Ef + mu, and ZEnergy is
             ! wrt. mu, we don't need to subtract mu again
             ZEnergy = ce(iEn)
             i_mean = 0._dp
             
             ! loop over the repeated cell...
             q_loop: do iqpt = 1 , nq

                H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
                S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
                H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
                S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)
                if ( nS == nuS ) then
                   ! instead of doing a copy afterward, we can
                   ! put it the correct place immediately
                   GS => Gq((    iqpt-1)*nS+1:      iqpt *nS)
                end if

                ! Calculate the surface Green function
                ! Zenergy is wrt. to the system Fermi-level
                if ( CalcDOS ) then
                   lDOS = 0._dp
                   call SSR_sGreen_DOS(nuo_E,ZEnergy,H00,S00,H01,S01, &
                        El%accu, GS, &
                        lDOS,i_mean,9*nS,zwork, &
                        iterations=iters(iqpt,iEn,ikpt,1), final_invert = reduce_size)
                   
                   ! We also average the k-points.
                   DOS(:,iEn,ispin) = DOS(:,iEn,ispin) + lDOS * wq * kweight(ikpt)
                   if ( CalcT ) T(iEn,ispin) = T(iEn,ispin) + i_mean

                else
                   call SSR_sGreen_NoDos(nuo_E,ZEnergy,H00,S00,H01,S01, &
                        El%accu, GS, &
                        8*nS,zwork, &
                        iterations=iters(iqpt,iEn,ikpt,1), final_invert = reduce_size)
                   
                end if
                  
                ! Copy over surface Green function
                i = (iqpt-1)*nuS
                if ( reduce_size ) then
                   if ( is_left ) then
                      ! Left, we use the last orbitals
                      off = nuo_E - nuou_E + 1
                      do jo = off - 1 , nuo_E - 1
                         do io = off , nuo_E
                            i = i + 1
                            Gq(i) = GS(nuo_E*jo+io)
                         end do           ! io
                      end do              ! jo
                   else
                      ! Right, the first orbitals
                      do jo = 0 , nuou_E-1
                         do io = 1 , nuou_E
                            i = i + 1
                            Gq(i) = GS(nuo_E*jo+io)
                         end do           ! io
                      end do              ! jo
                   end if

                   if ( nq == 1 ) then
                      ! We invert back here, instead of in
                      ! the SCF (this is important as the
                      ! decreased size of the surface-Green function
                      ! would otherwise yield a different result)
                      call mat_invert(Gq(1:nuS),zwork(1:nuS),&
                           nuou_E, &
                           MI_IN_PLACE_LAPACK)

                   end if

                end if
                
             end do q_loop

             if ( pre_expand ) then
                ! Expand this energy-point
                call update_UC_expansion_A(nuou_E,no_X,El,nq, &
                     El%na_used,El%lasto_used,Gq,X)
                if ( reduce_size ) then
                   call mat_invert(X(1:n_X),zwork(1:n_X),&
                        no_X, &
                        MI_IN_PLACE_LAPACK)
                end if
             end if
                
#ifdef MPI
             ! If not IONode we should send message
             ! This message parsing is directly connected to 
             ! a predefined size of the message, see right before
             ! spin loop.
             ! It communicates the Gq array to the Gq array
             if ( .not. IONode ) then
                call MPI_Start(req,MPIerror)
                call MPI_Wait(req,status,MPIerror)
             end if
             
          end if E_Nodes

#endif

          ! Save the surface Green function file
          call store_GS()

       end do Econtour_loop
          

       if ( itt_last(it2,2) ) then
#ifdef MPI
          call MPI_Reduce(iters(1,1,1,1), iters(1,1,1,2), nq*NEn*nkpnt, &
               MPI_Integer, MPI_Sum, 0, MPI_Comm_World, MPIerror)
!$OMP parallel default(shared), private(j,i,iqpt)
#else
!$OMP parallel default(shared), private(j,i,iqpt)

          iters(:,:,:,2) = iters(:,:,:,1)
#endif
          if ( IONode ) then
             i_mean = sum(iters(:,:,:,2)) / real(nq*NEn*nkpnt,dp)

!$OMP single
             i_std = 0._dp
!$OMP end single ! keep barrier

!$OMP do reduction(+:i_std)
             do j = 1 , nkpnt
             do i = 1 , NEn
             do iqpt = 1 , nq
                i_std = i_std + ( iters(iqpt,i,j,2) - i_mean ) ** 2
             end do
             end do
             end do
!$OMP end do

!$OMP master
             i_std = sqrt(i_std/real(NEn*nq*nkpnt,dp))
             ! TODO if new surface-Green function scheme is implemented, fix here
             write(*,'(1x,a,f10.4,'' / '',f10.4)') 'Lopez Sancho, Lopez Sancho & Rubio: &
                  &Mean/std iterations: ', i_mean             , i_std
             write(*,'(1x,a,i10,'' / '',i10)')     'Lopez Sancho, Lopez Sancho & Rubio: &
                  &Min/Max iterations : ', minval(iters(:,:,:,2)) , maxval(iters(:,:,:,2))
!$OMP end master
             
          end if
!$OMP end parallel
       end if

    end do
!*******************************************************************
!         Green function calculation is done
!*******************************************************************

    deallocate(iters)

#ifdef MPI
    ! Free requests made for the communications
    if ( IONode ) then
       do i = 1 , Nodes - 1 
          call MPI_Request_Free(reqs(i),MPIerror)
       end do
       call memory('D','I',Nodes-1,'create_green')
       deallocate(reqs)
    else
       call MPI_Request_Free(req,MPIerror)
    end if
#endif

    ! Close and finish
    call finish_TSGF()
    
    ! Clean up computational arrays
    if ( nS /= nuS ) then
       call memory('D','Z',size(GS),'create_green')
       deallocate(GS)
    end if

    if ( Gq_allocated ) then
       call memory('D','Z',size(Gq),'create_green')
       deallocate(Gq)
    end if

    ! Work-arrays
    call memory('D','Z',size(zwork),'create_green')
    deallocate(zwork)
    call memory('D','Z',size(zHS),'create_green')
    deallocate(zHS)

    if ( CalcDOS ) deallocate(lDOS)

    if ( pre_expand ) deallocate(X)

    call itt_destroy(it2)

    call clear_mat_inversion()

#ifdef MPI
    if ( CalcDOS ) then
       ! Sum the bulkdensity of states
       ! Here we can safely use the array as temporary (Gq)
       allocate(lDOS(nuo_E*NEn*nspin))
    else if ( CalcT ) then
       allocate(lDOS(NEn*nspin))
    end if
    if ( allocated(lDOS) ) then
       call memory('A','D',size(lDOS),'create_green')
    end if
    
    if ( CalcDOS ) then
       call MPI_AllReduce(DOS(1,1,1),lDOS(1),nuo_E*NEn*nspin, &
            MPI_double_precision, &
            MPI_Sum,MPI_Comm_World,MPIerror)
       i = 0
       do jo = 1 , nspin
          do io = 1 , NEn
             DOS(1:nuo_E,io,jo) = lDOS(i+1:i+nuo_E)
             i = i + nuo_E
          end do
       end do
       
    end if
    if ( CalcT ) then
       call MPI_AllReduce(T(1,1),lDOS(1),NEn*nspin, &
            MPI_double_precision, &
            MPI_Sum,MPI_Comm_World,MPIerror)
       i = 0
       do jo = 1 , nspin
          T(1:NEn,jo) = lDOS(i+1:i+NEn)
          i = i + NEn
       end do
       
    end if

    if ( allocated(lDOS) ) then
       call memory('D','D',size(lDOS),'create_green')
       deallocate(lDOS)
    end if
#endif

    call timer('TS_SE',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS create_Green' )
#endif

  contains


    subroutine init_TSGF()
      
      real(dp), allocatable :: kE(:,:)

      if ( .not. IONode ) return
      
      call io_assign(uGF)
      open(FILE=El%GFfile,UNIT=uGF,FORM='UNFORMATTED')
      
      ! Electrode information
      write(uGF) El%nspin, El%cell
      write(uGF) El%na_u, El%no_u
      write(uGF) El%na_used, El%no_used
      write(uGF) El%xa_used, El%lasto_used
      write(uGF) El%repeat, El%Bloch%B(:), El%pre_expand
      write(uGF) El%mu%mu
      
      ! Write out explicit information about this content
      write(uGF) nkpnt
      ! Notice that we write the k-points for the ELECTRODE
      ! They will be stored in units of the reciprocal lattice vector
      !   1/b 
      allocate(kE(3,nkpnt))
      do i = 1 , nkpnt
         ! Store the k-points in units of reciprocal lattice
         call Elec_kpt(El,ucell,kpoint(:,i),kE(:,i), opt = 2)
      end do
      write(uGF) kE, kweight
      deallocate(kE)
      
      ! write out the contour information
      write(uGF) NEn
      write(uGF) ce ! energy points
      
    end subroutine init_TSGF

    subroutine store_HS()
      
      if ( .not. IONode ) return
      ! k-point and energy-point is in front of Hamiltonian

      write(uGF) ikpt, 1, ce(1) ! k-point and energy point
      if ( reduce_size ) then
         if ( pre_expand .and. El%pre_expand > 1 ) then
            call update_UC_expansion_A(nuou_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,Hq,X)
            write(uGF) X
            call update_UC_expansion_A(nuou_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,Sq,X)
            write(uGF) X
         else
            write(uGF) Hq
            write(uGF) Sq
         end if
      else
         H00 => zHS(      1:nq*nS  )
         S00 => zHS(nq*nS+1:nq*nS*2)
         if ( pre_expand .and. El%pre_expand > 1 ) then
            call update_UC_expansion_A(nuo_E,no_X,El,nq, &
                 El%na_used,El%lasto_used,H00,X)
            write(uGF) X
            call update_UC_expansion_A(nuo_E,no_X,El,nq,&
                 El%na_used,El%lasto_used,S00,X)
            write(uGF) X
         else
            write(uGF) H00
            write(uGF) S00
         end if
      end if

    end subroutine store_HS

    subroutine store_GS()
      
      ! Save Surface Green function
      if ( .not. IONode ) return

#ifdef MPI
      if ( curNode /= Node ) then
         call MPI_Start(reqs(curNode),MPIerror)
         call MPI_Wait(reqs(curNode),status,MPIerror)
      end if
#endif
         
      ! Write out calculated information at E point
      if ( iEn /= 1 ) write(uGF) ikpt, iEn, ce(iEn)
      if ( pre_expand ) then
         write(uGF) X
      else
         write(uGF) Gq
      end if
      
    end subroutine store_GS

    subroutine finish_TSGF()
      
      ! Close file
      if ( .not. IONode ) return
      
      call io_close(uGF)
      write(*,'(a,/)') "Done creating '"//trim(El%GFfile)//"'."  

    end subroutine finish_TSGF
    
  end subroutine create_Green


  subroutine init_Electrode_HS(El)
    use fdf, only: fdf_get
#ifdef MPI
    use mpi_siesta
#endif
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    use m_ts_electype

    type(Elec), intent(inout) :: El
#ifdef MPI
    integer :: error
#endif
    logical :: neglect_conn
    
    ! Read-in and create the corresponding transfer-matrices
    call delete(El) ! ensure clean electrode
    call read_Elec(El,Bcast=.true.)

    if ( .not. associated(El%isc_off) ) then
       call die('An electrode file needs to be a non-Gamma calculation. &
            &Ensure good periodicity in the T-direction.')
    end if

    ! Create the default sparsity patterns in the sub-spaces needed
    ! for the self-energy calculations
    ! This is also important to create before running
    ! check_connectivity because of used-atoms possibly being set.
    call create_sp2sp01(El)

    ! print out the precision of the electrode (whether it extends
    ! beyond first principal layer)
    if ( check_connectivity(El) ) then
       neglect_conn = .true.
    else
       neglect_conn = fdf_get('TS.Elecs.Neglect.Principal', .false.)
#ifdef TBTRANS
       neglect_conn = fdf_get('TBT.Elecs.Neglect.Principal', neglect_conn)
#endif
    end if
    
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,error)
#endif
    if ( .not. neglect_conn ) then
       call die('Electrode connectivity is not perfect, &
            &refer to the manual for achieving a perfect electrode.')
    end if
    
    ! Clean-up, we will not need these!
    ! we should not be very memory hungry now, but just in case...
    call delete(El%H)
    call delete(El%S)
   
    ! We do not accept onlyS files
    if ( .not. initialized(El%H00) ) then
       call die('An electrode file must contain the Hamiltonian')
    end if

    call delete(El%sp)

  end subroutine init_Electrode_HS


!**********
! Create the Hamiltonian for the electrode as well
! as creating the transfer matrix.
!**********
  subroutine set_HS_Transfer_1d(ispin,El,n_s,sc_off,kq, &
       no,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
    use m_ts_electype
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)    :: ispin, no
    type(Elec), intent(inout) :: El
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3)   ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(no**2) :: Hk,Sk,Hk_T,Sk_T

    call set_HS_Transfer_2d(ispin,El,n_s,sc_off,kq, &
         no,Hk,Sk,Hk_T,Sk_T)
    
  end subroutine set_HS_Transfer_1d
  
  subroutine set_HS_Transfer_2d(ispin,El,n_s,sc_off,kq, &
       no,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
    use m_ts_electype
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in) :: ispin, no
    type(Elec), intent(inout) :: El
    integer, intent(in) :: n_s
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    real(dp), intent(in) :: kq(3)   ! k + q-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(no,no) :: Hk, Sk, Hk_T, Sk_T

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: Ef
    complex(dp) :: ph(0:n_s-1)
    integer :: i, j, io, jo, ind, is
    integer, pointer :: ncol00(:), l_ptr00(:), l_col00(:)
    integer, pointer :: ncol01(:), l_ptr01(:), l_col01(:)
    real(dp), pointer :: H00(:,:) , S00(:), H01(:,:), S01(:)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE elec_HS_Transfer' )
#endif

    ! we need to subtract as the below code shifts to Ef
    Ef = El%Ef - El%mu%mu

    if ( El%no_u /= no ) call die('Wrong size of the electrode array')

    ! retrieve values
    call attach(El%sp00,n_col=ncol00,list_ptr=l_ptr00,list_col=l_col00)
    call attach(El%sp01,n_col=ncol01,list_ptr=l_ptr01,list_col=l_col01)
    ! point to the data-segments...
    H00 => val(El%H00)
    H01 => val(El%H01)
    S00 => val(El%S00)
    S01 => val(El%S01)

    ! The algorithm outside should take care of the
    ! nullification of the k-point in the semi-infinite direction
    do i = 0 , n_s - 1
       ph(i) = cdexp(dcmplx(0._dp, sum(kq*sc_off(:,i))) )
    end do

    ! Initialize arrays
    Hk(:,:) = z_0
    Sk(:,:) = z_0
    Hk_T(:,:) = z_0
    Sk_T(:,:) = z_0

!$OMP parallel default(shared), private(i,j,io,jo,ind,is)

    ! We will not have any data-race condition here
!$OMP do 
    do io = 1 , no

       ! Create 00
       do j = 1 , ncol00(io)
          ind = l_ptr00(io) + j
          jo = ucorb(l_col00(ind),no)
          is = (l_col00(ind)-1) / no
          
          Hk(io,jo) = Hk(io,jo) + H00(ind,ispin) * ph(is)
          Sk(io,jo) = Sk(io,jo) + S00(ind)       * ph(is)
       enddo

       ! Create 01
       do j = 1 , ncol01(io)
          ind = l_ptr01(io) + j
          jo = ucorb(l_col01(ind),no)
          is = (l_col01(ind)-1) / no

          Hk_T(io,jo) = Hk_T(io,jo) + H01(ind,ispin) * ph(is)
          Sk_T(io,jo) = Sk_T(io,jo) + S01(ind)       * ph(is)
          
       end do

    end do
!$OMP end do

    ! Symmetrize 00 and make EF the energy-zero
    ! We will not have any data-race condition here
!$OMP do
    do io = 1 , no
       do jo = 1 , io - 1

          Sk(io,jo) = 0.5_dp*( Sk(io,jo) + dconjg(Sk(jo,io)) )
          Sk(jo,io) = dconjg(Sk(io,jo))

          Hk(io,jo) = 0.5_dp*( Hk(io,jo) + dconjg(Hk(jo,io)) ) - &
               Ef * Sk(io,jo)
          Hk(jo,io) = dconjg(Hk(io,jo))

          ! Transfer matrix is not symmetric, so do not symmetrize
          Hk_T(jo,io) = Hk_T(jo,io) - Ef * Sk_T(jo,io)
          Hk_T(io,jo) = Hk_T(io,jo) - Ef * Sk_T(io,jo)

       end do
       
       Sk(io,io) = real(Sk(io,io),dp)
       Hk(io,io) = real(Hk(io,io),dp) - Ef * Sk(io,io)
       
       ! Transfer matrix
       Hk_T(io,io) = Hk_T(io,io) - Ef * Sk_T(io,io)
       
    end do
!$OMP end do nowait

!$OMP end parallel

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS elec_HS_Transfer' )
#endif

  end subroutine set_HS_Transfer_2d

  subroutine calc_next_GS_Elec(El,ispin,bkpt,Z,nzwork,in_zwork,DOS,T)
    use precision,  only : dp

    use m_ts_electype
    use m_mat_invert

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use alloc, only : re_alloc, de_alloc

    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    type(Elec), intent(inout) :: El
    integer, intent(in) :: ispin
    ! the k-point in reciprocal units of the electrode
    ! also with / Rep
    real(dp), intent(in) :: bkpt(3)
    complex(dp), intent(in) :: Z
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: in_zwork(nzwork)
    ! Possibly the bulk density of states from the electrode
    ! If the DOS, also BULK transmission
    real(dp), intent(inout), optional :: DOS(:), T

    ! ***********************
    ! * LOCAL variables     *
    ! ***********************
    integer  :: iq
    real(dp) :: kpt(3), kq(3), rcell(3,3)
    
    ! Dimensions
    integer :: nq, nw
    integer :: nuo_E, nS, nuou_E, nuS, nuouT_E

    ! Electrode transfer and hamiltonian matrix
    complex(dp), pointer :: H00(:), H01(:), S00(:), S01(:)
    complex(dp), pointer :: zwork(:)
    complex(dp), pointer :: zHS(:) => null()
    real(dp), allocatable :: sc_off(:,:)

    ! Green function variables
    complex(dp), pointer :: GS(:)

    ! size requirement
    integer :: size_req(2)
    ! Counters
    integer :: i, ios, ioe, off, n_s
    logical :: is_left, reduce_size
    logical :: zHS_allocated
    logical :: same_k, calc_DOS

    ! Check input for what to do
    is_left = El%inf_dir == INF_NEGATIVE
    calc_DOS = present(DOS)
    if ( calc_DOS .and. .not. present(T) ) then
       call die('Need both DOS and T')
    end if

    zHS_allocated = .false.

    ! constants for this electrode
    nuo_E  = El%no_u
    nS     = nuo_E ** 2
    nuou_E = El%no_used
    nuS    = nuou_E ** 2
    ! create expansion k-points
    nq     = El%Bloch%size()
    ! We also need to invert to get the contribution in the
    ! reduced region
    reduce_size = nuo_E /= nuou_E
    nuouT_E = TotUsedOrbs(El)

    if ( calc_DOS ) then
       if ( nuo_E > size(DOS) ) &
            call die('Error in DOS size for calculation bulk DOS')
       ! Initialize density of states
       DOS(1:nuo_E) = 0._dp
       T = 0._dp
    end if

    n_s = size(El%isc_off,dim=2)
    allocate(sc_off(3,n_s))
    sc_off = matmul(El%cell,El%isc_off)
    
    ! whether we already have the H and S set correctly, 
    ! update accordingly, it will save a bit of time, but not much
    same_k = abs( bkpt(1) - El%bkpt_cur(1) ) < 1.e-8_dp
    same_k = same_k .and. abs( bkpt(2) - El%bkpt_cur(2) ) < 1.e-8_dp
    same_k = same_k .and. abs( bkpt(3) - El%bkpt_cur(3) ) < 1.e-8_dp
    if ( .not. same_k ) then
      El%bkpt_cur(:) = bkpt

      ! In case we do not need the hamiltonian
      ! This will be the case for non-bias points and when using bulk electrode
      same_k = .not. associated(El%HA)
    end if

    ! determine whether there is room enough
    if ( reduce_size ) then
       size_req(1) = 5 * nS
    else
       size_req(1) = 4 * nS
    end if
    if ( calc_DOS ) then
       size_req(2) = 9 * nS
    else
       size_req(2) = 8 * nS
    end if
    if ( size_req(1) + size_req(2) <= nzwork ) then

       ! we have enough room in the regular work-array for everything
       i = 0
       H00 => in_zwork(i+1:i+nS)
       i = i + nS
       S00 => in_zwork(i+1:i+nS)
       i = i + nS
       H01 => in_zwork(i+1:i+nS)
       i = i + nS
       S01 => in_zwork(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => in_zwork(i+1:i+nS)
          i = i + nS
       end if
       zwork => in_zwork(i+1:nzwork)

    else if ( size_req(2) <= nzwork ) then

       ! we will allocate H00,H01,S00,S01,GS arrays
       call re_alloc(zHS,1,size_req(1),routine='next_GS')
       zHS_allocated = .true.

       i = 0
       H00 => zHS(i+1:i+nS)
       i = i + nS
       S00 => zHS(i+1:i+nS)
       i = i + nS
       H01 => zHS(i+1:i+nS)
       i = i + nS
       S01 => zHS(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => zHS(i+1:i+nS)
       end if
       
       ! the work-array fits the input work-array
       zwork => in_zwork(1:nzwork)

    else if ( size_req(1) <= nzwork ) then
       ! we will allocate 8*nS work array

       i = 0
       H00 => in_zwork(i+1:i+nS)
       i = i + nS
       S00 => in_zwork(i+1:i+nS)
       i = i + nS
       H01 => in_zwork(i+1:i+nS)
       i = i + nS
       S01 => in_zwork(i+1:i+nS)
       i = i + nS
       if ( reduce_size ) then
          GS => in_zwork(i+1:i+nS)
       end if

       call re_alloc(zHS,1,size_req(2),routine='next_GS')
       zHS_allocated = .true.
       zwork => zHS(:)

    else

       call die('Your electrode is too large compared &
            &to your system in order to utilize the in-core &
            &calculation of the self-energies.')

    end if
    ! Get actual size of work-array
    nw = size(zwork)

    call init_mat_inversion(nuo_E)

    ! prepare the indices for the Gamma array
    ios = 1
    ioe = nuS 

    ! create the offset to be used for copying over elements
    off = nuo_E - nuou_E + 1

    call reclat(El%cell,rcell,1)

    ! loop over the repeated cell...
    q_loop: do iq = 1 , nq

       ! init qpoint in reciprocal lattice vectors
       kpt(:) = bkpt(:) + q_exp(El,iq)
      
       ! Convert to 1/Bohr
       call kpoint_convert(rcell,kpt,kq,-2)

       ! Calculate transfer matrices @Ef (including the chemical potential)
       call set_HS_Transfer(ispin, El, n_s,sc_off, kq, &
            nuo_E, H00,S00,H01,S01)
       
       if ( .not. same_k ) then
          ! we only need to copy over the data if we don't already have it calculated
!$OMP parallel default(shared)
          call copy_over(is_left,nuo_E,H00,nuou_E,El%HA(:,:,iq),off)
          call copy_over(is_left,nuo_E,S00,nuou_E,El%SA(:,:,iq),off)
!$OMP end parallel
       end if

       if ( .not. reduce_size ) then
          ! Instead of doing a copy, we store it directly
          GS => El%GA(ios:ioe)
       end if

       ! calculate the contribution for this q-point
       if ( calc_DOS ) then
          call SSR_sGreen_DOS(nuo_E,Z,H00,S00,H01,S01, &
               El%accu, GS, &
               DOS(1:nuo_E), T, &
               nw, zwork, &
               final_invert = reduce_size)
       else
          call SSR_sGreen_NoDOS(nuo_E,Z,H00,S00,H01,S01, &
               El%accu, GS, &
               nw, zwork, &
               final_invert = reduce_size)
       end if

       if ( reduce_size ) then
          ! Copy over surface Green function
          ! first we need to determine the correct placement
!$OMP parallel default(shared)
          call copy_over(is_left,nuo_E,GS,nuou_E,El%GA(ios:ioe),off)
!$OMP end parallel

          ! we need to invert back as we don't need to
          ! expand. And the algorithm expects it to be in correct format
          if ( nq == 1 ) then
             call mat_invert(El%GA(ios:ioe),zwork(1:nuS),&
                  nuou_E, &
                  MI_IN_PLACE_LAPACK)
          end if
       end if

       ! correct indices of Gamma-array
       ios = ios + nuS
       ioe = ioe + nuS

    end do q_loop

    ! We normalize DOS as this will be comparable to a bulk
    ! calculation.
    if ( calc_DOS .and. nq > 1 ) then
       DOS(1:nuo_E) = DOS(1:nuo_E) / nq
    end if
       
    if ( zHS_allocated ) then
       call de_alloc(zHS, routine='next_GS')
    end if

    deallocate(sc_off)
    call clear_mat_inversion()

  contains
    
    subroutine copy_over(is_left,fS,from,tS,to,off)
      logical, intent(in) :: is_left
      integer, intent(in) :: fS, tS, off
      complex(dp), intent(in) :: from(fS,fS)
      complex(dp), intent(inout) :: to(tS,tS)

      integer :: i, j, ioff

      if ( is_left ) then
         ! Left, we use the last orbitals
         ioff = 1 - off ! ioff is private in OMP orphaned routines
!$OMP do private(j,i)
         do j = off , fS
            do i = off , fS
               to(ioff+i,ioff+j) = from(i,j)
            end do
         end do
!$OMP end do nowait
      else
         ! Right, the first orbitals
!$OMP do private(j,i)
         do j = 1 , tS
            do i = 1 , tS
               to(i,j) = from(i,j)
            end do
         end do
!$OMP end do nowait
      end if

    end subroutine copy_over

  end subroutine calc_next_GS_Elec

end module m_ts_electrode
