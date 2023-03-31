module m_mat_invert
  
  use precision, only : dp
  use m_pivot_array, only : Npiv, ipiv
  use m_pivot_array, only : init_mat_inversion => init_pivot
  use m_pivot_array, only : clear_mat_inversion => clear_pivot
  
  implicit none

  private

  public :: mat_invert
  public :: init_mat_inversion
  public :: clear_mat_inversion

  ! This generic inversion module can perform 3 different kinds of inversions on square matrices
  ! 1) Invert in-place
  !  a) use LAPACK
  !  b) partition the matrix into 4 parts and do recursive inversion
  ! 2) Invert to work array

  integer, public, parameter :: MI_IN_PLACE_LAPACK = 1
  integer, public, parameter :: MI_IN_PLACE_RECURSIVE = 2
  integer, public, parameter :: MI_WORK = 3

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

contains

  subroutine mat_invert(M,zwork,no,method,ierr)
    use intrinsic_missing, only : EYE
    integer, intent(in) :: no ! Size of problem
    complex(dp), target :: M(no*no), zwork(no*no)
    integer, intent(in), optional :: method
    integer, intent(out), optional :: ierr

    integer :: lmethod, lierr
    
    if ( present(method) ) then
       lmethod = method
    else
       lmethod = MI_IN_PLACE_LAPACK
    end if

    if ( Npiv < no ) then
       write(0,'(a,i0)') 'Current pivoting array size: ',Npiv
       write(0,'(a,i0)') 'Wanted pivoting array size : ',no
       call die('Pivoting arrays not initialized')
    end if

    ! Initialize the error
    if ( present(ierr) ) ierr = 0

    select case ( lmethod )
    case ( MI_IN_PLACE_LAPACK ) 
       call zgetrf(no, no, M, no, ipiv, lierr )
       if ( lierr /= 0 ) then
          if ( present(ierr) ) then
             ierr = lierr
             return
          else
             call die('Error in LU-decomposition')
          end if
       end if
       call zgetri(no, M, no, ipiv, zwork, no**2, lierr)
       if ( lierr /= 0 ) then
          if ( present(ierr) ) then
             ierr = lierr
             return
          else
             call die('Error in inversion')
          end if
       end if
    case ( MI_IN_PLACE_RECURSIVE ) 
       call mat_invert_recursive(M,zwork,no,ierr=ierr)
    case ( MI_WORK )
       call EYE(no,zwork)
       call zgesv(no,no,M,no,ipiv,zwork,no,lierr)
       if ( lierr /= 0 ) then
          if ( present(ierr) ) then
             ierr = lierr
             return
          else
             call die('Error in inversion')
          end if
       end if
    case default
       call die('Unknown type of inversion')
    end select

  end subroutine mat_invert

#ifdef _NOT_USED
  recursive subroutine mat_invert_recursive(M, zwork, no, ierr)
    use intrinsic_missing, only : EYE
    integer, intent(in) :: no ! Size of problem
    complex(dp), target :: M(no*no), zwork(no*no)
    integer, intent(out), optional :: ierr

    ! The maximum dimensionality of the problem before we turn to a direct inversion algorithm
    integer, parameter :: N_MAX = 64

    complex(dp), pointer :: A1(:), C2(:), B1(:), A2(:)
    complex(dp), pointer :: t1(:), t2(:)
    integer :: sA1, eA1, sC2, eC2
    integer :: sB1, eB1, sA2, eA2
    integer :: n1, n2, i1, i2, idx, i, j

    ! M is the matrix we wish to invert using the tri-diagonal
    ! routine.
    ! This will enable us to half no in order to "faster" 
    ! invert the matrix

    if ( no <= N_MAX ) then
       call mat_invert(M,zwork,no, method = MI_IN_PLACE_LAPACK , ierr=ierr)
       return
    end if

    if ( present(ierr) ) ierr = 0

    ! Calculate the partition sizes of the matrix problem
    n1 = no / 2
    n2 = no - n1

    if ( Npiv < max(n1,n2) ) then
       call die('Error in initialization of the pivoting &
            &array.')
    end if

    ! Point to the partitions
    sA1 = 1
    eA1 = n1 ** 2
    sB1 = eA1 + 1
    eB1 = sB1 + n1 * n2 - 1
    sC2 = eB1 + 1
    eC2 = sC2 + n1 * n2 - 1
    sA2 = eC2 + 1
    eA2 = sA2 + n2 ** 2 - 1
    A1 => zwork(sA1:eA1)
    B1 => zwork(sB1:eB1)
    C2 => zwork(sC2:eC2)
    A2 => zwork(sA2:eA2)

    ! Copy over the parts to ready the inversion
    idx = 0
    i1 = 0
    i2 = 0
    do i = 1 , n1
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             A1(i1) = M(idx)
          else
             i2 = i2 + 1
             B1(i2) = M(idx)
          end if
       end do
    end do
    i1 = 0
    i2 = 0
    do i = n1+1 , no
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             C2(i1) = M(idx)
          else
             i2 = i2 + 1
             A2(i2) = M(idx)
          end if
       end do
    end do
    
    ! Now we have created the matrices in the correct space
    ! Calculate both X1 and Y2
    
    ! Copy over A2 array
    t1 => M(sA2:eA2)
    t2 => M(sB1:eB1)
    call zcopy(eA2-sA2+1, A2(1), 1, t1(1), 1)
    call zcopy(eB1-sB1+1, B1(1), 1, t2(1), 1)

    ! Calculate X1/C2 (store in B1 in original matrix)
    call zgesv(n2,n1,t1,n2,ipiv,t2,n2,i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting X1/C2')
       end if
    end if

    ! Copy over A1 array
    t1 => M(sA1:eA1)
    t2 => M(sC2:eC2)
    call zcopy(eA1-sA1+1, A1(1), 1, t1(1), 1)
    call zcopy(eC2-sC2+1, C2(1), 1, t2(1), 1)

    ! Calculate Y2/B1
    call zgesv(n1,n2,t1,n1,ipiv,t2,n1,i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting Y2/B1')
       end if
    end if

! <<<<<<<<<<<<< direct method
    ! Calculate the diagonal inverted matrix
    ! Calculate: A1 - X1
    ! Copy over A1 array
!    t1 => M(sA1:eA1)
!    t1(:) = A1(:)
!    t2 => M(sB1:eB1) ! this is X1/C2
!    call zgemm('N','N',n1,n1,n2, &
!         zm1, C2,n1, t2,n2,z1, t1,n1)

!    call EYE(n1,A1)
!    call zgesv(n1,n1,t1,n1,ipiv,A1,n1,i)
!    if ( i /= 0 ) call die('Error on inverting A1-X1')
! >>>>>>>>>>>>>>>>>>>> end of direct method

    ! Calculate the diagonal inverted matrix
    ! Calculate: A1 - X1
    ! Copy over A1 array
    t2 => M(sB1:eB1) ! this is X1/C2
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n1,n1,n2, &
         zm1, C2,n1, t2,n2,z1, A1,n1)

    t1 => M(sA1:eA1)
    call mat_invert(A1,t1,n1,ierr=ierr)
    if ( present(ierr) ) then
       if ( ierr /= 0 ) return
    end if

! <<<<<<<<<<<<< direct method
    ! Calculate: A2 - Y2
    ! Copy over A2 array
!    t1 => M(sA2:eA2)
!    t1(:) = A2(:)
!    t2 => M(sC2:eC2) ! this is Y2/B1
!    call zgemm('N','N',n2,n2,n1, &
!         zm1, B1,n2, t2,n1,z1, t1,n2)

!    call EYE(n2,A2)
!    call zgesv(n2,n2,t1,n2,ipiv,A2,n2,i)
!    if ( i /= 0 ) call die('Error on inverting A2-Y2')
! >>>>>>>>>>>>>>>>>>>> end of direct method

    ! Calculate: A2 - Y2
    ! Copy over A2 array
    t2 => M(sC2:eC2) ! this is Y2/B1
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n2,n2,n1, &
         zm1, B1,n2, t2,n1,z1, A2,n2)

    t1 => M(sA2:eA2)
    call mat_invert(A2,t1,n2,ierr=ierr)
    if ( present(ierr) ) then
       if ( ierr /= 0 ) return
    end if

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: X1/C2 * M11
    t1 => M(sB1:eB1) ! this is X1/C2
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n2,n1,n1, &
         zm1, t1,n2, A1,n1,z0, B1,n2)

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: Y2/B1 * M22
    t1 => M(sC2:eC2) ! this is Y2/B1
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n1,n2,n2, &
         zm1, t1,n1, A2,n2,z0, C2,n1)

    ! Copy back the result
    ! Notice that the inverted matrix is back in the same matrix
    idx = 0
    i1 = 0
    i2 = 0
    do i = 1 , n1
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             M(idx) = A1(i1)
          else
             i2 = i2 + 1
             M(idx) = B1(i2)
          end if
       end do
    end do
    i1 = 0
    i2 = 0
    do i = n1+1 , no
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             M(idx) = C2(i1)
          else
             i2 = i2 + 1
             M(idx) = A2(i2)
          end if
       end do
    end do

  end subroutine mat_invert_recursive
#endif

  recursive subroutine mat_invert_recursive(M, zwork, no, ierr)
    use intrinsic_missing, only : EYE
    integer, intent(in) :: no ! Size of problem
    complex(dp), target :: M(no*no), zwork(no*no)
    integer, intent(out), optional :: ierr

    ! The maximum dimensionality of the problem before we turn to a direct inversion algorithm
    integer, parameter :: N_MAX = 64
    integer, parameter :: sA1 = 1

    complex(dp), pointer :: A1(:), C2(:), B1(:), A2(:)
    integer :: sC2, sB1, sA2
    integer :: n1, n2, i, j

    ! M is the matrix we wish to invert using the tri-diagonal
    ! routine.
    ! This will enable us to half no in order to "faster" 
    ! invert the matrix

    if ( no <= N_MAX ) then
       call mat_invert(M,zwork,no, method = MI_IN_PLACE_LAPACK , ierr=ierr)
       return
    end if

    if ( present(ierr) ) ierr = 0

    ! Calculate the partition sizes of the matrix problem
    n1 = no / 2
    n2 = no - n1

    if ( Npiv < max(n1,n2) ) then
       call die('Error in initialization of the pivoting &
            &array.')
    end if

    ! Point to the partitions
    sB1 = sA1 + n1 ** 2
    sC2 = sB1 + n1 * n2
    sA2 = sC2 + n1 * n2
    A1 => zwork(sA1:)
    B1 => zwork(sB1:)
    C2 => zwork(sC2:)
    A2 => zwork(sA2:)

    ! Copy over everything to preserve values for inversion
!$OMP parallel do default(shared), private(i)
    do i = 1, no ** 2
       zwork(i) = M(i)
    end do
!$OMP end parallel do
    
    ! Calculate Y2/B1
    call zgesv(n1,n2,A1(1),no,ipiv,C2(1),no,i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting Y2/B1')
       end if
    end if

    ! Calculate X1/C2
    call zgesv(n2,n1,A2(1),no,ipiv,B1(1),no,i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting X1/C2')
       end if
    end if

    ! Calculate the diagonal inverted matrix
    ! Calculate: A1 - X1
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n1,n1,n2, &
         zm1,M(sC2),no,B1(1),no,z1,M(sA1),no)

    ! Calculate: A2 - Y2
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n2,n2,n1, &
         zm1,M(sB1),no,C2(1),no,z1,M(sA2),no)
    
    call zgetrf(n1, n1, M(sA1), no, ipiv, i )
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on LU-decomposition A1')
       end if
    end if

    call zgetrf(n2, n2, M(sA2), no, ipiv, i )
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on LU-decomposition A2')
       end if
    end if

    ! Now before we use A1 as work array we will copy
    !  B1 and C2
    ! because they are used to calculate the off-diagonal
    ! Note it has to be done in this order
    !   A1, B1, C2, A2 is the order of matrices in the array
    call copy(n1,n2,zwork(1),zwork(sB1),no)
    call copy(n1,n2,zwork(n1*n2+1),zwork(sC2),no)
    j = n1*n2*2 + 1

    call zgetri(n1, M(1), no, ipiv, zwork(j), no**2-j, i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting A1')
       end if
    end if

    call zgetri(n2, M(sA2), no, ipiv, zwork(j), no**2-j, i)
    if ( i /= 0 ) then
       if ( present(ierr) ) then
          ierr = i
          return
       else
          call die('Error on inverting A2')
       end if
    end if

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: X1/C2 * M11
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n2,n1,n1, &
         zm1, zwork(1),n2,M(sA1),no,z0, M(sB1),no)

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: Y2/B1 * M22
#ifdef USE_GEMM3M
    call zgemm3m( &
#else
    call zgemm( &
#endif
         'N','N',n1,n2,n2, &
         zm1, zwork(n1*n2+1),n1,M(sA2),no,z0, M(sC2),no)

  contains

    subroutine copy(n1,n2,A,B,LDB)
      integer, intent(in) :: n1, n2, LDB
      complex(dp), intent(inout) :: A(n1,n2)
      complex(dp), intent(in) :: B(LDB,n2)

      integer :: i, j

      do j = 1 , n2
         do i = 1 , n1
            A(i,j) = B(i,j)
         end do
      end do
      
    end subroutine copy

  end subroutine mat_invert_recursive

end module m_mat_invert
