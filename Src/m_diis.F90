! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_diis

  use precision, only: dp
  use parallel,  only: Node
  use class_dData1D
  use class_Pair_dData1D
  use class_Fstack_Pair_dData1D
  !
  implicit none
  !
  private

  public :: diis
  !
CONTAINS
  !
  subroutine diis(stack,scalar_product,coeff)
  ! Compute the optimal extrapolation coefficients in the DIIS
  ! sense for a set of residual vectors.
  ! The residual vectors are stored in a finite stack of pairs
  ! of vectors. The second component of the pair is the residual
  ! This is not completely general for this routine, but it allows
  ! an easy bookeeping for the extrapolation.
  ! A more general version might accept an array of pointers to the
  ! residual vectors.

#ifdef MPI
    use mpi_siesta
#endif
    use fdf, only: fdf_get
    use m_svd,      only : solve_with_svd
    !
    type(Fstack_Pair_dData1D), intent(in) :: stack
    !
    ! The scalar product function is abstracted and
    ! passed as an argument. This function is really
    ! all we need to know to carry out the procedure.

    interface
       function scalar_product(a,b) result(sp)
         use precision, only: dp
         real(dp), intent(in), dimension(:) :: a, b
         real(dp)                           :: sp
       end function scalar_product
    end interface
    
    ! The extrapolation coefficients
    real(dp), intent(out), dimension(:)     :: coeff

    logical :: debug_diis 

    type(Pair_dData1D), pointer     :: pairp
    type(dData1D), pointer           :: vp
    real(dp), dimension(:), pointer :: diff_i, diff_j

    integer :: nmix, i, j, info

    logical            :: use_svd_in_diis = .false.
    real(dp)           :: rcond_svd_diis  = 1.0e-8_dp
    integer            :: rank

#ifdef MPI
    integer  MPIerror
#endif
    !
    real(dp), dimension(:,:), allocatable ::  b , bi
    real(dp), dimension(:), allocatable   ::  rhs, beta, sigma
    real(dp), dimension(:), allocatable   ::  buffer
    !
    debug_diis = fdf_get("DebugDIIS",.false.)
    nmix= n_items(stack)
    if (size(coeff) < nmix) then
       call die("coeff array too small")
    endif

    use_svd_in_diis = fdf_get("SCF.RhoG.DIIS.UseSVD",.true.)
    ! Note that 1.0e-6 seems too conservative
    rcond_svd_diis = fdf_get("SCF.Rhog.DIIS.RcondSVD",1.0e-8_dp)

    ! Allocate local arrays
    !
    allocate(b(nmix+1,nmix+1), bi(nmix+1,nmix+1))
    allocate(rhs(nmix+1),beta(nmix+1),sigma(nmix+1))
    allocate(buffer(nmix))
    !
    !  calculate mixing coefficients
    !
    do i=1,nmix
       pairp => get_pointer(stack,i)
       call secondp(pairp,vp)
       diff_i => val(vp)
       b(i,i) = scalar_product(diff_i,diff_i)
       !
       do j=1,i-1
          pairp => get_pointer(stack,j)
          call secondp(pairp,vp)
          diff_j => val(vp)
          b(i,j) = scalar_product(diff_i,diff_j)
          b(j,i) = b(i,j)
       enddo
       ! Now extend the matrix with ones in an extra column
       ! and row ...
       b(i,nmix+1)=1.0_dp   ! This value is really arbitrary
       b(nmix+1,i)=1.0_dp   ! This represents the Sum_{Ci} = 1 constraint
    enddo
    !      ! ... except in the extra diagonal entry
    b(nmix+1,nmix+1)=0.0_dp
    !
#ifdef MPI
    ! Global operations, but only for the non-extended entries
    ! To save in latency, a contiguous auxiliary array could be
    ! used, and perform a single 'reduce' call.
    do i=1,nmix
       call MPI_AllReduce(b(1:nmix,i),buffer,nmix,     &
            MPI_double_precision,            &
            MPI_sum,MPI_Comm_World,MPIerror)
       do j=1,nmix
          b(j,i)=buffer(j)
       enddo
    enddo
#endif

    if ((Node == 0) .and. debug_diis) call print_mat(b,nmix+1)

    !
    ! Get coefficients for DIIS mixing
    ! (Last column of the inverse matrix, corresponding to solving a
    ! linear system with (0,0,0,...,0,1) in the right-hand side)

    if (use_svd_in_diis) then
       rhs(1:nmix) = 0.0_dp
       rhs(nmix+1) = 1.0_dp

       call solve_with_svd(b,rhs,beta,info,rcond=rcond_svd_diis, &
                                      rank_out=rank,sigma=sigma)

       if (Node == 0 .AND. debug_diis) then
          print "(a,i2,7g12.5)", "SVD rank, s(i):", rank, sigma(:)
       endif
       if (info == 0) then
          coeff(1:nmix) = beta(1:nmix)
       else
          if (Node == 0) then
             write(6,"(a,i5)") "SVD failed - fallback to linear mixing"
          endif
          coeff(1:nmix-1) = 0.0_dp
          coeff(nmix) = 1.0_dp
       endif
    else
       call inverse(b,bi,nmix+1,nmix+1,info,debug_diis)
       !
       ! If inver was successful, get coefficients for DIIS mixing
       if (info .eq. 0) then
          do i=1,nmix
             coeff(i)=bi(i,nmix+1)
          enddo
       else
          ! Otherwise, use only last step
          if (Node == 0) then
            write(6,"(a)")  &
            "Warning: unstable inversion in DIIS - fallback to linear mixing"
          endif
          coeff(1:nmix-1) = 0.0_dp
          coeff(nmix) = 1.0_dp
       endif

    endif ! SVD
    !
    if ((Node == 0) .and. debug_diis) then
       do i = 1, n_items(stack)
          print "(a,i2,f12.5)", "DIIS coeff - ", i, coeff(i)
       enddo
    endif

    ! Deallocate local arrays
    !
    deallocate(b,bi,buffer,beta,sigma,rhs)
    !
  CONTAINS

!---------------------------------------------------
    SUBROUTINE inverse(A,B,N,NDIM,INFO,debug_inverse)

      IMPLICIT NONE
      INTEGER, intent(in) ::  N,NDIM
      real(dp), intent(in)  ::  A(NDIM,NDIM)
      real(dp), intent(out) ::  B(NDIM,NDIM)
      integer, intent(out) :: info
      logical, intent(in)  :: debug_inverse

      real(dp)  ::  X

!!$C Routine to obtain the inverse of a general, nonsymmetric
!!$C square matrix, by calling the appropriate LAPACK routines
!!$C The matrix A is of order N, but is defined with a 
!!$C size NDIM >= N.
!!$C If the LAPACK routines fail, try the good old Numerical Recipes
!!$C algorithm
!!$C
!!$C P. Ordejon, June 2003
!!$C
!!$C **** INPUT ****
!!$C A(NDIM,NDIM)   real*8      Matrix to be inverted
!!$C N              integer     Size of the matrix
!!$C NDIM           integer     Defined size of matrix
!!$C **** OUTPUT ****
!!$C B(NDIM,NDIM)   real*8      Inverse of A
!!$C INFO           integer     If inversion was unsucessful, 
!!$C                            INFO <> 0 is returned
!!$C                            Is successfull, INFO = 0 returned
!!$C ***************


      real(dp) ::  WORK(N),C,ERROR,DELTA,TOL, pm(N,N)
      INTEGER IPIV(N)
      INTEGER I,J,K

      logical :: debug

      debug = debug_inverse .and. (node == 0)

      TOL=1.0D-4
      INFO = 0

      DO I=1,N
      DO J=1,N
        B(I,J)=A(I,J)
      ENDDO
      ENDDO

      CALL DGETRF(N,N,B,NDIM,IPIV,INFO)

      IF (INFO .NE. 0) THEN
       if (debug) print *,   &
           'inver:  ERROR: DGETRF exited with error message', INFO
        GOTO 100
      ENDIF

      CALL DGETRI(N,B,NDIM,IPIV,WORK,N,INFO)

      IF (INFO .NE. 0) THEN
       if (debug) print *,   &
          'inver:  ERROR: DGETRI exited with error message', INFO
        GOTO 100
      ENDIF

! CHECK THAT THE INVERSE WAS CORRECTLY CALCULATED

      pm = matmul(a(1:n,1:n),b(1:n,1:n))

      ERROR=0.0D0
      DO I=1,N
      DO J=1,N
        C=0.0D0
        DO K=1,N
          C=C+A(I,K)*B(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
        C=0.0D0
        DO K=1,N
          C=C+B(I,K)*A(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
      ENDDO
      ENDDO

      IF (ERROR/N .GT. TOL) THEN
       if (debug) then
          print *,   &
          'inver:  ERROR in lapack inverse. Error: ', error
          call print_mat(a,n)
          call print_mat(b,n)
          call print_mat(pm,n)
       endif
        GOTO 100
      ENDIF

      INFO = 0
      RETURN

100   CONTINUE

! Try simple, old algorithm:

      DO I=1,N
        DO J=1,N
          B(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,N
        IF (B(I,I) .EQ. 0.0D0) THEN
           if (debug) print *,   &
            'inver:  zero pivot in fallback algorithm'
          INFO = -777
          RETURN
        ENDIF
        X=B(I,I)
        B(I,I)=1.0d0
        DO J=1,N
          B(J,I)=B(J,I)/X
        ENDDO
        DO K=1,N
          IF ((K-I).NE.0) THEN 
            X=B(I,K)
            B(I,K)=0.0d0
            DO J=1,N
              B(J,K)=B(J,K)-B(J,I)*X
            ENDDO
          ENDIF
        ENDDO
      ENDDO

! CHECK THAT THE INVERSE WAS CORRECTLY CALCULATED

      pm = matmul(a(1:n,1:n),b(1:n,1:n))
      ERROR=0.0D0
      DO I=1,N
      DO J=1,N
        C=0.0D0
        DO K=1,N
          C=C+A(I,K)*B(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
        C=0.0D0
        DO K=1,N
          C=C+B(I,K)*A(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
      ENDDO
      ENDDO

      IF (ERROR/N .GT. TOL) THEN
       if (debug) then
          print *,   &
            'inver:  INVER unsuccessful, ERROR = ', ERROR
          call print_mat(a,n)
          call print_mat(b,n)
          call print_mat(pm,n)
       endif
        INFO = -888
        RETURN
      ENDIF

      INFO = 0
      RETURN

    end SUBROUTINE inverse

    subroutine print_mat(a,n)
      integer, intent(in)  :: n
      real(dp), intent(in) :: a(n,n)
      integer :: i

      print *, "mat:"
      do i = 1, n
         print "(6g15.7)", (a(i,j),j=1,n)
      enddo
      print *, "------"
    end subroutine print_mat

  end subroutine diis
  !
end module m_diis

