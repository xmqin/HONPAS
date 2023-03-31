! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_pulay

  use precision, only: dp
  use m_spin,    only: h_spin_dim, SpOrb
  implicit none

  private

  save

  real(dp), pointer :: auxpul(:,:)  => null()
  integer           :: n_records_saved = 0
  real(dp)          :: alpha_pulay
  real(dp)          :: alpha_before_pulay
  real(dp)          :: alpha_after_pulay
  real(dp)          :: max_dmax_for_pulay
  integer           :: pulay_minimum_history = 2
  logical           :: last_was_pulay = .false.
  integer           :: LM_before_pulay = 0
  logical           :: LM_after_pulay = .false.
  logical           :: use_svd_in_pulay = .false.
  logical           :: debug_svd_in_pulay = .true.
  real(dp)          :: rcond_svd_pulay = 1.0e-8_dp

  public :: pulayx, init_pulay_arrays, resetPulayArrays

CONTAINS

  subroutine init_pulay_arrays()
    use siesta_options, only: maxsav
    use alloc
    use atomlist, only: no_l
    use m_spin, only: nspin
    use sparse_matrices, only: numh
    use fdf

    integer :: ntmp, nauxpul

    n_records_saved = 0

    if (maxsav .le. 0) then
       ! No need for auxiliary arrays
    else
       nauxpul = sum(numh(1:no_l)) * h_spin_dim * maxsav
       !
       call re_alloc(auxpul,1,nauxpul,1,2,name="auxpul",        &
            routine="pulay")
    endif
  
    ! Read some operational parameters

    alpha_pulay = fdf_get("SCF.Pulay.Damping",-1.0_dp) ! Default will set it to alpha
    alpha_after_pulay = fdf_get("SCF.MixingWeightAfterPulay",0.5_dp)

    alpha_after_pulay = fdf_get("SCF.Pulay.MixingWeight.After",alpha_after_pulay)
    LM_after_pulay = fdf_get("SCF.LinearMixingAfterPulay",.false.)
    LM_after_pulay = fdf_get("SCF.Pulay.LinearMixing.After",LM_after_pulay)

    alpha_before_pulay = fdf_get("SCF.Pulay.MixingWeight.Before",0.01_dp)
    LM_before_pulay = fdf_get("SCF.Pulay.LinearMixing.Before",0)

    pulay_minimum_history = fdf_get("SCF.PulayMinimumHistory",2)
    max_dmax_for_pulay = fdf_get("SCF.PulayDmaxRegion",1000.0_dp) ! No effect by default
    use_svd_in_pulay = fdf_get("SCF.Pulay.UseSVD",.false.)
    debug_svd_in_pulay = fdf_get("SCF.Pulay.DebugSVD",.true.)
    ! Note that 1.0e-6 seems too conservative
    rcond_svd_pulay = fdf_get("SCF.Pulay.RcondSVD",1.0e-8_dp)

  end subroutine init_pulay_arrays
  !---------------------------------------------------------------------

  subroutine pulayx( iscf,mix1,no_l,maxnd,numd,             &
       listdptr,nspin,maxsav,alpha,nkick,alphakick,   &
       dmnew,dmold,dmax)

    ! Pulay mixing implemented to accelerate the self-consistency
    ! Mixes MAXSAV previous steps.  Linear mixing if MAXAV =< 0
    ! 
    ! Written by In-Ho Lee, Beckman Inst., Univ. of Illinois, Mar. 25 '97
    ! Modified and partly re-written by P. Ordejon, July'97
    ! Modified and optimized by P. Ordejon, November'97
    ! Modified by A. Garcia, April-June 2010
    !
    ! ************************** INPUT **************************************
    ! integer iscf               : Current SCF iteration
    ! logical mix1               : Mix first SCF step (T or F)
    ! integer no_l               : Number of atomic orbitals stored locally
    ! integer maxnd              : First dimension of D.M., and 
    !                              maximum number of nonzero elements of D.M.
    ! integer numd(no_l)         : Control vector of D.M.
    !                              (number of nonzero elements of each row)
    ! integer listdptr(no_l)     : Pointer to start of rows in listd
    ! integer nspin              : Spin polarization (1= no spin, 2=polarized)
    ! integer maxsav             : Number of history records for Pulay mixing.
    ! real*8 alpha               : Mixing parameter (for linear mixing)
    ! integer nkick              : Do a linear mixing (kick) each nkick cycles
    ! real*8 alphakick           : Mixing parameter for kick cycles
    !
    ! ********************* INPUT AND OUTPUT*********************************
    ! real*8 dmnew(maxnd)        : Density Matrix
    !                              Input: d.m. output in current SCF step
    !                              Output: d.m. input for next SCF iteration
    ! real*8 dmold(maxnd)        : Density matrix
    !                              Input: d.m. input in current SCF step
    !                              Output: d.m. input for next SCF iteration
    ! ************************** OUTPUT *************************************
    ! real*8 dmax                : Maximum change of a DM element between 
    !                              input and output
    ! ************************ BEHAVIOUR ************************************
    !
    ! Algorithm changed!! (2/19/99)
    ! Pulay mixing is now done EVERY iteration, except the first one
    ! (which is done with linear mixing or with no mixing, according
    ! to logical variable mix1)
    ! Mixes the last maxsav iterations (normal linear mixing if maxsav =< 1)
    !
    ! The mixing is done the following way (Anderson, Pulay):
    !
    ! Modified input and output matrices from the former step are obtained
    ! by mixing maxsav prefious steps:
    !
    !   D'_in (n) = Sum_i=1,maxsav  beta_i D_in (n-maxsav+i)
    !   D'_out(n) = Sum_i=1,maxsav  beta_i D_out(n-maxsav+i)
    !
    ! The beta coefficients are obtained by minimizing the norm between
    ! D'_in and D'_out.
    ! The input charge for step (n+1) is done by simple mixing the D's
    !
    !   D_in (n+1) = alpha D'_out(n) + (1-alpha) D_in(n)
    !              = Sum_i=1,maxsav D_in (n-maxsav+i) +
    !                alpha Sum_i=1,maxsav deltaD(n-maxsav+i)
    !
    ! where deltaD(n) is the residual of step n.
    !
    ! The density matrices of BOTH spins are mixed at the same
    ! time, with the same coefficients (to ensure conservation of
    ! total number of electrons).
    ! (and for NONCOL?)
    ! ***********************************************************************
    !
    !
    !  Modules
    !
    use precision,  only : dp
    use parallel,   only : Node
    use sys,        only : die
    use alloc
    use siesta_options, only: avoid_first_after_kick
    use fdf,        only : fdf_get
    use m_svd,      only : solve_with_svd
#ifdef MPI
    use mpi_siesta
    use m_mpi_utils, only: globalize_max
#endif

    integer , intent(in) ::  iscf,maxsav,maxnd, no_l,nkick,nspin

    integer , intent(in) :: numd(*),listdptr(*)
    logical, intent(in) ::  mix1

    real(dp), intent(in)    ::  alpha,alphakick
    real(dp), intent(inout) ::  dmnew(maxnd,nspin),dmold(maxnd,nspin)
    real(dp), intent(out)   :: dmax

    real(dp), pointer :: savedm(:), saveres(:)
    real(dp), pointer :: res1(:), res2(:)

    ! Internal variables ....................................................
    integer :: i0,i,ii,in,is,isite,j,jj,numel,ind,info,maxmix
    logical :: after_kick
    logical :: no_new_information
    logical :: debug_inverse = .false.

!    logical, save :: kick_due = .false.
    integer :: rank

#ifdef MPI
    integer  MPIerror
    real(dp) :: buffer1
#endif
    real(dp) :: ssum, cur_alpha
    real(dp), dimension(:,:), pointer ::  b, bi
    real(dp), dimension(:), pointer   ::  buffer
    real(dp), dimension(:), pointer   ::  coeff, sigma, rhs, beta

#ifdef DEBUG_PULAY_INVERSE
    ! Typically, inversion problems in the DIIS procedure
    ! occur when the error estimate for a step is very small
    ! compared to the previous ones. This is harmless, but
    ! this hook would allow a full diagnosis if compiled in.
    debug_inverse = .true.
#endif

    if (maxsav > 1) then
       savedm => auxpul(:,1)          ! DM in former iterations
       saveres => auxpul(:,2)         ! Residuals in former iterations
       !
       ! Check some input and dimensions 
       numel = sum(numd(1:no_l))
       call assert(numel == maxnd, "sum(numd) /= maxnd")
       numel = numel*nspin
       call assert(numel*maxsav == size(savedm), &
            "maxnd*spin*maxsav /= savedm")

       if (size(saveres) .lt. numel*maxsav) &
            call die('pulayx: dimensions are too small')
    endif

    ! For the purposes of history keeping, the very first iteration
    ! is treated as the first step after a kick: a possibly large
    ! jump in the DM

    ! The following code implements this in a safe way (avoid mod if nkick==0)
    ! after_kick = (iscf == 1 .OR.     &
    !             (nkick > 0 .AND. mod(iscf,nkick) == 1) ) 
    after_kick = .false.
    if (iscf == 1) then
       after_kick = .true.
    else if (nkick > 0) then
       if (mod(iscf,nkick)==1) after_kick = .true.
    endif

    !  Compute current maximum deviation ...........
    dmax = 0.0_dp
    do is = 1,nspin
       do ind = 1 , maxnd
          dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
       end do
    enddo

#ifdef MPI
    ! Ensure that dmax is the same on all nodes for convergence test/output
    call globalize_max(dmax,buffer1)
    dmax = buffer1
#endif

    if ( (after_kick .and. avoid_first_after_kick) .or. &
         (last_was_pulay .and. LM_after_pulay) .or. &
         (dmax > max_dmax_for_pulay) ) then

       ! Do not keep the residual of the very first iteration, or of 
       ! the first iteration after a kick, if so instructed
       ! Or if the deviation is too large
       no_new_information = .true.
    else
       if ( iscf <= LM_before_pulay ) then
          ! ensure that we overwrite the old
          ! history for bookkeeping
          n_records_saved = 0
       end if

       call WriteCurrentDinAndResidualInStore()
       no_new_information = .false.
       ! updates n_records_saved
    endif

    if ( last_was_pulay ) then
       last_was_pulay = .false.
       if ( LM_after_pulay ) then
          if (Node==0) write(6,*) "Linear mixing after Pulay step"
          call linear_mixing(alpha_after_pulay)
          RETURN
       endif
    endif

    ! Perform linear mixing if we do not have new information
    ! or not enough history
    ! (by default, 2, but it is settable by the user)
    !
    if ( no_new_information .or. &
         n_records_saved < pulay_minimum_history .or. &
         iscf <= LM_before_pulay ) then

       ! Could do simply:  
       ! if (iscf > 1 .or. mix1) DMnew = (1.0_dp-alpha)*DMold + alpha*DMnew
       ! DMold = DMnew

       if (iscf .gt. 1 .or. mix1) then
          if ( iscf <= LM_before_pulay ) then
             call linear_mixing(alpha_before_pulay)
          else
             call linear_mixing(alpha)
          end if
       endif

       RETURN

    endif

    ! Perform linear mixing if iscf = N x nkick
    !
    if ( nkick > 0 ) then
      if ( mod(iscf,nkick).eq.0) then
        ! Reset the history information
        if (Node == 0) then
          write(6,"(a,i4)") "Applying kick at iscf: ", iscf
          write(6,"(a,i4)") "Re-setting Pulay history"
        endif
        n_records_saved = 0
        
        call linear_mixing(alphakick)
        RETURN
      end if
    endif


    ! Perform Pulay mixing if we have reached this point
    !
    ! sanity check
    call assert(n_records_saved > 1, "N_SCF_records <= 1")
    ! Use only as many records as we have in the history store
    maxmix = min(n_records_saved, maxsav)

    ! Respect historical value for alpha_pulay if not specified
    if ( alpha_pulay < 0.0_dp ) then
       cur_alpha = alpha
    else
       cur_alpha = alpha_pulay
    endif

    ! Allocate local arrays
    !
    nullify( b )
    call re_alloc( b, 1, maxmix+1, 1, maxmix+1, name='b', &
         routine='pulayx' )
    nullify( bi )
    call re_alloc( bi, 1, maxmix+1, 1, maxmix+1, name='bi', &
         routine='pulayx' )
    nullify( coeff )
    call re_alloc( coeff, 1, maxmix, name='coeff', &
         routine='pulayx' )
    nullify( rhs )
    call re_alloc( rhs, 1, maxmix+1, name='rhs', &
         routine='pulayx' )
    nullify( beta )
    call re_alloc( beta, 1, maxmix+1, name='beta', &
         routine='pulayx' )
    nullify( sigma )
    call re_alloc( sigma, 1, maxmix+1, name='sigma', &
         routine='pulayx' )
    nullify( buffer )
    i = maxmix + 1
    call re_alloc( buffer, 1, i*maxmix, name='buffer', &
         routine='pulayx' )

    !  calculate mixing coefficients
    !
    do i = 1 , maxmix

       ! Get i'th residual array
       i0 = (i-1) * numel
       do is = 1,nspin
          do ii = 1,no_l
             do jj = 1,numd(ii)
                ind = listdptr(ii) + jj
                i0 = i0 + 1
                dmnew(ind,is) = saveres(i0)
             enddo
          enddo
       enddo
       !
       ! B(i,i) = dot_product(Res(i)*Res(i))
       b(i,i) = 0.0_dp
       ssum=0.0_dp
! CC RC Added for on-site SO
       if ( .not. SpOrb ) then
        do is=1,nspin
         do ii=1,no_l
          do jj=1,numd(ii)
           ind = listdptr(ii) + jj
           ssum=ssum+dmnew(ind,is)*dmnew(ind,is)
          enddo
         enddo
        enddo
       elseif ( SpOrb ) then
        do ii=1,no_l
         do jj=1,numd(ii)
          ind = listdptr(ii) + jj
          ssum=ssum+dmnew(ind,1)*dmnew(ind,1)
          ssum=ssum+dmnew(ind,2)*dmnew(ind,2)
          ssum=ssum+  &
          0.5*(dmnew(ind,3)*dmnew(ind,3)+dmnew(ind,7)*dmnew(ind,7))
          ssum=ssum+  &
          0.5*(dmnew(ind,4)*dmnew(ind,4)+dmnew(ind,8)*dmnew(ind,8))
         enddo
        enddo
       endif
! CC RC

       b(i,i)=ssum
       !
       do j=1,i-1

          i0 = (j-1) * numel
          do is = 1,nspin
             do ii = 1,no_l
                do jj = 1,numd(ii)
                   ind = listdptr(ii) + jj
                   i0 = i0 + 1
                   dmold(ind,is) = saveres(i0)
                enddo
             enddo
          enddo

          !          ! B(i,j) = B(j,i) = dot_product(Res(i)*Res(j))

          b(i,j)=0.0_dp
          ssum=0.0_dp
          if ( .not. SpOrb ) then
           do is=1,nspin
            do ii=1,no_l
             do jj=1,numd(ii)
              ind = listdptr(ii) + jj
              ssum=ssum+dmold(ind,is)*dmnew(ind,is)
             enddo
            enddo
           enddo
          elseif ( SpOrb ) then
           do ii=1,no_l
            do jj=1,numd(ii)
             ind = listdptr(ii) + jj
             ssum=ssum+dmold(ind,1)*dmnew(ind,1)
             ssum=ssum+dmold(ind,2)*dmnew(ind,2)
             ssum=ssum+  &
             0.5*(dmnew(ind,3)*dmold(ind,3)+dmold(ind,7)*dmnew(ind,7))
             ssum=ssum+  &
             0.5*(dmnew(ind,4)*dmold(ind,4)+dmold(ind,8)*dmnew(ind,8))
            enddo
           enddo
          endif
          b(i,j)=ssum
          b(j,i)=ssum
       enddo

       ! Now extend the matrix with ones in an extra colum
       ! and row ...
       b(i,maxmix+1) = 1.0_dp
       b(maxmix+1,i) = 1.0_dp

    end do
    
    ! ... except in the extra diagonal entry
    b(maxmix+1,maxmix+1) = 0.0_dp

#ifdef MPI
    ! Global operations, but only for the non-extended entries
    i = maxmix + 1
    call MPI_AllReduce(b(1,1),buffer(1),i*maxmix, &
         MPI_double_precision, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    i0 = 0
    do j = 1 , maxmix
       do i = 1 , maxmix
          i0 = i0 + 1
          b(i,j) = buffer(i0)
       end do
       i0 = i0 + 1 ! skipping non-extended entries
    end do
#endif

    if (use_svd_in_pulay) then
       rhs(1:maxmix) = 0.0_dp
       rhs(maxmix+1) = 1.0_dp

       call solve_with_svd(b,rhs,beta,info,rcond=rcond_svd_pulay, &
            rank_out=rank,sigma=sigma)

       if (Node == 0 .AND. debug_svd_in_pulay) then
          print "(a,i2,7g12.5)", "SVD rank, s(i):", rank, sigma(:)
       endif
       if (info == 0) then
          coeff(1:maxmix) = beta(1:maxmix)
       else
          if (Node == 0) then
             write(6,"(a,i5)") "SVD failed - fallback to linear mixing"
          endif
          coeff(1:maxmix-1) = 0.0_dp
          coeff(maxmix) = 1.0_dp
       endif
          
       ! Use a more sophisticated history of rank degradation...
       !kick_due =  (rank < maxmix + 1)
    else
       call inverse(b,bi,maxmix+1,maxmix+1,info,debug_inverse)
       !
       ! If inver was successful, get coefficients for Pulay mixing
       if (info .eq. 0) then
          do i= 1 , maxmix
             coeff(i) = bi(i,maxmix+1)
          enddo
       else
          ! Otherwise, use only last step
          if (Node == 0) then
            write(6,"(a)")  &
            "Warning: unstable inversion in Pulayx - fallback to linear mixing"
          endif
          do i=1,maxmix
             coeff(i)=0.0_dp
          enddo
          j=mod(iscf,maxmix)
          if(j.eq.0) j=maxmix
          coeff(j) = 1.0_dp
       endif

    endif ! SVD

    ! Read former matrices for mixing .........
    !
    dmnew(1:maxnd,1:nspin) = 0._dp
    do i = 1 , maxmix
       i0 = (i-1) * numel
       do is = 1 , nspin
          res1 => savedm (i0+1:i0+maxnd)
          res2 => saveres(i0+1:i0+maxnd)
          i0 = i0 + maxnd
!$OMP parallel do default(shared), &
!$OMP&private(ind)
          do ind = 1 , maxnd
             dmnew(ind,is) = dmnew(ind,is) + coeff(i) * &
                  ( res1(ind) + cur_alpha * res2(ind) )
          end do
!$OMP end parallel do
       end do
    end do
    dmold(1:maxnd,1:nspin) = dmnew(1:maxnd,1:nspin)

    if (LM_after_pulay) last_was_pulay = .true.

    ! Deallocate local arrays
    call de_alloc( b, name='b', routine='pulayx' )
    call de_alloc( bi, name='bi', routine='pulayx' )
    call de_alloc( coeff, name='coeff', routine='pulayx' )
    call de_alloc( rhs, name='rhs', routine='pulayx' )
    call de_alloc( beta, name='beta', routine='pulayx' )
    call de_alloc( sigma, name='sigma', routine='pulayx' )
    call de_alloc( buffer, name='buffer', routine='pulayx' )

  CONTAINS

    subroutine linear_mixing(alp)
      real(dp), intent(in) :: alp
      real(dp) :: alp1m

      alp1m = 1._dp - alp
      do is = 1 , nspin
!$OMP parallel do default(shared), private(ind)
         do ind = 1 , maxnd
            dmnew(ind,is) = &
                 alp1m * dmold(ind,is) + alp * dmnew(ind,is)
            dmold(ind,is) = dmnew(ind,is)
         end do
!$OMP end parallel do
      end do

    end subroutine linear_mixing
    
    subroutine WriteCurrentDinAndResidualInStore()

      ! The store is really a circular array of size maxsav. It is
      ! implemented inline here, but it could be done more cleanly
      ! with a finite circular stack.

      if (maxsav > 1) then
         n_records_saved = n_records_saved + 1
         ! isite marks the point to write 
         isite = mod(n_records_saved,maxsav)
         if (isite .eq. 0) isite = maxsav
         i0 = (isite-1) * numel
         if (Node == 0) then
            !write(6,"(a,i4)") "Saving scf history record. iscf: ", iscf
         endif

         ! Both spins are written, one block after the other,
         ! but the position on the array will not be conmmensurate
         ! for the different processors...
         ! ... that is why it is difficult to parallelize the on-file version
   
         do is = 1 , nspin
            res1 => savedm(i0+1:i0+maxnd)
            res2 => saveres(i0+1:i0+maxnd)
            i0 = i0 + maxnd
!$OMP parallel do default(shared), private(ind)
            do ind = 1 , maxnd
               res1(ind) = dmold(ind,is)
               res2(ind) = dmnew(ind,is) - dmold(ind,is)
            end do
!$OMP end parallel do
         end do
      endif

    end subroutine WriteCurrentDinAndResidualInStore

    subroutine assert(condition, message)
    ! A simple assertion checker
      logical, intent(in)          :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         if (Node == 0) then
            write(6,"(2a)") "** Assertion failed **: ", message
            write(0,"(2a)") "** Assertion failed **: ", message
         endif
      endif
    end subroutine assert
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
  end subroutine pulayx
  !
  subroutine resetPulayArrays( )
    use alloc, only: de_alloc
    implicit none
    if (associated(auxpul))      &
      call de_alloc( auxpul, 'auxpul', 'pulay' )
  end subroutine resetPulayArrays


end module m_pulay

