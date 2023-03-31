! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.

! Module for performing diagonalization of both real symmetric
! and complex hermitian matrices.
!
! Its current implementation has been developed by:
!   Nick R. Papior, 2017.
!
! Its current structure is heavily based on the original
!   rdiag and cdiag
! routines and has been optimized and shortened.
!
! The basis of this routine is to provide diagonalization
! and it differs from the original routines in these respects:
!  1. It does NOT implement the generalized algorithms.
!     Every diagonalization is performed in a 2-step manner,
!     a) convert to generalized form (*potrf, *gst),
!     b) solve standard problem
!     This leverages the requirement of the external library
!     to implement the generalized forms (which does
!     exactly the same).
!  2. All pre-rotate and save-eigenvector options has been removed.
!     There is no point in having them if they were not used.
!     Using them should require changes to this.
!  3. In 2D parallel the original implementations *always*
!     allocated the equivalent 2D distribution arrays (H2D, S2D, Z2D).
!     However, when the node's distributed 2D number of elements
!     is <= size(H) then we do not need to allocate room for these
!     nodes. I.e. we can save 3 * size(H2D) elements in memory!
!     This will typically be so, and it may even happen that only a subset
!     of the nodes actually require to allocate additional memory.
!  4. Change the default symmetric part from 'Upper' to 'Lower'.
!     This is because the pzhengst (in ScaLAPACK) performs better
!     for 'Lower' (see pdyengst/pzhengst documentations).
!  5. Added more solvers
!     - MRRR
!     - 2stage solvers (only in LAPACK, and currently only for jobz='N')
!     - no-expert drivers (regular ones)
!  6. Allowed only calculating a subset of the eigenvalues.
!     There are 3 options for neig:
!      neig < 0:
!       calculate -neig eigenvalues (jobz == 'N')
!      neig > 0:
!       calculate neig eigenvalues and eigenvectors (jobz == 'V')
!      neig == 0:
!       same as calculating all eigenvalues (jobz == 'N')
!      
module m_diag

#ifndef MPI
! Make sure ELPA is only used for MPI
# undef SIESTA__ELPA
#endif

  use precision
  use parallel, only : Node, Nodes, IONode
  
  implicit none

  private

  save

#ifdef MPI
  ! The BLACS context
  integer :: iCTXT = -1
  integer :: iCTXT2D = -1
#endif
  
#ifdef _DIAG_WORK
  logical :: diag_work_r(2) = .true.
  logical :: diag_work_c(2) = .true.
#endif

  ! Initialization and exit routines for the blacs communicators
  ! Note that the real and complex algorithms share the blacs communicators.
  public :: diag_init, diag_exit

  ! Diagonalization routines (real and complex)
  public :: diag_r, diag_c

contains

  subroutine diag_init()

    use m_diag_option, only: Serial, ParallelOverK

#ifdef MPI
    use m_diag_option, only: Use2D, ProcessorY
    use mpi_siesta, only: MPI_Comm_World
# ifdef SIESTA__ELPA
    use m_diag_option, only: algorithm, ELPA_1stage, ELPA_2stage
    use elpa, only: elpa_initialized, elpa_init, ELPA_OK
# endif
#endif

#ifdef MPI
    integer :: nr, nc
#endif

    call diag_exit()

#ifdef _DIAG_WORK
    diag_work_r = .true.
    diag_work_c = .true.
#endif

    if ( Serial ) return
    if ( ParallelOverK ) return
    
#ifdef MPI
    
    ! Create a new context
    iCTXT = MPI_Comm_World
    nr = 1
    nc = Nodes
    call blacs_gridinit( iCTXT, 'C', nr, nc )

    if ( Use2D ) then

       ! Setup secondary grid for the 2D distribution
       nr = ProcessorY
       nc = max(1, Nodes/ProcessorY)
       call blacs_get(iCTXT, 10, iCTXT2D)
       call blacs_gridinit(iCTXT2D, 'R', nr, nc)

    end if
#endif

#ifdef SIESTA__ELPA
    if ( algorithm == ELPA_1stage .or. &
         algorithm == ELPA_2stage ) then

       ! Check whether ELPA is initialized
       if ( elpa_initialized() /= ELPA_OK ) then
          if ( elpa_init(20170403) /= ELPA_OK ) then
             call die('diag: ELPA initialization could not use API 20170403')
          end if
       end if

    end if
#endif
    
  end subroutine diag_init

  subroutine diag_exit()
#ifdef SIESTA__ELPA
    use elpa, only: elpa_initialized, elpa_uninit, ELPA_OK
#endif

#ifdef SIESTA__ELPA
    integer :: info
    ! Check whether ELPA is initialized, uninitialize it if it is.
    if ( elpa_initialized() == ELPA_OK ) then
       call elpa_uninit(info)
       call elpa_check(info, 'uninit')
    end if
#endif

#ifdef MPI
    if ( iCTXT >= 0 ) then
       call blacs_gridexit(iCTXT)
       iCTXT = -1
    end if

    if ( iCTXT2D >= 0 ) then
       call blacs_gridexit(iCTXT2D)
       iCTXT2D = -1
    end if
#endif

  end subroutine diag_exit

  subroutine diag_correct_input(algo, jobz, range, uplo, trans, neig, n)

    use m_diag_option
    
    integer, intent(inout) :: algo
    character, intent(out) :: jobz, range, uplo, trans
    integer, intent(inout) :: neig
    integer, intent(in) :: n

    ! Set general LAPACK/ScaLAPACK parameters
    if ( neig > 0 ) then
       ! both eigenvalues and eigenvectors
       jobz = 'V'
       
    else if ( neig < 0 ) then
       jobz = 'N'
       ! only a subset of eigenvalues
       neig = -neig
       
    else if ( neig == 0 ) then
       jobz = 'N'
       ! correct to all eigenvalues
       neig = n
       
    end if

    ! Set range to 'A' for all values
    if ( neig == n ) then
       range = 'A'
    else
       range = 'I'
    end if

    ! Use user option
    uplo = UpperLower
    if ( uplo == 'U' ) then
       trans = 'N'
    else
       trans = 'C'
    end if

    
    ! Correct for special routines
    if ( Serial ) then

       if ( jobz == 'V' ) then
          select case ( algo )
          case ( DivideConquer_2stage )
             algo = DivideConquer
          case ( MRRR_2stage )
             algo = MRRR
          case ( Expert_2stage )
             algo = Expert
          case ( QR_2stage )
             algo = QR
          end select
       end if

#ifdef MPI
    else

       ! These checks will only happen where it starts serial
       ! but then goes to parallel
       select case ( algo )
       case ( DivideConquer_2stage )
          algo = DivideConquer
       case ( MRRR_2stage )
          algo = MRRR
       case ( Expert_2stage )
          algo = Expert
       case ( QR_2stage )
          algo = QR
       end select

       if ( algo == DivideConquer ) then
          ! regardless of neig
          jobz = 'V'
          range = 'A'
          neig = n
       end if
       
#endif
    end if

  end subroutine diag_correct_input


  subroutine diag_c( H, S, n, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
! ***************************************************************************
! Subroutine to solve all eigenvalues and eigenvectors of the
! complex general eigenvalue problem  H z = w S z,  with H and S
! complex hermitian matrices.
! Written by G.Fabricius and J.Soler, March 1998
! Rewritten by Julian Gale, August 2004
! Rewritten by Nick R. Papior, July 2017
! ************************** INPUT ******************************************
! complex*16 H(nml,nm)             : Hermitian H matrix
! complex*16 S(nml,nm)             : Hermitian S matrix
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neig                     : No. of eigenvalues(<0)/vectors(>0) to calculate
! integer iscf                     : SCF cycle
! integer BlockSize                : Effective parallel block size
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! complex*16 Z(nml,nm)             : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ************************* PARALLEL ****************************************
! When running in parallel this routine now uses Scalapack to perform a
! parallel matrix diagonalisation. This requires Scalapack and Blacs to
! be installed first. Although globally a 1-D block cyclic data distribution
! is employed, locally 1 or 2-D distributions are allowed for.
! The blocksize is now explicitly passed to the routine (A. Garcia, 2016)      
! The routine allows access to all the phases of diagonalisation for fuller
! control, and allows for parallel divide and conquer with reduced memory.
! The presence of eigenvalue clusters is checked for and the memory adjusted
! accordingly to try to guarantee convergence.
! Note that the blacs grid is only initialised on the first call to avoid
! exceeding MPI limits for split/group creation.
!
! When running in parallel and using a 2D distribution it is crucial that the
! Z array is also allocated.
! When using a 2D distribution the algorithm tries to reuse memory without
! double allocation. This is allowed if the 2D distribution # of elements
! is less than or equal to the number of elements in the local arrays.
! (NOTE this is typically so).
! 
!
! It allows for the following routines:
!   LAPACK (Serial or ParallelOverK):
!    - zheevd, zheevd_2stage
!    - zheevr, zheevr_2stage
!    - zheevx, zheevx_2stage
!    - zheev, zheev_2stage
!   ScaLAPACK (Parallel):
!    - pzheevd
!    - pzheevr
!    - pzheevx
!    - pzheev
! To enable the *vr* routines this file needs to be compiled with:
!   -DSIESTA__MRRR
! To enable the *v*_2stage routines this file needs to be compiled with:
!   -DSIESTA__DIAG_2STAGE
! Note that the LAPACK implementation of the 2stage solvers does (in 3.7.1) not
! implement the jobz='V' stage.
! ***************************************************************************

    ! Modules
    use m_diag_option ! just all
    
    use alloc
    use sys, only : die
#ifdef SIESTA__ELPA
    use mpi_siesta, only: MPI_Comm_World
    use elpa
#endif

    ! Passed variables
    integer, intent(out) :: ierror
    integer, intent(in) :: iscf
    integer, intent(in) :: n
    integer, intent(in) :: neig
    integer, intent(in) :: nm, nml
    integer, intent(in) :: BlockSize
    real(dp), intent(inout) :: w(nml)
    complex(dp), intent(inout), target :: H(nml,nm)
    complex(dp), intent(inout), target :: S(nml,nm)
    complex(dp), intent(inout), target :: Z(nml,nm)

    ! Local variables
    type(allocDefaults) :: oldDefaults
    integer :: algo
#ifdef MPI

    ! Scale when transforming to generalized eigenvalue problem
    real(dp) :: scale

    ! Expert and MRRR drivers
    integer :: nz, mclustr
    integer,  pointer :: iclustr(:) => null()
    real(dp), pointer :: gap(:) => null()

    ! BLACS descriptors
    integer, target :: desc_1d(9), desc_2d(9)
    integer, pointer :: desc(:) => null()

    ! Additional variables for a 2D parallel grid
    integer :: np2d(2), my2d(2), mat_2d(2)

    ! Pointers of H, S and Z
    complex(dp), pointer :: Hp(:,:) => null()
    complex(dp), pointer :: Sp(:,:) => null()
    complex(dp), pointer :: Zp(:,:) => null()

# ifdef SIESTA__ELPA
    class(elpa_t), pointer :: ELPAt => null()
# endif

#endif

    ! Arguments to routines
    character :: jobz, range, uplo, trans

    ! Lapack info
    integer :: info

    ! Ok-eigenvalues
    integer :: neigok

    integer, pointer :: ifail(:) => null()
    integer, pointer :: isuppz(:) => null()
    integer :: il, iu
    real(dp) :: vl, vu
    
    ! Work sizes
    integer :: liwork, lrwork, lwork
    integer, save :: lrwork_add = 0

    complex(dp), pointer :: work(:) => null()
    real(dp), pointer :: rwork(:) => null()
    integer, pointer :: iwork(:) => null()

    ! Local variables for loops etc.
    integer :: i

    integer, external :: numroc

    ! Start time count
    call timer('cdiag',1)

#ifdef MPI
    ! Only re-initialize if the routine hasn't been setup.
    if ( iCTXT < 0 .and. .not. Serial ) call diag_init()
#endif

!*******************************************************************************
! Setup                                                                        *
!*******************************************************************************
      
    ! Initialise error flag
    ierror = 0

    ! Trap n=1 case, which is not handled correctly otherwise (JMS 2011/07/19)
    if ( n == 1 ) then

       w(:) = 0._dp
       w(1) = real(H(1,1), dp) / real(S(1,1), dp)
       Z(:,:) = 0._dp
       Z(1,1) = 1._dp / sqrt( real(S(1,1), dp) )

       call timer('cdiag', 2)

       return

    end if

    ! Get old allocation defaults and set new ones
    call alloc_default( old=oldDefaults, &
         copy=.false., shrink=.true., &
         imin=1, routine='cdiag' )

    ! vl/il and vu/iu are not currently used, set them to not
    ! confuse a new programmer
    vl = -huge(0._dp)
    vu = huge(0._dp)
    il = 1
    iu = neig

    
    ! Correct the input according to the diagonalization
    ! routines and queries
    algo = algorithm
    call diag_correct_input(algo, jobz, range, uplo, trans, iu, n)

    
#ifdef MPI
    if ( .not. Serial) then

       ! Set up blacs descriptors for 1D case
       call descinit( desc_1d, n, n, BlockSize, BlockSize, 0, 0, &
            iCTXT, n, info)
       if ( info /= 0 ) then
          call die('cdiag: Blacs setup has failed!')
       end if

       if ( Use2D ) then

          ! Retrieve information about the BLACS 2D distribution
          call blacs_gridinfo(iCTXT2D, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of 2D matrices
          mat_2d(1) = numroc(n, diag_BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, diag_BlockSize, my2d(2), 0, np2d(2))

          ! Set up blacs descriptors for 2D case
          call descinit(desc_2d, n, n, diag_BlockSize, diag_BlockSize, 0, 0, &
               iCTXT2D, mat_2d(1), info)
          if ( info /= 0 ) then
             call die('cdiag: Blacs setup has failed!')
          end if

          desc => desc_2d(:)

       else

          ! Retrieve information about the BLACS 1D distribution
          call blacs_gridinfo(iCTXT, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of the 1D matrices
          mat_2d(1) = numroc(n, BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, BlockSize, my2d(2), 0, np2d(2))

          desc => desc_1d(:)
          
       end if
       
# ifdef SIESTA__ELPA
       ! Initialize the elpa type
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then

          ! ELPA setup
          ELPAt => elpa_allocate(info)
          call elpa_check(info, 'allocate')
          call ELPAt%set('mpi_comm_parent', MPI_Comm_World, info)
          call elpa_check(info, 'mpi_comm_parent')
          
          call ELPAt%set('process_row', my2d(1), info)
          call elpa_check(info, 'process_row')
          call ELPAt%set('process_col', my2d(2), info)
          call elpa_check(info, 'process_col')
                    
          if ( algorithm == ELPA_1stage ) then
             call ELPAt%set('solver', ELPA_SOLVER_1STAGE, info)
          else if ( algorithm == ELPA_2stage ) then
             call ELPAt%set('solver', ELPA_SOLVER_2STAGE, info)
          end if
          call elpa_check(info, 'solver')
          
          call ELPAt%set('na', n, info)
          call elpa_check(info, 'na')
          call ELPAt%set('local_nrows', mat_2d(1), info)
          call elpa_check(info, 'local_nrows')
          call ELPAt%set('local_ncols', mat_2d(2), info)
          call elpa_check(info, 'local_ncols')
          if ( Use2D ) then
             call ELPAt%set('nblk', diag_BlockSize, info)
          else
             call ELPAt%set('nblk', BlockSize, info)
          end if
          call elpa_check(info, 'nblk')
          ! Set the number of calculated eigenvalues/vectors
          call ELPAt%set('nev', iu, info)
          call elpa_check(info, 'nev')

          ! Now ELPA should be ready to be setup
          info = ELPAt%setup()
          call elpa_check(info, 'setup')

          if (elpa_use_gpu) then
             call ELPAt%set('gpu', 1, info)
             call elpa_check(info, 'gpu set')
             if ( algorithm == ELPA_2stage ) then
                call ELPAt%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,info) 
                call elpa_check(info, 'complex kernel gpu')
             endif
          endif


          ! There is no need to check the info parameter
          ! because ELPA will fail if one sets a value that
          ! has already been set.
          info = 0

       end if
# endif
       
    end if
#endif

    ! Initialize the variables for the different routines
    if ( Serial ) then

#ifdef SIESTA__MRRR
       if ( algo == MRRR .or. &
            algo == MRRR_2stage ) then
          call re_alloc(isuppz, 1, 2*n, name='isuppz')
       end if
#endif

#ifdef MPI
    else

       if ( algo == Expert .or. &
            algo == Expert_2stage ) then

          ! We will add a max-cluster-size of 12 to the
          ! lrwork array
          lrwork_add = max(lrwork_add, (12-1) * n)

          call re_alloc(gap, 1, Nodes, name='gap')
          call re_alloc(iclustr, 1, 2*Nodes, name='iclustr')

       end if

       if ( Use2D ) then

          ! When requesting a 2D redistribution of data
          ! we can in some cases re-use the data-arrays to
          ! limit the required used memory.
          ! This may be so IFF the 2D distributed # of matrix elements
          ! is less than or equal to the 1D (intrinsic Siesta)
          ! distributed # of matrix elements.

          if ( product(mat_2d) > nml*nm ) then
             ! I.e, here we allocate more memory

             call re_alloc(Hp, 1, mat_2d(1), 1, mat_2d(2), name='H2D')
             call re_alloc(Sp, 1, mat_2d(1), 1, mat_2d(2), name='S2D')
             call re_alloc(Zp, 1, mat_2d(1), 1, mat_2d(2), name='Z2D')

          else
             ! I.e, here we re-use memory to substantially reduce the
             ! required memory of Siesta

             ! This order means that nothing gets overwritten
             ! when we distribute to the 2D distribution.
             Hp => S
             Sp => Z
             Zp => H
             
          end if
          
       else
          
          Hp => H
          Sp => S
          Zp => Z

       end if

#endif
    end if

    if ( algo == Expert .or. &
         algo == Expert_2stage ) then
       call re_alloc(ifail, 1, n, name='ifail')
    end if

    ! Perform work-size query
    ! ScaLAPACK typically uses a bit more of the work-size elements
    lwork = 10
    lrwork = 10
    liwork = 10
    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(rwork, 1, lrwork, name='rwork')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Get memory requirements
    call work_query()
    if ( info /= 0 ) then
       write(*, *) 'cdiag: work-query info ', info
       call die('cdiag: work-query error')
    end if

    ! Add lrwork_add
    lrwork = lrwork + lrwork_add
    
#ifdef _DIAG_WORK
    if ( jobz == 'N' ) then
       if ( diag_work_c(1) ) then
          write(*,'(3a,i2,a,i5,3(a,i12))') &
               'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', rwork=', lrwork, ', iwork=', liwork
          diag_work_c(1) = .false.
       end if
    else
       if ( diag_work_c(2) ) then
          write(*,'(3a,i2,a,i5,3(a,i12))') &
               'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', rwork=', lrwork, ', iwork=', liwork
          diag_work_c(2) = .false.
       end if
    end if
#endif

    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(rwork, 1, lrwork, name='rwork')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Begin calculation
    ! Default OK eigenvalues to the queried entries
    neigok = iu


#ifdef MPI
    if ( Use2D .and. .not. Serial ) then
       ! Redistribute to 2D layout
       ! Note that the call sequence HAS to be in this order (see pointers
       ! above).
       call pzgemr2d(n, n, S, 1, 1, desc_1d, Sp, 1, 1, desc_2d, iCTXT)
       call pzgemr2d(n, n, H, 1, 1, desc_1d, Hp, 1, 1, desc_2d, iCTXT)
    end if
#endif
    
!*******************************************************************************
! Factorise overlap matrix                                                     *
!*******************************************************************************
    call timer('cdiag1',1)
    if ( Serial ) then
       call zpotrf(uplo,n,S,n,info)
#ifdef MPI
    else
# ifdef SIESTA__ELPA
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then
          ! ELPA only calculates the Cholesky on the upper
          ! half of the matrix.
          call ELPAt%cholesky(Sp, info)
          call elpa_check(info, 'cdiag: Cholesky factorisation')
          info = 0
       else
          call pzpotrf(uplo,n,Sp,1,1,desc,info)
       end if
# else
       call pzpotrf(uplo,n,Sp,1,1,desc,info)
# endif
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('cdiag: Error in Cholesky factorisation')
    end if
    call timer('cdiag1',2)

!*******************************************************************************
! Transform problem to standard eigenvalue problem                             *
!*******************************************************************************
    call timer('cdiag2',1)
    if ( Serial ) then
       call zhegst(1,uplo,n,H,n,S,n,info)
#ifdef MPI
    else
# ifdef SIESTA__ELPA
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then

          ! Set the scale to be 1 (ELPA does not use the scale)
          scale = 1._dp

          ! The following routine requires the use of the upper
          ! routines, due to the hermitian multiply function
          ! of ELPA.

          ! Invert the upper triangular Cholesky decomposition.
          call ELPAt%invert_triangular(Sp, info)
          call elpa_check(info, 'cdiag: Triangular inversion')
          
          ! Left hand of the inverse
          ! We do require the full matrix calculated
          ! due to the ELPA solver.
          call zcopy(size(Hp),Hp,1,Zp,1)
          call ELPAt%hermitian_multiply('U','F',n, &
               Sp,Zp,mat_2d(1),mat_2d(2), &
               Hp,mat_2d(1),mat_2d(2),info)
          call elpa_check(info, 'cdiag: Hermitian multiply Left')
          
          ! Right hand of the inverse.
          ! Note the ELPA Hermitian multiply always takes the Hermitian
          ! conjugate of the left operator, hence we cannot use it here
          ! We could, possibly use pztranc
          ! and then use hermitian_multiply... (?)
          call pztrmm('R','U','N','N',n,n,dcmplx(1._dp,0._dp), &
               Sp,1,1,desc,Hp,1,1,desc)
          
          info = 0
       else
          call pzhengst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
               desc,scale,work,lwork,info)
       end if
# else
       call pzhengst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
            desc,scale,work,lwork,info)
# endif
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('cdiag: Error in forward transformation')
    end if
    call timer('cdiag2',2)

!*******************************************************************************
! Solve standard eigenvalue problem                                            *
!*******************************************************************************
    call timer('cdiag3',1)
    if ( Serial ) then

       select case ( algo )
       case ( DivideConquer ) 
          call zheevd(jobz,uplo,n,H,n,&
               w,work,lwork,rwork,lrwork,iwork,liwork, &
               info)
          if ( neig > 0 ) then
             call zcopy(n*neig,H,1,Z,1)
          end if
          
#ifdef SIESTA__MRRR
       case ( MRRR )
          call zheevr(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)
#endif
          
       case ( Expert )
          call zheevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,rwork,iwork,ifail, &
               info)
          
       case ( QR )
          call zheev(jobz,uplo,n,H,n,w, &
               work,lwork,rwork, &
               info)
          if ( neig > 0 ) then
             call zcopy(n*neig,H,1,Z,1)
          end if
          
#ifdef SIESTA__DIAG_2STAGE
       case ( DivideConquer_2stage ) 
          call zheevd_2stage(jobz,uplo,n,H,n,&
               w,work,lwork,rwork,lrwork,iwork,liwork, &
               info)
          if ( neig > 0 ) then
             call zcopy(n*neig,H,1,Z,1)
          end if
          
# ifdef SIESTA__MRRR
       case ( MRRR_2stage )
          call zheevr_2stage(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)
# endif

       case ( Expert_2stage )
          call zheevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,rwork,iwork,ifail, &
               info)

       case ( QR_2stage )
          call zheev_2stage(jobz,uplo,n,H,n,w, &
               work,lwork,rwork, &
               info)
          if ( neig > 0 ) then
             call zcopy(n*neig,H,1,Z,1)
          end if
#endif

       end select
       
#ifdef MPI
    else

       select case ( algo )
       case ( DivideConquer )
          call pzheevd(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)

#ifdef SIESTA__MRRR
       case ( MRRR )
          call pzheevr(jobz,range,uplo,n,Hp,1,1,desc, &
               vl,vu,il,iu,neigok,nz,w, &
               Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)
#endif

#ifdef SIESTA__ELPA
       case ( ELPA_1stage, ELPA_2stage )

          ! Calculate the eigenvector or eigenvalues
          if ( neig > 0 ) then
             call ELPAt%eigenvectors(Hp, w, Zp, info)
          else
             call ELPAt%eigenvalues(Hp, w, info)
          end if
          call elpa_check(info, 'cdiag: Solver')
          info = 0
#endif

       case ( Expert ) 
          call pzheevx(jobz,range,uplo,n,Hp,1,1,desc,vl,vu,il,iu, &
               abstol,neigok,nz,w,orfac,Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               ifail,iclustr,gap, &
               info)

          mclustr = 0
          do i = 1, Nodes
             mclustr = max(iclustr(2*i) - iclustr(2*i-1), mclustr)
          end do

          if ( info == -25 ) then
             ! LRWORK is too small to compute all the eigenvectors
             ! However, I do not know by how much... ???
             call die('cdiag: Requires bigger rwork')

          else if ( mod(info,2) /= 0 .or. mod(info/8,2) /= 0 ) then
             
             ! One or more than one eigenvector failed to
             ! converge, we should warn the user to decrease
             ! the tolerance.
             if ( IONode ) then
                write(*,*) "cdiag: Decrease the absolute tolerance "//&
                     "due to insufficient eigenvector convergence..."
             end if
             call die('cdiag: Decrease the absolute tolerance!')
             
          else if ( mod(info/2, 2) /= 0 ) then

             ! We need to signal an increase in workspace
             if ( IONode ) then
                write(*,*) "cdiag: Increasing memory and trying diagonalization again"
             end if
             ierror = -1

             i = (mclustr-1) * n
             if ( lrwork_add < i ) then
                
                ! Try to increase the work-size
                lrwork_add = i
                call clean_memory()
                
                call timer('cdiag3', 2)
                call timer('cdiag', 2)

                return
             end if
             
          end if

       case ( QR )
          call pzheev(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,rwork,lrwork, &
               info)

       end select

#endif
    end if

    ! Check error flag
    if ( info /= 0 ) then
       ierror = 1
       if ( info < 0 ) then
          call die('cdiag: Illegal argument to standard eigensolver')
       else
          call die('cdiag: Failure to converge standard eigenproblem')
       end if
       if ( neigok < iu ) then
          call die('cdiag: Insufficient eigenvalues converged')
       end if
    end if
    ! Ensure that the eigenvalues that haven't been calculated
    ! are "extreme" and hence not applicable
    ! This is merely a precaution
    if ( neigok < n ) then
       do i = neigok + 1 , n
          w(i) = huge(1._dp)
       end do
    end if
    call timer('cdiag3',2)

    
!*******************************************************************************
! Back transformation of eigenvectors                                          *
!*******************************************************************************
    if ( neig > 0 ) then
       call timer('cdiag4',1)
       if ( Serial ) then
          call ztrsm('L',uplo,trans,'N',n,neig,dcmplx(1._dp,0._dp),S,n,Z,n)
#ifdef MPI
       else
# ifdef SIESTA__ELPA
          if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then
             ! Back-transform the eigenvectors
             ! Note the ELPA Hermitian multiply always takes the Hermitian
             ! conjugate of the left operator, so we cannot use it.
             ! Also, the ELPA routines have already calculated
             ! the inverse of Cholesky(S), hence we only need
             ! a matrix-multiplication
             call pztrmm('L','U','N','N',n,neig,dcmplx(1._dp,0._dp), &
                  Sp,1,1,desc,Zp,1,1,desc)
          else
             call pztrsm('L',uplo,trans,'N',n,neig,dcmplx(1._dp,0._dp), &
                  Sp,1,1,desc,Zp,1,1,desc)
          end if
# else
          call pztrsm('L',uplo,trans,'N',n,neig,dcmplx(1._dp,0._dp), &
               Sp,1,1,desc,Zp,1,1,desc)
# endif
          if ( Use2D ) then
             call pzgemr2d(n,neig,Zp,1,1,desc_2d,Z,1,1,desc_1d,iCTXT)
          end if
#endif
       end if
       call timer('cdiag4',2)
    end if
#ifdef MPI
    ! Rescale the eigenvalues
    if ( scale /= 1.0_dp .and. .not. Serial ) then
       call dscal(neigok,scale,w,1)
    end if
#endif
    if ( info /= 0 ) then
       call die('cdiag: Error in back transformation')
    end if

    call clean_memory()

    ! Stop time count
    call timer('cdiag',2)
    
  contains

    subroutine clean_memory()
      
!*******************************************************************************
! Clean up                                                                     *
!*******************************************************************************
      
#ifdef SIESTA__ELPA
      integer :: info
#endif
      ! Deallocate workspace arrays
      if ( Serial ) then
#ifdef SIESTA__MRRR
         if ( algo == MRRR .or. &
              algo == MRRR_2stage ) then
            call de_alloc(isuppz, name='isuppz')
         end if
#endif
#ifdef MPI
      else

         if ( algo == Expert .or. &
              algo == Expert_2stage ) then
            call de_alloc(gap, name='gap')
            call de_alloc(iclustr, name='iclustr')
         end if

         if ( Use2D .and. product(mat_2d) > nml*nm ) then
            call de_alloc(Hp, name='H2D')
            call de_alloc(Sp, name='S2D')
            call de_alloc(Zp, name='Z2D')
         end if

# ifdef SIESTA__ELPA
         if ( associated(ELPAt) ) then
            call elpa_deallocate(ELPAt,info)
            call elpa_check(info, 'deallocate')
            nullify(ELPAt)
         end if
# endif
#endif
      end if

      if ( algo == Expert .or. &
           algo == Expert_2stage ) then
         call de_alloc(ifail, name='ifail')
      end if

      call de_alloc(work, name='work')
      call de_alloc(iwork, name='iwork')
      call de_alloc(rwork, name='rwork')


      !  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

    end subroutine clean_memory
    
    subroutine work_query()
      integer :: l_lwork

      ! Initialize
      l_lwork = 0
      
      work(1) = 1
      rwork(1) = 1
      iwork(1) = 1

      lwork = -1
      lrwork = -1
      liwork = -1

      ! Get memory requirements
      if ( Serial ) then
         
         select case ( algo )
         case ( DivideConquer ) 
            call zheevd(jobz,uplo,n,H,n,&
                 w,work,lwork,rwork,lrwork,iwork,liwork, &
                 info)

#ifdef SIESTA__MRRR
         case ( MRRR )
            call zheevr(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,rwork,lrwork,iwork,liwork, &
                 info)
#endif

         case ( Expert )
            call zheevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,rwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork or iwork
            rwork(1) = max(nint(rwork(1)), 7*n)
            iwork(1) = max(iwork(1), 5*n)

         case ( QR )
            call zheev(jobz,uplo,n,H,n,w, &
                 work,lwork,rwork, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork
            rwork(1) = max(nint(rwork(1)), 3*n)

#ifdef SIESTA__DIAG_2STAGE
         case ( DivideConquer_2stage ) 
            call zheevd_2stage(jobz,uplo,n,H,n,&
                 w,work,lwork,rwork,lrwork,iwork,liwork, &
                 info)

# ifdef SIESTA__MRRR
         case ( MRRR_2stage )
            call zheevr_2stage(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,rwork,lrwork,iwork,liwork, &
                 info)
# endif

         case ( Expert_2stage )
            call zheevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,rwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork or iwork
            rwork(1) = max(nint(rwork(1)), 7*n)
            iwork(1) = max(iwork(1), 5*n)

         case ( QR_2stage )
            call zheev_2stage(jobz,uplo,n,H,n,w, &
                 work,lwork,rwork, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork
            rwork(1) = max(nint(rwork(1)), 3*n)
#endif

         case default
                        
            call die('cdiag: error in work_query')

         end select

#ifdef MPI
      else

         ! They all require pzhengst
         ! Well pzheevx exists in the generalized form, but
         ! we currently do not use it
         call pzhengst(1, uplo, &
              n, Hp, 1, 1, desc, Sp, 1, 1, desc, &
              scale, &
              work, lwork, &
              info)
         l_lwork = work(1)
         work(1) = 1
         lwork = -1

         select case ( algo )
         case ( DivideConquer )
            call pzheevd(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 info)
            
# ifdef SIESTA__MRRR
         case ( MRRR )
            call pzheevr(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, neigok, nz, w, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 info)
# endif

# ifdef SIESTA__ELPA
         case ( ELPA_1stage, ELPA_2stage )
            ! skip the l_lwork for pzhengst routine.
            l_lwork = 1
# endif
            
         case ( Expert )
            call pzheevx(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, abstol, neigok, nz, w, &
                 orfac, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 ifail, iclustr, gap, &
                 info)

         case ( QR )
            call pzheev(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, &
                 info)
            ! Possible bug in scalapack
            ! At least this makes it work!
#ifdef _DIAG_WORK
            if ( jobz == 'N' ) then
               if ( diag_work_c(1) ) then
                  if ( nint(rwork(1)) < 2*n ) then
                     write(*,'(3a,i2,a,i5,a)') &
                          'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
                          ' BUG in ScaLAPACK work-query'
                  end if
               end if
            else
               if ( diag_work_c(2) ) then
                  if ( nint(rwork(1)) < 4*n-2 ) then
                     write(*,'(3a,i2,a,i5,a)') &
                          'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
                          ' BUG in ScaLAPACK work-query'
                  end if
               end if
            end if
#endif
            if ( jobz == 'N' ) then
               rwork(1) = max(nint(rwork(1)), 2*n)
            else
               rwork(1) = max(nint(rwork(1)), 4*n-2)
            end if

         case default

            call die('cdiag: error in work_query')

         end select
#endif         
      end if

      lwork = nint(max(nint(real(work(1), dp)), l_lwork) * mem_factor)
      lrwork = nint(rwork(1) * mem_factor)
      liwork = iwork(1)

      lwork = max(1, lwork)
      lrwork = max(1, lrwork)
      liwork = max(1, liwork)

    end subroutine work_query
    
  end subroutine diag_c

  
  subroutine diag_r( H, S, n, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
! ***************************************************************************
! Subroutine to solve all eigenvalues and eigenvectors of the
! real general eigenvalue problem  H z = w S z,  with H and S
! real symmetry matrices.
! Written by G.Fabricius and J.Soler, March 1998
! Rewritten by Julian Gale, August 2004
! Rewritten by Nick R. Papior, July 2017
! ************************** INPUT ******************************************
! real*8 H(nml,nm)                 : Symmetric H matrix
! real*8 S(nml,nm)                 : Symmetric S matrix
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neig                     : No. of eigenvalues(<0)/vectors(>0) to calculate
! integer iscf                     : SCF cycle
! integer BlockSize                : Effective parallel block size
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! real*8 Z(nml,nm)                 : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ************************* PARALLEL ****************************************
! When running in parallel this routine now uses Scalapack to perform a
! parallel matrix diagonalisation. This requires Scalapack and Blacs to
! be installed first. Although globally a 1-D block cyclic data distribution
! is employed, locally 1 or 2-D distributions are allowed for.
! The blocksize is now explicitly passed to the routine (A. Garcia, 2016)      
! The routine allows access to all the phases of diagonalisation for fuller
! control, and allows for parallel divide and conquer with reduced memory.
! The presence of eigenvalue clusters is checked for and the memory adjusted
! accordingly to try to guarantee convergence.
! Note that the blacs grid is only initialised on the first call to avoid
! exceeding MPI limits for split/group creation.
!
! When running in parallel and using a 2D distribution it is crucial that the
! Z array is also allocated.
! When using a 2D distribution the algorithm tries to reuse memory without
! double allocation. This is allowed if the 2D distribution # of elements
! is less than or equal to the number of elements in the local arrays.
! (NOTE this is typically so).
!
! It allows for the following routines:
!   LAPACK (Serial or ParallelOverK):
!    - dsyevd, dsyevd_2stage
!    - dsyevr, dsyevr_2stage
!    - dsyevx, dsyevx_2stage
!    - dsyev, dsyev_2stage
!   ScaLAPACK (Parallel):
!    - pdsyevd
!    - pdsyevr
!    - pdsyevx
!    - pdsyev
! To enable the *vr* routines this file needs to be compiled with:
!   -DSIESTA__MRRR
! To enable the *v*_2stage routines this file needs to be compiled with:
!   -DSIESTA__DIAG_2STAGE
! Note that the LAPACK implementation of the 2stage solvers does (in 3.7.1) not
! implement the jobz='V' stage.
! ***************************************************************************

    ! Modules
    use m_diag_option ! just all
    
    use alloc
    use sys, only : die

#ifdef SIESTA__ELPA
    use mpi_siesta, only: MPI_Comm_World
    use elpa
#endif

    ! Passed variables
    integer, intent(out) :: ierror
    integer, intent(in) :: iscf
    integer, intent(in) :: n
    integer, intent(in) :: neig
    integer, intent(in) :: nm, nml
    integer, intent(in) :: BlockSize
    real(dp), intent(inout) :: w(nml)
    real(dp), intent(inout), target :: H(nml,nm)
    real(dp), intent(inout), target :: S(nml,nm)
    real(dp), intent(inout), target :: Z(nml,nm)

    ! Local variables
    type(allocDefaults) :: oldDefaults
    integer :: algo
#ifdef MPI

    ! Scale when transforming to generalized eigenvalue problem
    real(dp) :: scale

    ! Expert and MRRR drivers
    integer :: nz, mclustr
    integer,  pointer :: iclustr(:) => null()
    real(dp), pointer :: gap(:) => null()

    ! BLACS descriptors
    integer, target :: desc_1d(9), desc_2d(9)
    integer, pointer :: desc(:) => null()

    ! Additional variables for a 2D parallel grid
    integer :: np2d(2), my2d(2), mat_2d(2)

    ! Pointers of H, S and Z
    real(dp), pointer :: Hp(:,:) => null()
    real(dp), pointer :: Sp(:,:) => null()
    real(dp), pointer :: Zp(:,:) => null()

# ifdef SIESTA__ELPA
    class(elpa_t), pointer :: ELPAt => null()
# endif

#endif
    ! Arguments to routines
    character :: jobz, range, uplo, trans

    ! Lapack info
    integer :: info

    ! Ok-eigenvalues
    integer :: neigok

    integer, pointer :: ifail(:) => null()
    integer, pointer :: isuppz(:) => null()
    integer :: il, iu
    real(dp) :: vl, vu
    
    ! Work sizes
    integer :: liwork, lwork
    integer, save :: lwork_add = 0

    real(dp), pointer :: work(:) => null()
    integer, pointer :: iwork(:) => null()

    ! Local variables for loops etc.
    integer :: i

    integer, external :: numroc

    ! Start time count
    call timer('rdiag',1)

#ifdef MPI
    ! Only re-initialize if the routine hasn't been setup.
    if ( iCTXT < 0 .and. .not. Serial ) call diag_init()
#endif

!*******************************************************************************
! Setup                                                                        *
!*******************************************************************************
      
    ! Initialise error flag
    ierror = 0

    ! Trap n=1 case, which is not handled correctly otherwise (JMS 2011/07/19)
    if ( n == 1 ) then

       w(:) = 0._dp
       w(1) = H(1,1) / S(1,1)
       Z(:,:) = 0._dp
       Z(1,1) = 1._dp / sqrt( S(1,1) )

       call timer('rdiag', 2)

       return

    end if

    ! Get old allocation defaults and set new ones
    call alloc_default( old=oldDefaults, &
         copy=.false., shrink=.true., &
         imin=1, routine='rdiag' )

    ! vl/il and vu/iu are not currently used, set them to not
    ! confuse a new programmer
    vl = -huge(0._dp)
    vu = huge(0._dp)
    il = 1
    iu = neig
    

    ! Correct the input according to the diagonalization
    ! routines and queries
    algo = algorithm
    call diag_correct_input(algo, jobz, range, uplo, trans, iu, n)

    
#ifdef MPI
    if ( .not. Serial) then

       ! Set up blacs descriptors for 1D case
       call descinit( desc_1d, n, n, BlockSize, BlockSize, 0, 0, &
            iCTXT, n, info)
       if ( info /= 0 ) then
          call die('rdiag: Blacs setup has failed!')
       end if

       if ( Use2D ) then

          ! Retrieve information about the BLACS 2D distribution
          call blacs_gridinfo(iCTXT2D, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of 2D matrices
          mat_2d(1) = numroc(n, diag_BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, diag_BlockSize, my2d(2), 0, np2d(2))

          ! Set up blacs descriptors for 2D case
          call descinit(desc_2d, n, n, diag_BlockSize, diag_BlockSize, 0, 0, &
               iCTXT2D, mat_2d(1), info)
          if ( info /= 0 ) then
             call die('rdiag: Blacs setup has failed!')
          end if

          desc => desc_2d(:)

       else

          ! Retrieve information about the BLACS 1D distribution
          call blacs_gridinfo(iCTXT, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of the 1D matrices
          mat_2d(1) = numroc(n, BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, BlockSize, my2d(2), 0, np2d(2))

          desc => desc_1d(:)
          
       end if

# ifdef SIESTA__ELPA
       ! Initialize the elpa type
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then
          
          ! ELPA setup
          ELPAt => elpa_allocate(info)
          call elpa_check(info, 'allocate')
          call ELPAt%set('mpi_comm_parent', MPI_Comm_World, info)
          call elpa_check(info, 'mpi_comm_parent')
          
          call ELPAt%set('process_row', my2d(1), info)
          call elpa_check(info, 'process_row')
          call ELPAt%set('process_col', my2d(2), info)
          call elpa_check(info, 'process_col')
                    
          if ( algorithm == ELPA_1stage ) then
             call ELPAt%set('solver', ELPA_SOLVER_1STAGE, info)
          else if ( algorithm == ELPA_2stage ) then
             call ELPAt%set('solver', ELPA_SOLVER_2STAGE, info)
          end if
          call elpa_check(info, 'solver')
          
          call ELPAt%set('na', n, info)
          call elpa_check(info, 'na')
          call ELPAt%set('local_nrows', mat_2d(1), info)
          call elpa_check(info, 'local_nrows')
          call ELPAt%set('local_ncols', mat_2d(2), info)
          call elpa_check(info, 'local_ncols')
          if ( Use2D ) then
             call ELPAt%set('nblk', diag_BlockSize, info)
          else
             call ELPAt%set('nblk', BlockSize, info)
          end if
          call elpa_check(info, 'nblk')
          ! Set the number of calculated eigenvalues/vectors
          call ELPAt%set('nev', iu, info)
          call elpa_check(info, 'nev')

          ! Now ELPA should be ready to be setup
          info = ELPAt%setup()
          call elpa_check(info, 'setup')

          if (elpa_use_gpu) then
             call ELPAt%set('gpu', 1, info)
             call elpa_check(info, 'gpu set')
             if ( algorithm == ELPA_2stage ) then
                call ELPAt%set("real_kernel",ELPA_2STAGE_REAL_GPU,info) 
                call elpa_check(info, 'real kernel gpu')
             endif
          endif

          ! There is no need to check the info parameter
          ! because ELPA will fail if one sets a value that
          ! has already been set.
          info = 0

       end if
# endif

    end if
#endif

    
    ! Initialize the variables for the different routines
    if ( Serial ) then

#ifdef SIESTA__MRRR
       if ( algo == MRRR .or. &
            algo == MRRR_2stage ) then
          call re_alloc(isuppz, 1, 2*n, name='isuppz')
       end if
#endif

#ifdef MPI
    else

       if ( algo == Expert .or. &
            algo == Expert_2stage ) then

          ! We will add a max-cluster-size of 12 to the
          ! lwork array
          lwork_add = max(lwork_add, (12-1) * n)

          call re_alloc(gap, 1, Nodes, name='gap')
          call re_alloc(iclustr, 1, 2*Nodes, name='iclustr')

       end if

       if ( Use2D ) then

          ! When requesting a 2D redistribution of data
          ! we can in some cases re-use the data-arrays to
          ! limit the required used memory.
          ! This may be so IFF the 2D distributed # of matrix elements
          ! is less than or equal to the 1D (intrinsic Siesta)
          ! distributed # of matrix elements.

          if ( product(mat_2d) > nml*nm ) then
             ! I.e, here we allocate more memory

             call re_alloc(Hp, 1, mat_2d(1), 1, mat_2d(2), name='H2D')
             call re_alloc(Sp, 1, mat_2d(1), 1, mat_2d(2), name='S2D')
             call re_alloc(Zp, 1, mat_2d(1), 1, mat_2d(2), name='Z2D')

          else
             ! I.e, here we re-use memory to substantially reduce the
             ! required memory of Siesta

             ! This order means that nothing gets overwritten
             ! when we distribute to the 2D distribution.
             Hp => S
             Sp => Z
             Zp => H
             
          end if
          
       else
          
          Hp => H
          Sp => S
          Zp => Z

       end if

#endif
    end if

    if ( algo == Expert .or. &
         algo == Expert_2stage ) then
       call re_alloc(ifail, 1, n, name='ifail')
    end if

    ! Perform work-size query
    ! ScaLAPACK typically uses a bit more of the work-size elements
    lwork = 10
    liwork = 10
    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Get memory requirements
    call work_query()
    if ( info /= 0 ) then
       write(*, *) 'rdiag: work-query info ', info
       call die('rdiag: work-query error')
    end if

    ! Add lwork_add
    lwork = lwork + lwork_add

#ifdef _DIAG_WORK
    if ( jobz == 'N' ) then
       if ( diag_work_r(1) ) then
          write(*,'(3a,i2,a,i5,2(a,i12))') &
               'rdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', iwork=', liwork
          diag_work_r(1) = .false.
       end if
    else
       if ( diag_work_r(2) ) then
          write(*,'(3a,i2,a,i5,2(a,i12))') &
               'rdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', iwork=', liwork
          diag_work_r(2) = .false.
       end if
    end if
#endif

    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Begin calculation
    ! Default OK eigenvalues to the queried entries
    neigok = iu


#ifdef MPI
    if ( Use2D .and. .not. Serial ) then
       ! Redistribute to 2D layout
       ! Note that the call sequence HAS to be in this order (see pointers
       ! above).
       call pdgemr2d(n, n, S, 1, 1, desc_1d, Sp, 1, 1, desc_2d, iCTXT)
       call pdgemr2d(n, n, H, 1, 1, desc_1d, Hp, 1, 1, desc_2d, iCTXT)
    end if
#endif
    
!*******************************************************************************
! Factorise overlap matrix                                                     *
!*******************************************************************************
    call timer('rdiag1',1)
    if ( Serial ) then
       call dpotrf(uplo,n,S,n,info)
#ifdef MPI
    else
# ifdef SIESTA__ELPA
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then
          ! ELPA only calculates the Cholesky on the upper
          ! half of the matrix.
          call ELPAt%cholesky(Sp, info)
          call elpa_check(info, 'rdiag: Cholesky factorisation')
          info = 0
       else
          call pdpotrf(uplo,n,Sp,1,1,desc,info)
       end if
# else
       call pdpotrf(uplo,n,Sp,1,1,desc,info)
# endif
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('rdiag: Error in Cholesky factorisation')
    end if
    call timer('rdiag1',2)

!*******************************************************************************
! Transform problem to standard eigenvalue problem                             *
!*******************************************************************************
    call timer('rdiag2',1)
    if ( Serial ) then
       call dsygst(1,uplo,n,H,n,S,n,info)
#ifdef MPI
    else
# ifdef SIESTA__ELPA
       if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then

          ! Set the scale to be 1 (ELPA does not use the scale)
          scale = 1._dp

          ! The following routine requires the use of the upper
          ! routines, due to the hermitian multiply function
          ! of ELPA.

          ! Invert the upper triangular Cholesky decomposition.
          call ELPAt%invert_triangular(Sp, info)
          call elpa_check(info, 'rdiag: Triangular inversion')
          
          ! Left hand of the inverse
          ! We do require the full matrix calculated
          ! due to the ELPA solver.
          call dcopy(size(Hp),Hp,1,Zp,1)
          call ELPAt%hermitian_multiply('U','F',n, &
               Sp,Zp,mat_2d(1),mat_2d(2), &
               Hp,mat_2d(1),mat_2d(2),info)
          call elpa_check(info, 'rdiag: Hermitian multiply Left')

          ! Right hand of the inverse.
          ! Note the ELPA Hermitian multiply always takes the Hermitian
          ! conjugate of the left operator, hence we cannot use it here
          ! We could, possibly use pdtran
          ! and then use hermitian_multiply... (?)
          call pdtrmm('R','U','N','N',n,n,1._dp, &
               Sp,1,1,desc,Hp,1,1,desc)
          
          info = 0
       else
          call pdsyngst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
               desc,scale,work,lwork,info)
       end if
# else
       call pdsyngst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
            desc,scale,work,lwork,info)
# endif
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('rdiag: Error in forward transformation')
    end if
    call timer('rdiag2',2)

!*******************************************************************************
! Solve standard eigenvalue problem                                            *
!*******************************************************************************
    call timer('rdiag3',1)
    if ( Serial ) then

       select case ( algo )
       case ( DivideConquer ) 
          call dsyevd(jobz,uplo,n,H,n,&
               w,work,lwork,iwork,liwork, &
               info)
          if ( neig > 0 ) then
             call dcopy(n*neig,H,1,Z,1)
          end if
          
#ifdef SIESTA__MRRR
       case ( MRRR )
          call dsyevr(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,iwork,liwork, &
               info)
#endif
          
       case ( Expert )
          call dsyevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,iwork,ifail, &
               info)
          
       case ( QR )
          call dsyev(jobz,uplo,n,H,n,w, &
               work,lwork, &
               info)
          if ( neig > 0 ) then
             call dcopy(n*neig,H,1,Z,1)
          end if
          
#ifdef SIESTA__DIAG_2STAGE
       case ( DivideConquer_2stage ) 
          call dsyevd_2stage(jobz,uplo,n,H,n,&
               w,work,lwork,iwork,liwork, &
               info)
          if ( neig > 0 ) then
             call dcopy(n*neig,H,1,Z,1)
          end if
          
# ifdef SIESTA__MRRR
       case ( MRRR_2stage )
          call dsyevr_2stage(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,iwork,liwork, &
               info)
# endif

       case ( Expert_2stage )
          call dsyevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,iwork,ifail, &
               info)

       case ( QR_2stage )
          call dsyev_2stage(jobz,uplo,n,H,n,w, &
               work,lwork, &
               info)
          if ( neig > 0 ) then
             call dcopy(n*neig,H,1,Z,1)
          end if
#endif

       end select
       
#ifdef MPI
    else

       select case ( algo )
       case ( DivideConquer )
          call pdsyevd(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               info)

#ifdef SIESTA__MRRR
       case ( MRRR )
          call pdsyevr(jobz,range,uplo,n,Hp,1,1,desc, &
               vl,vu,il,iu,neigok,nz,w, &
               Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               info)
#endif

#ifdef SIESTA__ELPA
       case ( ELPA_1stage, ELPA_2stage )

          ! Calculate the eigenvector or eigenvalues
          if ( neig > 0 ) then
             call ELPAt%eigenvectors(Hp, w, Zp, info)
          else
             call ELPAt%eigenvalues(Hp, w, info)
          end if
          call elpa_check(info, 'rdiag: Solver')
          info = 0
#endif

       case ( Expert ) 
          call pdsyevx(jobz,range,uplo,n,Hp,1,1,desc,vl,vu,il,iu, &
               abstol,neigok,nz,w,orfac,Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               ifail,iclustr,gap, &
               info)

          mclustr = 0
          do i = 1, Nodes
             mclustr = max(iclustr(2*i) - iclustr(2*i-1), mclustr)
          end do

          if ( info == -25 ) then
             ! However, I do not know by how much... ???
             call die('rdiag: Requires bigger work')

          else if ( mod(info,2) /= 0 .or. mod(info/8,2) /= 0 ) then
             
             ! One or more than one eigenvector failed to
             ! converge, we should warn the user to decrease
             ! the tolerance.
             if ( IONode ) then
                write(*,*) "rdiag: Decrease the absolute tolerance "//&
                     "due to insufficient eigenvector convergence..."
             end if
             call die('rdiag: Decrease the absolute tolerance!')
             
          else if ( mod(info/2, 2) /= 0 ) then

             ! We need to signal an increase in workspace
             if ( IONode ) then
                write(*,*) "rdiag: Increasing memory and trying diagonalization again"
             end if
             ierror = -1

             i = (mclustr-1) * n
             if ( lwork_add < i ) then
                
                ! Try to increase the work-size
                lwork_add = i
                call clean_memory()
                
                call timer('rdiag3', 2)
                call timer('rdiag', 2)
                
                return
             end if
             
          end if

       case ( QR )
          call pdsyev(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork, &
               info)

       end select

#endif
    end if

    ! Check error flag
    if ( info /= 0 ) then
       ierror = 1
       if ( info < 0 ) then
          call die('rdiag: Illegal argument to standard eigensolver')
       else
          call die('rdiag: Failure to converge standard eigenproblem')
       end if
       if ( neigok < iu ) then
          call die('rdiag: Insufficient eigenvalues converged')
       end if
    end if
    ! Ensure that the eigenvalues that haven't been calculated
    ! are "extreme" and hence not applicable
    if ( neigok < n ) then
       do i = neigok + 1 , n
          w(i) = huge(1._dp)
       end do
    end if
    call timer('rdiag3',2)

    
!*******************************************************************************
! Back transformation of eigenvectors                                          *
!*******************************************************************************
    if ( neig > 0 ) then
       call timer('rdiag4',1)
       if ( Serial ) then
          call dtrsm('L',uplo,trans,'N',n,neig,1._dp,S,n,Z,n)
#ifdef MPI
       else
# ifdef SIESTA__ELPA
          if ( algo == ELPA_1stage .or. algo == ELPA_2stage ) then
             ! Back-transform the eigenvectors
             ! Note the ELPA Hermitian multiply always takes the Hermitian
             ! conjugate of the left operator
             call pdtrmm('L','U','N','N',n,neig,1._dp, &
                  Sp,1,1,desc,Zp,1,1,desc)
          else
             call pdtrsm('L',uplo,trans,'N',n,neig,1._dp, &
                  Sp,1,1,desc,Zp,1,1,desc)
          end if
# else
          call pdtrsm('L',uplo,trans,'N',n,neig,1._dp, &
               Sp,1,1,desc,Zp,1,1,desc)
# endif
          if ( Use2D ) then
             call pdgemr2d(n,neig,Zp,1,1,desc_2d,Z,1,1,desc_1d,iCTXT)
          end if
#endif
       end if
       call timer('rdiag4',2)
    end if
#ifdef MPI
    ! Rescale the eigenvalues
    if ( scale /= 1.0_dp .and. .not. Serial ) then
       call dscal(neigok,scale,w,1)
    end if
#endif
    if ( info /= 0 ) then
       call die('rdiag: Error in back transformation')
    end if

    call clean_memory()

    ! Stop time count
    call timer('rdiag',2)
    
  contains

    subroutine clean_memory()
      
!*******************************************************************************
! Clean up                                                                     *
!*******************************************************************************
      
#ifdef SIESTA__ELPA
      integer :: info
#endif
      ! Deallocate workspace arrays
      if ( Serial ) then
#ifdef SIESTA__MRRR
         if ( algo == MRRR .or. &
              algo == MRRR_2stage ) then
            call de_alloc(isuppz, name='isuppz')
         end if
#endif
#ifdef MPI
      else

         if ( algo == Expert .or. &
              algo == Expert_2stage ) then
            call de_alloc(gap, name='gap')
            call de_alloc(iclustr, name='iclustr')
         end if

         if ( Use2D .and. product(mat_2d) > nml*nm ) then
            call de_alloc(Hp, name='H2D')
            call de_alloc(Sp, name='S2D')
            call de_alloc(Zp, name='Z2D')
         end if

# ifdef SIESTA__ELPA
         if ( associated(ELPAt) ) then
            call elpa_deallocate(ELPAt,info)
            call elpa_check(info, 'deallocate')
            nullify(ELPAt)
         end if
# endif
#endif
      end if

      if ( algo == Expert .or. &
           algo == Expert_2stage ) then
         call de_alloc(ifail, name='ifail')
      end if

      call de_alloc(work, name='work')
      call de_alloc(iwork, name='iwork')


      !  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

    end subroutine clean_memory
    
    subroutine work_query()
      integer :: l_lwork

      ! Initialize
      l_lwork = 0
      
      work(1) = 1
      iwork(1) = 1

      lwork = -1
      liwork = -1

      ! Get memory requirements
      if ( Serial ) then
         
         select case ( algo )
         case ( DivideConquer ) 
            call dsyevd(jobz,uplo,n,H,n,&
                 w,work,lwork,iwork,liwork, &
                 info)

#ifdef SIESTA__MRRR
         case ( MRRR )
            call dsyevr(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,iwork,liwork, &
                 info)
#endif

         case ( Expert )
            call dsyevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! iwork
            iwork(1) = max(iwork(1), 5*n)

         case ( QR )
            call dsyev(jobz,uplo,n,H,n,w, &
                 work,lwork, &
                 info)

#ifdef SIESTA__DIAG_2STAGE
         case ( DivideConquer_2stage ) 
            call dsyevd_2stage(jobz,uplo,n,H,n,&
                 w,work,lwork,iwork,liwork, &
                 info)

# ifdef SIESTA__MRRR
         case ( MRRR_2stage )
            call dsyevr_2stage(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,iwork,liwork, &
                 info)
# endif

         case ( Expert_2stage )
            call dsyevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! iwork
            iwork(1) = max(iwork(1), 5*n)

         case ( QR_2stage )
            call dsyev_2stage(jobz,uplo,n,H,n,w, &
                 work,lwork, &
                 info)
#endif

         case default
                        
            call die('rdiag: error in work_query')

         end select

#ifdef MPI
      else

         ! They all require pdsyngst
         ! Well pdsyngst exists in the generalized form, but
         ! we currently do not use it
         call pdsyngst(1, uplo, &
              n, Hp, 1, 1, desc, Sp, 1, 1, desc, &
              scale, &
              work, lwork, &
              info)
         l_lwork = work(1)
         work(1) = 1
         lwork = -1

         select case ( algo )
         case ( DivideConquer )
            call pdsyevd(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 info)
            
# ifdef SIESTA__MRRR
         case ( MRRR )
            call pdsyevr(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, neigok, nz, w, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 info)
# endif

# ifdef SIESTA__ELPA
         case ( ELPA_1stage, ELPA_2stage )
            l_lwork = 1
# endif
            
         case ( Expert ) 
            call pdsyevx(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, abstol, neigok, nz, w, &
                 orfac, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 ifail, iclustr, gap, &
                 info)

         case ( QR )
            call pdsyev(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, Zp, 1, 1, desc, &
                 work, lwork, &
                 info)

         case default

            call die('rdiag: error in work_query')

         end select
#endif         
      end if

      lwork = nint(max(nint(work(1)), l_lwork) * mem_factor)
      liwork = iwork(1)

      lwork = max(1, lwork)
      liwork = max(1, liwork)

    end subroutine work_query
    
  end subroutine diag_r


#ifdef SIESTA__ELPA
  subroutine elpa_check(err, name)
    use elpa, only: elpa_strerr, ELPA_OK
    integer, intent(in) :: err
    character(len=*), intent(in) :: name
    
    if ( err == ELPA_OK ) return
    
    write(*,'(a)') 'diag: ELPA error on ' //trim(name)
    write(*,'(a)') elpa_strerr(err)
    call die('diag: ELPA error, see output')
    
  end subroutine elpa_check
#endif
  
end module m_diag


subroutine cdiag( H, S, n, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
  use precision, only: dp
  use m_diag, only: diag_c

  implicit none
  
  integer, intent(in) :: nml, nm
  complex(dp), intent(inout), target :: H(nml,nm)
  complex(dp), intent(inout), target :: S(nml,nm)
  real(dp), intent(inout) :: w(nml)
  complex(dp), intent(inout), target :: Z(nml,nm)
  integer, intent(in) :: n
  integer, intent(in) :: neig
  integer, intent(in) :: iscf
  integer, intent(out) :: ierror
  integer, intent(in) :: BlockSize
  
  call diag_c(H, S, N, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
  
end subroutine cdiag

subroutine rdiag( H, S, n, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
  use precision, only: dp
  use m_diag, only: diag_r
  
  implicit none
  
  integer, intent(in) :: nml, nm
  real(dp), intent(inout), target :: H(nml,nm)
  real(dp), intent(inout), target :: S(nml,nm)
  real(dp), intent(inout) :: w(nml)
  real(dp), intent(inout), target :: Z(nml,nm)
  integer, intent(in) :: n
  integer, intent(in) :: neig
  integer, intent(in) :: iscf
  integer, intent(out) :: ierror
  integer, intent(in) :: BlockSize
  
  call diag_r(H, S, N, nm, nml, w, Z, neig, iscf, ierror, BlockSize)
  
end subroutine rdiag
