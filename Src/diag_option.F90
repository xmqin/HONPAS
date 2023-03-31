! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
#ifndef MPI
! Make sure ELPA is only used for MPI
# undef SIESTA__ELPA
#endif
module m_diag_option

  use precision, only: dp
  
  implicit none
  
!  public :: read_diag, print_diag
!  public :: diag_recognizes_neigwanted

  public
  
  save

  !> Whether diagonalization calls are made via LAPACK (true) or ScaLAPACK (false)
  logical :: Serial = .true.

  !> Whether ScaLAPACK uses a 2D distribution
  logical :: Use2D = .true.
  !> Number of processors on the columns for the 2D distribution
  integer :: ProcessorY = 1
  !> The block-size
  integer :: diag_BlockSize = 24
  
  !> Whether we should use the upper or lower part of the Hermitian/symmetric matrices
  character :: UpperLower = 'L'

  ! Different choices of algorithms.

  ! Choose the algorithm
  ! Note that if any of these are true, the others should be false

  ! The 2stage solvers in LAPACK are (as of 3.7.1) not able to
  ! calculate the eigenvectors.
  ! Hence, they will not be documented.

  !> Use the divide-and-conquer algorithm
  integer, parameter :: DivideConquer = 1
  !> Use the 2-stage divide-and-conquer algorithm (only LAPACK)
  integer, parameter :: DivideConquer_2stage = 2
  !> Use the MRRR algorithm (RRR for LAPACK)
  integer, parameter :: MRRR = 3
  !> Use the 2-stage MRRR algorithm (only LAPACK)
  integer, parameter :: MRRR_2stage = 4
  !> Use the expert driver (if the above are all false)
  integer, parameter :: Expert = 5
  !> Use the 2-stage expert driver (only LAPACK ;\)
  integer, parameter :: Expert_2stage = 6
  !> Use the QR driver
  integer, parameter :: QR = 7
  !> Use the 2-stage QR driver (only LAPACK)
  integer, parameter :: QR_2stage = 8
  !> Use the ELPA driver
  integer, parameter :: ELPA_1stage = 9
  !> Use the 2-stage ELPA driver
  integer, parameter :: ELPA_2stage = 10

  integer :: algorithm = DivideConquer

  !> Tolerance for MRRR (LAPACK) and expert drivers
  real(dp) :: abstol = 1.e-8_dp
  !> Tolerance for expert ScaLAPACK driver
  real(dp) :: orfac = 1.e-3_dp

  !> Memory factor for the real work arrays
  real(dp) :: mem_factor = 1._dp

  logical :: ParallelOverK = .false.
#ifdef SIESTA__ELPA
  logical :: elpa_use_gpu = .false.
#endif
  
contains

  subroutine read_diag(Gamma, nspin)

    use parallel, only: IONode, Nodes, BlockSize
    use fdf, only: fdf_get, leqi

    logical, intent(in) :: Gamma
    integer, intent(in) :: nspin

    integer :: nr

    character(len=32) :: algo

    ! No matter what we always re-read them.
    ! In case one will change the algorithm
    ! we will allow that
    ! (However, note that rdiag and cdiag haven't implemented this yet)

#ifdef MPI
    if ( Nodes > 1 .and. .not. Gamma ) then
       ParallelOverK = fdf_get( 'Diag.ParallelOverK', .false.)
    end if

    if ( Nodes == 1 ) then
       Serial = .true.
       ParallelOverK = .false.
    else if ( ParallelOverK ) then
       Serial = .true.
    else
       Serial = .false.
    end if
#endif

    ! Determine whether we should default a 2D block-cyclic
    ! A 2D block-cyclic for small # of nodes may not be the best
    ! choice, say if it a 1 x Y distribution, it makes no sense
    ! to use a 2D block-cyclic one.

    ! The first step is to calculate the number
    ! of processors in the row
    !  row X col
    do nr = nint(sqrt(real(Nodes))), 1, -1
       if ( mod(Nodes, nr) == 0 ) exit
    end do
    ! Ensure it is minimally 1
    nr = max(1, nr)

    ! Query the requested number of processors in the rows
    ProcessorY = fdf_get('Diag.ProcessorY', nr)
    ProcessorY = max(1, ProcessorY)
    ! Assert that it is a valid distribution
    ! If not, correct it by reducing it until a common
    ! multiple is reached. This will prevent
    ! SIESTA from crashing if a wrong input is provided.
    if ( mod(Nodes, ProcessorY) /= 0 ) then
       nr = ProcessorY
       do ProcessorY = nr, 1, -1
          if ( mod(Nodes, ProcessorY) == 0 ) exit
       end do
       ProcessorY = max(1, ProcessorY)
    end if

    ! Retrieve the blocksize
    diag_BlockSize = fdf_get('Diag.BlockSize', BlockSize)

    ! If there are very few processors, it makes
    ! no sense to make them 2D distributed
    Use2D = (ProcessorY > 1) .and. (Nodes / ProcessorY > 1)
    Use2D = Use2D .or. (BlockSize /= diag_BlockSize)
    Use2D = fdf_get('Diag.Use2D', Use2D)
    
    ! Fall back to original BlockSize when not requesting 2D
    ! distribution (the 1D distribution is not implemented
    ! for different blocksize)
    if ( .not. Use2D ) then
       diag_BlockSize = BlockSize
    end if

    algo = fdf_get('Diag.UpperLower', 'lower')
    if ( leqi(algo, 'lower') .or. leqi(algo, 'l') ) then
       UpperLower = 'L'
    else if ( leqi(algo, 'upper') .or. leqi(algo, 'u') ) then
       UpperLower = 'U'
    else
       call die('diag: Unknown argument to Diag.UpperLower U|L')
    end if

    
    ! Decide the default algorithm
    algo = ' '
    
    if ( fdf_get('Diag.DivideAndConquer',.true.) ) then
       algo = 'Divide-and-Conquer'
    end if
    
#ifdef SIESTA__MRRR
    if ( fdf_get('Diag.MRRR',.false.) ) then
       algo = 'MRRR'
    end if
#endif
    
#ifdef SIESTA__ELPA
    if ( fdf_get('Diag.ELPA',.false.) ) then
       algo = 'ELPA'
    end if
#endif
    
    if ( fdf_get('Diag.NoExpert',.false.) ) then
       algo = 'QR'
    end if

    
    ! Assert that it has been set, or default to the
    ! expert driver
    if ( len_trim(algo) == 0 ) then
       algo = 'Expert'
    end if


    ! Get the requested algorithm by using the above default
    algo = fdf_get('Diag.Algorithm', trim(algo))
    
    
    ! Determine the global method
    if ( leqi(algo, 'D&C') .or. leqi(algo, 'divide-and-conquer') .or. &
         leqi(algo, 'DandC') .or. leqi(algo, 'vd') ) then
       algorithm = DivideConquer

    else if ( leqi(algo, 'D&C-2') .or. leqi(algo, 'D&C-2stage') .or. &
         leqi(algo, 'divide-and-conquer-2stage') .or. leqi(algo, 'DandC-2stage') .or. &
         leqi(algo, 'DandC-2') .or. leqi(algo, 'vd_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = DivideConquer_2stage
       else
          algorithm = DivideConquer
       end if
#else
       algorithm = DivideConquer
#endif

#ifdef SIESTA__ELPA
    else if ( leqi(algo, 'elpa-1') .or. leqi(algo, 'elpa-1stage') ) then
       algorithm = ELPA_1stage
       
       ! The current ELPA implementation requires non-serial
       Serial = .false.
       ParallelOverK = .false.

    else if ( leqi(algo, 'elpa') .or. &
         leqi(algo, 'elpa-2stage') .or. leqi(algo, 'elpa-2') ) then
       algorithm = ELPA_2stage

       ! The current ELPA implementation requires non-serial
       Serial = .false.
       ParallelOverK = .false.

#endif
       
#ifdef SIESTA__MRRR
    else if ( leqi(algo, 'MRRR') .or. leqi(algo, 'RRR') .or. &
         leqi(algo, 'vr') ) then
       algorithm = MRRR

    else if ( leqi(algo, 'MRRR-2stage') .or. leqi(algo, 'RRR-2stage') .or. &
         leqi(algo, 'MRRR-2') .or. leqi(algo, 'RRR-2') .or. &
         leqi(algo, 'vr_2stage') ) then
# ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = MRRR_2stage
       else
          algorithm = MRRR
       end if
# else
       algorithm = MRRR
# endif
#endif

    else if ( leqi(algo, 'expert') .or. leqi(algo, 'vx') ) then
       algorithm = Expert
       
    else if ( leqi(algo, 'expert-2stage') .or. leqi(algo, 'expert-2') .or. &
         leqi(algo, 'vx_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = Expert_2stage
       else
          algorithm = Expert
       end if
#else
       algorithm = Expert
#endif
       
    else if ( leqi(algo, 'noexpert') .or. leqi(algo, 'qr') .or. &
         leqi(algo, 'v') ) then
       algorithm = QR

    else if ( leqi(algo, 'noexpert-2stage') .or. leqi(algo, 'noexpert-2') .or. &
         leqi(algo, 'qr-2stage') .or. leqi(algo, 'qr-2') .or. &
         leqi(algo, 'v_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = QR_2stage
       else
          algorithm = QR
       end if
#else
       algorithm = QR
#endif

    else

       write(*,'(a)') 'diag: Queried algorithm: '//trim(algo)
      
#ifndef SIESTA__MRRR
       write(*,'(a)') 'diag: Algorithm cannot be MRRR '// &
            '(not compiled with -DSIESTA__MRRR)'
#endif
#ifndef SIESTA__ELPA
       write(*,'(a)') 'diag: Algorithm cannot be ELPA '// &
            '(not compiled with -DSIESTA__ELPA)'
#endif

       call die('diag: Unknown routine requested for the diagonalization')

    end if

#ifdef SIESTA__ELPA
    elpa_use_gpu = fdf_get('Diag.ELPA.UseGPU',.false.)
#endif

    ! Retrieve tolerances for the expert drivers
    abstol = fdf_get('Diag.AbsTol', 1.e-16_dp)
    orfac = fdf_get('Diag.OrFac', 1.e-6_dp)

    ! Currently this is not used (it shouldn't be needed)
    mem_factor = fdf_get('Diag.Memory', 1.0_dp)
    mem_factor = max(mem_factor, 1.0_dp)

  end subroutine read_diag

  !> A convenience function to keep track of which solvers in 'diag' are
  !> capable of working with a reduced number of eigenvectors/eigenvalues
  
  function diag_recognizes_neigwanted() result (neig_capable)
    logical :: neig_capable
    
    neig_capable =  ( algorithm == MRRR .or. &
                     algorithm == MRRR_2stage .or. &
                     algorithm == Expert .or.  &
                     algorithm == Expert_2stage .or. &
                     algorithm == ELPA_1stage .or. &
                     algorithm == ELPA_2stage )
  end function diag_recognizes_neigwanted

  subroutine print_diag()
    use parallel, only: IONode, Nodes

    if ( .not. IONode ) return

    write(*,*) ! new-line
    
    select case ( algorithm )
    case ( DivideConquer ) 
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'D&C'
    case ( DivideConquer_2stage ) 
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'D&C-2stage'
    case ( MRRR )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'MRRR'
    case ( MRRR_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'MRRR-2stage'
    case ( ELPA_1stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'ELPA-1stage'
    case ( ELPA_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'ELPA-2stage'
    case ( Expert )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'Expert'
    case ( Expert_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'Expert-2stage'
    case ( QR )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'QR'
    case ( QR_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'QR-2stage'
    end select


#ifdef MPI
    write(*,'(a,t53,''= '',tr2,l1)') 'diag: Parallel over k', ParallelOverK
    write(*,'(a,t53,''= '',tr2,l1)') 'diag: Use parallel 2D distribution', Use2D
    write(*,'(a,t53,''= '',i0)') 'diag: Parallel block-size', diag_BlockSize
    if ( Use2D ) then
       write(*,'(a,t53,''= '',i5,'' x '',i5)') 'diag: Parallel distribution', &
            ProcessorY, max(1,Nodes / ProcessorY)
    else
       write(*,'(a,t53,''= '',i5,'' x '',i5)') 'diag: Parallel distribution', &
            1, Nodes
    end if
#endif

    if ( UpperLower == 'L' ) then
       write(*,'(a,t53,''= '',a)') 'diag: Used triangular part', 'Lower'
    else
       write(*,'(a,t53,''= '',a)') 'diag: Used triangular part', 'Upper'
    end if
    
    write(*,'(a,t53,''= '', e10.3)') 'diag: Absolute tolerance', abstol
    write(*,'(a,t53,''= '', e10.3)') 'diag: Orthogonalization factor', orfac

    write(*,'(a,t53,''= '',f7.4)') 'diag: Memory factor', mem_factor

  end subroutine print_diag
  
end module m_diag_option
