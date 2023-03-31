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

! Module for converting a sparsity pattern to an
! "optimal" tri-diagonal matrix...

! Routine for converting a TranSiesta sparsity pattern to a tri-diagonal form
! It requires that the sparsity pattern is fully contained in the current
! processor.

! This module has been fully developed by Nick Papior Andersen, 2013
! nickpapior@gmail.com

! Please contact the author before utilization in other routines.

module m_ts_rgn2trimat

  ! Use regions...
  use precision, only : dp, i8b
  use m_region

  use m_ts_electype
  use m_ts_tri_common, only : GFGGF_needed_worksize
  use m_ts_tri_common, only : nnzs_tri_i8b

  ! method of BTD matrix
  use m_ts_method, only: TS_BTD_A_PROPAGATION, TS_BTD_A_COLUMN
  use m_ts_method, only: ts_A_method

  implicit none

  private

  public :: ts_rgn2trimat

  integer, parameter :: VALID = 0
  integer, parameter :: NONVALID_SIZE = 1
  integer, parameter :: NONVALID_ELEMENT_CONTAIN = 2
  integer, parameter :: NONVALID_TS_ELECTRODE = 3
  
contains

  ! IF parts == 0 will create new partition
  subroutine ts_rgn2TriMat(N_Elec, Elecs, IsVolt, &
       dit, sp, r, parts, n_part, method, last_eq, par)

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use parallel, only : IONode, Node, Nodes
    use fdf, only : fdf_get
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc, de_alloc
    use geom_helper, only: ucorb

    ! electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! Whether an entire column should be calculated
    logical, intent(in) :: IsVolt
    ! the distribution
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The region that we will create a tri-diagonal matrix on.
    type(tRgn), intent(in) :: r
    ! The sizes of the parts in the tri-diagonal matrix
    integer, intent(out) :: parts
    integer, pointer :: n_part(:)
    ! Which kind of method should be used to create the tri-diagonal
    integer, intent(in) :: method
    ! Whether we should retain the last partition to a fixed size.
    integer, intent(in) :: last_eq
    ! Whether the search should be performed in parallel or not
    logical, intent(in), optional :: par

    ! Local variables
    integer, pointer :: guess_part(:) => null()
    integer, pointer :: mm_col(:,:) => null()
    integer :: i, no, guess_parts, max_block
    ! In case of parallel
    integer :: guess_start, guess_step
    logical :: copy_first, lpar
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    character(len=64) :: fname
    integer :: io, jo, jr, j, ind, no_u, iu
#ifdef MPI
    integer :: MPIerror
#endif

    call timer('TS-rgn2tri',1)

    lpar = .true.
    if ( present(par) ) lpar = par
    if ( Nodes == 1 ) lpar = .false.

    ! This is the size of the regional 3-diagonal matrix
    no = r%n

    if ( no == 1 ) then
       
       ! Simple case, return immediately with default block-size
       
       call re_alloc(n_part , 1, 1, &
            routine='tsR2TM', name='n_part')
       n_part(1) = 1
       parts = 1

       call timer('TS-rgn2tri', 2)
       return
       
    end if

    ! Establish a guess on the partition of the tri-diagonal 
    ! matrix...
    call re_alloc(guess_part, 1, no, &
         routine='tsR2TM', name='guess_part')
    call re_alloc(n_part    , 1, no, &
         routine='tsR2TM', name='n_part')
    guess_part(:) = 0

    ! create array containing max-min for each ts-orbital
    call re_alloc(mm_col, 1, 2, 1, no, &
         routine='tsR2TM', name='mm_col')

    ! Set the min/max column indices in the pivoted matrix
    call set_minmax_col(sp, r, mm_col)

    parts = 2
    n_part(1) = no / 2
    n_part(2) = no / 2 + mod(no,2)
    if ( last_eq > 0 ) then
       parts = 2
       n_part(2) = last_eq
       n_part(1) = no - last_eq
       ! Initialize the guess for the 
       call guess_TriMat_last(no,mm_col,guess_parts,guess_part,last_eq)
       if ( valid_tri(no,mm_col,guess_parts, guess_part,last_eq) == VALID ) then
          parts = guess_parts
          n_part(1:parts) = guess_part(1:parts)
       end if
    end if

    if ( lpar ) guess_step = Nodes

    ! If the first one happens to be the best partition, 
    ! but non-valid, we need to make sure to overwrite it
    copy_first = .false. ! currently TODO THIS COULD BE A PROBLEM

    ! If the blocks are known by the user to not exceed a certain
    ! size, then we can greatly reduce the guessing step
    ! for huge systems
    ind = mm_col(2,1)
    ! Now figure out the min/max connections.
    ! We do this by cheking the first connection-block and the
    ! min/max ranges

    ! Here we find the orbital in the first range
    ! that connects to the fewest amount of parent orbitals.
    ! This means that the block *could* potentially be as small
    ! as the one found.
    io = no
    jo = 0
    do i = 1, ind
      ! bandwidth at row i
      j = mm_col(2, i) - mm_col(1, i)
      io = min(io, j)
      jo = max(jo, j)
    end do
    ! Allow 5% of the minimum block size +/- to search
    ! This *is* too much but to be on the safe side...
    ! On this note we also increase the step to 1%
    i = max(int(io * 0.05), 2)
    ! We allow the splitting of blocks to:
    !   [1 , 4]
    ! but this does not allow any
    !   [1 , >4]
    ! splittings.
    guess_start = max(1, io / 4 - i)
    guess_start = fdf_get('TS.BTD.Guess1.Min',guess_start)
#ifdef TBTRANS
    guess_start = fdf_get('TBT.BTD.Guess1.Min',guess_start)
#endif
    ! Define the stepping
    guess_step = max(int(io * 0.01), 1)
    max_block = max(min( no / 4, jo + i), 1)
    max_block = fdf_get('TS.BTD.Guess1.Max',max_block)
#ifdef TBTRANS
    max_block = fdf_get('TBT.BTD.Guess1.Max',max_block)
#endif
    ! In case the orbitals of this region is much smaller than
    ! max-block, then use the half 'no'
    max_block = max(max_block , guess_start + guess_step)
    max_block = min(max_block , no / 2)
    guess_start = max(min(guess_start, max_block), 1)

    ! Correct starting guess for the node
    if ( lpar ) guess_start = min(guess_start + Node, max_block)
    
    ! We loop over all possibilities from the first part having size
    ! 2 up to and including total number of orbitals in the 
    ! In cases of MPI we do it distributed (however, the collection routine
    ! below could be optimized)
    do i = guess_start , max_block , guess_step

      ! Make new guess...
      call guess_TriMat(no,mm_col,i,guess_parts,guess_part,last_eq)

      ! Quick escape if the first part is zero (signal from guess_trimat)
      if ( guess_part(1) == 0 ) cycle

      ! If not valid tri-pattern, simply jump...
      if ( valid_tri(no,mm_col,guess_parts, guess_part,last_eq) /= VALID ) then
        cycle
      end if

      ! Try and even out the different parts
      call full_even_out_parts(N_Elec,Elecs,IsVolt,method, &
          no,mm_col,guess_parts,guess_part,last_eq)

      if ( copy_first ) then
        ! ensure to copy it over (the initial one was not valid)
        copy_first = .false.
        parts = guess_parts
        n_part(1:parts) = guess_part(1:parts)
      else
        call select_better(method, parts,n_part, guess_parts, guess_part)
      end if

    end do

#ifdef MPI
    if ( lpar ) then
       ! Select the most optimal partition scheme...
       ! Only check up-till the largest block that is actually searched
       do i = 0 , Nodes - 1
          if ( i == Node ) then
             call MPI_Bcast(parts, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(n_part(1), parts, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
          else
             call MPI_Bcast(guess_parts, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(guess_part(1), guess_parts, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             ! Only all the other nodes are allowed to check...
             ! TODO this could be made to a communication tree to limit communication
             call select_better(method, &
                  parts,n_part, guess_parts, guess_part)
          end if
       end do
    end if
#endif

    call de_alloc(guess_part,routine='tsR2TM',name='guess_part')
    ! Shrink to the found parts
    call re_alloc(n_part,1, parts, copy=.true., shrink=.true., &
         routine='tsR2TM', name='n_part')

    if ( parts < 2 ) then
       
       if ( IONode ) then 
          write(*,'(a)') 'Could not determine an optimal tri-diagonalization &
               &partition'
          write(*,'(2a)') 'Running on region: ',trim(r%name)
          if ( parts > 0 ) then
             write(*,'(a,i0)') 'Found: ',parts
             write(*,'(1000000(tr1,i0))') n_part
          else
             write(*,'(a)') 'None found...'
          end if
       end if
       call re_alloc(n_part, 1, 3, routine='tsR2TM',name='n_part')
       call die('Not yet implemented')

    end if

    ! The parts now have almost the same size and we will check that it
    ! is a valid thing, if not, we will revert to the other method of
    ! creating the tri-diagonal sparsity pattern
    ! We do not expect this to fail. After all we check that we can
    ! even out the partitions.
    ! The most probable thing is that the electrodes are not
    ! contained in the first two parts.
    i = valid_tri(no,mm_col,parts, n_part,last_eq)
    if ( i /= VALID ) then
       write(*,'(2a)') 'Running on region: ',trim(r%name)
       write(*,'(a,i0)') 'TranSiesta system size: ',no
       write(*,'(a,i0)') 'Current parts: ',parts
       write(*,'(10000000(tr1,i0))') n_part
       write(*,'(a,i0)') 'Current part size: ',sum(n_part(:))
       select case ( i )
       case ( NONVALID_SIZE )
          write(*,'(a)') 'The size is not valid.'
       case ( NONVALID_ELEMENT_CONTAIN ) 
          write(*,'(a)') 'Some elements are not contained.'
       case ( NONVALID_TS_ELECTRODE ) 
          write(*,'(a)') 'The electrode is not fully encompassed.'
       case default
          write(*,'(a,i0,a)') 'Row ',-i,' not encompassed in the tri-matrix'
       end select
       call die('Contact the developers. (missing implementation). &
            &You appear to have a special form of electrode.')
    end if

    call de_alloc(mm_col,routine='tsR2TM',name='mm_col')

    call timer('TS-rgn2tri',2)

    if ( .not. IONode ) return

    fname = fdf_get('TS.BTD.Output',' ')
    if ( len_trim(fname) == 0 ) return
    
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows_g=no_u)

    ! Write out the BTD format in a file to easily be processed
    ! by python, this is the pivoted sparsity pattern
    call io_assign(iu)
    open(iu, file=trim(fname)//'.sp',action='write')
    write(iu,'(i0)') no
    do i = 1 , no
       io = r%r(i)
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          jo = UCORB(l_col(ind),no_u)
          jr = rgn_pivot(r,jo)
          if ( jr > i ) cycle ! only print lower half
          if ( jr <= 0 ) cycle
          write(iu,'(2(i0,tr1),i1)') i, jr, 1
       end do
    end do
    call io_close(iu)

    call io_assign(iu)
    open(iu, file=trim(fname)//'.pvt',action='write')
    write(iu,'(i0)') no
    do i = 1 , no
       ! store the pivoting table such that sp[i,j] = orig[ pvt[i] , pvt[j] ]
       write(iu,'(i0,tr1,i0)') i, r%r(i)
    end do
    call io_close(iu)

    call io_assign(iu)
    open(iu, file=trim(fname)//'.btd',action='write')
    ! write the tri-mat blocks
    write(iu,'(i0)') parts
    do i = 1 , parts
       write(iu,'(i0)') n_part(i)
    end do
    call io_close(iu)

  contains 

    subroutine select_better(method, cur_parts, cur_part, guess_parts, guess_part)
      
      integer, intent(in)    :: method
      integer, intent(inout) :: cur_parts
      integer, intent(in)    :: guess_parts
      integer, intent(inout) :: cur_part(max(cur_parts,guess_parts))
      integer, intent(in)    :: guess_part(guess_parts)
      logical :: copy
      integer :: cur_work, guess_work
      integer :: cur_pad, guess_pad

      ! We check whether the number of elements is smaller
      ! or that the number of parts is greater (however, this should
      ! in principle always go together)
      ! If the method of optimization is memory:
      if ( method == 0 ) then
        
        copy = faster_parts(cur_parts,cur_part,guess_parts,guess_part)
        
      else if ( method == 1 ) then
        
        ! We optimize for memory, i.e. we check for number of elements
        ! in this regard we also check whether we should allocate
        ! a work-array in case of bias calculations.
        if ( IsVolt .and. ts_A_method == TS_BTD_A_COLUMN ) then
          call GFGGF_needed_worksize(guess_parts, guess_part, &
              N_Elec, Elecs, guess_pad, guess_work)
          call GFGGF_needed_worksize(cur_parts, cur_part, &
              N_Elec, Elecs, cur_pad, cur_work)

          ! Update values
          guess_work = guess_work + guess_pad
          cur_work = cur_work + cur_pad

          ! Calculate difference between old and new
          cur_pad = cur_work - guess_work
        else
          ! No default values
          cur_pad = 0
        end if

        ! Difference between the two BTD matrices
        guess_pad = nnzs_tri_i8b(cur_parts, cur_part) - &
            nnzs_tri_i8b(guess_parts, guess_part)

        ! total difference in number of elements
        ! If this is positive the guessed BTD matrix has fewer
        ! memory elements allocated than the currently selected one
        cur_pad = guess_pad + cur_pad

        copy = cur_pad > 0
        if ( .not. copy ) then
          ! in case the work-size is the same we fall-back to the fastests method
          if ( cur_pad == 0 ) then
            copy = faster_parts(cur_parts,cur_part,guess_parts,guess_part)
          end if
        end if

      else
        call die('Unknown optimization scheme for the tri-mat')
      end if

      if ( copy ) then
        cur_parts = guess_parts
        cur_part(1:cur_parts) = guess_part(1:cur_parts)
      end if

    end subroutine select_better

  end subroutine ts_rgn2TriMat

  subroutine guess_TriMat(no,mm_col,first_part,n_part,parts,last_eq)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(in) :: first_part
    integer, intent(inout) :: n_part
    integer, intent(inout) :: parts(no)
    integer, intent(in) :: last_eq

    ! Local variables
    integer :: N

    if ( first_part > no ) &
        call die('Not allowed to do 1 tri-diagonal part')

    n_part = 1
    parts(1) = first_part
    N = parts(1)
    
    do while ( N < no )

      ! Step the currently searched part
      n_part = n_part + 1
      
      if ( n_part > no ) then
        print *,'Error',n_part,no
        call die('Size error when guessing the tri-mat size')
      end if
      
      call guess_next_part_size(no, mm_col, N, n_part, parts)
      
      N = N + parts(n_part)
      ! if a last-part was "forced" we do this here...
      if ( N + last_eq > no ) then
        ! We need to add the former part with the "too many"
        ! orbitals
        parts(n_part) = parts(n_part) + no - N
        N = no
      end if

    end do

    if ( last_eq > 0 ) then
      
      ! Correct so that we actually do contain the last_eq
      ! in the last one
      N = parts(n_part) - last_eq
      parts(n_part) = last_eq
      parts(n_part-1) = parts(n_part-1) + N
      
      ! Signal that this is a faulty TRIMAT
      if ( N < 0 ) parts(1) = 0
      
    end if

  end subroutine guess_TriMat

  subroutine guess_TriMat_last(no,mm_col,parts,n_part,last_eq)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(out) :: parts
    integer, intent(out) :: n_part(:)
    integer, intent(in) :: last_eq

    ! Local variables
    integer :: N

    parts = 1
    n_part(1) = last_eq
    N = n_part(1)
    do while ( N < no )
       parts = parts + 1
       if ( parts > size(n_part) ) then
          print *,'Error',parts,size(n_part)
          call die('Size error when guessing the tri-mat size')
       end if
       call guess_prev_part_size(no, mm_col, parts, parts, n_part)
       N = N + n_part(parts)
    end do

    ! Reverse
    n_part(1:parts) = n_part(parts:1:-1)

  end subroutine guess_TriMat_last

  function faster_parts(np,n_part,ng,g_part) result(faster)
    integer, intent(in) :: np, n_part(np)
    integer, intent(in) :: ng, g_part(ng)
    logical :: faster

    integer :: i, n
    real(dp) :: p_N, g_N, diff
    real(dp), parameter :: off_diag = 1.2_dp

    ! We estimate the fastest algorithm
    ! by the number of operations the matrices make

    ! 2 * 5 / 3 + 4 ~~ 1 + 1.2
    p_N = R(n_part(1)) ** 2 + off_diag * R(n_part(2))
    p_N = p_N * R(n_part(1))
    
    g_N = R(g_part(1)) ** 2 + off_diag * R(g_part(2))
    g_N = g_N * R(g_part(1))
    
    diff = p_N - g_N
    
    n = min(np, ng)
!$OMP parallel do default(shared), private(i,p_N,g_N), reduction(+:diff), if(n>1000)
    do i = 2, n

      p_N = R(n_part(i)) ** 2 + off_diag * R(n_part(i-1))
      p_N = p_N * R(n_part(i))

      g_N = R(g_part(i)) ** 2 + off_diag * R(g_part(i-1))
      g_N = g_N * R(g_part(i))
      
      diff = diff + p_N - g_N
      
    end do
!$OMP end parallel do

    if ( np > ng ) then
      faster = .true.

      do i = ng + 1, np
        
        p_N = R(n_part(i)) ** 2 + off_diag * R(n_part(i-1))
        p_N = p_N * R(n_part(i))

        diff = diff + p_N

        if ( diff > 0._dp ) return

      end do
      
    else
      faster = .false.
      
      do i = np + 1, ng
        
        g_N = R(g_part(i)) ** 2 + off_diag * R(g_part(i-1))
        g_N = - g_N * R(g_part(i))
        
        diff = diff + g_N

        if ( diff < 0._dp ) return

      end do
      
    end if

    faster = diff > 0._dp

  contains

    elemental function R(i) result(o)
      integer, intent(in) :: i
      real(dp) :: o
      o = real(i, dp)
    end function R
      
  end function faster_parts


  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_next_part_size(no,mm_col,no_cur,n_part,parts)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    ! no_cur is sum(parts(:-1))
    integer, intent(in) :: no_cur, n_part
    integer, intent(inout) :: parts(n_part)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    eRow = no_cur
    sRow = eRow - parts(n_part-1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    mcol = 0
    do i = sRow, eRow
      ! this is the # of elements from the RHS of the 'part-1'
      ! part of the tridiagonal matrix and out to the last element of
      ! this row...
      mcol = max(mcol, mm_col(2, i))
    end do

    ! In case there is no connection, we should
    ! force the next-part to be 1!
    parts(n_part) = max(1, mcol - eRow)

  end subroutine guess_next_part_size

  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_prev_part_size(no,mm_col,part,parts,n_part)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    eRow = no
    if ( part > 2 ) then
       do i = 1 , part - 2
          eRow = eRow - n_part(i)
       end do
    end if
    sRow = eRow - n_part(part-1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'part-1'
       ! part of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = sRow - mm_col(1,i)
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do

    ! In case there is actually no connection, we should
    ! force the next-part to be 1!
    if ( n_part(part) == 0 ) then
       n_part(part) = 1
    end if

  end subroutine guess_prev_part_size

  subroutine full_even_out_parts(N_Elec,Elecs,IsVolt, &
      method,no,mm_col,n_part,parts, last_eq)

    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    logical, intent(in) :: IsVolt
    integer, intent(in) :: method ! the method used for creating the parts
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: n_part
    integer, intent(inout) :: parts(n_part)
    integer, intent(in) :: last_eq
    ! Local variables
    integer, allocatable :: mem_parts(:), cum_parts(:)
    integer :: n, d
    logical :: changed

    if ( n_part == 1 ) return

    allocate(mem_parts(n_part), cum_parts(n_part))

    ! Setup the cumultative parts
    cum_parts(1) = parts(1)
    mem_parts(1) = parts(1)
    do n = 2, n_part
      cum_parts(n) = parts(n) + cum_parts(n-1)
      mem_parts(n) = parts(n)
    end do

    select case ( method )
    case ( 0 ) ! speed

      do
        changed = .false.
        do n = 2 , n_part - 1
          call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
          call diff_perf(n, n_part, parts, mem_parts, d)
          call select()
        end do
        n = 1
        call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
        call diff_perf(n, n_part, parts, mem_parts, d)
        call select()
        n = n_part
        call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
        call diff_perf(n, n_part, parts, mem_parts, d)
        call select()
        
        if ( .not. changed ) exit
      end do

    case ( 1 ) ! memory
      
      do
        changed = .false.
        do n = 2 , n_part - 1
          call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
          call diff_mem(n, n_part, parts, mem_parts, d)
          call select()
        end do
        n = 1
        call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
        call diff_mem(n, n_part, parts, mem_parts, d)
        call select()
        n = n_part
        call even_out_parts(no, mm_col, n_part, parts, cum_parts, n, last_eq)
        call diff_mem(n, n_part, parts, mem_parts, d)
        call select()
        
        if ( .not. changed ) exit
      end do
      
    case default

      call die('not implemented')

    end select

    deallocate(mem_parts, cum_parts)

  contains

    subroutine select()
      if ( d == 0 ) then
        ! Simply store. It could be that we swapped a few things
        call store_part(n_part, parts, mem_parts, n)
        call store_cum_part(n_part, parts, cum_parts, n)
      else if ( d > 0 ) then
        ! Copy back
        call store_part(n_part, mem_parts, parts, n)
        ! the cumultative parts are not changed since we restore them
      else
        call store_part(n_part, parts, mem_parts, n)
        call store_cum_part(n_part, parts, cum_parts, n)
        changed = .true.
      end if
    end subroutine select

    subroutine store_part(n_part, parts, store_parts, n)
      integer, intent(in) :: n, n_part
      integer, intent(in) :: parts(n_part)
      integer, intent(inout) :: store_parts(n_part)
      
      store_parts(n) = parts(n)
      if ( n == 1 ) then
        store_parts(n+1) = parts(n+1)
      else if ( n == n_part ) then
        store_parts(n-1) = parts(n-1)
      else
        store_parts(n-1) = parts(n-1)
        store_parts(n+1) = parts(n+1)
      end if
      
    end subroutine store_part

    subroutine store_cum_part(n_part, parts, cum_parts, n)
      integer, intent(in) :: n, n_part
      integer, intent(in) :: parts(n_part)
      integer, intent(inout) :: cum_parts(n_part)

      if ( n == 1 ) then
        cum_parts(n) = parts(n)
        cum_parts(n+1) = cum_parts(n) + parts(n+1)
      else if ( n == n_part ) then
        if ( n_part == 2 ) then
          cum_parts(n-1) = parts(n-1)
        else
          cum_parts(n-1) = cum_parts(n-2) + parts(n-1)
        end if
        cum_parts(n) = cum_parts(n-1) + parts(n)
      else if ( n == 2 ) then
        cum_parts(n-1) = parts(n-1)
        cum_parts(n) = cum_parts(n-1) + parts(n)
        cum_parts(n+1) = cum_parts(n) + parts(n+1)
      else
        cum_parts(n-1) = cum_parts(n-2) + parts(n-1)
        cum_parts(n) = cum_parts(n-1) + parts(n)
        cum_parts(n+1) = cum_parts(n) + parts(n+1)
      end if

    end subroutine store_cum_part
    
    function changed_part(n_part, parts, store_parts, n) result(changed)
      integer, intent(in) :: n, n_part
      integer, intent(in) :: parts(n_part)
      integer, intent(inout) :: store_parts(n_part)
      logical :: changed
      
      changed = store_parts(n) /= parts(n)
      if ( n == 1 ) then
        changed = changed .or. store_parts(n+1) /= parts(n+1)
      else if ( n == n_part ) then
        changed = changed .or. store_parts(n-1) /= parts(n-1)
      else
        changed = changed .or. store_parts(n-1) /= parts(n-1)
        changed = changed .or. store_parts(n+1) /= parts(n+1)
      end if
      
    end function changed_part

  end subroutine full_even_out_parts

  subroutine even_out_parts(no,mm_col,n_part,parts, cum_parts, n, last_eq)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: n_part, cum_parts(n_part)
    integer, intent(inout) :: parts(n_part)
    integer, intent(in) :: n, last_eq
    
    ! Local variables
    integer :: copy_part
    integer :: sRow, eRow

    if ( last_eq > 0 ) then
      
      ! We need the last one to be of a certain size.
      if ( n >= n_part - 1 ) return
      
    end if

    select case ( n_part )
    case ( 1 )
      call die('You cannot use tri-diagonalization &
          &without having at least 2 parts')
    case ( 2 )
      ! For only two blocks there is no gain/loss regardless of matrix
      ! So set them equal
      parts(1) = no / 2
      parts(2) = parts(1) + mod(no, 2)
      return
    end select

    ! Initialize
    copy_part = 0

    ! We do not allow to diminish the edges
    ! This is because we need additional checks of their
    ! regions, hence, we let the central parts partition out
    ! the regions
    if ( n == 1 ) then
      
      ! We can always align the two blocks (sign is pointless here)
      sRow = 0
      do while ( parts(1) /= copy_part )
        copy_part = parts(1)
        call even_if_larger(sRow,parts(1),parts(2),sign=-1)
      end do

      return

    else if ( n == n_part ) then

      ! We can always align the two blocks (sign is pointless here)
      sRow = 0
      do while ( parts(n_part) /= copy_part )
        copy_part = parts(n_part)
        call even_if_larger(sRow,parts(n_part),parts(n_part-1),sign=1)
      end do

      return

    end if

    ! Find min/max columns checked
    sRow = cum_parts(n-1) + 1
    eRow = cum_parts(n)

    ! We will continue to shift columns around
    ! until we do not shift columns any more...
    copy_part = 0
    do while ( parts(n) /= copy_part )
        
       ! Copy the current partition so that we can check in the
       ! next iteration...
       copy_part = parts(n)

       ! TODO - consider adding a flag for memory reduced utilization of TRI
       ! this will require the two electrodes parts to be larger
       ! than the other parts...

       ! 1. if you wish to shrink it left, then:
       !    the first row must not have any elements
       !    extending into the right part
       if ( mm_col(2,sRow) <= eRow ) then
          call even_if_larger(sRow,parts(n),parts(n-1),sign= 1)
       end if

       ! 2. if you wish to shrink it right, then:
       !    the last row must not have any elements
       !    extending into the left part
       if ( sRow <= mm_col(1,eRow) ) then
          call even_if_larger(eRow,parts(n),parts(n+1),sign=-1)
       end if

    end do
     
  contains

    pure subroutine even_if_larger(Row,p1,p2,sign)
      integer , intent(inout) :: Row, p1, p2
      integer, intent(in) :: sign
      ! If we don't have p2 + 1 we could end up
      ! in 12 -- 11 swaps
      if ( p1 > p2 + 1) then
         p1 = p1 - 1
         p2 = p2 + 1
         Row = Row + sign
      end if
    end subroutine even_if_larger

  end subroutine even_out_parts

! Min and max column requires that the sparsity pattern
! supplied has already stripped off the buffer orbitals.
! Otherwise this will fail
  subroutine set_minmax_col(sp, r, mm_col)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    type(tRgn), intent(in) :: r
    integer, intent(out) :: mm_col(2,r%n)
    ! The results
    type(tRgn) :: pvt
    integer :: ir, row, ptr, nr, j
    integer, pointer :: l_col(:), l_ptr(:), ncol(:)

    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col,nrows_g=nr)

    ! Using a pivoting table reduces overhead
    ! of performing rgn_pivot on a non-sorted
    ! region! SUBSTANTIALLY!
    call rgn_init(pvt, nr, val=0)

!$OMP parallel default(shared), private(ir,row,ptr,j)

    ! Create the back-pivoting array
!$OMP do
    do ir = 1 , r%n
      pvt%r(r%r(ir)) = ir
    end do
!$OMP end do
    
!$OMP do
    do ir = 1 , r%n

       ! Get original sparse matrix row
       row = r%r(ir)
       
       ! initialize to region row
       mm_col(1,ir) = ir
       mm_col(2,ir) = ir

       ! Loop on the sparse entries
       do ptr = l_ptr(row) + 1 , l_ptr(row) + ncol(row)
          j = pvt%r( ucorb(l_col(ptr),nr) )
          if ( j > 0 ) then
             if ( j < mm_col(1,ir) ) mm_col(1,ir) = j
             if ( j > mm_col(2,ir) ) mm_col(2,ir) = j
          end if
       end do

    end do
!$OMP end do nowait
    
!$OMP end parallel

    call rgn_delete(pvt)

  end subroutine set_minmax_col

  !< Calculate the difference in memory requirement for two BTD
  !< part1 - part2 == mem
  recursive subroutine diff_mem(part, n, part1, part2, mem, final)
    integer, intent(in) :: part
    integer, intent(in) :: n
    integer, intent(in) :: part1(n), part2(n)
    integer, intent(out) :: mem
    logical, intent(in), optional :: final
    integer :: next
    logical :: lfinal

    mem = 0
    if ( part < 1 ) return
    if ( part > n ) return
    lfinal = .false.
    if ( present(final) ) lfinal = final

    if ( part == 1 ) then
      
      mem = part1(1) * ( part1(1) + part1(2) * 2 ) - &
          part2(1) * ( part2(1) + part2(2) * 2 )

      if ( .not. lfinal ) then
        call diff_mem(part+1, n, part1, part2, next, .true.)
        mem = mem + next
      end if
      
    else if ( part == n ) then
      
      mem = part1(n) * ( part1(n) + part1(n-1) * 2 ) - &
          part2(n) * ( part2(n) + part2(n-1) * 2 )

      if ( .not. lfinal ) then
        call diff_mem(part-1, n, part1, part2, next, .true.)
        mem = mem + next
      end if
      
    else
      
      mem = part1(part) * ( part1(part) + part1(part-1) * 2 + &
          part1(part+1) * 2 ) + part1(part-1) ** 2 + part1(part+1) ** 2 &
          - ( &
          part2(part) * ( part2(part) + part2(part-1) * 2 + &
          part2(part+1) * 2 ) + part2(part-1) ** 2 + part2(part+1) ** 2)

      if ( .not. lfinal ) then
        call diff_mem(part-1, n, part1, part2, next, .true.)
        mem = mem + next
        call diff_mem(part+1, n, part1, part2, next, .true.)
        mem = mem + next
      end if
      
    end if

  end subroutine diff_mem

  !< Calculate the difference in memory requirement for two BTD
  !< part1 - part2 == mem
  recursive subroutine diff_perf(part, n, part1, part2, perf, final)
    use precision, only: dp
    integer, intent(in) :: part
    integer, intent(in) :: n
    integer, intent(in) :: part1(n), part2(n)
    integer, intent(out) :: perf
    logical, intent(in), optional :: final
    real(dp) :: one_third = 0.3333333333333333333333_dp
    real(dp) :: p
    integer :: next
    logical :: lfinal

    perf = 0
    if ( part < 1 ) return
    if ( part > n ) return
    lfinal = .false.
    if ( present(final) ) lfinal = final

    if ( part == 1 ) then

      p = p_inv(part1(1)) - p_inv(part2(1)) + & ! inv
          (p_mm(part1(1), part1(2)) - p_mm(part2(1), part2(2))) ! mm
      
      if ( .not. lfinal ) then
        ! This takes inv of part + 1 and mm 
        call diff_perf(part+1, n, part1, part2, next, .true.)
        p = p + real(next, dp) ** 3
      end if
      
    else if ( part == n ) then

      p = p_inv(part1(n)) - p_inv(part2(n)) + & ! inv
          (p_mm(part1(n), part1(n-1)) - p_mm(part2(n), part2(n-1))) ! mm

      if ( .not. lfinal ) then
        ! This takes inv of part - 1 and mm 
        call diff_perf(part-1, n, part1, part2, next, .true.)
        p = p + real(next, dp) ** 3
      end if
      
    else

      p = p_inv(part1(part)) - p_inv(part2(part)) + & ! inv
          (p_mm(part1(part), part1(part-1)) - p_mm(part2(part), part2(part-1))) + & ! mm
          (p_mm(part1(part), part1(part+1)) - p_mm(part2(part), part2(part+1))) ! mm
      
      if ( .not. lfinal ) then
        call diff_perf(part-1, n, part1, part2, next, .true.)
        p = p + real(next, dp) ** 3
        call diff_perf(part+1, n, part1, part2, next, .true.)
        p = p + real(next, dp) ** 3
      end if
      
    end if
    
    ! This should remove possible overflows
    ! We could essentially also just return +1/0/-1
    ! But perhaps we can use the actual value to something useful?
    if ( p < 0 ) then
      perf = - int( (-p) ** one_third )
    else
      perf = int(p ** one_third)
    end if
      
  contains

    pure function p_inv(s) result(p)
      use precision, only: dp
      integer, intent(in) :: s
      real(dp) :: p
      p = real(s, dp) ** 3
    end function p_inv

    pure function p_mm(m, n) result(p)
      use precision, only: dp
      integer, intent(in) :: m, n
      real(dp) :: p
      p = real(m, dp) ** 2 * real(n, dp)
    end function p_mm

  end subroutine diff_perf
  
  function valid_tri(no,mm_col,n_part,parts,last_eq) result(val)
    integer, intent(in) :: no, mm_col(2,no)
    integer, intent(in) :: n_part, parts(n_part), last_eq
    integer :: val
    ! Local variables
    integer :: i, N, Nm1, Np1, col_min, col_max
    logical :: first

    if ( last_eq > 0 ) then
      
      ! We do not allow something to not end in the requested number of orbitals.
      if ( parts(n_part) /= last_eq ) then
        val = NONVALID_SIZE
        return
      end if
      
      if ( parts(1) == 0 ) then
        val = NONVALID_SIZE
        return
      end if
      
    end if

    ! Default to a valid BTD
    val = VALID

    ! Check that every element is contained in the 
    ! tri-diagonal matrix...
    first = .true.
    Nm1 = 1
    Np1 = parts(1) + parts(2)

    ! Find min/max for first block
    call find_min_max(no, mm_col, 1, parts(1), col_min, col_max)
    if ( col_min < Nm1 ) then
      ! If this ever occur it suggests that the 
      ! sparsity pattern is not fully symmetric !
      i = 1
      call print_error(first)
      val = - 1
    else if ( Np1 < col_max ) then
      i = 1
      call print_error(first)
      val = - 1
    end if

    ! Initialize loop counters
    N = parts(1) + 1
    do i = 2 , n_part - 1

      ! Find min/max for this block
      call find_min_max(no, mm_col, N, Np1, col_min, col_max)

      ! Update the size of the part after this
      Np1 = Np1 + parts(i+1)

      if ( col_min < Nm1 ) then
        ! If this ever occur it suggests that the 
        ! sparsity pattern is not fully symmetric !
        call print_error(first)
        val = - i
        
      else if ( Np1 < col_max ) then
        call print_error(first)
        val = - i
        
      end if

      ! Update loop
      N = N + parts(i)
      
      ! Update the previous part
      Nm1 = Nm1 + parts(i-1)

    end do

    ! Find min/max for the last block
    call find_min_max(no, mm_col, N, Np1, col_min, col_max)
    if ( col_min < Nm1 ) then
      ! If this ever occur it suggests that the 
      ! sparsity pattern is not fully symmetric !
      i = n_part
      call print_error(first)
      val = - n_part
      
    else if ( Np1 < col_max ) then
      i = n_part
      call print_error(first)
      val = - n_part
      
    end if
    if ( val /= valid ) return
    
    ! Update total number of orbitals
    N = N + parts(n_part) - 1

    ! Size of tri-matrix is already calculated
    ! Easy check, if the number of rows
    ! does not sum up to the total number of rows.
    ! Then it must be invalid...
    if ( N /= no ) then
      val = NONVALID_SIZE
    end if
    
  contains

    pure subroutine find_min_max(no, mm_col, N1, N2, col_min, col_max)
      integer, intent(in) :: no, mm_col(2,no), N1, N2
      integer, intent(out) :: col_min, col_max
      integer :: ir
      col_min = mm_col(1,N1)
      col_max = mm_col(2,N1)
      do ir = N1 + 1, N2
        col_min = min(col_min, mm_col(1,ir))
        col_max = max(col_max, mm_col(2,ir))
      end do
    end subroutine find_min_max

    subroutine print_error(first)
      use parallel, only : IONode
      logical, intent(inout) :: first
      if ( IONode ) then
        if ( first ) then
          write(*,'(a)') 'BTD: Found non-symmetric matrix!'
          write(*,'(a)') '     If you are not using delta methods you are probably doing something wrong!'
          write(*,'(a)') ' block: block row is located in'
          write(*,'(a)') ' N    : size of block'
          write(*,'(a)') ' min_C: minimum element row connects to'
          write(*,'(a)') ' min_B: minimum element in the block'
          write(*,'(a)') ' max_B: maximum element in the block'
          write(*,'(a)') ' max_C: maximum element row connects to'
          write(*,'(6a8)') 'block', 'N', 'min_C', 'min_B', 'max_B', 'max_C'
          first = .false.
        end if
        write(*,'(6i8)') i, parts(i), col_min, Nm1, Np1, col_max
     end if
    end subroutine print_error

  end function valid_tri

end module m_ts_rgn2trimat

