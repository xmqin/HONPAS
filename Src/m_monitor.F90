!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
!

! When dealing with high bias one often realizes that 
! the integral becomes rather murky.
! One can then do a monitor run and check which energy ranges
! account for the most noise in the signal. 
! This allows the user to customize the energy contour
! segments to optimize the convergence
! It will probably only be really useful for large bias, but should be rather
! easy to use.

module m_monitor

  use precision, only : dp

  implicit none

contains

  subroutine read_monitor(bName,block_dist,sparse,N,orb_list)

    use class_OrbitalDistribution
    use class_Sparsity
    use parallel, only : IONode, Node, Nodes
    use fdf
    use alloc
    use geom_helper, only : UCORB
#ifdef MPI
    use mpi_siesta
#endif

    ! the block for reading in the monitor list
    character(len=*), intent(in) :: bName    
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! the sparse pattern we wish to monitor
    type(Sparsity), intent(inout) :: sparse
    integer, intent(inout) :: N
    integer, pointer :: orb_list(:,:)

    integer, pointer :: tmp(:,:) => null()
#ifdef MPI
    integer, pointer :: tmpMPI(:,:) => null()
#endif

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: iM
    integer :: no_u, no_local
    integer :: io, j, ptr, ic, jc

    character(len=50) :: c
#ifdef MPI
    integer :: MPIerror
#endif

    ! This will read in the monitored orbitals
    if ( N /= 0 ) return

    ! if the block does not exist simply return
    if ( .not. fdf_block(trim(bName),bfdf) ) return

    ! The block we read is formatted like this:
    
    ! %block <name>
    !   # Onsite orbital
    !   200
    !   # Connection orbital (if it exists)
    !   # if it doesn't, nothing will happen
    !   200 210
    ! %endblock <name>

    ! First we read the number of entries in the block
    
    N = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnintegers(pline) > 0 ) then
          N = N + 1
       end if
    end do

    ! if it is empty, return
    if ( N == 0 ) return

    ! allocate space
    call re_alloc(tmp,1,3,1,N, routine='monitor')
    tmp(3,:) = -1

    ! read in the actual numbers
    ! first rewind...
    call fdf_brewind(bfdf)
    
    N = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnintegers(pline) > 0 ) then
          N = N + 1

          tmp(1,N) = fdf_bintegers(pline,1)
          if ( fdf_bnintegers(pline) > 1 ) then
             tmp(2,N) = fdf_bintegers(pline,2)
          else
             ! onsite...
             tmp(2,N) = tmp(1,N)
          end if
       end if
    end do

    ! done reading

    ! find the indices in the sparse pattern
    call attach(sparse,nrows=no_local,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    do io = 1 , no_local

       ! Get the global position
       ic = index_local_to_global(block_dist,io,Node)

       ! quick loop if not requested
       if ( all(tmp(1,:) /= ic) ) cycle

       ! Loop number of entries in the row...
       do j = 1 , l_ncol(io)

          ! The column
          jc = ucorb(l_col(l_ptr(io) + j),no_u)
          
          do iM = 1 , N
             if ( tmp(1,iM) == ic .and. &
                  tmp(2,iM) == jc ) then
                tmp(3,iM) = l_ptr(io) + j
             end if
          end do

       end do
    end do

    ! remove any dead references
    iM = count(tmp(3,:) > 0 )
    if ( iM == 0 ) then
       ! everything the user has requested is non-existing
       call de_alloc(tmp,routine='monitor')
       N = 0
       return
    end if
    if ( associated(orb_list) ) &
         call de_alloc(orb_list,routine='monitor')
    nullify(orb_list)
    call re_alloc(orb_list,1,iM,1,4, &
         routine='monitor')
    iM = 0
    do io = 1 , N
       if ( tmp(3,io) == -1 ) cycle
       iM = iM + 1
       ! copy over orb_list
       orb_list(iM,1) = tmp(1,io)
       orb_list(iM,2) = tmp(2,io)
       orb_list(iM,3) = 1
       orb_list(iM,4) = tmp(3,io)
    end do
    call de_alloc(tmp,routine='monitor')
    nullify(tmp)

    return

    ! TODO finalize this to be able to deal with a distributed 
    ! orbital
#ifdef MPI
    ! retrieve the number of orbitals on each node
    !call re_alloc(tmpMPI,0,Nodes-1)
    !tmpMPI = 0
    !tmpMPI(Node) = N
    ! The IONode should have the globalized array
    !call MPI_Reduce(N,iM,1,MPI_Integer,MPI_Sum,0, MPI_Comm_World, MPIerror)
    !call MPI_Reduce(MPI_InPlace,tmpMPI,Nodes,MPI_Integer, &
    !     MPI_Sum,0, MPI_Comm_World, MPIerror)
    
    !if ( IONode ) then
    !   call re_alloc(tmp,1,2,1,iM,routine='monitor')
    !end if
    
#endif
    
  end subroutine read_monitor

  
  function fname_monitor(iB,iC,slabel,suffix,basename) result(fname)
    integer, intent(in) :: iB, iC
    character(len=*), intent(in), optional :: slabel, suffix, basename
    character(len=200) :: fname
    ! Initialize the files
    if ( present(basename) ) then
       if ( iB == iC ) then
          write(fname,'(2a,i0)') trim(basename), &
               '_',iB
       else
          write(fname,'(a,2(a,i0))') trim(basename), &
               '_',iB,'_',iC
       end if
    else
       if ( .not. &
            ( present(slabel) .and. present(suffix) ) ) &
            call die('Error in filename input')
       if ( iB == iC ) then
          write(fname,'(4a,i0)') trim(slabel),'.', &
               trim(suffix),'_',iB
       else
          write(fname,'(3a,2(a,i0))') trim(slabel),'.', &
               trim(suffix),'_',iB,'_',iC
       end if
    end if
  end function fname_monitor
  
end module m_monitor

    
    
    

    
    
    


