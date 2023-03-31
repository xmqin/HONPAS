!     
!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996- .
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     
      module m_matio

      implicit none

      private

      PUBLIC :: write_mat

CONTAINS

!------------------------------------------------------
  subroutine write_mat (maxnd, no_l, nspin, &
       numd, listdptr, listd, mat, userfile, compatible)

#ifdef MPI
    use mpi_siesta
    use parallel, only: blocksize, SIESTA_comm
    use alloc
#endif

    integer, parameter :: dp = selected_real_kind(10,100)


    integer, intent(in) :: maxnd
    integer, intent(in) :: no_l
    integer, intent(in) :: nspin
    integer, intent(in) :: numd(1:no_l)
    integer, intent(in) :: listdptr(1:no_l)
    integer, intent(in) :: listd(maxnd)
    real(dp), intent(in) :: mat(maxnd,nspin)

    character(len=*), intent(in), optional :: userfile
    ! Whether to create files with one-orbital data blocks,
    ! compatible with the current reading routines
    ! The default is to do so
    logical, intent(in), optional          :: compatible

    integer :: no_u, m, ml, im, ndmaxg, lun, is
    integer :: n_l, norbs_l, n_g, norbs_g, node, myrank, nprocs
    integer :: nnzbs, base, maxnnzbs, nblocks, norbs, nsize
    integer :: base_l, nnzs_bl, nnzs_bg
    integer :: i, j, ptr, nnzsi
#ifdef MPI
    integer  :: MPIerror, stat(MPI_STATUS_SIZE)
    real(dp), dimension(:), pointer :: buffer => null()
    integer,  dimension(:), pointer :: ibuffer => null()
    integer, dimension(:), pointer  :: numdg => null()
#endif

    character(len=256) :: filename
    logical            :: bck_compat

    call timer("WriteMat",1)

#ifndef MPI

    ! Serial mode.
       call io_assign(lun)
       open( lun, file=trim(filename), form="unformatted", status='unknown' )
       rewind(lun)
       no_u = no_l
       write(lun) no_u, nspin       
       write(lun) (numd(m),m=1,no_u)
       do i = 1, no_u
          write(lun) (listd(j),j=listdptr(i)+1,listdptr(i)+numd(i))
       enddo
       do is=1,nspin
          do i = 1, no_u
             write(lun) (mat(j,is),j=listdptr(i)+1,listdptr(i)+numd(i))
          enddo
       enddo
       call io_close(lun)

#else  
    call MPI_Comm_Size( SIESTA_Comm, nprocs, MPIerror )
    call MPI_Comm_Rank( SIESTA_Comm, myrank, MPIerror )

!     Find total number of orbitals over all Nodes
    call MPI_Reduce(no_l,no_u,1,MPI_integer,MPI_sum,0,SIESTA_Comm,MPIerror)

    if (myrank.eq.0) then
       if (.not. present(userfile)) then
          filename = "SPMAT"
       else
          filename = userfile
       endif
       if (.not. present(compatible)) then
          bck_compat = .true.
       else
          bck_compat = compatible
       endif

       call io_assign(lun)
       open( lun, file=trim(filename), form="unformatted", status='unknown' )
       rewind(lun)
       if (bck_compat) then
          write(lun) no_u, nspin
       else
          write(lun) no_u, nspin, blocksize
       endif
       
       call re_alloc( numdg, 1, no_u, 'numdg', 'write_mat' )
    endif

!     Get info about numd
    n_g = 0
    n_l = 0
    node = -1
    DO

       node = node + 1
       if (node == nprocs) node = 0

       !print *, " node: ", node, " myrank: ", myrank

       if (myrank == node) then
          norbs_l = min(blocksize,no_l-n_l)
          !print *, "myrank: ", myrank, " processing size:", norbs_l
          if (node==0) then
             norbs_g = min(blocksize,no_u-n_g)
             !print *, "myrank: ", myrank, " will just copy: (l,g)", &
             !                      norbs_l, norbs_g
             numdg(n_g+1:n_g+norbs_g) = numd(n_l+1:n_l+norbs_l)
             n_g = n_g + norbs_g
             !print *, "root has received so far: ", n_g
          else
             !print *, "myrank: ", myrank, " will send to 0: ", norbs_l
             call MPI_Send(numd(n_l+1),norbs_l,MPI_integer, &
                  0,1,SIESTA_Comm,MPIerror)
             !print *, "myrank: ", myrank, " completed send ", norbs_l
          endif
          n_l = n_l + norbs_l

       else if (myrank == 0) then
          norbs_g = min(blocksize,no_u-n_g)
          !print *, "root will receive from ", node, " size: ", norbs_g
          !print *, "will put it starting at: ", n_g+1

          call MPI_Recv(numdg(n_g+1:),norbs_g,MPI_integer, &
                node,1,SIESTA_Comm,stat,MPIerror)
          n_g = n_g + norbs_g
          !print *, "root has received so far: ", n_g
       endif

       if (myrank == 0) then
          if (n_g == no_u) EXIT
       else
          if (n_l == no_l) then
             !print *, "rank ", myrank, " exiting loop"
             EXIT
          endif
       endif
          
    enddo
          
!     Write out numd array
      if (myrank.eq.0) then
         write(lun) (numdg(m),m=1,no_u)
      endif

!     Find out how big the buffer has to be
      if (myrank.eq.0) then
         maxnnzbs = 0
         nblocks = 0
         norbs = 0
         do 
            base = nblocks*blocksize
            nsize = min(blocksize,no_u-norbs)
            nnzbs = sum(numdg(base+1:base+nsize))
            if (nnzbs > maxnnzbs) maxnnzbs = nnzbs
            norbs = norbs + nsize
            if (norbs == no_u) EXIT
            nblocks = nblocks + 1
         enddo
         !print *, "Maznnzbs = ", maxnnzbs

         call re_alloc( buffer,  1, maxnnzbs, 'buffer',  'write_mat' )
         call re_alloc( ibuffer, 1, maxnnzbs, 'ibuffer', 'write_mat' )
      endif

!     Get listh

      call mpi_barrier(SIESTA_comm, mpierror)


    n_g = 0
    n_l = 0
    node = -1
    DO

       node = node + 1
       if (node == nprocs) node = 0

       !print *, " node: ", node, " myrank: ", myrank

       if (myrank == node) then

          base_l = listdptr(n_l+1) 

          norbs_l = min(blocksize,no_l-n_l)
          nnzs_bl = sum(numd(n_l+1:n_l+norbs_l))
          !print *, "myrank: ", myrank, " processing size:", norbs_l, nnzs_bl

          if (node==0) then
             norbs_g = min(blocksize,no_u-n_g)
             nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
             !print *, "myrank: ", myrank, " will just copy: (l,g)", &
             !                      nnzs_bl, nnzs_bg
             ibuffer(1:nnzs_bg) = listd(base_l+1:base_l+nnzs_bl)
          else
             !print *, "myrank: ", myrank, " will send to 0: ", norbs_l
             !print *, "myrank: ", myrank, " will start at: ", base_l+1
             call MPI_Send(listd(base_l+1:),nnzs_bl,MPI_integer, &
                  0,1,SIESTA_Comm,MPIerror)
             !print *, "myrank: ", myrank, " completed send ", nnzs_bl
          endif
          n_l = n_l + norbs_l

       else if (myrank == 0) then
          norbs_g = min(blocksize,no_u-n_g)
          nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
          !print *, "root will receive from ", node, " size: ", nnzs_bg

          call MPI_Recv(ibuffer,nnzs_bg,MPI_integer, &
                node,1,SIESTA_Comm,stat,MPIerror)
       endif

       if (myrank == 0) then
          !           if (old_style) then
	  if (bck_compat) then
             ptr = 0
             do i = 1, norbs_g
                nnzsi = numdg(n_g+i)
                write(lun) (ibuffer(j),j=ptr+1,ptr+nnzsi)
                ptr = ptr + nnzsi
             enddo
          else	
             write(lun) (ibuffer(j),j=1,nnzs_bg)
          endif

             n_g = n_g + norbs_g
             !print *, "root has received so far: ", n_g
                
          if (n_g == no_u) EXIT
       else
          if (n_l == no_l) then
             !print *, "rank ", myrank, " exiting loop"
             EXIT
          endif
       endif
          
    enddo

!     Get values

 do is = 1, nspin

     call mpi_barrier(SIESTA_comm, mpierror)


    n_g = 0
    n_l = 0
    node = -1
    DO

       node = node + 1
       if (node == nprocs) node = 0

       !print *, " node: ", node, " myrank: ", myrank

       if (myrank == node) then

          base_l = listdptr(n_l+1) 

          norbs_l = min(blocksize,no_l-n_l)
          nnzs_bl = sum(numd(n_l+1:n_l+norbs_l))
          !print *, "myrank: ", myrank, " processing size:", norbs_l, nnzs_bl

          if (node==0) then
             norbs_g = min(blocksize,no_u-n_g)
             nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
             !print *, "myrank: ", myrank, " will just copy: (l,g)", &
             !                      nnzs_bl, nnzs_bg
             buffer(1:nnzs_bg) = mat(base_l+1:base_l+nnzs_bl,is)
          else
             !print *, "myrank: ", myrank, " will send to 0: ", norbs_l
             !print *, "myrank: ", myrank, " will start at: ", base_l+1
             call MPI_Send(mat(base_l+1:,is),nnzs_bl,MPI_double_precision, &
                  0,1,SIESTA_Comm,MPIerror)
             !print *, "myrank: ", myrank, " completed send ", nnzs_bl
          endif
          n_l = n_l + norbs_l

       else if (myrank == 0) then
          norbs_g = min(blocksize,no_u-n_g)
          nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
          !print *, "root will receive from ", node, " size: ", nnzs_bg

          call MPI_Recv(buffer,nnzs_bg,MPI_double_precision, &
                node,1,SIESTA_Comm,stat,MPIerror)
       endif

       if (myrank == 0) then
          !           if (old_style) then
          if (bck_compat) then
             ptr = 0
             do i = 1, norbs_g
                nnzsi = numdg(n_g+i)
                write(lun) (buffer(j),j=ptr+1,ptr+nnzsi)
                ptr = ptr + nnzsi
             enddo
          else
             write(lun) (buffer(j),j=1,nnzs_bg)
          endif



             n_g = n_g + norbs_g
             !print *, "root has received so far: ", n_g
                
          if (n_g == no_u) EXIT
       else
          if (n_l == no_l) then
             !print *, "rank ", myrank, " exiting loop"
             EXIT
          endif
       endif
          
    enddo

 enddo

 if (myrank.eq.0) then
    call io_close(lun)
    call de_alloc( numdg, 'numdg', 'write_mat')
    call de_alloc( buffer, 'buffer',  'write_mat' )
    call de_alloc( ibuffer, 'ibuffer', 'write_mat' )
 endif

#endif 

 call timer("WriteMat",2)

end subroutine write_mat

  end module m_matio
