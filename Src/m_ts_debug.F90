module m_ts_debug

  use precision, only :dp

  implicit none

  private :: dp

#ifdef TRANSIESTA_DEBUG
contains

  subroutine write_TriMat(iu,tri)
    use class_zTriMat
    use parallel, only : Nodes
    integer, intent(inout) :: iu
    type(zTriMat), intent(inout) :: tri
    complex(dp), pointer :: z(:), zf(:)
    integer :: i,j, p, np, n,idx, ntmp

    ! Number of parts
    np = parts(tri)
    n = 0
    zf => val(tri)
    
    do p = 1 , np - 1

       z => val(tri,p,p)
       if ( size(z) /= nrows_g(tri,p)**2 ) call die('SIZE1')
       do j = 1 , nrows_g(tri,p)
          idx = (j-1)*nrows_g(tri,p)
          do i = 1 , nrows_g(tri,p)
             call out_write(iu,n+i,n+j,z(idx+i))
             call test(z(idx+i),n+i,n+j)
          end do
       end do

       z => val(tri,p+1,p)
       if ( size(z) /= nrows_g(tri,p)*nrows_g(tri,p+1) ) &
            call die('SIZE2')
       ntmp = n + nrows_g(tri,p)
       do j = 1 , nrows_g(tri,p)
          idx = (j-1)*nrows_g(tri,p+1)
          do i = 1 , nrows_g(tri,p+1)
             call out_write(iu,ntmp+i,n+j,z(idx+i))
             call test(z(idx+i),ntmp+i,n+j)
          end do
       end do

       z => val(tri,p,p+1)
       if ( size(z) /= nrows_g(tri,p)*nrows_g(tri,p+1) ) &
            call die('SIZE3')
       ntmp = n + nrows_g(tri,p)
       do j = 1 , nrows_g(tri,p+1)
          idx = (j-1)*nrows_g(tri,p)
          do i = 1 , nrows_g(tri,p)
             call out_write(iu,n+i,ntmp+j,z(idx+i))
             call test(z(idx+i),n+i,ntmp+j)
          end do
       end do
       
       n = n + nrows_g(tri,p)
    end do
    
    z => val(tri,np,np)
    if ( size(z) /= nrows_g(tri,np)**2 ) call die('SIZE1')
    do j = 1 , nrows_g(tri,np)
       idx = (j-1)*nrows_g(tri,np)
       do i = 1 , nrows_g(tri,np)
          call test(z(idx+i),n+i,n+j)
          call out_write(iu,n+i,n+j,z(idx+i))
       end do
    end do

    iu = iu + Nodes

  contains
    
    subroutine test(z1,i,j)
      complex(dp), intent(in) :: z1
      integer, intent(in) :: i,j
      if ( cdabs(z1-zf(index(tri,i,j))) > 1e-10_dp ) then
         write(*,*) i,j,z1,zf(index(tri,i,j))
         call die('Not same element')
      end if
    end subroutine test

  end subroutine write_TriMat

  subroutine write_Full(iu,no,GF)
    use parallel, only : Nodes
    integer, intent(inout) :: iu
    integer, intent(in) :: no
    complex(dp), intent(in) :: GF(no,no)
    integer :: i,j 
    
    do j = 1 , no
       do i = 1 , no
          call out_write(iu,i,j,GF(i,j))
       end do
    end do

    iu = iu + Nodes

  end subroutine write_Full

  subroutine out_write(iu,i,j,z)
    integer, intent(in) :: iu,i,j
    complex(dp), intent(in) :: z
    write(iu,'(2(tr1,i5),2(tr1,e20.13))') i,j,real(z),aimag(z)
  end subroutine out_write
  
  subroutine sp_to_file(u,sp)
    use class_Sparsity
    use geom_helper, only : UCORB
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
#endif
    integer, intent(in) :: u
    type(Sparsity), intent(inout) :: sp
    integer :: io,jo,j,ind, no, no_u
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    character(len=20) :: f

    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no,nrows_g=no_u)

    write(f,'(a,i0)') 'SP.',u
    open(unit=5555,file=trim(f),form='formatted')

    write(5555,'(i5)') no

    do io = 1 , no
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          jo = UCORB(l_col(ind),no_u)
          write(5555,'(2(i7,tr1),i1)') io,jo,1
       end do
    end do
    if ( ind /= nnzs(sp) ) then
       call die('Have not looped through all things')
    end if

    close(5555)

#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,io)
#endif

  end subroutine sp_to_file

#endif

end module m_ts_debug
