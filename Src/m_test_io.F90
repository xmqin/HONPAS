module m_test_io

  implicit none

contains

  ! Module to test various IO-routines

  subroutine time_io(nspin,nsc,D_2D)

    use precision, only : dp
    use parallel, only : Node
    use class_Sparsity
    use class_dSpData2D
    use m_matio
    use m_iodm

    integer, intent(in) :: nspin
    integer, intent(in) :: nsc(3)
    type(dSpData2D), intent(inout) :: D_2D

    type(Sparsity), pointer :: sp
    integer :: n_nzs, no_l, i
    integer, parameter :: N_WRITES = 50
    real(dp), pointer :: a_2D(:,:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    sp => spar(D_2D)
    call attach(sp, nrows=no_l, nnzs=n_nzs, &
         n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    a_2d => val(D_2D)

    if ( Node == 0 ) then
       write(*,*)' Starting the mat-io non-compatible'
    end if
    
    call timer('MATIO-F',1)
    do i = 1 , N_WRITES
       call write_mat(n_nzs,no_l,nspin,ncol,l_ptr,l_col,a_2d, &
            userfile='TESTIO.DM',compatible=.false.)
    end do
    call timer('MATIO-F',2)

    if ( Node == 0 ) then
       write(*,*)' Starting the mat-io compatible'
    end if

    call timer('MATIO-T',1)
    do i = 1 , N_WRITES
       call write_mat(n_nzs,no_l,nspin,ncol,l_ptr,l_col,a_2d, &
            userfile='TESTIO.DM',compatible=.true.)
    end do
    call timer('MATIO-T',2)

    if ( Node == 0 ) then
       write(*,*)' Starting the Nick-io'
    end if

    call timer('NICK-io',1)
    do i = 1 , N_WRITES
       call write_dm('TESTIO.DM',nsc,D_2D)
    end do
    call timer('NICK-io',2)

    call timer('MATIO-F',3)
    call timer('MATIO-T',3)
    call timer('NICK-io',3)

    call die('Stopping...')

  end subroutine time_io

end module m_test_io
