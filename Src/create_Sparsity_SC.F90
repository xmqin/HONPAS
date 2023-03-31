module create_Sparsity_SC

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)
  integer, parameter :: TM_ALL = -999999

  private

  public :: TM_ALL
  public :: crtSparsity_SC

contains

  ! Subroutine will create a new pattern which matches that of the
  ! 'in' sparse pattern. Furthermore, the pattern will only contain
  ! the UC reduced matrix. 
  ! That is all supercells are "downfolded" to this region.
  ! The 'out' sparsity pattern should not be allocated in any way.
  ! It will be allocated upon exit...
  ! The user, however, has to setup the matrix values afterwards by calling
  ! the appropriate routine
  ! Then routine name reflects its use
  ! crt for "create"
  ! Sparsity for "Sparsity" and
  ! SC for operation on sparsity patterns which has to do with super-cell sparsity patterns
  subroutine crtSparsity_SC(in,out, &
       DUMMY, &
       MASK, & ! 1. option
       UC, &   ! 2. option
       TM,ucell,lasto,xa,xij,isc_off) ! 3. option
    use class_Sparsity

    ! The matrix which contains the supercell.
    type(Sparsity), intent(in out) :: in
    ! Maybe (in) is not needed, however, for this, it does't make a
    ! difference, everything is overwritten
    type(Sparsity), intent(in out) :: out
    
    ! this IS a DUMMY argument, supplying it will cause the program to DIE!
    logical, intent(in), optional :: DUMMY 

    ! A mask of same size af list_col(in)
    ! Every true value contributes to the sparse matrix
    logical, intent(in out), optional :: MASK(:) ! IT HAS TO HAVE THE SAME SIZE
    ! AS nnzs(in)

    ! Determines whether we wish to refer the sparsity
    ! to the unit-cell columns...
    logical, intent(in), optional :: UC

    ! TM stands for Transfer Matrix
    ! It enables to retrieve the coupling between the UC and the requested
    ! it tells which TM we retieve (x,y,z) (can be negative as well)
    integer, intent(in), optional :: TM(3)

    real(dp), intent(in), optional :: ucell(3,3)

    ! The last orbital of each atom (needed, for determining the 
    ! orbital on a supercell orbital)
    integer, intent(in), optional :: lasto(:) ! lasto(0:na_u)
    ! The atomic coordinates in the unit-cell.
    real(dp),intent(in), optional :: xa(:,:) ! (xa(3,na_u)

    ! Supply the xij matrix
    ! *NOTE* This is necessary if one requests a TM matrix
    real(dp), intent(in), optional :: xij(:,:)

    ! Supply the super-cell offset array
    integer, intent(in), optional :: isc_off(:,:)

    ! We need space to create the new sparsity pattern:
    integer :: n_rows, n_rows_g, n_nzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer, allocatable :: entries(:)
    ! This will accomodate a scheme up to 2 decimals in the integer region,
    ! comma seperated:    
    character(len=11) :: cTM
    integer :: ir,i,iu

    ! If DUMMY is present, kill...
    if (present(DUMMY)) call die('Create sparsity: name arguments!')

    if ( present(MASK) ) then
       if ( ubound(MASK,dim=1) /= nnzs(in) ) then
          call die('Could not recognize the MASK attributed. &
               &Please supply a MASK of same size as in-sparsity &
               &pattern.')
       end if
    end if

    if ( present(xij) ) then
       if ( ubound(xij,dim=2) /= nnzs(in) ) then
          call die('Could not recognize the xij attributed. &
               &Please supply a xij of same size as in-sparsity &
               &pattern.')
       end if
    end if

    ! Save the rows ( this is the same for all cases)
    ! Even for TM which typically have 0 entries
    ! in some rows. However, it provides the full information.
    n_rows   = nrows  (in)
    n_rows_g = nrows_g(in)

    ! Prepare creation of num and listptr arrays
    allocate(num(n_rows))
    allocate(listptr(n_rows))

    ! The list pointer for the first entry is always the "first"
    ! element, hence we can already initialize it here...
    listptr(1) = 0

    ! We initialize the sparsity creation via the first row
    do ir = 1 , n_rows
       call sparsity_row_entries(in,ir,n=num(ir), &
            MASK=MASK, &
            UC=UC,TM=TM,ucell=ucell,lasto=lasto,xa=xa,xij=xij, &
            isc_off=isc_off)
       ! Update list pointer
       if ( ir > 1 ) &
            listptr(ir) = listptr(ir-1) + num(ir-1)
    end do

    ! The number of non-zero elements in the requested array
    n_nzs = listptr(n_rows) + num(n_rows)

    ! Now we will create the list array...
    allocate(list(n_nzs))

    ! allocate a temporary array which holds
    ! the maximum number of entries per row...
    iu = num(1)
    do ir = 2, n_rows
       if ( num(ir) > iu ) iu = num(ir)
    end do       
    allocate(entries(iu))

    do ir = 1, n_rows
       if ( num(ir) > 0 ) then
          call sparsity_row_entries(in,ir,n=num(ir),entries=entries, &
               MASK=MASK, &
               UC=UC,TM=TM,ucell=ucell,lasto=lasto,xa=xa,xij=xij, &
               isc_off=isc_off)
          ! Step to the element corresponding to the pointer
          iu = listptr(ir)
          ! Update all the entries
          do i = 1 , num(ir)
             list(iu+i) = entries(i)
          end do
       end if
    end do

    ! Deallocate un-needed space before we create the actual sparsity
    deallocate(entries)

    ! Down here we are sure that the options passed has
    ! been sufficient to create a sparsity pattern.
    ! Thus we simply check for existence of options
    if ( present(MASK) ) then
       call newSparsity(out,n_rows,n_rows_g,n_nzs,num,listptr,list, &
            name='(M of: '//name(in)//')', &
            ncols=ncols(in),ncols_g=ncols_g(in))
    else if ( present(TM) ) then
       do i = 1 , 3
          if ( TM(i) == TM_ALL ) then
             cTM((i-1)*3+1:i*3) = '--,'
          else
             write(cTM((i-1)*3+1:i*3),'(i2,'','')') TM(i)
          end if
       end do
       cTM(9:) = ' '
       call newSparsity(out,n_rows,n_rows_g,n_nzs,num,listptr,list, &
            name='(TM ['//trim(cTM)//'] of: '//name(in)//')', &
            ncols=ncols(in),ncols_g=nrows_g(in))
    else if ( present(UC) ) then
       if ( UC ) then
          call newSparsity(out,n_rows,n_rows_g,n_nzs,num,listptr,list, &
               name='(UC of: '//name(in)//')', &
               ncols=ncols(in),ncols_g=nrows_g(in))
       else
          ! THIS SHOULD NEVER HAPPEN (the dependency check
          ! should be performed in the initialization of this routine)
          call die('INITIALIZATION GONE WRONG IN SPARSITY_UC')
       end if
    end if

    ! We need to deallocate the arrays used, remember, that they
    ! are allocated.
    ! The newSparsity copies the values...
    deallocate(num,listptr,list)

  end subroutine crtSparsity_SC

  subroutine sparsity_row_entries(sp,row,n,entries, &
       DUMMY, &
       MASK, &
       UC,TM,ucell,lasto,xa,xij, &
       isc_off)
    use class_Sparsity

    type(Sparsity), intent(in out) :: sp
    integer, intent(in)            :: row
    integer, intent(in out), optional :: n
    integer, intent(out), optional :: entries(:)

    ! this IS a DUMMY argument, supplying it will cause the program to DIE!
    logical, intent(in), optional :: DUMMY 

    ! A MASK array which directly determines which segments
    ! should be included
    logical, intent(in), optional :: MASK(:)

    ! Determines whether we wish to refer the sparsity
    ! to the unit-cell columns...
    logical, intent(in), optional :: UC

    ! TM stands for Transfer Matrix
    ! It enables to retrieve the coupling between the UC and the requested
    ! it tells which TM we retieve (x,y,z) (can be negative as well)
    integer, intent(in), optional :: TM(3)

    real(dp), intent(in), optional :: ucell(3,3)

    ! The last orbital of each atom (needed, for determining the 
    ! orbital on a supercell orbital)
    integer, intent(in), optional :: lasto(:) ! lasto(0:na_u)
    ! The atomic coordinates in the unit-cell.
    real(dp),intent(in), optional :: xa(:,:) ! (xa(3,na_u)

    ! Supply the xij matrix
    real(dp), intent(in), optional :: xij(:,:)

    ! Supply the offset matrix
    integer, intent(in), optional :: isc_off(:,:)

    ! Check the input
    if ( present(UC) ) then
       if ( UC ) then
          call sparsity_row_entries_UC(sp,row,n=n,entries=entries)
          ! We return as the end will call a DIE
          return
       end if
    end if
    
    if ( present(TM) ) then
       ! Check consistency in call
       if ( present(isc_off) ) then
          call sparsity_row_entries_TM_off(sp,row,TM,isc_off, &
               n=n,entries=entries)
       else
          if ( .not. present(ucell) ) call die('TM creation needs ucell')
          if ( .not. present(lasto) ) call die('TM creation needs lasto')
          if ( .not. present(xa) ) call die('TM creation needs xa')
          if ( .not. present(xij) ) call die('TM creation needs xij')
          call sparsity_row_entries_TM(sp,row,TM,ucell,lasto,xa,xij, &
               n=n,entries=entries)
       end if
       return
    end if

    if ( present(MASK) ) then
       call sparsity_row_entries_MASK(sp,row,MASK, &
            n=n,entries=entries)
       return
    end if

    call die('Could not determine what you requested? &
         &Please state UC, TM or MASK!')

  end subroutine sparsity_row_entries


  ! A unit-cell row entry counter
  subroutine sparsity_row_entries_UC(sp,row, &
       n,entries)
    use class_Sparsity
    use intrinsic_missing, only : UNIQC, UNIQ
    use intrinsic_missing, only : SORT_QUICK
    use geom_helper

    type(Sparsity), intent(in out) :: sp
    integer, intent(in)            :: row
    integer, intent(in out), optional :: n
    integer, intent(out), optional :: entries(:)

    ! Local variables...
    integer, pointer :: l_col(:)
    ! We probably need this for safety reasons
    integer, allocatable :: vals(:)
    integer :: ncol, ptr, i, no_e, nr

    ! Retrieve the pointer providing the index of the columns
    ncol  =  n_col   (sp,row)
    ptr   =  list_ptr(sp,row)
    l_col => list_col(sp)
    nr    =  nrows_g (sp)

    allocate(vals(ncol))
    do i = 1 , ncol
       vals(i) = UCORB(l_col(ptr+i),nr)
    end do

    ! If the user requests the entries
    ! Then a previous call to this routine must have been
    ! performed (where entries was NOT present)
    if ( present(n) .and. present(entries) ) then
       no_e = n
    end if

    ! To retrieve the number of elements in
    ! a transfer matrix element then
    ! we do the following
    if ( .not. present(entries) ) then
       ! We wish to count the number of entries in each column
       if ( ncol > 0 ) then
          no_e = UNIQC(vals)
       else
          no_e = 0
       end if

       if ( present(n) ) then
          n = no_e
       end if

    end if
    
    if ( present(entries) ) then
       ! We wish to count the number of entries in each column
       if ( ncol > 0 ) then
          ! We need to do it in steps (gnu has troubles with non-allocated
          ! arrays in nesting constructs)...
          entries(1:no_e) = UNIQ(vals)
          call sort_quick(no_e, entries)
       end if
       ! When it returns the programmer should already know that no_e is
       ! zero
    end if

    deallocate(vals)

  end subroutine sparsity_row_entries_UC


  ! Handles sparsity row entry counts when
  ! requesting TM
  subroutine sparsity_row_entries_TM(sp,row, & 
       TM,ucell,lasto,xa,xij,&
       n,entries)
    use class_Sparsity
    use intrinsic_missing, only : SORT_QUICK
    use geom_helper

    type(Sparsity), intent(in out) :: sp
    integer, intent(in)            :: row

    ! TM stands for Transfer Matrix
    ! It enables to retrieve the coupling between the UC and the requested
    ! it tells which TM we retieve (x,y,z) (can be negative as well)
    integer, intent(in) :: TM(3)

    real(dp), intent(in) :: ucell(3,3)

    ! The last orbital of each atom (needed, for determining the 
    ! orbital on a supercell orbital)
    integer, intent(in) :: lasto(:) ! lasto(0:na_u)
    ! The atomic coordinates in the unit-cell.
    real(dp),intent(in) :: xa(:,:) ! (xa(3,na_u)

    ! Supply the xij matrix
    ! *NOTE* This is necessary if one requests a TM matrix
    real(dp), intent(in) :: xij(:,:)

    integer, intent(in out), optional :: n
    integer, intent(out), optional :: entries(:)

    ! Local variables...
    real(dp) :: recell(3,3)
    integer, pointer :: l_col(:) => null()
    integer :: iar
    integer :: ncol, ptr, i,j, no_e, t(3)

    call attach(sp,nrows=i,nrows_g=j)
    if ( i /= j ) then
       call die('Creating TM sparsity pattern requires correct &
            &atomic placement. You must not use a distributed &
            &sparsity pattern.')
    end if

    ! Retrieve the pointer providing the index of the columns
    ncol  =  n_col   (sp,row)
    ptr   =  list_ptr(sp,row)
    l_col => list_col(sp)
    iar   =  iaorb(row ,lasto)

    call reclat(ucell,recell,0) ! without 2Pi
    
    ! If the user requests the entries
    ! Then a previous call to this routine must have been
    ! performed (where entries was NOT present)
    if ( present(n) .and. present(entries) ) then
       no_e = n
    end if

    ! To retrieve the number of elements in
    ! a transfer matrix element then
    ! we do the following
    if ( .not. present(entries) ) then
       no_e = 0
       ! Retrieve the data pointer value.
       ! By doing this, we explicitly assume that the
       ! sparsity pattern of xij and *in* are the same! 
       ! TODO, check this in the beginning of the routine!
       do i = ptr+1 , ptr+ncol
          t = cell_abc(recell, &
               xa(:,iar)                  , & ! xa_i
               xa(:,iaorb(l_col(i),lasto)), & ! xa_j
               xij(:,i))
          if ( (TM(1) == TM_ALL .or. TM(1) == t(1)) .and. &
               (TM(2) == TM_ALL .or. TM(2) == t(2)) .and. &
               (TM(3) == TM_ALL .or. TM(3) == t(3)) ) then
             no_e = no_e + 1
          end if
       end do
    end if
    
    if ( present(n) .and. .not. present(entries) ) then
       n = no_e
    end if

    if ( present(entries) ) then
       j = 0
       do i = ptr+1 , ptr+ncol
          t = cell_abc(recell, &
               xa(:,iar)                  , & ! xa_i
               xa(:,iaorb(l_col(i),lasto)), & ! xa_j
               xij(:,i))
          if ( (TM(1) == TM_ALL .or. TM(1) == t(1)) .and. &
               (TM(2) == TM_ALL .or. TM(2) == t(2)) .and. &
               (TM(3) == TM_ALL .or. TM(3) == t(3)) ) then

             j = j + 1
             entries(j) = l_col(i)
          end if
       end do
       ! lets sort the entries !
       call sort_quick(j, entries)

       ! We need to check that the entries in fact does
       ! match the requested number of entries.
       if ( j /= no_e ) then
          call die('TM count of the entries does not match')
       end if
    end if
    
  end subroutine sparsity_row_entries_TM

  ! Handles sparsity row entry counts when
  ! requesting TM
  subroutine sparsity_row_entries_TM_off(sp,row, & 
       TM,isc_off,&
       n,entries)
    use class_Sparsity
    use intrinsic_missing, only : SORT_QUICK
    use geom_helper

    type(Sparsity), intent(in out) :: sp
    integer, intent(in)            :: row

    ! TM stands for Transfer Matrix
    ! It enables to retrieve the coupling between the UC and the requested
    ! it tells which TM we retieve (x,y,z) (can be negative as well)
    integer, intent(in) :: TM(3)

    ! Supply the sc_off 
    integer, intent(in) :: isc_off(:,:)

    integer, intent(in out), optional :: n
    integer, intent(out), optional :: entries(:)

    ! Local variables...
    integer, pointer :: l_col(:) => null()
    integer :: is
    integer :: ncol, ptr, i,j, no_e, nr, t(3)

    call attach(sp,nrows=i,nrows_g=nr)
    if ( i /= nr ) then
       call die('Creating TM sparsity pattern requires correct &
            &atomic placement. You must not use a distributed &
            &sparsity pattern.')
    end if

    ! Retrieve the pointer providing the index of the columns
    ncol  =  n_col   (sp,row)
    ptr   =  list_ptr(sp,row)
    l_col => list_col(sp)

    ! If the user requests the entries
    ! Then a previous call to this routine must have been
    ! performed (where entries was NOT present)
    if ( present(n) .and. present(entries) ) then
       no_e = n
    end if

    ! To retrieve the number of elements in
    ! a transfer matrix element then
    ! we do the following
    if ( .not. present(entries) ) then
       no_e = 0
       ! Retrieve the data pointer value.
       do i = ptr+1 , ptr+ncol
          is = (l_col(i)-1)/nr + 1
          t = isc_off(:,is)
          if ( (TM(1) == TM_ALL .or. TM(1) == t(1)) .and. &
               (TM(2) == TM_ALL .or. TM(2) == t(2)) .and. &
               (TM(3) == TM_ALL .or. TM(3) == t(3)) ) then
             no_e = no_e + 1
          end if
       end do
    end if
    
    if ( present(n) .and. .not. present(entries) ) then
       n = no_e
    end if

    if ( present(entries) ) then
       j = 0
       do i = ptr+1 , ptr+ncol
          is = (l_col(i)-1)/nr + 1
          t = isc_off(:,is)
          if ( (TM(1) == TM_ALL .or. TM(1) == t(1)) .and. &
               (TM(2) == TM_ALL .or. TM(2) == t(2)) .and. &
               (TM(3) == TM_ALL .or. TM(3) == t(3)) ) then

             j = j + 1
             entries(j) = l_col(i)
          end if
       end do
       ! lets sort the entries !
       call sort_quick(j, entries)

       ! We need to check that the entries in fact does
       ! match the requested number of entries.
       if ( j /= no_e ) then
          call die('TM count of the entries does not match')
       end if
    end if
    
  end subroutine sparsity_row_entries_TM_off


  ! This subroutine handles the desimation of the sparsity
  ! pattern in regards of a MASK
  subroutine sparsity_row_entries_MASK(sp,row, &
       MASK, &
       n,entries)
    use class_Sparsity
    use intrinsic_missing, only : UNIQC, UNIQ
    use intrinsic_missing, only : SORT_QUICK
    use geom_helper

    type(Sparsity), intent(in out) :: sp
    integer, intent(in)            :: row
    ! A MASK array which directly determines which segments
    ! should be included
    logical, intent(in)            :: MASK(:)
    integer, intent(in out), optional :: n
    integer, intent(out), optional :: entries(:)

    integer, pointer :: l_col(:) => null()
    integer :: ncol, ptr, i, j, no_e

    ! Retrieve the pointer providing the index of the columns
    ncol  =  n_col   (sp,row)
    ptr   =  list_ptr(sp,row)

    ! If the user requests the entries
    ! Then a previous call to this routine must have been
    ! performed (where entries was NOT present)
    if ( present(n) .and. present(entries) ) then
       no_e = n
    end if

    ! To retrieve the number of elements in
    ! a transfer matrix element then
    ! we do the following
    if ( .not. present(entries) ) then
       no_e = 0
       ! Retrieve the data pointer value.
       ! By doing this, we explicitly assume that the
       ! sparsity pattern of xij and *in* are the same! 
       ! TODO, check this in the beginning of the routine!
       do i = 1 , ncol
          if ( MASK(ptr+i) ) then ! MASK-check
             no_e = no_e + 1
          end if
       end do

    end if
    
    ! We need to return the number of elements...
    if ( present(n) .and. .not. present(entries) ) then
       n = no_e
    end if

    if ( present(entries) ) then
       ! We only need the column pointer here.
       l_col => list_col(sp)
       j = 0
       do i = 1 , ncol
          if ( MASK(ptr+i) ) then
             j = j + 1
             entries(j) = l_col(ptr+i)
          end if
       end do
       ! lets sort the entries !
       call sort_quick(j, entries)
       
       ! We need to check that the entries in fact does
       ! match the requested number of entries.
       if ( j /= no_e ) then
          call die('MASK provided inconsistency results...')
       end if
       
    end if

  end subroutine sparsity_row_entries_MASK

end module create_Sparsity_SC
