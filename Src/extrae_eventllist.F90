#ifdef TRACING
! linked list for extrae events, holding the event number and the section name
! functionalities: 
!   - add entry to list  -> addToList(list, name)
!   - find ID for name   -> getNumber(list, name)
!   - write list to file -> writeList(list)
!   - delete list        -> deleteList(list)
module extrae_eventllist
  implicit none

  integer, parameter :: extrae_MAXNAMELEN = 30
  integer            :: extrae_maxEventNumber
  integer, parameter :: NOT_FOUND = -1

  type extrae_event
    integer                           :: eventnumber
    character (len=extrae_MAXNAMELEN) :: sectionname
    type(extrae_event), pointer       :: next => null()
  end type extrae_event

  type(extrae_event), pointer :: eventlist => null()

  contains
  
  ! finds the event number associated with 'name' 
  ! returns NOT_FOUND if the entry is not found
  recursive function getNumber(list, name) result(getNumber_result)
!use parallel,only: node
    implicit none

    type(extrae_event), pointer, intent(in)  :: list
    character (len=*),           intent(in)  :: name  
    integer :: getNumber_result

    if (associated(list)) then
      if (trim(list%sectionname) == trim(name)) then
        getNumber_result = list%eventnumber
!if (Node == 0) write (*,*) 'extraeLIST found                ', list%eventnumber, name
      else
!if (Node == 0) write (*,*) 'extraeLIST continue looking for ', name
        getNumber_result = getNumber(list%next, name)
      end if
    else
      getNumber_result = NOT_FOUND
    end if

  end function getNumber


  ! adds a new element to the HEAD of the list, 
  ! resulting in descending numbering in the list
  ! this order reduces searching time for repeated calls
  function addToList(list, name)
!use parallel,only: node
    implicit none

    type(extrae_event), pointer, intent(inout)  :: list
    character (len=*),           intent(in)     :: name  
    integer    :: addToList
    
    type(extrae_event), pointer :: newElement

    extrae_maxEventNumber = extrae_maxEventNumber + 1
    allocate(newElement)
    newElement%sectionname = trim(name)
    newElement%eventnumber = extrae_maxEventNumber
    newElement%next        => list
    list                   => newElement
    
!if (Node == 0) write (*,*) 'extraeLIST adding               ', newElement%eventnumber, name
    addToList = newElement%eventnumber
  end function addToList


  ! writes the content of the list in pcf format to the file 'siesta_user_labels.pcf'
  ! its content has to be appended to the trace's pcf file for labeling the user-functions traced
  subroutine writeList(list)
    implicit none
    type(extrae_event), pointer, intent(in)  :: list

    character(len=*),parameter :: userfuncFile = 'siesta_user_labels.pcf'
    integer                    :: iu

    call io_assign(iu)
    open(unit=iu, file=userfuncFile, form='formatted', status='unknown')

    write (iu,*) "EVENT_TYPE"
    write (iu,*) "9   1000    PEXSI user functions"
    write (iu,*) "VALUES"
    write (iu,*) "    0   End"

    call writeElements(list, iu)

    call io_close(iu)

  end subroutine writeList


  recursive subroutine writeElements(list, iu)
    implicit none
    type(extrae_event), pointer, intent(in)  :: list
    integer,                     intent(in)  :: iu

    ! elements are in reverse order in list, so reverse them again  
    if (associated(list)) then
      call writeElements(list%next, iu)
      write(iu,'(i6, a)') list%eventnumber, '   ' // list%sectionname
    end if
   
  end subroutine writeElements

  
  recursive subroutine deleteList(list)
    implicit none
    type(extrae_event), pointer, intent(inout)  :: list

    if (associated(list)) then
      if (associated(list%next)) then
        call deleteList(list%next)
        deallocate(list)
      end if
    end if

  end subroutine deleteList

end module extrae_eventllist
#endif /* TRACING */
