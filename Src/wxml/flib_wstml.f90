MODULE flib_wstml

  use flib_wxml, only: xmlf_t
  use flib_wxml, only: str, xml_AddArray
  use flib_wxml, only: xml_NewElement, xml_AddPcData, xml_AddAttribute
  use flib_wxml, only: xml_EndElement, xml_AddComment
  use m_wxml_overloads

  implicit none

  private


  integer, private, parameter ::  sp = selected_real_kind(6,30)
  integer, private, parameter ::  dp = selected_real_kind(14,100)


  PUBLIC :: stmAddScalar
  PUBLIC :: stmAddArray
  PUBLIC :: stmAddMatrix
  PUBLIC :: stmAddTriangle
  PUBLIC :: stmAddStartTag

  INTERFACE stmAddScalar
     MODULE PROCEDURE stmAddString
     MODULE PROCEDURE stmAddInteger
     MODULE PROCEDURE stmAddFloatSP
     MODULE PROCEDURE stmAddFloatDP
     MODULE PROCEDURE stmAddLogical
  END INTERFACE

  INTERFACE stmAddArray
     MODULE PROCEDURE stmAddStringArray
     MODULE PROCEDURE stmAddIntegerArray
     MODULE PROCEDURE stmAddFloatArraySP
     MODULE PROCEDURE stmAddFloatArrayDP
     MODULE PROCEDURE stmAddLogicalArray
  END INTERFACE

  INTERFACE stmAddMatrix
     MODULE PROCEDURE stmAddStringMatrix
     MODULE PROCEDURE stmAddIntegerMatrix
     MODULE PROCEDURE stmAddFloatMatrixSP
     MODULE PROCEDURE stmAddFloatMatrixDP
  END INTERFACE

  INTERFACE stmAddTriangle
     MODULE PROCEDURE stmAddTriangleSP
     MODULE PROCEDURE stmAddTriangleDP
  END INTERFACE

 
CONTAINS
  
  
  ! =================================================
  ! STMML convenience routines
  ! =================================================
  
  ! -------------------------------------------------
  ! create STMML start tag in xml channel
  ! -------------------------------------------------
  
  SUBROUTINE stmAddStartTag(xf, name, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

    implicit none
    type(xmlf_t) :: xf
    character(len=*), intent(in) :: name                   ! the element name
    character(len=*), intent(in), optional :: id           ! the element id; if whitespace, is omitted
    character(len=*), intent(in), optional :: title        ! the title; if whitespace, is omitted
    character(len=*), intent(in), optional :: dictref      ! the dictionary reference; if whitespace, is omitted
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, name)
    if (present(id))         call xml_AddAttribute(xf, 'id', id)
    if (present(title))      call xml_AddAttribute(xf, 'title', title)
    if (present(dictref))    call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(dataType))   call xml_AddAttribute(xf, 'dataType', dataType)
    if (present(convention)) call xml_AddAttribute(xf, 'convention', convention)
    if (present(errorValue)) call xml_AddAttribute(xf, 'errorValue', errorValue)
    if (present(errorBasis)) call xml_AddAttribute(xf, 'errorBasis', errorBasis)
    if (present(min))        call xml_AddAttribute(xf, 'min', min)
    if (present(max))        call xml_AddAttribute(xf, 'max', max)
    if (present(units))      call xml_AddAttribute(xf, 'units', units)

  END SUBROUTINE stmAddStartTag


  ! -------------------------------------------------
  ! outputs STMML scalar in xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddString(xf, value, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    character(len=*), intent(in)           :: value         ! the value to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'scalar')
    if (present(id))         call xml_AddAttribute(xf, 'id', id)
    if (present(title))      call xml_AddAttribute(xf, 'title', title)
    if (present(dictref))    call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(dataType))   call xml_AddAttribute(xf, 'dataType', dataType)
    if (present(convention)) call xml_AddAttribute(xf, 'convention', convention)
    if (present(errorValue)) call xml_AddAttribute(xf, 'errorValue', errorValue)
    if (present(errorBasis)) call xml_AddAttribute(xf, 'errorBasis', errorBasis)
    if (present(min))        call xml_AddAttribute(xf, 'min', min)
    if (present(max))        call xml_AddAttribute(xf, 'max', max)
    if (present(units))      call xml_AddAttribute(xf, 'units', units)
    call xml_AddPcdata(xf, value, line_feed=.false.)
    call xml_EndElement(xf, 'scalar')

  END SUBROUTINE stmAddString

  ! -------------------------------------------------
  ! outputs STMML logical in xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddLogical(xf, value, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    logical,          intent(in)           :: value        ! the value to be output
    character(len=*), intent(in), optional :: id           ! the id
    character(len=*), intent(in), optional :: title        ! the title
    character(len=*), intent(in), optional :: dictref      ! the dictionary reference
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: units        ! units (default = none)

    call stmAddString(xf, str(value), id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

  END SUBROUTINE stmAddLogical


  ! -------------------------------------------------
  ! outputs STMML integer in xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddInteger(xf, value, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    integer,          intent(in)           :: value        ! the value to be output
    character(len=*), intent(in), optional :: id           ! the id
    character(len=*), intent(in), optional :: title        ! the title
    character(len=*), intent(in), optional :: dictref      ! the dictionary reference
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: units        ! units (default = none)

    call stmAddString(xf, str(value), id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units)

  END SUBROUTINE stmAddInteger


  ! -------------------------------------------------
  ! 1. create an STMML <scalar> DP float in xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatDP(xf, value, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units, fmt)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    real(kind=dp),     intent(in)           :: value        ! the value to be output
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! the title
    character(len=*), intent(in), optional :: dictref      ! the dictionary reference
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: units        ! units
    character(len=*), intent(in), optional :: fmt          ! the format 


       call stmAddString(xf, str(value,fmt), id, title, dictref, dataType, &
         convention, errorValue, errorBasis, min, max, units)

  END SUBROUTINE stmAddFloatDP

  ! -------------------------------------------------
  ! 2. create an STMML <scalar> SP float in xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatSP(xf, value, id, title, dictref, dataType, &
       convention, errorValue, errorBasis, min, max, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    real(kind=sp), intent(in)               :: value        ! the value to be output
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! the title
    character(len=*), intent(in), optional :: dictref      ! the dictionary reference
    character(len=*), intent(in), optional :: units        ! units (' ' = none)
    character(len=*), intent(in), optional :: dataType  
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: errorValue
    character(len=*), intent(in), optional :: errorBasis
    character(len=*), intent(in), optional :: min
    character(len=*), intent(in), optional :: max
    character(len=*), intent(in), optional :: fmt
    
       call stmAddString(xf, str(value,fmt), id, title, dictref, dataType, &
         convention, errorValue, errorBasis, min, max, units)

  END SUBROUTINE stmAddFloatSP


  ! -------------------------------------------------
  ! outputs string array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddStringArray(xf, nvalue, array, id, title, dictref, type, delim, ref)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue        ! number of values to be output
    character(len=*), intent(in)           :: array(*)      ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: type          ! the dataType
    character(len=*), intent(in), optional :: delim         ! delimiter
    character(len=*), intent(in), optional :: ref           ! delimiter

    ! splits data into lines whenever it overflows workspace/linelength
    character(len=1) :: delim1
    integer          :: i


    if (present(delim)) then
       delim1 = delim
    else
       delim1 = ' '
    endif

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(type))    call xml_AddAttribute(xf, 'type', type)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    ! We must not trim the value of the delimiter - it might be a single space
    call xml_AddAttribute(xf, 'delimiter', delim1)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_AddAttribute(xf, 'dataType', 'xsd:string')

    call xml_AddPcdata(xf, array(1))
    do i = 2, nvalue
      call xml_AddPcdata(xf, delim1//array(i))
    enddo
    call xml_EndElement(xf, 'array')

  END SUBROUTINE stmAddStringArray


  ! -------------------------------------------------
  ! outputs string array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddLogicalArray(xf, nvalue, array, id, title, dictref, type, delim, ref)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    integer,          intent(in)           :: nvalue        ! number of values to be output
    logical,          intent(in)           :: array(*)      ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: type          ! the dataType
    character(len=*), intent(in), optional :: delim         ! delimiter
    character(len=*), intent(in), optional :: ref           ! delimiter

    ! splits data into lines whenever it overflows workspace/linelength
    character(len=1) :: delim1
    integer          :: i


    if (present(delim)) then
       delim1 = delim
    else
       delim1 = ' '
    endif

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(type))    call xml_AddAttribute(xf, 'type', type)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call xml_AddAttribute(xf, 'delimiter', delim1)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_AddAttribute(xf, 'dataType', 'xsd:boolean')

    call xml_AddPcdata(xf, array(1))
    do i = 2, nvalue
       if (delim1 .eq. ' ') then
          call xml_AddPcdata(xf, ' ')
          call xml_AddPcdata(xf, array(i))
       else
          call xml_AddPcdata(xf, delim1)
          call xml_AddPcdata(xf, array(i))
       endif
    enddo
    call xml_EndElement(xf, 'array')

  END SUBROUTINE stmAddLogicalArray


  ! -------------------------------------------------
  ! outputs integer array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddIntegerArray(xf, nvalue, array, id, title, dictref, ref, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    integer,          intent(in)           :: nvalue        ! the number of values to be output
    integer,          intent(in)           :: array(*)      ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')
    character(len=*), intent(in), optional :: ref           ! scienitific units (default ' ')

    ! splits data into lines wherever it overflows the workspace
    integer          :: i

    ! Flush on entry and exit

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_AddAttribute(xf, 'dataType', 'xsd:integer')

    call xml_AddArray(xf,array(1:nvalue))

    call xml_EndElement(xf, 'array')
    
  END SUBROUTINE stmAddIntegerArray


  ! -------------------------------------------------
  ! 1. outputs DP float array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatArrayDP(xf, nvalue, array, id, title, dictref, units, ref, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue        ! number of values to be output
    real(kind=dp), intent(in)              :: array(*)      ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')
    character(len=*), intent(in), optional :: ref           ! 
    character(len=*), intent(in), optional :: fmt           ! the output format

    ! splits data into lines whenever it overflows workspace/linelength
    ! Flush on entry and exit

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_AddAttribute(xf, 'dataType', 'xsd:double')

    call xml_AddArray(xf,array(1:nvalue),fmt)
    call xml_EndElement(xf, 'array')

  END SUBROUTINE stmAddFloatArrayDP

  ! -------------------------------------------------
  ! 2. outputs SP float array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatArraySP(xf, nvalue, array, id, title, dictref, units, ref, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue        ! number of values to be output
    real(kind=sp), intent(in)               :: array(*)       ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')
    character(len=*), intent(in), optional :: fmt           ! the output format
    character(len=*), intent(in), optional :: ref           ! the output format


    ! splits data into lines whenever it overflows workspace/linelength
    ! Flush on entry and exit

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_AddAttribute(xf, 'dataType', 'xsd:float')

    call xml_AddArray(xf,array(1:nvalue),fmt)
    call xml_EndElement(xf, 'array')

  END SUBROUTINE stmAddFloatArraySP


  ! -------------------------------------------------
  ! outputs integer matrix to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddStringMatrix(xf, nrows, ncols, matrix, id, title, dictref, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    integer,          intent(in)           :: nrows         ! the number of rows to be output
    integer,          intent(in)           :: ncols         ! the number of rows to be output
    character(len=*), intent(in)           :: matrix(nrows,ncols) ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')

    ! splits data into lines wherever it overflows the workspace
    integer ::  i, j

    call xml_AddComment(xf,"In matrix, row (first) index is fastest")
    call xml_NewElement(xf, 'matrix')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'columns', ncols)
    call xml_AddAttribute(xf, 'rows', nrows)
    call xml_AddAttribute(xf, 'dataType', 'xsd:string')

    do i = 1, ncols
       do j = 1, nrows
          call xml_AddPcdata(xf, matrix(j, i), space=.true.)
       enddo
    enddo
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddstringMatrix



  ! -------------------------------------------------
  ! outputs integer matrix to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddIntegerMatrix(xf, nrows, ncols, matrix, id, title, dictref, units)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nrows         ! the number of rows to be output
    integer, intent(in)                    :: ncols         ! the number of rows to be output
    integer, intent(in)                    :: matrix(nrows,ncols) ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')

    ! splits data into lines wherever it overflows the workspace
    ! Flush on entry and exit
    integer ::  i, j



    call xml_AddComment(xf,"In matrix, row (first) index is fastest")
    call xml_NewElement(xf, 'matrix')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'columns', ncols)
    call xml_AddAttribute(xf, 'rows', nrows)
    call xml_AddAttribute(xf, 'dataType', 'xsd:integer')

    do i = 1, ncols
       call xml_AddArray(xf,matrix(1:nrows,i))
    enddo
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddIntegerMatrix

  ! -------------------------------------------------

  SUBROUTINE stmAddLogicalMatrix(xf, nrows, ncols, matrix, id, title, dictref, units)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nrows         ! the number of rows to be output
    integer, intent(in)                    :: ncols         ! the number of rows to be output
    logical, intent(in)                    :: matrix(nrows,ncols) ! the values to be output
    character(len=*), intent(in), optional :: id            ! the id
    character(len=*), intent(in), optional :: title         ! the title
    character(len=*), intent(in), optional :: dictref       ! the dictionary reference
    character(len=*), intent(in), optional :: units         ! scienitific units (default ' ')
    
    ! splits data into lines wherever it overflows the workspace
    ! Flush on entry and exit
    integer ::  i, j 

    call xml_AddComment(xf,"In matrix, row (first) index is fastest")
    call xml_NewElement(xf, 'matrix')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(dictref))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'columns', ncols)
    call xml_AddAttribute(xf, 'rows', nrows)
    call xml_AddAttribute(xf, 'dataType', 'xsd:boolean')
    
    do i = 1, ncols
       do j = 1, nrows
          call xml_AddPcdata(xf, matrix(j, i), space=.true.)
       enddo
    enddo
    call xml_EndElement(xf, 'matrix')
    
  END SUBROUTINE stmAddLogicalMatrix



  ! -------------------------------------------------
  ! 1. outputs DP float matrix to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatMatrixDP(xf, ncols, nrows, matrix, id, title, dictref, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: ncols                ! the number of cols to be output
    integer, intent(in)                    :: nrows                ! the number of rows to be output
    real(kind=dp), intent(in)               :: matrix(nrows,ncols)  ! the values to be output
    character(len=*), intent(in), optional :: id                   ! the id
    character(len=*), intent(in), optional :: title                ! the title
    character(len=*), intent(in), optional :: dictref              ! the dictionary reference
    character(len=*), intent(in), optional :: units                ! scienitific units (default ' ')
    character(len=*), intent(in), optional :: fmt                  ! format

    integer ::  i, j

    ! splits data into lines wherever it overflows the workspace
    ! Flush on entry and exit      
    !-------------
    call xml_AddComment(xf,"In matrix, row (first) index is fastest")
    call xml_NewElement(xf, 'matrix')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'columns', ncols)
    call xml_AddAttribute(xf, 'rows', nrows)
    call xml_AddAttribute(xf, 'dataType', 'xsd:double')
    !-------------
    do i = 1, ncols
       call xml_AddArray(xf,matrix(1:nrows,i),fmt)
    enddo
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddFloatMatrixDP

  ! -------------------------------------------------
  ! 2. outputs SP float matrix to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddFloatMatrixSP(xf, ncols, nrows, matrix, id, title, dictref, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: ncols               ! the number of cols to be output
    integer, intent(in)                    :: nrows               ! the number of rows to be output
    real(kind=sp), intent(in)              :: matrix(nrows,ncols) ! the values to be output
    character(len=*), intent(in), optional :: id                  ! the id
    character(len=*), intent(in), optional :: title               ! the title
    character(len=*), intent(in), optional :: dictref             ! the dictionary reference
    character(len=*), intent(in), optional :: units               ! scienitific units (default ' ')
    character(len=*), intent(in), optional :: fmt                 ! format

    ! internal variable
    integer ::  i, j

    ! splits data into lines wherever it overflows the workspace
    ! Flush on entry and exit      
    !
    call xml_AddComment(xf,"In matrix, row (first) index is fastest")
    call xml_NewElement(xf, 'matrix')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'columns', ncols)
    call xml_AddAttribute(xf, 'rows', nrows)
    call xml_AddAttribute(xf, 'dataType', 'xsd:float')
    do i = 1, ncols
       call xml_AddArray(xf,matrix(1:nrows,i),fmt)
    enddo
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddFloatMatrixSP


  ! -------------------------------------------------
  ! 1. outputs DP lower triangle array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddTriangleDP(xf, nvalue, array, id, title, dictref, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue         ! number of values to be output
    real(kind=dp), intent(in)              :: array(*)       ! the values to be output
    character(len=*), intent(in), optional :: id             ! the id
    character(len=*), intent(in), optional :: title          ! the title
    character(len=*), intent(in), optional :: dictref        ! the dictionary reference
    character(len=*), intent(in), optional :: units          ! units (' ' = none)
    character(len=*), intent(in), optional :: fmt            ! the output format

    ! splits data into lines whenever it overflows workspace/linelength
    ! Flush on entry and exit
    integer           :: size
    character(len=200) :: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    size = (nvalue*(nvalue+1))/2
    call xml_NewElement(xf, 'array')
    call xml_AddAttribute(xf, 'size', size)
    call xml_AddAttribute(xf, 'rows', nvalue)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call STMARCF9DP(xf, size, array, formt)
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddTriangleDP

  ! -------------------------------------------------
  ! 2. outputs SP lower triangle array to xml channel
  ! -------------------------------------------------

  SUBROUTINE stmAddTriangleSP(xf, nvalue, array, id, title, dictref, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue         ! number of values to be output
    real(kind=sp), intent(in)              :: array(*)       ! the values to be output
    character(len=*), intent(in), optional :: id             ! the id
    character(len=*), intent(in), optional :: title          ! the title
    character(len=*), intent(in), optional :: dictref        ! the dictionary reference
    character(len=*), intent(in), optional :: units          ! units (' ' = none)
    character(len=*), intent(in), optional :: fmt            ! the output format

    ! splits data into lines whenever it overflows workspace/linelength
    ! Flush on entry and exit
    integer           :: size
    character(len=200) :: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    size = (nvalue*(nvalue+1))/2
    call xml_NewElement(xf, 'array')
    call xml_AddAttribute(xf, 'size', size)
    call xml_AddAttribute(xf, 'rows', nvalue)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call STMARCF9SP(xf, size, array, formt)
    call xml_EndElement(xf, 'matrix')

  END SUBROUTINE stmAddTriangleSP


  ! -------------------------------------------------
  ! outputs fatal error message
  ! -------------------------------------------------

  SUBROUTINE stmErrorMessage(xf, msg, id, title, dictref)

    implicit none
    type(xmlf_t) :: xf
    character(len=*), intent(in)           :: msg            ! the message
    character(len=*), intent(in), optional :: id             ! the id
    character(len=*), intent(in), optional :: title          ! the title
    character(len=*), intent(in), optional :: dictref        ! the dictionary reference

    call xml_NewElement(xf, 'message')
    call xml_AddAttribute(xf, 'severity', 'fatal')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    call xml_AddPcdata(xf, msg)
    call xml_EndElement(xf, 'message')

  END SUBROUTINE stmErrorMessage


  ! -------------------------------------------------
  ! outputs informational message
  ! -------------------------------------------------

  SUBROUTINE stmInfoMessage(xf, msg, id, title, dictref)

    implicit none
    type(xmlf_t) :: xf
    character(len=*), intent(in)           :: msg            ! the message
    character(len=*), intent(in), optional :: id             ! the id
    character(len=*), intent(in), optional :: title          ! the title
    character(len=*), intent(in), optional :: dictref        ! the dictionary reference

    call xml_NewElement(xf, 'message')
    call xml_AddAttribute(xf, 'severity', 'warning')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    call xml_AddPcdata(xf, msg)
    call xml_EndElement(xf, 'message')

  END SUBROUTINE stmInfoMessage


  ! -------------------------------------------------
  ! outputs warning message
  ! -------------------------------------------------

  SUBROUTINE stmWarningMessage(xf, msg, id, title, dictref)

    implicit none
    type(xmlf_t) :: xf
    character(len=*), intent(in)           :: msg            ! the message
    character(len=*), intent(in), optional :: id             ! the id
    character(len=*), intent(in), optional :: title          ! the title
    character(len=*), intent(in), optional :: dictref        ! the dictionary reference

    call xml_NewElement(xf, 'message')
    call xml_AddAttribute(xf, 'severity', 'info')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    call xml_AddPcdata(xf, msg)
    call xml_EndElement(xf, 'message')

  END SUBROUTINE stmWarningMessage



  ! =================================================
  ! basic STMML routines
  ! =================================================


  ! -------------------------------------------------
  ! creates STMML <scalar> string
  ! -------------------------------------------------

  SUBROUTINE STMSCAS9(xf, value, id, title, dictref, type)

    implicit none
    type(xmlf_t) :: xf
    character(len=*), intent(in) :: value                 ! the value to be output
    character(len=*), intent(in), optional :: id          ! the id
    character(len=*), intent(in), optional :: title       ! the title
    character(len=*), intent(in), optional :: dictref     ! the dictionary reference
    character(len=*), intent(in), optional :: type        ! the data type (default 'xsd:string')

    call xml_NewElement(xf, 'scalar')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(dictref)) call xml_AddAttribute(xf, 'title', title)
    if (present(type))    call xml_AddAttribute(xf, 'dataType', type)

!    if (XMLCHKS9(value)) then
    call xml_AddPcdata(xf, value)
    call xml_EndElement(xf, 'scalar')

  END SUBROUTINE STMSCAS9



  ! -------------------------------------------------
  ! output start tag for an STMML array
  ! -------------------------------------------------

  SUBROUTINE STMARST9(xf, nvalue, id, title, dictref, tuval, delim)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue      ! the number of values to be output
    character(len=*), intent(in), optional :: id          ! the id
    character(len=*), intent(in), optional :: title       ! the title
    character(len=*), intent(in), optional :: dictref     ! the dictionary reference
    character(len=*), intent(in), optional :: tuval       ! data type (default 'xsd:string') or units (' ' = none)
    character(len=*), intent(in), optional :: delim       ! the delimiter (default ' ')

    ! Internal Variables
    character(len=1) :: delim1

    if (present(delim)) then
       delim1 = delim
    else
       delim1 = ' '
    endif

    call xml_NewElement(xf, 'array')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(dictref)) call xml_AddAttribute(xf, 'title', title)
    if (present(tuval))   call xml_AddAttribute(xf, 'type', tuval)
    call xml_AddAttribute(xf, 'delimiter', delim1)
    call xml_AddAttribute(xf, 'size', nvalue)
    call xml_EndElement(xf, 'array')

  END SUBROUTINE STMARST9



  ! -------------------------------------------------
  ! 2. outputs SP float array to channel
  ! -------------------------------------------------

  SUBROUTINE STMARF9SP(xf, nvalue, arrf, id, title, dictref, units, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)                    :: nvalue      ! the number of values to be output
    real(kind=sp), intent(in)               :: arrf(*)     ! the values to be output
    character(len=*), intent(in), optional :: id          ! the id
    character(len=*), intent(in), optional :: title       ! the title
    character(len=*), intent(in), optional :: dictref     ! the dictionary reference
    character(len=*), intent(in), optional :: units       ! units (' ' = none)
    character(len=*)                       :: fmt         ! the output format

    call xml_NewElement(xf, 'scalar')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(dictref)) call xml_AddAttribute(xf, 'title', title)
    if (present(units))   call xml_AddAttribute(xf, 'units', units)
    call xml_AddAttribute(xf, 'size', nvalue)
    call STMARCF9SP(xf, nvalue, arrf, fmt)
    call xml_NewElement(xf, 'scalar')

  END SUBROUTINE STMARF9SP


  ! -------------------------------------------------
  ! 1. outputs content of DP float array to channel
  ! -------------------------------------------------

  SUBROUTINE STMARCF9DP(xf, nvalue, arrf, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer, intent(in)         :: nvalue   
    real(kind=dp), intent(in)   :: arrf(*) 
    character(len=*), optional :: fmt

    integer :: i
    character(len=30) :: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    call xml_AddPcdata(xf, arrf(1), formt)
    do i = 2, nvalue
       call xml_AddPcdata(xf, arrf(i), formt, space=.true.)
    enddo
  END SUBROUTINE STMARCF9DP


  ! -------------------------------------------------
  ! 2. outputs content of SP float array to channel
  ! -------------------------------------------------

  SUBROUTINE STMARCF9SP(xf, nvalue, arrf, fmt)

    implicit none
    type(xmlf_t) :: xf
    integer,       intent(in)   :: nvalue     ! the number of values to be output
    real(kind=sp), intent(in)   :: arrf(*)    ! the values to be output
    character(len=*), optional :: fmt

    integer :: i
    character(len=30) :: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    call xml_AddPcdata(xf, arrf(1), formt)
    do i = 2, nvalue
       call xml_AddPcdata(xf, arrf(i), formt, space=.true.)
    enddo
  END SUBROUTINE STMARCF9SP

end module flib_wstml
