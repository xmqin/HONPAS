module m_wcml_coma
  ! Implements routines relating to the (currently unfinished)
  ! CML Condensed Matter schema

  use flib_wxml, only: xmlf_t
  use flib_wxml, only: xml_NewElement, xml_AddAttribute
  use flib_wxml, only: xml_EndElement
  use flib_wstml, only: stmAddArray
  use m_wxml_error
  
  Implicit None

  Private

  integer, private, parameter ::  sp = selected_real_kind(6,30)
  integer, private, parameter ::  dp = selected_real_kind(14,100)

  public :: cmlAddBand
  public :: cmlAddKpoint
  public :: cmlStartBandList
  public :: cmlEndBandList

Contains

  subroutine cmlAddKPoint(xf, coords, weight, kptfmt, wtfmt &
       ,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)             :: xf
    real(kind=dp), dimension(3), intent(in) :: coords
    real(kind=dp), intent(in), optional     :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label
    call cmlStartKpoint(xf, coords, weight, kptfmt, wtfmt &
         ,dictRef,convention,title,id,ref,label)
    call cmlEndKpoint(xf)
  end subroutine cmlAddKPoint

  subroutine cmlStartKPoint(xf, coords, weight, kptfmt, wtfmt &
       ,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)              :: xf
    real(kind=dp), dimension(3), intent(in)  :: coords
    real(kind=dp), intent(in), optional      :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label

    character(len=60) :: aux_str
    
    call xml_NewElement(xf, "kpoint")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(ref)) call xml_addAttribute(xf, "ref", ref)
    if (present(label)) call xml_addAttribute(xf, "label", label)

    if (present(kptfmt)) then
       write(aux_str,fmt=kptfmt) coords
    else
       write(aux_str,fmt="(3f12.8)") coords
    endif
    call xml_AddAttribute(xf, "coords", trim(aux_str))

       
    if (present(weight)) then
       if (present(wtfmt)) then
          write(aux_str,fmt=wtfmt) weight
       else
          write(aux_str,fmt="(f12.8)") weight
       endif
       
       call xml_AddAttribute(xf, "weight", trim(aux_str))
    endif
    
  end subroutine cmlStartKPoint

  subroutine cmlEndKpoint(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf, "kpoint")
  end subroutine cmlEndKpoint
  
  subroutine cmlAddBand(xf, kpoint, kweight, bands, kptfmt, eigfmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(dp), intent(in) :: kpoint(3)
    real(dp), intent(in) :: kweight
    real(dp), intent(in) :: bands(:)
    character(len=*), intent(in), optional :: kptfmt
    character(len=*), intent(in), optional :: eigfmt
    
    integer :: i, n
    character(len=20) :: k_fmt
    character(len=20) :: kp_c(3)

    if (present(kptfmt)) then
      if (len_trim(kptfmt) > 20) stop 'k-point format too long' 
      k_fmt=kptfmt
    else
      k_fmt='(f10.7)'
    endif

    n = size(bands)

    do i = 1, 3
      write(kp_c(i), k_fmt, err=100) kpoint(i)
    enddo

    call xml_NewElement(xf, 'band')
    call xml_AddAttribute(xf, 'kpoint', trim(kp_c(1))//' '//trim(kp_c(2))//' '//trim(kp_c(3)))
    call xml_AddAttribute(xf, 'weight', kweight)
    call stmAddArray(xf, array=bands, nvalue=n, fmt=eigfmt)
    call xml_EndElement(xf, 'band')

    return
    
100 call wxml_fatal('internal write failed in cmlAddBand')
  end subroutine cmlAddBand
    

  subroutine cmlStartBandList(xf, id, title, conv, dictref, ref, role)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: role
    
    call xml_NewElement(xf, 'bandList')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    if (present(role)) call xml_AddAttribute(xf, 'role', role)
    
  end subroutine cmlStartBandList


  subroutine cmlEndBandList(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'bandList')
    
  end subroutine cmlEndBandList


end module m_wcml_coma
