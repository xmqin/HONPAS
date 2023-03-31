module flib_wcml

  use flib_wxml, only: xmlf_t, str, xml_OpenFile, xml_Close
  use flib_wxml, only: xml_NewElement, xml_AddPcData, xml_AddAttribute
  use flib_wxml, only: xml_EndElement, xml_AddArray
  use flib_wxml, only: xml_AddXMLDeclaration, xml_AddComment
  use flib_wstml, only: stmAddScalar
  use flib_wstml, only: stmAddMatrix
  use flib_wstml, only: stmAddArray
  use flib_wstml, only: stmAddStartTag
  
  use m_wcml_coma


  PRIVATE

  integer, private, parameter ::  sp = selected_real_kind(6,30)
  integer, private, parameter ::  dp = selected_real_kind(14,100)
  
! CMLUnits
  character(len=*), parameter :: U_ANGSTR = 'units:angstrom'
  character(len=*), parameter :: U_PMETER = 'units:pm'
  character(len=*), parameter :: U_DEGREE = 'units:degree'
  character(len=*), parameter :: U_RADIAN = 'units:radian'
  character(len=*), parameter :: U_INVCM  = 'units:cm-1'
  character(len=*), parameter :: U_KCALMO = 'units:kcal-mole'
  character(len=*), parameter :: U_EVOLT  = 'units:ev'
  character(len=*), parameter :: U_SECOND = 'units:second'
  character(len=*), parameter :: U_VOLT   = 'units:volt'


! CMLCore
  PUBLIC :: cmlAddCoordinates
  PUBLIC :: cmlAddLattice
  PUBLIC :: cmlAddCrystal
  PUBLIC :: cmlAddAngle
  PUBLIC :: cmlAddLength
  PUBLIC :: cmlAddEigenvalue
  PUBLIC :: cmlAddMolecule
  PUBLIC :: cmlAddMetadata
  PUBLIC :: cmlAddProperty
  PUBLIC :: cmlStartMetadataList
  PUBLIC :: cmlEndMetadataList
  PUBLIC :: cmlStartModule
  PUBLIC :: cmlEndModule
  PUBLIC :: cmlStartParameterList
  PUBLIC :: cmlEndParameterList
  PUBLIC :: cmlStartPropertyList
  PUBLIC :: cmlEndPropertyList

! CMLComp
  PUBLIC :: cmlAddParameter
  PUBLIC :: cmlStartStep
  PUBLIC :: cmlEndStep

! CMLComa
  PUBLIC :: cmlStartBandList
  PUBLIC :: cmlEndBandList
  PUBLIC :: cmlAddBand
  public :: cmlAddKpoint

! CMLCore
  INTERFACE cmlAddCoordinates
     MODULE PROCEDURE cmlAddCoordinatesSP
     MODULE PROCEDURE cmlAddCoordinatesDP
  END INTERFACE

  INTERFACE cmlAddLattice
     MODULE PROCEDURE cmlAddLatticeSP
     MODULE PROCEDURE cmlAddLatticeDP
  END INTERFACE

  INTERFACE cmlAddCrystal
     MODULE PROCEDURE cmlAddCrystalSP
     MODULE PROCEDURE cmlAddCrystalDP
  END INTERFACE

  INTERFACE cmlAddAngle
     MODULE PROCEDURE cmlAddAngleSP
     MODULE PROCEDURE cmlAddAngleDP
  END INTERFACE

  INTERFACE cmlAddLength
     MODULE PROCEDURE cmlAddLengthSP
     MODULE PROCEDURE cmlAddLengthDP
  END INTERFACE

  INTERFACE cmlAddEigenvalue
     MODULE PROCEDURE cmlAddEigenvalueSP
     MODULE PROCEDURE cmlAddEigenvalueDP
  END INTERFACE

  INTERFACE cmlAddMolecule
     MODULE PROCEDURE cmlAddMoleculeSP
     MODULE PROCEDURE cmlAddMoleculeDP
     MODULE PROCEDURE cmlAddMolecule3SP
     MODULE PROCEDURE cmlAddMolecule3DP
  END INTERFACE

  INTERFACE cmlAddProperty
     MODULE PROCEDURE cmlAddPropScalarDP
     MODULE PROCEDURE cmlAddPropScalarSP
     MODULE PROCEDURE cmlAddPropScalarI
     MODULE PROCEDURE cmlAddPropScalarCh
     MODULE PROCEDURE cmlAddPropScalarLG  ! only for scalar
     MODULE PROCEDURE cmlAddPropArrayDPSi
     MODULE PROCEDURE cmlAddPropArrayDPSh
     MODULE PROCEDURE cmlAddPropArraySPSi
     MODULE PROCEDURE cmlAddPropArraySPSh
     MODULE PROCEDURE cmlAddPropArrayISi
     MODULE PROCEDURE cmlAddPropArrayISh
     MODULE PROCEDURE cmlAddPropArrayChSi
     MODULE PROCEDURE cmlAddPropArrayChSh
     MODULE PROCEDURE cmlAddPropMatrixDPSi
     MODULE PROCEDURE cmlAddPropMatrixDPSh
     MODULE PROCEDURE cmlAddPropMatrixSPSi
     MODULE PROCEDURE cmlAddPropMatrixSPSh
     MODULE PROCEDURE cmlAddPropMatrixISi
     MODULE PROCEDURE cmlAddPropMatrixISh
     MODULE PROCEDURE cmlAddPropMatrixChSi
     MODULE PROCEDURE cmlAddPropMatrixChSh
  END INTERFACE

  INTERFACE cmlAddParameter
     MODULE PROCEDURE cmlAddParameterCH
     MODULE PROCEDURE cmlAddParameterI
     MODULE PROCEDURE cmlAddParameterSP
     MODULE PROCEDURE cmlAddParameterDP
     MODULE PROCEDURE cmlAddParameterLG
  END INTERFACE

  INTERFACE cmlAddMetadata
     MODULE PROCEDURE cmlAddMetaDataCh
     MODULE PROCEDURE cmlAddMetaDataI
  END INTERFACE

  public :: cmlBeginFile, cmlFinishFile
  public :: cmlStartCML, cmlEndCML
  public :: cmlNamespaceAttribute, cmlAddComment
  
CONTAINS

    subroutine cmlAddComment(xf, comment)
      type(xmlf_t), intent(inout) :: xf
      character(len=*) :: comment

      call xml_AddComment(xf,comment)
    end subroutine cmlAddComment

    subroutine cmlBeginFile(xf, filename, unit, replace)
      type(xmlf_t), intent(out) :: xf
      character(len=*), intent(in) :: filename
      integer, intent(in) :: unit
      logical, intent(in), optional :: replace

      call xml_OpenFile(filename, xf, indent=.true.)
      call xml_AddXMLDeclaration(xf,encoding="UTF-8")

    end subroutine cmlBeginFile

    subroutine cmlFinishFile(xf)
      type(xmlf_t), intent(inout) :: xf

      call xml_Close(xf)

    end subroutine cmlFinishFile
    
    subroutine cmlNamespaceAttribute(xf, prefix, URI)
      ! A fake namespace declaration
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: URI

      call xml_AddAttribute(xf,prefix,URI)

    end subroutine cmlNamespaceAttribute
    
    subroutine cmlStartCml(xf, id, title, convention, dictref, fileId, version)
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: id
      character(len=*), intent(in), optional :: title
      character(len=*), intent(in), optional :: convention
      character(len=*), intent(in), optional :: dictref
      character(len=*), intent(in), optional :: fileId
      character(len=*), intent(in), optional :: version

      call xml_NewElement(xf, 'cml')
      if (present(id)) call xml_AddAttribute(xf, 'id', id)
      if (present(title)) call xml_AddAttribute(xf, 'title', title)
      if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

      call xml_AddAttribute(xf, 'convention', 'CMLComp')

      if (present(fileId)) then
         call xml_AddAttribute(xf, 'fileId', fileId)
      endif
      if (present(version)) then
         call xml_AddAttribute(xf, 'version', version)
      endif

    end subroutine cmlStartCml
    !--------------------------------------------------------------------
    
    subroutine cmlEndCml(xf)
      type(xmlf_t), intent(inout) :: xf

      call cmlAddMetadata(xf, name='dc:contributor', content='Siesta-CML')
      call xml_EndElement(xf, 'cml')

    end subroutine cmlEndCml

  ! =================================================
  ! convenience CML routines
  ! =================================================
  ! -------------------------------------------------
  ! writes a metadataList start/end Tag to xml channel
  ! -------------------------------------------------
  SUBROUTINE cmlStartMetadataList(xf, id, title, conv, dictref, ref, role)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: role
    
    call xml_NewElement(xf, 'metadataList')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    if (present(role)) call xml_AddAttribute(xf, 'role', role)
    
  END SUBROUTINE cmlStartMetadataList

  SUBROUTINE cmlEndMetadataList(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'metadataList')
    
  END SUBROUTINE cmlEndMetadataList

  ! -------------------------------------------------
  ! writes a Module start/end Tag to xml channel
  ! -------------------------------------------------
  SUBROUTINE cmlStartModule(xf, id, title, conv, dictref, ref, role, serial)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: role
    character(len=*), intent(in), optional :: serial
    
    call xml_NewElement(xf, 'module')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    if (present(role)) call xml_AddAttribute(xf, 'role', role)
    if (present(serial)) call xml_AddAttribute(xf, 'serial', serial)
    
  END SUBROUTINE cmlStartModule

  SUBROUTINE cmlEndModule(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'module')
    
  END SUBROUTINE cmlEndModule


  ! -------------------------------------------------
  ! writes a propertyList start/end Tag to xml channel
  ! -------------------------------------------------
  SUBROUTINE cmlStartPropertyList(xf, id, title, conv, dictref, ref, role)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: role
    
    call xml_NewElement(xf, 'propertyList')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    if (present(role)) call xml_AddAttribute(xf, 'role', role)
    
  END SUBROUTINE cmlStartPropertyList

  SUBROUTINE cmlEndPropertyList(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'propertyList')
    
  END SUBROUTINE cmlEndPropertyList


  ! -------------------------------------------------
  ! writes a parameterList start/end Tag to xml channel
  ! -------------------------------------------------
  SUBROUTINE cmlStartParameterList(xf, id, title, conv, dictref, ref, role)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: role
    
    call xml_NewElement(xf, 'parameterList')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    if (present(role)) call xml_AddAttribute(xf, 'role', role)
    
  END SUBROUTINE cmlStartParameterList

  SUBROUTINE cmlEndParameterList(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'parameterList')
    
  END SUBROUTINE cmlEndParameterList

  ! -------------------------------------------------
  ! writes a step start Tag to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlStartStep(xf, type, index, id, title, conv, ref)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: type
    character(len=*), intent(in), optional :: id
    integer, intent(in), optional :: index
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref

    if (present(index)) then
      call cmlStartModule(xf, id, title, conv, type, ref, 'step', str(index))
    else
      call cmlStartModule(xf, id, title, conv, type, ref, 'step')
    endif
    
  END SUBROUTINE cmlStartStep

  SUBROUTINE cmlEndStep(xf)

    implicit none
    type(xmlf_t), intent(inout) :: xf

    Call xml_EndElement(xf, 'module')
    
  END SUBROUTINE cmlEndStep


  
 
  ! -------------------------------------------------
  ! 1. writes a DP property to xml channel
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddPropScalarDP(xf, value, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)               :: value
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))    call xml_AddAttribute(xf, 'ref', ref)
    call stmAddScalar(xf=xf, value=value, datatype='xsd:double', units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')

  END SUBROUTINE cmlAddPropScalarDP

  ! -------------------------------------------------
  ! 2. writes a Scalar SP property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropScalarSP(xf, property, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    real(kind=sp),    intent(in)           :: property
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    call stmAddScalar(xf=xf, value=property, datatype='xsd:float', units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropScalarSP
  
  ! -------------------------------------------------
  ! 3. writes a Scalar integer property to xml channel
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddPropScalarI(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in) :: value
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddScalar(xf=xf, value=value, datatype='xsd:integer', units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropScalarI

  ! -------------------------------------------------
  ! 4. writes a DP Float matrix property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropMatrixDPSi(xf, value, nrows, ncols, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: nrows
    integer, intent(in)                    :: ncols
    real(kind=dp), Intent(in)              :: value(ncols, nrows)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixDPSi

  SUBROUTINE cmlAddPropMatrixDPSh(xf, value, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: value(:,:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units
    
    integer :: nrows, ncols

    ncols=size(value, 2)
    nrows=size(value, 1)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixDPSh

  ! -------------------------------------------------
  ! 5. writes an SP Float matrix property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropMatrixSPSi(xf, property, nrows, ncols, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: nrows
    integer, intent(in)                    :: ncols
    real(kind=sp), intent(in)              :: property(ncols,nrows)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf,matrix=property, ncols=ncols, nrows=nrows, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixSPSi

  SUBROUTINE cmlAddPropMatrixSPSh(xf, property, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: property(:,:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    integer :: nrows, ncols

    ncols=size(property, 2)
    nrows=size(property, 1)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf,matrix=property, ncols=ncols, nrows=nrows, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixSPSh

  ! -------------------------------------------------
  ! 6. writes an Integer matrix property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropMatrixISi(xf, value, nrows, ncols, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: nrows
    integer, intent(in)                    :: ncols
    integer, intent(in)                    :: value(nrows,ncols)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(conv))    call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixISi

  SUBROUTINE cmlAddPropMatrixISh(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: value(:,:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    integer :: nrows, ncols

    ncols=size(value, 2)
    nrows=size(value, 1)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(conv))    call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixISh


  ! -------------------------------------------------
  ! 7. writes an Array DP property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropArrayDPSi(xf, value, nvalue, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: value(*)
    integer, intent(in)                    :: nvalue
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf=xf, array=value, nvalue=nvalue, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayDPSi

  SUBROUTINE cmlAddPropArrayDPSh(xf, value, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: value(:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    integer :: nvalue
    
    nvalue=size(value)

    call xml_NewElement(xf, 'property')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref)) call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf=xf, array=value, nvalue=nvalue, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayDPSh

  ! -------------------------------------------------
  ! 8. writes an Array SP property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropArraySPSi(xf, property, nvalue, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: property(*)
    integer, intent(in)                    :: nvalue
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf=xf, array=property, nvalue=nvalue, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArraySPSi

  SUBROUTINE cmlAddPropArraySPSh(xf, property, id, title, conv, dictref, ref, units, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)              :: property(:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: fmt
    character(len=*), intent(in), optional :: units

    integer :: nvalue
    
    nvalue=size(property)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf=xf, array=property, nvalue=nvalue, units=units, fmt=fmt)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArraySPSh

  ! -------------------------------------------------
  ! 9. writes an Array integer property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropArrayISi(xf, value, nvalue, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: value(*)
    integer, intent(in)                    :: nvalue
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf, array=value, nvalue=nvalue, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayISi

  SUBROUTINE cmlAddPropArrayISh(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: value(:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    integer :: nvalue

    nvalue=size(value)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf, array=value, nvalue=nvalue, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayISh


  ! -------------------------------------------------
  ! 10. writes a character property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropScalarCH(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in)           :: value
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddScalar(xf=xf, value=value, units=units)
    call xml_EndElement(xf, 'property')

  END SUBROUTINE cmlAddPropScalarCH


  ! -------------------------------------------------
  ! 11. writes an character matrix property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropMatrixCHSi(xf, value, nrows, ncols, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    integer,          intent(in)           :: nrows
    integer,          intent(in)           :: ncols
    character(len=*), intent(in)           :: value(ncols,nrows)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixCHSi


  SUBROUTINE cmlAddPropMatrixCHSh(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    character(len=*), intent(in)           :: value(:,:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    integer :: nrows, ncols

    ncols=size(value, 2)
    nrows=size(value, 1)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddMatrix(xf=xf, matrix=value, ncols=ncols, nrows=nrows, units=units)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropMatrixCHSh

  ! -------------------------------------------------
  ! 12. writes an character array value to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropArrayCHSi(xf, value, nvalue, id, title, conv, dictref, ref)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    character(len=*), intent(in)           :: value(*)
    integer,          intent(in)           :: nvalue
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf, array=value, nvalue=nvalue)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayCHSi

  SUBROUTINE cmlAddPropArrayCHSh(xf, value, id, title, conv, dictref, ref)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    character(len=*), intent(in)           :: value(:)
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref

    integer :: nvalue

    nvalue=size(value)

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddArray(xf, array=value, nvalue=nvalue)
    call xml_EndElement(xf, 'property')
  END SUBROUTINE cmlAddPropArrayCHSh


  ! -------------------------------------------------
  ! 13. writes a logical property to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddPropScalarLG(xf, value, id, title, conv, dictref, ref, units)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    logical,          intent(in)           :: value
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: units

    call xml_NewElement(xf, 'property')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    call stmAddScalar(xf=xf, value=value, units=units)
    call xml_EndElement(xf, 'property')

  END SUBROUTINE cmlAddPropScalarLG

  !------------------------------------------------------------
  ! END OF PROPERTIES 
  !------------------------------------------------------------



  ! -------------------------------------------------
  ! 1. writes complete DP molecule to xml channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddMoleculeDP(xf, natoms, elements, atomRefs, coords, style, id, title, dictref, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in)                    :: natoms             ! number of atoms
    real(kind=dp), intent(in)              :: coords(3, natoms)  ! atomic coordinates
    character(len=*), intent(in)           :: elements(natoms)   ! chemical element types
    character(len=*), intent(in), optional :: atomRefs(natoms)       ! id
    character(len=*), intent(in), optional :: id                 ! id
    character(len=*), intent(in), optional :: title              ! the title
    character(len=*), intent(in), optional :: dictref            ! the dictionary reference
    character(len=*), intent(in), optional :: fmt                ! format for coords
    character(len=*), intent(in), optional :: style              ! type of coordinates 

    ! 'x3' for Cartesians, 
    ! 'xFrac' for fractionals
    ! default => cartesians

    ! Internal Variables
    character(len=6) :: id1, id0
    character(len=20):: formt, stylei
    integer          :: i

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f14.8)'
    endif
    if (present(style)) then
       stylei = style
    else
       stylei = 'x3'
    endif

    call stmAddStartTag(xf, 'molecule', id, title, dictref)
    call xml_NewElement(xf, 'atomArray')
    do i = 1, natoms
       write(id0, '(i4)') i
       id0 = adjustl(id0)
       id1 = 'a'
       id1(2:) = id0
       call cmlAddAtom(xf=xf, elem=elements(i), id=trim(id1))
       if (present(atomRefs)) call xml_AddAttribute(xf, 'ref', trim(atomRefs(i)))
       if (stylei .eq. 'x3') then
          call CMLATX39DP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xFrac') then
          call CMLATXF9DP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xyz3') then
          call CMLATXYZ39DP(xf, coords(1, i), coords (2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xyzFrac') then
          call CMLATXYZFRACT9DP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       endif
       call xml_EndElement(xf, 'atom')
    enddo

    call xml_EndElement(xf, 'atomArray')
    call xml_EndElement(xf, 'molecule')
    
  END SUBROUTINE cmlAddMoleculeDP

  
  ! -------------------------------------------------
  ! 2. writes complete SP molecule to xml channel
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddMoleculeSP(xf, natoms, elements, atomRefs, coords, style, id, title, dictref, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in)                    :: natoms          ! number of atoms
    character(len=*), intent(in)           :: elements(*)     ! chemical element types
    real(kind=sp), intent(in)              :: coords(3, *)    ! atomic coordinates
    character(len=*), intent(in), optional :: atomRefs(natoms)    ! id
    character(len=*), intent(in), optional :: id              ! id
    character(len=*), intent(in), optional :: title           ! the title
    character(len=*), intent(in), optional :: dictref         ! the dictionary reference
    character(len=*), intent(in), optional :: fmt             ! format for coords
    character(len=*), intent(in), optional :: style           ! type of coordinates ('x3'for Cartesians, 'xFrac' 
    ! for fractionals; ' ' = default => cartesians)
    ! Flush on entry and exit
    character(len=6) :: id1, id0
    integer          :: i
    character(len=10):: formt, stylei

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(style)) then
       stylei = style
    else
       stylei = 'x3'
    endif

    call stmAddStartTag(xf, 'molecule', id, title, dictref)
    call xml_NewElement(xf, 'atomArray')
    do i = 1, natoms
       write(id0, '(i4)') i
       id0 = adjustl(id0)
       id1 = 'a'
       id1(2:) = id0
       call cmlAddAtom(xf=xf, elem=elements(i), id=trim(id1))
       if (present(atomRefs)) call xml_AddAttribute(xf, 'ref', trim(atomRefs(i)))
       if (stylei .eq. 'x3') then
          call CMLATX39SP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xFrac') then
          call CMLATXF9SP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xyz3') then
          call CMLATXYZ39SP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       elseif (stylei .eq. 'xyzFrac') then
          call CMLATXYZFRACT9SP(xf, coords(1, i), coords(2, i), coords(3, i), formt)
       endif
       call xml_EndElement(xf, 'atom')
    enddo

    call xml_EndElement(xf, 'atomArray')
    call xml_EndElement(xf, 'molecule')
    
    
  END SUBROUTINE cmlAddMoleculeSP
  
  
  ! -------------------------------------------------
  ! 1. writes complete DP molecule to xml channel (No. 2)
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddMolecule3DP(xf, natoms, elements, x, y, z, style, id, title, dictref, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in)                    :: natoms          ! number of atoms
    real(kind=dp), intent(in)               :: x(*)
    real(kind=dp), intent(in)               :: y(*)
    real(kind=dp), intent(in)               :: z(*)
    character(len=*), intent(in)           :: elements(*)     ! chemical element types
    character(len=*), intent(in), optional :: id              ! id
    character(len=*), intent(in), optional :: title           ! the title
    character(len=*), intent(in), optional :: dictref         ! the dictionary reference
    character(len=*), intent(in), optional :: fmt             ! format for coords
    character(len=*), intent(in), optional :: style           ! type of coordinates ('x3' for Cartesians, 'xFrac' 
    ! for fractionals; ' ' = default => cartesians)
    character(len=6)  :: id1, id0
    integer           :: i
    character(len=10) :: formt, stylei

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(style)) then
       stylei = trim(style)
    else
       stylei = 'x3'
    endif

    call stmAddStartTag(xf=xf, name='molecule', id=id, title=title, dictref=dictref)
    call xml_NewElement(xf, 'atomArray')

    do i = 1, natoms
       write(id0, '(i4)') i
       id0 = adjustl(id0)
       id1 = 'a'
       id1(2:) = id0
       call cmlAddAtom(xf=xf, elem=elements(i), id=trim(id1))
       if (trim(stylei) .eq. 'x3') then
          call CMLATX39DP(xf, x(i), y(i), z(i), formt)
       elseif (stylei .eq. 'xFrac') then
          call CMLATXF9DP(xf, x(i), y(i), z(i), formt)
       elseif (stylei .eq. 'xyz3') then
          call CMLATXYZ39DP(xf, x(i), y(i), z(i), formt)
       elseif (stylei .eq. 'xyzFrac') then
          call CMLATXYZFRACT9DP(xf, x(i), y(i), z(i), formt)
       endif
       call xml_EndElement(xf, 'atom')
    enddo

    call xml_EndElement(xf, 'atomArray')
    call xml_EndElement(xf, 'molecule')
    
  END SUBROUTINE cmlAddMolecule3DP
  
  
  ! -------------------------------------------------
  ! 2. writes complete SP molecule to xml channel (No. 2)
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddMolecule3SP(xf, natoms, elements, x, y, z, style, id, title, dictref, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    ! 10 Arguments
    integer, intent(in)                    :: natoms          ! number of atoms
    real(kind=sp), intent(in)               :: x(*)
    real(kind=sp), intent(in)               :: y(*)
    real(kind=sp), intent(in)               :: z(*)
    character(len=*), intent(in)           :: elements(*)      ! chemical element types
    character(len=*), intent(in), optional :: id               ! id
    character(len=*), intent(in), optional :: title            ! the title
    character(len=*), intent(in), optional :: dictref          ! the dictionary reference
    character(len=*), intent(in), optional :: fmt              ! format for coords
    character(len=*), intent(in), optional :: style            ! type of coordinates ('x3' for Cartesians, 'xFrac' 
    ! for fractionals; ' ' = default => cartesians)
    ! Internal variables
    character(len=6)  :: id1, id0
    integer           :: i
    character(len=10) :: formt, stylei

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(style)) then
       stylei = style
    else
       stylei = 'x3'
    endif

    call xml_NewElement(xf, 'molecule')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'title', title)
    call xml_AddAttribute(xf, 'dictref', dictref)
    call xml_NewElement(xf, 'atomArray')
    do i = 1, natoms
       write(id0, '(i4)') i
       id0 = adjustl(id0)
       id1 = 'a'
       id1(2:) = id0
       call cmlAddAtom(xf=xf, elem=elements(i), id=trim(id1))
       if (stylei .eq. 'x3') then
          call CMLATX39SP(xf, x(i), y(i), z(i), formt)
       else if (stylei .eq. 'xFrac') then
          call CMLATXF9SP(xf, x(i), y(i), z(i), formt)
       else if (stylei .eq. 'xyz3') then
          call CMLATXYZ39SP(xf, x(i), y(i), z(i), formt)
       else if (stylei .eq. 'xyzFrac') then
          call CMLATXYZFRACT9SP(xf, x(i), y(i), z(i), formt)
       endif
           call xml_EndElement(xf, 'atom')
    enddo

    call xml_EndElement(xf, 'atomArray')
    call xml_EndElement(xf, 'molecule')

  END SUBROUTINE cmlAddMolecule3SP
  
  ! -------------------------------------------------
  ! writes an <atom> start tag
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddAtom(xf, elem, id, charge, hCount, occupancy, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in), optional           :: charge     ! formalCharge
    integer, intent(in), optional           :: hCount     ! hydrogenCount
    real(kind=sp), intent(in), optional     :: occupancy  ! hydrogenCount
    character(len=*), intent(in), optional  :: elem       ! chemical element name
    character(len=*), intent(in), optional  :: id         ! atom id
    character(len=*), intent(in), optional  :: fmt        ! format

    ! internal Variable
    character(len=10):: formt
    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    call xml_NewElement(xf, 'atom')
    if (present(elem))      call xml_AddAttribute(xf, 'elementType', trim(elem))
    if (present(id))        call xml_AddAttribute(xf, 'id', id)
    if (present(charge))    call xml_AddAttribute(xf, 'formalCharge', charge)
    if (present(hCount))    call xml_AddAttribute(xf, 'hydrogenCount', hCount)
    if (present(occupancy)) call xml_AddAttribute(xf, 'occupancy', occupancy, formt)

  END SUBROUTINE cmlAddAtom
  
  
  ! -------------------------------------------------
  ! 1. append SP coordinates to atom tag
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddCoordinatesSP(xf, x, y, z, style, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)               :: x, y
    real(kind=sp), intent(in), optional     :: z
    character(len=*), intent(in), optional :: style
    character(len=*), intent(in), optional :: fmt

    ! Internal variable
    character(len=10):: formt
    character(len=10):: stylei
    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(style)) then
       stylei = style
    else
       stylei = 'x3'
    endif

    if (present(z) .and. stylei .eq. 'x3') then
       call CMLATX39SP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xFrac') then
       call CMLATXF9SP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xyz3') then
       call CMLATXYZ39SP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xyzFrac') then
       call CMLATXYZFRACT9SP(xf, x, y, z, formt)
    elseif (.not. present(z) .and. stylei .eq. 'xy2') then
       call CMLATXY9SP(xf, x, y, formt)           
    endif

  END SUBROUTINE cmlAddCoordinatesSP

  ! -------------------------------------------------
  ! 2. append DP coordinates to atom tag
  ! -------------------------------------------------

  SUBROUTINE cmlAddCoordinatesDP(xf, x, y, z, style, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf 
    real(kind=dp), intent(in)               :: x, y
    real(kind=dp), intent(in), optional     :: z
    character(len=*), intent(in), optional :: style
    character(len=*), intent(in), optional :: fmt
    
    ! Internal variable
    character(len=10):: formt
    character(len=10):: stylei
    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(style)) then
       stylei = style
    else
       stylei = 'x3'
    endif
    
    if (present(z) .and. stylei .eq. 'x3') then
       call CMLATX39DP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xFrac') then
       call CMLATXF9DP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xyz3') then
       call CMLATXYZ39DP(xf, x, y, z, formt)
    else if (present(z) .and. stylei .eq. 'xyzFrac') then
       call CMLATXYZFRACT9DP(xf, x, y, z, formt)
    else if (.not. present(z) .and. stylei .eq. 'xy2') then
       call CMLATXY9DP(xf, x, y, formt)           
    endif
    
  END SUBROUTINE cmlAddCoordinatesDP

  
  ! -------------------------------------------------
  ! 1. writes a DP <length> element to output channel
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddLengthDP(xf, length, id, atomRef1, atomRef2, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)     :: length     ! length
    character(len=*), intent(in) :: id         ! length id
    character(len=*), intent(in) :: atomRef1   ! ref to first atom
    character(len=*), intent(in) :: atomRef2   ! ref to second atom
    character(len=*), intent(in) :: fmt        ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLLEN9DP(xf, id, atomRef1, atomRef2, length, formt)
  END SUBROUTINE cmlAddLengthDP

  ! -------------------------------------------------
  ! 2. writes a SP <length> element to output channel
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddLengthSP(xf, length, id, atomRef1, atomRef2, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: length     ! the length
    character(len=*), intent(in) :: id         ! length id
    character(len=*), intent(in) :: atomRef1   ! ref to first atom
    character(len=*), intent(in) :: atomRef2   ! ref to second atom
    character(len=*), intent(in) :: fmt        ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLLEN9SP(xf, id, atomRef1, atomRef2, length, formt)
  END SUBROUTINE cmlAddLengthSP


  ! -------------------------------------------------
  ! 1. writes an DP <angle> element to output channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddAngleDP(xf, angle, id, atomRef1, atomRef2, atomRef3, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)     :: angle        ! the angle
    character(len=*), intent(in) :: id           ! angle id
    character(len=*), intent(in) :: atomRef1     ! ref to first atom
    character(len=*), intent(in) :: atomRef2     ! ref to second atom
    character(len=*), intent(in) :: atomRef3     ! ref to third atom
    character(len=*), intent(in) :: fmt          ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLANG9DP(xf, id, atomRef1, atomRef2, atomRef3, angle, formt)
  END SUBROUTINE cmlAddAngleDP

  ! -------------------------------------------------
  ! 2. writes an SP <angle> element to output channel
  ! -------------------------------------------------

  SUBROUTINE cmlAddAngleSP(xf, angle, id, atomRef1, atomRef2, atomRef3, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: angle        ! the angle
    character(len=*), intent(in) :: id           ! angle id
    character(len=*), intent(in) :: atomRef1     ! ref to first atom
    character(len=*), intent(in) :: atomRef2     ! ref to second atom
    character(len=*), intent(in) :: atomRef3     ! ref to third atom
    character(len=*), intent(in) :: fmt          ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLANG9SP(xf, id, atomRef1, atomRef2, atomRef3, angle, formt)
  END SUBROUTINE cmlAddAngleSP


  ! -------------------------------------------------
  ! 1. creates and writes a DP <torsion> element
  ! -------------------------------------------------

  SUBROUTINE cmlAddTorsionDP(xf, torsion, id, atomRef1, atomRef2, atomRef3, atomRef4, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)     :: torsion         ! the torsion
    character(len=*), intent(in) :: id              ! torsion id
    character(len=*), intent(in) :: atomRef1        ! ref to first atom
    character(len=*), intent(in) :: atomRef2        ! ref to second atom
    character(len=*), intent(in) :: atomRef3        ! ref to third atom
    character(len=*), intent(in) :: atomRef4        ! ref to fourth atom
    character(len=*), intent(in) :: fmt             ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLTOR9DP(xf, id, atomRef1, atomRef2, atomRef3, atomRef4, torsion, formt)
  END SUBROUTINE cmlAddTorsionDP
  
  ! -------------------------------------------------
  ! 2. creates and writes a SP <torsion> element
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddTorsionSP(xf, torsion, id, atomRef1, atomRef2, atomRef3, atomRef4, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: torsion         ! the torsion
    character(len=*), intent(in) :: id              ! torsion id
    character(len=*), intent(in) :: atomRef1        ! ref to first atom
    character(len=*), intent(in) :: atomRef2        ! ref to second atom
    character(len=*), intent(in) :: atomRef3        ! ref to third atom
    character(len=*), intent(in) :: atomRef4        ! ref to fourth atom
    character(len=*), intent(in) :: fmt             ! format

    optional         :: fmt
    character(len=10):: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif

    ! Flush on entry and exit
    call CMLTOR9SP(xf, id, atomRef1, atomRef2, atomRef3, atomRef4, torsion, formt)
  END SUBROUTINE cmlAddTorsionSP


  ! -------------------------------------------------
  ! 1. creates and writes an SP Lattice element
  ! -------------------------------------------------

  SUBROUTINE cmlAddLatticeSP(xf, cell, units, title, id, dictref, conv, lattType, spaceType, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)               :: cell(3,3)
    character(len=*), intent(in), optional :: units       
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! title
    character(len=*), intent(in), optional :: dictref      ! dictref
    character(len=*), intent(in), optional :: conv         ! convention
    character(len=*), intent(in), optional :: lattType 
    character(len=*), intent(in), optional :: spaceType    !
    character(len=*), intent(in), optional :: fmt         

    integer :: i
    character(len=10) :: formt
    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'   
    endif

    call xml_NewElement(xf, 'lattice')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(lattType)) call xml_AddAttribute(xf, 'latticeType', lattType)
    if (present(spaceType)) call xml_AddAttribute(xf, 'spaceType', spaceType)

    do i = 1,3
       call xml_NewElement(xf, 'latticeVector')
       if (present(units)) call xml_AddAttribute(xf, 'units', units)
       call xml_AddAttribute(xf, 'dictRef', 'cml:latticeVector')
       call xml_AddPcdata(xf, cell(1,i), formt)
       call xml_AddPcdata(xf, cell(2,i), formt, space=.true.)
       call xml_AddPcdata(xf, cell(3,i), formt, space=.true.)
       call xml_EndElement(xf, 'latticeVector')
    enddo
    call xml_EndElement(xf, 'lattice')
    
  END SUBROUTINE cmlAddLatticeSP


  ! -------------------------------------------------
  ! 2. creates and writes DP Lattice element
  ! -------------------------------------------------

  SUBROUTINE cmlAddLatticeDP(xf, cell, units, title, id, dictref, conv, lattType, spaceType, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)               :: cell(3,3)
    character(len=*), intent(in), optional :: units       
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! title
    character(len=*), intent(in), optional :: dictref      ! dictref
    character(len=*), intent(in), optional :: conv         ! 
    character(len=*), intent(in), optional :: lattType     ! 
    character(len=*), intent(in), optional :: spaceType    !
    character(len=*), intent(in), optional :: fmt         

    integer :: i
    character(len=10) :: formt
    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'   
    endif

    call xml_NewElement(xf, 'lattice')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    if (present(lattType)) call xml_AddAttribute(xf, 'latticeType', lattType)
    if (present(spaceType)) call xml_AddAttribute(xf, 'spaceType', spaceType)

    do i = 1,3
       call xml_NewElement(xf, 'latticeVector')
       if (present(units)) call xml_AddAttribute(xf, 'units', units)
       call xml_AddAttribute(xf, 'dictRef', 'cml:latticeVector')
       call xml_AddArray(xf,cell(1:,i),fmt)
!       call xml_AddPcdata(xf, cell(1,i), formt)
!       call xml_AddPcdata(xf, cell(2,i), formt, space=.true.)
!       call xml_AddPcdata(xf, cell(3,i), formt, space=.true.)
       call xml_EndElement(xf, 'latticeVector')
    enddo
    call xml_EndElement(xf, 'lattice')

  END SUBROUTINE cmlAddLatticeDP


  ! -------------------------------------------------
  ! 1. creates and writes a DP <cell> element
  ! -------------------------------------------------

  SUBROUTINE cmlAddCrystalDP(xf, a, b, c, alpha, beta, gamma, id, title, dictref, conv, lenunits, angunits, spaceType, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)               :: a, b, c      ! cell parameters
    real(kind=dp), intent(in)               :: alpha        ! alpha cell parameter
    real(kind=dp), intent(in)               :: beta         ! beta cell parameter
    real(kind=dp), intent(in)               :: gamma        ! gamma cell parameter
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! title
    character(len=*), intent(in), optional :: dictref      ! dictref
    character(len=*), intent(in), optional :: conv         ! convention
    character(len=*), intent(in), optional :: lenunits     ! units for length (default = angstrom)
    character(len=*), intent(in), optional :: angunits     ! units for angles (default = degree)
    character(len=*), intent(in), optional :: spaceType    ! spacegroup
    character(len=*), intent(in), optional :: fmt          ! format

    ! Flush on entry and exit
    character(len=30) ::  lunits, aunits

    if (present(lenunits)) then
       lunits = lenunits
    else
       lunits = U_ANGSTR
    endif
    if (present(angunits)) then
       aunits = angunits
    else
       aunits = U_DEGREE
    endif

    call xml_NewElement(xf=xf, name='crystal')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictRef)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    call stmAddScalar(xf=xf, value=a, title='a', dictref='cml:a', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=b, title='b', dictref='cml:b', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=c, title='c', dictref='cml:c', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=alpha, title='alpha', dictref='cml:alpha', units=aunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=beta,  title='beta',  dictref='cml:beta',  units=aunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=gamma, title='gamma', dictref='cml:gamma', units=aunits, fmt=fmt)
    if (present(spaceType)) then
      call xml_NewElement(xf, 'symmetry')
      call xml_AddAttribute(xf, 'spaceGroup', spaceType)
      call xml_EndElement(xf, 'symmetry')
    endif
    call xml_EndElement(xf, 'crystal')

  END SUBROUTINE cmlAddCrystalDP

  ! -------------------------------------------------
  ! 2. creates and writes a SP <cell> element
  ! -------------------------------------------------

  SUBROUTINE cmlAddCrystalSP(xf, a, b, c, alpha, beta, gamma, id, title, dictref, conv, lenunits, angunits, spaceType, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: a, b, c      ! cell parameters
    real(kind=sp), intent(in)     :: alpha        ! alpha cell parameter
    real(kind=sp), intent(in)     :: beta         ! beta cell parameter
    real(kind=sp), intent(in)     :: gamma        ! gamma cell parameter
    character(len=*), intent(in), optional :: id           ! id
    character(len=*), intent(in), optional :: title        ! title
    character(len=*), intent(in), optional :: dictref      ! dictref
    character(len=*), intent(in), optional :: conv         ! convention
    character(len=*), intent(in), optional :: lenunits     ! units for length (' ' = angstrom)
    character(len=*), intent(in), optional :: angunits     ! units for angles (' ' = degree)
    character(len=*), intent(in), optional :: spaceType    ! spacegroup
    character(len=*), intent(in), optional :: fmt          ! format

    ! Flush on entry and exit
    character(len=30) :: lunits, aunits
    character(len=10) :: formt

    if (present(fmt)) then
       formt = fmt
    else
       formt = '(f8.3)'
    endif
    if (present(lenunits)) then
       lunits = lenunits
    else
       lunits = U_ANGSTR
    endif
    if (present(angunits)) then
       aunits = angunits
    else
       aunits = U_DEGREE
    endif

    call xml_NewElement(xf=xf, name='crystal')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictRef)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    call stmAddScalar(xf=xf, value=a, title='a', dictref='cml:a', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=b, title='b', dictref='cml:b', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=c, title='c', dictref='cml:c', units=lunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=alpha, title='alpha', dictref='cml:alpha', units=aunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=beta,  title='beta',  dictref='cml:beta', units=aunits, fmt=fmt)
    call stmAddScalar(xf=xf, value=gamma, title='gamma', dictref='cml:gamma', units=aunits, fmt=fmt)
    if (present(spaceType)) then
      call xml_NewElement(xf, 'symmetry')
      call xml_AddAttribute(xf, 'spaceGroup', spaceType)
      call xml_EndElement(xf, 'symmetry')
    endif
    call xml_EndElement(xf, 'crystal')


  END SUBROUTINE cmlAddCrystalSP
  
  
  ! -------------------------------------------------
  ! 1. creates and writes an DP <eigen> element
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddEigenvalueDP(xf, n, eigvec, eigval, id, title, dictref, fmt)


    implicit none
    type(xmlf_t), intent(inout)            :: xf
    integer, intent(in)                    :: n              ! number of elements
    real(kind=dp), intent(in)              :: eigvec(n, *)   ! eigenvectors
    real(kind=dp), intent(in)              :: eigval(*)      ! eigenvalues
    character(len=*), intent(in), optional :: id             ! id
    character(len=*), intent(in), optional :: title          ! title
    character(len=*), intent(in), optional :: dictref        ! dictionary reference
    character(len=*), intent(in), optional :: fmt            ! format

    ! Flush on entry and exit
    call xml_NewElement(xf, 'eigen')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictRef)
    call stmAddArray(xf=xf, nvalue=n, array=eigval, title='eigenvalues', dictref=dictRef, fmt=fmt)
    call stmAddMatrix(xf=xf, ncols=n, nrows=n, matrix=eigvec, title='eigenvectors', fmt=fmt)
    call xml_EndElement(xf, 'eigen')

  END SUBROUTINE cmlAddEigenvalueDP



  ! -------------------------------------------------
  ! 2. creates and writes an SP <eigen> element
  ! -------------------------------------------------
  
  SUBROUTINE cmlAddEigenvalueSP(xf, n, eigvec, eigval, id, title, dictref, fmt)


    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in)          :: n              ! number of elements
    real(kind=sp), intent(in)     :: eigvec(n, *) ! eigenvectors
    real(kind=sp), intent(in)     :: eigval(*)      ! eigenvalues
    character(len=*), intent(in), optional :: id             ! id
    character(len=*), intent(in), optional :: title          ! title
    character(len=*), intent(in), optional :: dictref        ! dictionary reference
    character(len=*), intent(in), optional :: fmt            ! format


    ! Flush on entry and exit
    call xml_NewElement(xf, 'eigen')
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictRef))   call xml_AddAttribute(xf, 'dictRef', dictref)
    call stmAddArray(xf=xf, nvalue=n, array=eigval, title='eigenvalues', dictref=dictRef, fmt=fmt)
    call stmAddMatrix(xf=xf, ncols=n, nrows=n, matrix=eigvec, title='eigenvectors', fmt=fmt)
    call xml_EndElement(xf, 'eigen')

  END SUBROUTINE cmlAddEigenvalueSP


  SUBROUTINE cmlAddMetadataCh(xf, name, content, conv)
    
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: content
    character(len=*), optional, intent(in) :: conv
    
    call xml_NewElement(xf, 'metadata')
    call xml_AddAttribute(xf, 'name', name)
    call xml_AddAttribute(xf, 'content', content)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    call xml_EndElement(xf, 'metadata')

  END SUBROUTINE cmlAddMetadataCh

  SUBROUTINE cmlAddMetadataI(xf, name, content, conv)
    
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: content
    character(len=*), intent(in), Optional :: conv
    
    call xml_NewElement(xf, 'metadata')
    call xml_AddAttribute(xf, 'name', name)
    call xml_AddAttribute(xf, 'content', content)
    if (present(conv)) call xml_AddAttribute(xf, 'convention', conv)
    call xml_EndElement(xf, 'metadata')

  END SUBROUTINE cmlAddMetadataI


  ! -------------------------------------------------
  ! 1. creates and writes an Char <parameter> element
  ! -------------------------------------------------


  SUBROUTINE cmlAddParameterCh(xf, value, ref, id, title, conv, cons, units, name, role, dictref)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: value 
    character(len=*), optional :: ref 
    character(len=*), optional :: title
    character(len=*), optional :: id
    character(len=*), optional :: conv
    character(len=*), optional :: cons
    character(len=*), optional :: units
    character(len=*), optional :: name
    character(len=*), optional :: role
    character(len=*), intent(in), optional :: dictref

    call xml_NewElement(xf, 'parameter')
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(cons))    call xml_AddAttribute(xf, 'constraint', cons)
    if (present(name))    call xml_AddAttribute(xf, 'name', name)
    if (present(role))    call xml_AddAttribute(xf, 'role', role)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

    call stmAddScalar(xf,value,dataType="xsd:string",units=units)

    call xml_EndElement(xf, 'parameter')

  END SUBROUTINE CMLADDPARAMETERCH


  ! -------------------------------------------------
  ! 2. creates and writes an SP <parameter> element
  ! -------------------------------------------------


  SUBROUTINE cmlAddParameterSP(xf, value, ref, title, id, conv, cons, units, name, role, dictref, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp) :: value 
    character(len=*), optional :: ref 
    character(len=*), optional :: title
    character(len=*), optional :: id
    character(len=*), optional :: conv
    character(len=*), optional :: cons
    character(len=*), optional :: units
    character(len=*), optional :: name
    character(len=*), optional :: role
    character(len=*), optional :: fmt
    character(len=*), intent(in), optional :: dictref

    call xml_NewElement(xf, 'parameter')
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(cons))    call xml_AddAttribute(xf, 'constraint', cons)
    if (present(name))    call xml_AddAttribute(xf, 'name', name)
    if (present(role))    call xml_AddAttribute(xf, 'role', role)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

    call stmAddScalar(xf,value,dataType="xsd:float",fmt=fmt,units=units)

    call xml_EndElement(xf, 'parameter')

  END SUBROUTINE CMLADDPARAMETERSP


  ! -------------------------------------------------
  ! 3. creates and writes an DP <parameter> element
  ! -------------------------------------------------


  SUBROUTINE cmlAddParameterDP(xf, value, ref, title, id, conv, cons, units, name, role, dictref, fmt)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp) :: value 
    character(len=*), optional :: ref 
    character(len=*), optional :: title
    character(len=*), optional :: id
    character(len=*), optional :: conv
    character(len=*), optional :: cons
    character(len=*), optional :: units
    character(len=*), optional :: name
    character(len=*), optional :: role
    character(len=*), intent(in), optional :: dictref
    character(len=*), optional :: fmt    

    call xml_NewElement(xf, 'parameter')
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(cons))    call xml_AddAttribute(xf, 'constraint', cons)
    if (present(name))    call xml_AddAttribute(xf, 'name', name)
    if (present(role))    call xml_AddAttribute(xf, 'role', role)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

    call stmAddScalar(xf,value,dataType="xsd:double",fmt=fmt,units=units)

    call xml_EndElement(xf, 'parameter')

  END SUBROUTINE CMLADDPARAMETERDP


  ! -------------------------------------------------
  ! 4. creates and writes an Integer <parameter> element
  ! -------------------------------------------------


  SUBROUTINE cmlAddParameterI(xf, value, ref, id, title, conv, cons, units, name, role, dictref)

    implicit none
    type(xmlf_t), intent(inout) :: xf
    integer :: value 
    character(len=*), optional :: ref 
    character(len=*), optional :: title
    character(len=*), optional :: id
    character(len=*), optional :: conv
    character(len=*), optional :: cons
    character(len=*), optional :: units
    character(len=*), optional :: name
    character(len=*), optional :: role
    character(len=*), intent(in), optional :: dictref

    call xml_NewElement(xf, 'parameter')
    if (present(ref))   call xml_AddAttribute(xf, 'ref', ref)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(id))    call xml_AddAttribute(xf, 'id', id)
    if (present(conv))  call xml_AddAttribute(xf, 'convention', conv)
    if (present(cons))  call xml_AddAttribute(xf, 'constraint', cons)
    if (present(name))  call xml_AddAttribute(xf, 'name', name)
    if (present(role))  call xml_AddAttribute(xf, 'role', role)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

    call stmAddScalar(xf,value,dataType="xsd:integer",units=units)

    call xml_EndElement(xf, 'parameter')

  END SUBROUTINE CMLADDPARAMETERI


  SUBROUTINE cmlAddParameterLG(xf, value, ref, id, title, conv, cons, units, name, role, dictref)

    implicit none
    type(xmlf_t),     intent(inout)        :: xf
    logical,          intent(in)           :: value 
    character(len=*), intent(in), optional :: ref 
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: conv
    character(len=*), intent(in), optional :: cons
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: role
    character(len=*), intent(in), optional :: dictref

    call xml_NewElement(xf, 'parameter')
    if (present(ref))     call xml_AddAttribute(xf, 'ref', ref)
    if (present(title))   call xml_AddAttribute(xf, 'title', title)
    if (present(id))      call xml_AddAttribute(xf, 'id', id)
    if (present(conv))    call xml_AddAttribute(xf, 'convention', conv)
    if (present(cons))    call xml_AddAttribute(xf, 'constraint', cons)
    if (present(name))    call xml_AddAttribute(xf, 'name', name)
    if (present(role))    call xml_AddAttribute(xf, 'role', role)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)

    call stmAddScalar(xf,value,dataType="xsd:boolean",units=units)

    call xml_EndElement(xf, 'parameter')

  END SUBROUTINE CMLADDPARAMETERLG


! =================================================
! basic CML routines
! =================================================

  
  ! -------------------------------------------------
  ! 1. adds DP xyz3 to start tag
  ! -------------------------------------------------
  
  SUBROUTINE CMLATXYZ39DP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp)      :: x3, y3, z3 ! coordinates
    character(len=*)  :: fmt        ! format (default '(f8.3)')
    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'xyz3', trim(adjustl(x))//' '//trim(adjustl(y))//' '//trim(adjustl(z)) )

  END SUBROUTINE CMLATXYZ39DP


  ! -------------------------------------------------
  ! 2. adds SP xyz3 to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATXYZ39SP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp)      :: x3, y3, z3 ! coordinates
    character(len=*)  :: fmt        ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'xyz3', trim(adjustl(x))//' '//trim(adjustl(y))//' '//trim(adjustl(z)) )

  END SUBROUTINE CMLATXYZ39SP
  
  ! -------------------------------------------------
  ! 1. adds DP xyzFrac to start tag
  ! -------------------------------------------------
  
  SUBROUTINE CMLATXYZFRACT9DP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp)      :: x3, y3, z3 ! coordinates
    character(len=*)  :: fmt        ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'xyzFrac', trim(adjustl(x))//' '//trim(adjustl(y))//' '//trim(adjustl(z)) )

  END SUBROUTINE CMLATXYZFRACT9DP

  ! -------------------------------------------------
  ! 2. adds SP xyzFrac to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATXYZFRACT9SP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: x3, y3, z3 ! coordinates
    character(len=*), intent(in) :: fmt        ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'xyzFrac', trim(adjustl(x))//' '//trim(y)//' '//trim(z))

  END SUBROUTINE CMLATXYZFRACT9SP


  ! -------------------------------------------------
  ! 1. adds DP x3, y3, z3 to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATX39DP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)     :: x3, y3, z3 ! coordinates
    character(len=*), intent(in) :: fmt        ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'x3', x)
    call xml_AddAttribute(xf, 'y3', y)
    call xml_AddAttribute(xf, 'z3', z)

  END SUBROUTINE CMLATX39DP

  ! -------------------------------------------------
  ! 2. adds SP x3, y3, z3 to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATX39SP(xf, x3, y3, z3, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp), intent(in)     :: x3, y3, z3 ! coordinates
    character(len=*), intent(in) :: fmt        ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) x3
    write(y,fmt) y3
    write(z,fmt) z3

    call xml_AddAttribute(xf, 'x3', x)
    call xml_AddAttribute(xf, 'y3', y)
    call xml_AddAttribute(xf, 'z3', z)

  END SUBROUTINE CMLATX39SP


  ! -------------------------------------------------
  ! 1. adds DP xFract, yFract, zFract to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATXF9DP(xf, xFract, yFract, zFract, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp), intent(in)     :: xFract, yFract, zFract ! coordinates
    character(len=*), intent(in) :: fmt                    ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) xFract
    write(y,fmt) yFract
    write(z,fmt) zFract

    call xml_AddAttribute(xf, 'xFract', x)
    call xml_AddAttribute(xf, 'yFract', y)
    call xml_AddAttribute(xf, 'zFract', z)

  END SUBROUTINE CMLATXF9DP
  
  ! -------------------------------------------------
  ! 2. adds SP xfrac, yFractractrac, zFractrac to start tag
  ! -------------------------------------------------
  
  SUBROUTINE CMLATXF9SP(xf, xFract, yFract, zFract, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp)     :: xFract, yFract, zFract   ! fractional coordinates
    character(len=*) :: fmt                      ! format (default '(f8.3)')

    character(len=45) :: x, y, z

    write(x,fmt) xFract
    write(y,fmt) yFract
    write(z,fmt) zFract

    call xml_AddAttribute(xf, 'xFract', x)
    call xml_AddAttribute(xf, 'yFract', y)
    call xml_AddAttribute(xf, 'zFract', z)

  END SUBROUTINE CMLATXF9SP


  ! -------------------------------------------------
  ! 1. adds DP x2, y2 to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATXY9DP(xf, x2, y2, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=dp)     :: x2, y2   ! coordinates
    character(len=*) :: fmt      ! format (default f8.3)

    character(len=45) :: x, y

    write(x,fmt) x2
    write(y,fmt) y2

    call xml_AddAttribute(xf, 'x2', x)
    call xml_AddAttribute(xf, 'y2', y)
    call xml_AddPcdata(xf, '>')

  END SUBROUTINE CMLATXY9DP

  ! -------------------------------------------------
  ! 2. adds SP x2, y2 to start tag
  ! -------------------------------------------------

  SUBROUTINE CMLATXY9SP(xf, x2, y2, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    real(kind=sp)     :: x2, y2   ! coordinates
    character(len=*) :: fmt      ! format (default f8.3)

    character(len=45) :: x, y

    write(x,fmt) x2
    write(y,fmt) y2

    call xml_AddAttribute(xf, 'x2', x)
    call xml_AddAttribute(xf, 'y2', y)
    call xml_AddPcdata(xf, '>')

  END SUBROUTINE CMLATXY9SP


  ! -------------------------------------------------
  ! 1. creates a DP <length> element
  ! -------------------------------------------------

  SUBROUTINE CMLLEN9DP(xf, id, atomRef1, atomRef2, length, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id           ! length id
    character(len=*) :: atomRef1     ! ref to first atom
    character(len=*) :: atomRef2     ! ref to second atom
    real(kind=dp)     :: length       ! the length
    character(len=*) :: fmt          ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)

    call xml_NewElement(xf, 'length')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs2', temp)
    call xml_AddPcdata(xf, length, fmt)
    call xml_EndElement(xf, 'length')

  END SUBROUTINE CMLLEN9DP
  
  ! -------------------------------------------------
  ! 2. creates a SP <length> element
  ! -------------------------------------------------
  
  SUBROUTINE CMLLEN9SP(xf, id, atomRef1, atomRef2, length, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id           ! length id
    character(len=*) :: atomRef1     ! ref to first atom
    character(len=*) :: atomRef2     ! ref to second atom
    real(kind=sp)     :: length       ! the length
    character(len=*) :: fmt          ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)

    call xml_NewElement(xf, 'length')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs2', temp)
    call xml_AddPcdata(xf, length, fmt)
    call xml_EndElement(xf, 'length')

  END SUBROUTINE CMLLEN9SP


  ! -------------------------------------------------
  ! 1. creates a DP <angle> element
  ! -------------------------------------------------

  SUBROUTINE CMLANG9DP(xf, id, atomRef1, atomRef2, atomRef3, angle, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id              ! angle id
    character(len=*) :: atomRef1        ! ref to first atom
    character(len=*) :: atomRef2        ! ref to second atom
    character(len=*) :: atomRef3        ! ref to third atom
    real(kind=dp)     :: angle           ! the angle
    character(len=*) :: fmt             ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)//' '//adjustl(atomRef3)

    call xml_NewElement(xf, 'angle')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs3', temp)
    call xml_AddPcdata(xf, angle, fmt)
    call xml_EndElement(xf, 'angle')

  END SUBROUTINE CMLANG9DP

  ! -------------------------------------------------
  ! 2. creates a SP <angle> element
  ! -------------------------------------------------

  SUBROUTINE CMLANG9SP(xf, id, atomRef1, atomRef2, atomRef3, angle, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id              ! angle id
    character(len=*) :: atomRef1        ! ref to first atom
    character(len=*) :: atomRef2        ! ref to second atom
    character(len=*) :: atomRef3        ! ref to third atom
    real(kind=sp)     :: angle           ! the angle
    character(len=*) :: fmt             ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)//' '//adjustl(atomRef3)

    call xml_NewElement(xf, 'angle')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs3', temp)
    call xml_AddPcdata(xf, angle, fmt)
    call xml_EndElement(xf, 'angle')

  END SUBROUTINE CMLANG9SP


  ! -------------------------------------------------
  ! 1. creates a DP <torsion> element
  ! -------------------------------------------------
  
  SUBROUTINE CMLTOR9DP(xf, id, atomRef1, atomRef2, atomRef3, atomRef4, torsion, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id              ! torsion id
    character(len=*) :: atomRef1        ! ref to first atom
    character(len=*) :: atomRef2        ! ref to second atom
    character(len=*) :: atomRef3        ! ref to third atom
    character(len=*) :: atomRef4        ! ref to fourth atom
    real(kind=dp)    :: torsion         ! the torsion
    character(len=*) :: fmt             ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)//' '//adjustl(atomRef3)//' '//adjustl(atomRef4)

    call xml_NewElement(xf, 'torsion')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs4', temp)
    call xml_AddPcdata(xf, torsion, fmt)
    call xml_EndElement(xf, 'torsion')

  END SUBROUTINE CMLTOR9DP

  ! -------------------------------------------------
  ! 2. creates a SP <torsion> element
  ! -------------------------------------------------

  SUBROUTINE CMLTOR9SP(xf, id, atomRef1, atomRef2, atomRef3, atomRef4, torsion, fmt)
    implicit none
    type(xmlf_t), intent(inout) :: xf
    character(len=*) :: id              ! torsion id
    character(len=*) :: atomRef1        ! ref to first atom
    character(len=*) :: atomRef2        ! ref to second atom
    character(len=*) :: atomRef3        ! ref to third atom
    character(len=*) :: atomRef4        ! ref to fourth atom
    real(kind=sp)     :: torsion         ! the torsion
    character(len=*) :: fmt             ! format
    character(len=20) :: temp

    temp = atomRef1//' '//adjustl(atomRef2)//' '//adjustl(atomRef3)//' '//adjustl(atomRef4)

    call xml_NewElement(xf, 'torsion')
    call xml_AddAttribute(xf, 'id', id)
    call xml_AddAttribute(xf, 'atomRefs4', temp)
    call xml_AddPcdata(xf, torsion, fmt)
    call xml_EndElement(xf, 'torsion')

  END SUBROUTINE CMLTOR9SP

end module flib_wcml
