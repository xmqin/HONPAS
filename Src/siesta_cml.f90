! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
Module siesta_cml

  Use FoX_wcml, only: cmlStartModule, cmlEndModule
  Use FoX_wcml, only: cmlStartStep, cmlEndStep
  Use FoX_wcml, only: cmlStartPropertyList, cmlEndPropertyList
  Use FoX_wcml, only: cmlStartParameterList, cmlEndParameterList
  Use FoX_wcml, only: cmlAddProperty, cmlAddLattice, cmlAddKPoint
  Use FoX_wcml, only: cmlAddMolecule, cmlAddParameter, cmlAddCrystal
  Use FoX_common, only: str_fox => str
  Use FoX_wxml, only: xmlf_t      ! help pgf95...
  Use FoX_wcml, only: cmlBeginFile, cmlAddNamespace, cmlStartCml
  Use FoX_wcml, only: cmlStartMetadataList, cmlAddMetadata
  Use FoX_wcml, only: cmlEndMetadataList, cmlEndCml, cmlFinishFile
  Use FoX_common, only: FoX_set_fatal_warnings, FoX_set_fatal_errors

  Implicit None
  Logical, public      :: cml_p = .False.
  Type(xmlf_t), public, save :: mainXML

  Public :: str_fox
  public :: cmlStartModule, cmlEndModule
  public :: cmlStartStep, cmlEndStep
  public :: cmlStartPropertyList, cmlEndPropertyList
  public :: cmlStartParameterList, cmlEndParameterList
  public :: cmlAddProperty, cmlAddLattice, cmlAddKPoint
  public :: cmlAddMolecule, cmlAddParameter, cmlAddCrystal
  public :: cmlBeginFile, cmlAddNamespace, cmlStartCml
  public :: cmlStartMetadataList, cmlAddMetadata
  public :: cmlEndMetadataList, cmlEndCml, cmlFinishFile
  public :: FoX_set_fatal_warnings, FoX_set_fatal_errors

  private
  
End Module siesta_cml
