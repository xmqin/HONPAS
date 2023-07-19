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
Module siesta_cmlsubs

  use siesta_cml, only:   cml_p, mainXML

  Use siesta_cml, only: cmlBeginFile, cmlAddNamespace, cmlStartCml
  Use siesta_cml, only: cmlStartMetadataList, cmlAddMetadata
  Use siesta_cml, only: cmlEndMetadataList, cmlEndCml, cmlFinishFile
  Use siesta_cml, only: FoX_set_fatal_warnings, FoX_set_fatal_errors

  public :: siesta_cml_init, siesta_cml_exit

  Private

  Contains

    Subroutine siesta_cml_init( )
      Use fdf,   Only : fdf_boolean, fdf_string
      Use files, only : slabel, label_length
      Use parallel, only : nodes, ionode
      Use version_info
      Use m_timestamp, only: datestring

      Character(len=label_length+4) :: fname
      Character(len=10)             :: nodes_str 

      fname = ' '

      If (IOnode) Then
         cml_p = fdf_boolean( 'XML.Write', .True. )
         cml_p = fdf_boolean( 'Write.XML', cml_p ) ! For legacy symbol
         call FoX_set_fatal_errors(fdf_boolean('XML.AbortOnErrors', .false.))
         call FoX_set_fatal_warnings(fdf_boolean('XML.AbortOnWarnings', .false.))
      Else
         cml_p = .False.
      Endif !IOnode

      If (cml_p) Then
         Write(fname,'(a)') Trim(slabel)//'.xml'
         Call cmlBeginFile(mainXML, trim(fname), unit=-1)
         Call cmlAddNamespace(mainXML, 'siesta', 'http://www.uam.es/siesta/namespace')
         Call cmlAddNamespace(mainXML, 'siestaUnits', 'http://www.uam.es/siesta/namespace/units')
         Call cmlStartCml(mainXML, convention="CMLComp")
         Call cmlStartMetadataList(mainXML)
         Call cmlAddMetadata(mainXML, name='siesta:Program', content='Siesta')
         Call cmlAddMetadata(mainXML, name='siesta:Version', content=version_str)
         Call cmlAddMetadata(mainXML, name='siesta:Arch',    content=siesta_arch)
         Call cmlAddMetadata(mainXML, name='siesta:Flags',   content=fflags)
         Call cmlAddMetadata(mainXML, name='siesta:StartTime',content=datestring()) 
         If (nodes>1) Then
           Call cmlAddMetadata(mainXML, name='siesta:Mode', content='Parallel')
         Else
           Call cmlAddMetadata(mainXML, name='siesta:Mode', content='Serial')
         Endif
         ! Avoid using 'str' function, which causes trouble with PGI compilers
         write(nodes_str, "(i10)") nodes
         Call cmlAddMetadata(mainXML, name='siesta:Nodes', content=nodes_str)
#ifdef CDF
         Call cmlAddMetadata(mainXML, name='siesta:NetCDF',  content='true')
#else
         Call cmlAddMetadata(mainXML, name='siesta:NetCDF',  content='false')
#endif
         Call cmlEndMetadataList(mainXML)
      Endif !cml_p
      
    End Subroutine siesta_cml_init
       
    Subroutine siesta_cml_exit

      use m_timestamp, only : datestring


      If (cml_p) Then
        Call cmlAddMetadata(mainXML, name='siesta:EndTime',content=datestring())
        Call cmlEndCml(mainXML)
        Call cmlFinishFile(mainXML)
      Endif !cml_p

    End Subroutine siesta_cml_exit

End Module siesta_cmlsubs
