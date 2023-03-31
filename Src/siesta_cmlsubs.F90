! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
Module siesta_cmlsubs

  use siesta_cml, only:   cml_p, mainXML

  Use siesta_cml, only: cmlBeginFile, cmlStartCml, cmlNamespaceAttribute
  Use siesta_cml, only: cmlStartMetadataList, cmlAddMetadata
  Use siesta_cml, only: cmlEndMetadataList, cmlEndCml, cmlFinishFile

  public :: siesta_cml_init, siesta_cml_exit

  Private

  Contains

    Subroutine siesta_cml_init( )
      Use fdf,   Only : fdf_boolean, fdf_string
      Use files, only : slabel, label_length
      Use parallel, only : nodes, ionode
      Use version_info
      use m_uuid, only: generate_uuid
      Use m_timestamp, only: datestring

      Character(len=label_length+4) :: fname
      Character(len=10)             :: nodes_str 

      fname = ' '

      If (IOnode) Then
         cml_p = fdf_boolean( 'XML.Write', .false. )
      Else
         cml_p = .False.
      Endif !IOnode

      If (cml_p) Then
         Write(fname,'(a)') Trim(slabel)//'.xml'
         Call cmlBeginFile(mainXML, trim(fname), unit=-1)

         Call cmlStartCml(mainXML)

         ! Namespace declarations are faked as attributes
         !
         call cmlNamespaceAttribute(mainXML, 'xmlns', 'http://www.xml-cml.org/schema')
         Call cmlNamespaceAttribute(mainXML, 'xmlns:siesta', 'http://www.uam.es/siesta/namespace')
         Call cmlNamespaceAttribute(mainXML, 'xmlns:siestaUnits', 'http://www.uam.es/siesta/namespace/units')
         call cmlNamespaceAttribute(mainXML, 'xmlns:xsd', 'http://www.w3.org/2001/XMLSchema')
         call cmlNamespaceAttribute(mainXML, 'xmlns:fpx', 'http://www.uszla.me.uk/fpx')
         call cmlNamespaceAttribute(mainXML, 'xmlns:dc', 'http://purl.org/dc/elements/1.1/')
         call cmlNamespaceAttribute(mainXML, 'xmlns:units', 'http://www.uszla.me.uk/FoX/units')
         call cmlNamespaceAttribute(mainXML, 'xmlns:cmlUnits', 'http://www.xml-cml.org/units/units')
         call cmlNamespaceAttribute(mainXML, 'xmlns:siUnits', 'http://www.xml-cml.org/units/siUnits')
         call cmlNamespaceAttribute(mainXML, 'xmlns:atomicUnits', 'http://www.xml-cml.org/units/atomic')

         Call cmlStartMetadataList(mainXML)
         Call cmlAddMetadata(mainXML, name='siesta:Program', content='Siesta')
         Call cmlAddMetadata(mainXML, name='siesta:Version', content=version_str)
         Call cmlAddMetadata(mainXML, name='siesta:Arch',    content=siesta_arch)
         Call cmlAddMetadata(mainXML, name='siesta:Flags',   content=fflags)
         Call cmlAddMetadata(mainXML, name='siesta:PPFlags',   content=fppflags)
         Call cmlAddMetadata(mainXML, name='siesta:StartTime',content=datestring())
         !
         ! Generate the uuid at some Siesta top-level and pass it here
         !
         Call cmlAddMetadata(mainXML, name='siesta:run_UUID',content=generate_uuid(1))
         
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
