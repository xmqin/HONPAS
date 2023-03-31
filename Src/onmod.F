! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module on_main

        use precision, only: dp, on_p

        implicit none

        integer,  pointer, save :: numc(:)
        integer,  pointer, save :: numcold(:)
        integer,  pointer, save :: listc(:,:)
        integer,  pointer, save :: listcold(:,:)
        integer,  pointer, save :: ncG2L(:)
        integer,  pointer, save :: ncL2G(:)
        integer,  pointer, save :: ncP2T(:)
        integer,  pointer, save :: ncT2P(:)
C
C  Band pointers
C
C  nbandsloc = no. of bands to which local orbitals belong
C  nbL2G     = pointer from local band to global band index
C
        integer,           save :: nbandsloc = 1
        integer,  pointer, save :: nbL2G(:)
        integer,  pointer, save :: nbG2L(:)

        real(dp), pointer, save :: c(:,:,:)
        real(dp), pointer, save :: cold(:,:,:)
        real(dp), pointer, save :: xi(:,:,:)

        real(on_p),   pointer, save :: g(:,:,:)
        real(on_p),   pointer, save :: hg(:,:,:)

      end module on_main

      module on_core

        use precision, only: dp

        implicit none

        integer,  pointer, save :: numct(:)
        integer,  pointer, save :: numf(:)
        integer,  pointer, save :: numft(:)
        integer,  pointer, save :: numft2(:)
        integer,  pointer, save :: numhij(:)
        integer,  pointer, save :: indon(:)
        integer,  pointer, save :: nindv(:)
        integer,  pointer, save :: cttoc(:,:)
        integer,  pointer, save :: fttof(:,:)
        integer,  pointer, save :: listct(:,:)
        integer,  pointer, save :: listf(:,:)
        integer,  pointer, save :: listft(:,:)
        integer,  pointer, save :: listft2(:,:)
        integer,  pointer, save :: listhij(:,:)

        real(dp), pointer, save :: f(:,:)
        real(dp), pointer, save :: fs(:,:)

      end module on_core
