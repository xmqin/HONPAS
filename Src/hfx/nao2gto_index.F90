! *** Module: nao2gto_index ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Supercell indexing utilities for NAO2GTO
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_index

  private 

  public :: nao2gto_index_init, nao2gto_index_free, indexsc 

  integer, pointer, save :: indsc1(:) => null(), indsc2(:) => null()

contains

! *****************************************************************************
!> \brief Constructor for the internal indices of the module
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nsc: ...
!! \param[in] nuo: local number of orbitals within the unit cell in this node
! *****************************************************************************
  subroutine nao2gto_index_init(nsc, nuo)

    use alloc, only: re_alloc

    implicit none

    ! Arguments
    integer, intent(in) :: nsc(3)
    integer, intent(in) :: nuo

    ! Local variables
    integer :: i1, i2, i3, ic, j1, j2, j3, jc, kuo, lastio, lastjo, ncells, no

    ! -------------------------------------------------------------------------

    ncells = nsc(1) * nsc(2) * nsc(3)
    no = nuo * ncells

    call re_alloc(indsc1, 1, 8*no, &
&     name="indsc1", routine="nao2gto_index_init")
    call re_alloc(indsc2, 1, no, &
&     name="indsc2", routine="nao2gto_index_init")

    do j3 = 0,2*nsc(3)-1
      do j2 = 0,2*nsc(2)-1
        do j1 = 0,2*nsc(1)-1

          i1 = mod(j1,nsc(1))
          i2 = mod(j2,nsc(2))
          i3 = mod(j3,nsc(3))

          ic = i1 +   nsc(1)*i2 +   nsc(1)*nsc(2)*i3
          jc = j1 + 2*nsc(1)*j2 + 4*nsc(1)*nsc(2)*j3

          lastio = ic * nuo
          lastjo = jc * nuo

          do kuo = 1,nuo
            indsc1(lastjo+kuo) = lastio + kuo
          enddo

          if (i1+nsc(1).eq.j1 .and. i2+nsc(2).eq.j2 .and. i3+nsc(3).eq.j3) then
            do kuo = 1,nuo
              indsc2(lastio+kuo) = lastjo + kuo
            enddo
          endif

        enddo
      enddo
    enddo

  end subroutine nao2gto_index_init

! *****************************************************************************
!> \brief Destructor for the internal indices of the module
!!
!! \par History
!!      - 01.2018 Created for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
  subroutine nao2gto_index_free()

    use alloc, only: de_alloc

    if ( associated(indsc1) ) then
      call de_alloc(indsc1, name='indsc1', routine='nao2gto_index_free')
      nullify(indsc1)
    endif
    if ( associated(indsc2) ) then
      call de_alloc(indsc2, name='indsc2', routine='nao2gto_index_free')
      nullify(indsc2)
    endif

  end subroutine nao2gto_index_free

! *****************************************************************************
!> \brief Returns the index of the specified orbital
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] io: ...
!! \param[in] iuo: local index of the orbital within the unit cell in this node
!! \param[in] jo: ...
! *****************************************************************************
  function indexsc(io, iuo, jo)

    implicit none

    ! Arguments
    integer, intent(in) :: io, iuo, jo

    ! Local variables
    integer :: indexsc

    indexsc = indsc1(indsc2(jo) - indsc2(io) + indsc2(iuo))

  end function indexsc

end module nao2gto_index
