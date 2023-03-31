! *** Module: nao2gto_debug ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Utilities to debug the NAO2GTO implementation
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2018 Created [Yann Pouillon]
!!
! *****************************************************************************
module nao2gto_debug

  implicit none

  private

  public :: dump_hfx

contains

  ! ***************************************************************************
  !> \brief Dumps the contents of the NAO2GTO data structures
  !!
  !! \par History
  !!      - 10.2018 Created [Yann Pouillon]
  !!
  !! \param[in,out] libint_data: initialized Libint data structure
  !! \param[in] hfx_optdata: initialized Hartree-Fock options data structure
  !! \param[in] hfx_sysdata: initialized Hartree-Fock system data structure
  ! ***************************************************************************
  subroutine dump_hfx(libint_data, hfx_optdata, hfx_sysdata)

    use precision,       only: dp
    use parallel,        only: IOnode, Node, Nodes

    use nao2gto_libint, only: Libint_t
    use nao2gto_types , only: hfx_options_type, hfx_system_type

    implicit none

    ! Arguments
    type(Libint_t), intent(inout)      :: libint_data
    type(hfx_options_type), intent(in) :: hfx_optdata
    type(hfx_system_type), intent(in)  :: hfx_sysdata

    ! -------------------------------------------------------------------------

    ! Dump relevant parameters
    if ( IOnode ) then
      write(*, '(/,a,/)') &
&       "dump_hfx: input parameters begin --------------------------------------"
      write(*,'(a/a/a)') "# *** YAML START ***", "---", "dump_hfx:"
      write(*,'(2x,"nspin: ",I10)') hfx_sysdata%nspin
      write(*,'(2x,"norb: ",I10)') hfx_sysdata%norb
      write(*,'(2x,"iaorb_dims: [",I10,"]")') shape(hfx_sysdata%iaorb)
      write(*,'(2x,"iphorb_dims: [",I10,"]")') shape(hfx_sysdata%iphorb)
      write(*,'(2x,"nuo: ",I10)') hfx_sysdata%nuo
      write(*,'(2x,"nuotot: ",I10)') hfx_sysdata%nuotot
      write(*,'(2x,"nua: ",I10)') hfx_sysdata%nua
      write(*,'(2x,"na: ",I10)') hfx_sysdata%na
      write(*,'(2x,"isa_dims: [",I10,"]")') shape(hfx_sysdata%isa)
      write(*,'(2x,"xa_dims: [",I10,",",1X,I10,"]")') shape(hfx_sysdata%xa)
      write(*,'(2x,"indxua_dims: [",1X,I10,"]")') shape(hfx_sysdata%indxua)
      write(*,'(2x,"cell: [",E10.3,8(",",1X,E10.3),"]")') hfx_sysdata%cell(:,:)
      write(*,'(2x,"nsc: [",I10,2(",",1X,I10),"]")') hfx_sysdata%nsc(:)
      write(*,'(2x,"maxnh: ",I10)') hfx_sysdata%maxnh
      write(*,'(2x,"numh_dims: [",I10,"]")') shape(hfx_sysdata%numh)
      write(*,'(2x,"listhptr_dims: [",I10,"]")') shape(hfx_sysdata%listhptr)
      write(*,'(2x,"listh_dims: [",I10,"]")') shape(hfx_sysdata%listh)
      write(*,'(a/a)') "...", "# *** YAML STOP ***"
      write(*, '(/,a,/)') &
&       "dump_hfx: input parameters end ----------------------------------------"
    endif

  end subroutine dump_hfx

end module nao2gto_debug
