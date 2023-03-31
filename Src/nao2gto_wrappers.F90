! *** Module: nao2gto_wrappers ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Low-level wrappers for the communication between SIESTA and Libint
!!
!! This module takes care of managing the initialization, tracking, and
!! destruction, of NAO2GTO data.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_wrappers

  implicit none

  private

  public :: &
&   nao2gto_libint_init, &
&   nao2gto_libint_free, &
&   nao2gto_libint_dump, &
&   get_eris, &
&   nao2gto_libderiv_init, &
&   nao2gto_libderiv_free, &
&   nao2gto_libderiv_dump, &
&   get_derivs, &
&   nao2gto_system_init, &
&   nao2gto_system_free, &
&   nao2gto_system_update_cell

contains

  ! ***************************************************************************
  ! *** Library: libint                                                     ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Initializes the NAO2GTO data structures associated to Libint
  !!
  !! \author Yann Pouillon
  !!
  !! \param[out] libint_data: Libint data structure
  !! \param[in] l_max: maximum angular momentum to consider
  ! ***************************************************************************
  subroutine nao2gto_libint_init(libint_data, l_max)

    use, intrinsic :: iso_c_binding, only: c_int

    use nao2gto_libint

    implicit none

    ! Arguments
    integer           , intent(in)  :: l_max
    type(Libint_t)    , intent(inout) :: libint_data

    ! Local variables
    integer(kind=c_int) :: lib_storage, max_am, max_prim

    ! -------------------------------------------------------------------------

    ! Actual Libint initialisation
    max_am = l_max
    max_prim = 1
    call init_libint_base()
    lib_storage = init_libint(libint_data, max_am, max_prim)
    if ( lib_storage < 0 ) then
      call die("Maximum angular momentum too high for Libint")
    endif

  end subroutine nao2gto_libint_init

  ! ***************************************************************************
  !> \brief Terminates the NAO2GTO data structures associated to Libint
  !!
  !! \author Yann Pouillon
  !!
  !! \note: Calling \ref free_libint causes the program to abort because of
  !!        a double free, probably because the routine is automatically
  !!        called by the C++ garbage collector (to be checked).
  ! ***************************************************************************
  subroutine nao2gto_libint_free(libint_data)

    use nao2gto_data , only: subshell
    use nao2gto_index, only: nao2gto_index_free
    use nao2gto_libint
    use nao2gto_utils, only: deallocate_md_ftable

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data

    ! -------------------------------------------------------------------------

    ! Get rid of the Libint data structure
    !call free_libint(libint_data)

  end subroutine nao2gto_libint_free

  ! ***************************************************************************
  !> \brief Displays the current status of NAO2GTO data structures associated
  !!        to Libint
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] libint_data: Libint data structure
  ! ***************************************************************************
  subroutine nao2gto_libint_dump(libint_data)

    use, intrinsic :: iso_c_binding

    use nao2gto_libint

    implicit none

    ! Arguments
    type(Libint_t), intent(in) :: libint_data

    ! -------------------------------------------------------------------------

!!   For debugging
!    write(*, '(/,a,/)') "nao2gto_libint_status: Libint data structure ----------------------------------"
!    write(*,'(a,/,"---",/,a,":")') "# *** YAML START ***", "libint_data"
!    write(*,'(2x,a,":",1x,a)') "int_stack", status_string(libint_data%int_stack)
!    write(*,'(2x,a,":",1x,a)') "primquartet", status_string(libint_data%primquartet)
!    write(*,'(2x,a,":",1x,a)') "vrr_classes", status_string(libint_data%vrr_classes(1,1))
!    write(*,'(2x,a,":",1x,a)') "vrr_stack", status_string(libint_data%vrr_stack)
!    write(*,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "ab", dble(libint_data%ab(:))
!    write(*,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "cd", dble(libint_data%cd(:))
!    write(*,'(a)') "# *** YAML STOP ***"
!    write(*, '(/,a,/)') "nao2gto_libint_status: END ----------------------------------------------------"
!!   End debugging

  contains

    function status_string(var)

      use, intrinsic :: iso_c_binding, only: c_associated, c_ptr

      implicit none

      type(c_ptr), intent(in) :: var
      character(len=4) :: status_string

      if ( c_associated(var) ) then
        status_string = "addr"
      else
        status_string = "null"
      endif

    end function status_string

  end subroutine nao2gto_libint_dump

  ! ***************************************************************************
  !> \brief ...
  !!
  !! \param n_d ...
  !! \param n_c ...
  !! \param n_b ...
  !! \param n_a ...
  !! \param libint_data ...
  !! \param prim ...
  !! \param p_work ...
  !! \param[in] a_mysize: array shape to pass to c_f_pointer (must be a vector)
  ! ***************************************************************************
  subroutine get_eris(n_d, n_c, n_b, n_a, libint_data, prim, p_work, a_mysize)

    use, intrinsic :: iso_c_binding

    use precision, only: dp
    use nao2gto_libint

    implicit none

    ! Arguments
    integer, intent(in) :: n_d, n_c, n_b, n_a
    type(Libint_t), intent(inout) :: libint_data
    type(prim_data), target, intent(inout) :: prim
    real(dp), dimension(:), pointer, intent(out) :: p_work
    integer, intent(in) :: a_mysize(1)

    ! Local interfaces
    interface
      function build(libint_data, max_num_prim_comb) bind(c)
        import
        type(c_ptr)                :: build
        type(Libint_t)             :: libint_data
        integer(kind=c_int), value :: max_num_prim_comb
      end function build
    end interface

    ! Local variables
    integer :: i
    procedure(build), pointer :: pbuild
    type(c_ptr) :: pc_result = C_NULL_PTR

    ! -------------------------------------------------------------------------

    p_work => null()
    libint_data%primquartet = c_loc(prim)
    call c_f_procpointer(build_eri(n_d, n_c, n_b, n_a), pbuild)
    pc_result = pbuild(libint_data, 1)
    call c_f_pointer(pc_result, p_work, a_mysize)

  end subroutine get_eris

  ! ***************************************************************************
  ! *** Library: libderiv                                                   ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Initializes a Libderiv data structure
  !!
  !! \param[in,out] libderiv_data: Libderiv data structure
  !! \param[in] max_am: maximum angular momentum to consider
  ! ***************************************************************************
  subroutine nao2gto_libderiv_init(libderiv_data, max_am)

    use sys, only: die
    use nao2gto_data, only: nco
    use nao2gto_libint

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: libderiv_data
    integer, intent(in) :: max_am

    ! Local variables
    integer(kind=c_int) :: deriv_storage, max_am_local, max_classes, max_prim

    ! -------------------------------------------------------------------------

    max_am_local = max_am
    max_prim = 1
    max_classes = nco(max_am)**4

    call init_libderiv_base()

    deriv_storage = init_libderiv1(libderiv_data, max_am_local, max_prim, max_classes)
    if ( deriv_storage < 0 ) then
      call die("The angular momentum needed exceeds the value configured in libint")
    endif

  end subroutine nao2gto_libderiv_init

  ! ***************************************************************************
  !> \brief Properly terminates a Libderiv data structure
  !!
  !! \param[in,out] deriv: data structure to terminate
  ! ***************************************************************************
  subroutine nao2gto_libderiv_free(libderiv_data)

    use nao2gto_libint

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: libderiv_data

    ! -------------------------------------------------------------------------

    !call free_libderiv(libderiv_data)

  end subroutine nao2gto_libderiv_free

  ! ***************************************************************************
  !> \brief Displays the current status of NAO2GTO data structures associated
  !!        to Libderiv
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] libderiv_data: Libderiv data structure
  ! ***************************************************************************
  subroutine nao2gto_libderiv_dump(libderiv_data)

    use, intrinsic :: iso_c_binding

    use nao2gto_libint

    implicit none

    ! Arguments
    type(Libderiv_t), intent(in) :: libderiv_data

    ! -------------------------------------------------------------------------

!!   For debugging
!    write(*, '(/,a,/)') "nao2gto_libderiv_status: Libderiv data structure --------------------------------"
!    write(*,'(a,/,"---",/,a,":")') "# *** YAML START ***", "libderiv_data"
!    write(*,'(2x,a,":",1x,a)') "int_stack", status_string(libderiv_data%int_stack)
!    write(*,'(2x,a,":",1x,a)') "primquartet", status_string(libderiv_data%primquartet)
!    write(*,'(2x,a,":",1x,a)') "abcd(1)", status_string(libderiv_data%abcd(1))
!    write(*,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "ab", dble(libderiv_data%ab(:))
!    write(*,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "cd", dble(libderiv_data%cd(:))
!    write(*,'(2x,a,":",1x,a)') "deriv_classes", status_string(libderiv_data%deriv_classes(1,1,1))
!    write(*,'(2x,a,":",1x,a)') "deriv2_classes", status_string(libderiv_data%deriv2_classes(1,1,1))
!    write(*,'(2x,a,":",1x,a)') "dvrr_classes", status_string(libderiv_data%dvrr_classes(1,1))
!    write(*,'(2x,a,":",1x,a)') "dvrr_stack", status_string(libderiv_data%dvrr_stack)
!    write(*,'(a)') "# *** YAML STOP ***"
!    write(*, '(/,a,/)') "nao2gto_libderiv_status: END --------------------------------------------------"
!!   End debugging

  contains

    function status_string(var)

      use, intrinsic :: iso_c_binding, only: c_associated, c_ptr

      implicit none

      type(c_ptr), intent(in) :: var
      character(len=4) :: status_string

      if ( c_associated(var) ) then
        status_string = "addr"
      else
        status_string = "null"
      endif

    end function status_string

  end subroutine nao2gto_libderiv_dump

  ! ***************************************************************************
  !> \brief ...
  !!
  !! \param n_d ...
  !! \param n_c ...
  !! \param n_b ...
  !! \param n_a ...
  !! \param libderiv_data ...
  !! \param prim ...
  !! \param work_forces ...
  !! \param a_mysize: array shape to pass to c_f_pointer (must be a vector)
  ! ***************************************************************************
  subroutine get_derivs(n_d, n_c, n_b, n_a, libderiv_data, prim, &
&              work_forces, a_mysize)

    use precision, only: dp
    use nao2gto_data, only: nco
    use nao2gto_libint

    implicit none

    ! Arguments
    integer, intent(in) :: n_d, n_c, n_b, n_a
    type(Libderiv_t), intent(inout) :: libderiv_data
    type(prim_data), target, intent(inout) :: prim
    real(dp), dimension(nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d),12), &
&     intent(out) :: work_forces
    integer, intent(in) :: a_mysize(1)

    ! Local interfaces
    interface
      subroutine build_deriv1(libderiv_data, max_num_prim_comb) bind(c)
        import
        type(Libderiv_t)           :: libderiv_data
        integer(kind=c_int), value :: max_num_prim_comb
      end subroutine build_deriv1
    end interface

    ! Local variables
    integer :: i, k
    type(c_ptr) :: pc_result
    procedure(build_deriv1), pointer :: pbuild_deriv1 => null()
    real(c_double), dimension(:), pointer :: tmp_data

    ! -------------------------------------------------------------------------

    i = 1
    libderiv_data%primquartet = c_loc(prim)
    call c_f_procpointer(build_deriv1_eri(n_d, n_c, n_b, n_a), pbuild_deriv1)
    call pbuild_deriv1(libderiv_data, i)

    do k=1,12
      tmp_data => null()
      if ( (k == 4) .or. (k == 5) .or. (k == 6) ) cycle
      pc_result = libderiv_data%abcd(k)
      call c_f_pointer(pc_result, tmp_data, a_mysize)
      do i=1,a_mysize(1)
        work_forces(i,k) = tmp_data(i)
      enddo
    end do

  end subroutine get_derivs

  ! ***************************************************************************
  ! *** SIESTA-specific routines                                            ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Initializes the NAO2GTO data structures associated to SIESTA
  !!
  !! \author Yann Pouillon
  !!
  !! \param[out] hfx_sys: information about the Hartree-Fock system data
  !! \param[in] maxnh: maximum number of non-zero H matric elements
  !! \param[in] na: number of atoms in the supercell
  !! \param[in] norb: number of orbitals in the supercell
  !! \param[in] nspin: number of spìn degrees of freedom
  !! \param[in] nua: number of atoms in the unit cell
  !! \param[in] nuo: number of orbitals local to node
  !! \param[in] nuotot: number of orbitals in the unit cell
  !! \param[in] iaorb: atomic index of each orbital
  !! \param[in] indxua: ...
  !! \param[in] iphorb: orbital index of each orbital in its atom
  !! \param[in] isa: ...
  !! \param[in] listh: ...
  !! \param[in] listhptr: ...
  !! \param[in] nsc: number of supercells in each direction (diagonal of mscell)
  !! \param[in] numh: ...
  !! \param[in] ucell: unit cell in real space
  !! \param[in] xa: positions of the atoms in the unit cell
  ! ***************************************************************************
  subroutine nao2gto_system_init(hfx_sys, maxnh, na, norb, nspin, nua, &
&                nuo, nuotot, iaorb, indxua, iphorb, isa, &
&                listh, listhptr, nsc, numh, ucell, xa)

    use alloc         , only: re_alloc
    use atmfuncs      , only: lofio, mofio
    use precision, only: dp
    use nao2gto_data  , only: subshell    
                                          ! It has the dimension of the
                                          !  number of orbitals in the supercell
                                          !  and points to the index of the
                                          !  "different" angular momentum shell.
                                          !  For instance, bulk Si
                                          !  (two atom per unit cell)
                                          !  and a DZP basis, we have
                                          !  ten different shells in the 
                                          !  unit cell, 
                                          !  and (ten*number of repetitions on 
                                          !  the unit cell in the supercell)
                                          !  different shells in the supercell.
                                          !  The value of the subshell for each
                                          !  orbital would be
                                          !  Si #1, 3s, first  zeta, subshell=1
                                          !  Si #1, 3s, second zeta, subshell=2
                                          !  Si #1, 3p, first  zeta, subshell=3
                                          !  Si #1, 3p, second zeta, subshell=4
                                          !  Si #1, 3d, first  zeta, subshell=5
                                          !  Si #2, 3s, first  zeta, subshell=6
                                          !  Si #2, 3s, second zeta, subshell=7
                                          !  Si #2, 3p, first  zeta, subshell=8
                                          !  Si #2, 3p, second zeta, subshell=9
                                          !  Si #2, 3d, first  zeta, subshell=10
                                          !  Si #3, 3s, first  zeta, subshell=11
                                          !  ...
    use nao2gto_index , only: nao2gto_index_init
    use nao2gto_types , only: hfx_system_type, l_max
    use nao2gto_utils , only: init_md_ftable
!   For debugging
    use parallel      , only: Node, Nodes
!   End debugging

    implicit none

    ! Arguments
    integer, target, intent(in)  :: maxnh, na, norb, nspin, nua, nuo, nuotot
    integer, target, intent(in)  :: iaorb(norb), indxua(na), iphorb(norb), &
&     isa(na), listh(maxnh), listhptr(nuo), nsc(3), numh(nuo)
    real(dp), target, intent(in) :: ucell(3,3), xa(3,na)
    type(hfx_system_type), intent(out) :: hfx_sys

    ! Local variables
    integer :: i, io, ioa, is, l, m, ncells
    integer :: ic, jc

    ! -------------------------------------------------------------------------

    ! Make system information available to Hartree-Fock exchange routines,
    ! avoiding data copy
    hfx_sys%maxnh => maxnh
    hfx_sys%na => na
    hfx_sys%norb => norb
    hfx_sys%nspin => nspin
    hfx_sys%nua => nua
    hfx_sys%nuo => nuo
    hfx_sys%nuotot => nuotot
    hfx_sys%iaorb => iaorb
    hfx_sys%indxua => indxua
    hfx_sys%iphorb => iphorb
    hfx_sys%isa => isa
    hfx_sys%listh => listh
    hfx_sys%listhptr => listhptr
    hfx_sys%nsc => nsc
    hfx_sys%numh => numh
    hfx_sys%xa => xa

    ! Synchronize with the current SIESTA cell
    call nao2gto_system_update_cell(hfx_sys, nsc, ucell)

    ! Gamma function table initialisation
    call init_md_ftable(4*l_max)

    ! Supercell orbital initialisation
    call nao2gto_index_init(hfx_sys%nsc, hfx_sys%nuotot)

!   Orbital shell initialisation
    call re_alloc(subshell, 1, hfx_sys%norb, &
&     name='subshell', routine='nao2gto_libint_init')
    subshell(:) = 0

!   2l+1 orbitals form a shell
    i = 0
    do io=1,hfx_sys%norb
      is = hfx_sys%isa(hfx_sys%iaorb(io))
      ioa = hfx_sys%iphorb(io)
      l = lofio(is,ioa)
      m = mofio(is,ioa)
      if (m .eq. -l ) then
        i = i + 1
      endif
      subshell(io) = i
!!     For debugging
!      write(6,'(a,4i5)') &
! &      'nao2gto_system_init: Node, Nodes, io, subshell(io) = ', &
! &       Node, Nodes, io, subshell(io) 
!!     End debugging
    enddo

  end subroutine nao2gto_system_init

  ! ***************************************************************************
  !> \brief Terminates the NAO2GTO data structures associated to SIESTA
  !!
  !! \author Yann Pouillon
  ! ***************************************************************************
  subroutine nao2gto_system_free(hfx_sys)

    use precision, only: dp
    use nao2gto_data, only: subshell, D2Sindx, eri_prescreen
    use nao2gto_index , only: nao2gto_index_free
    use nao2gto_types, only: hfx_system_type
    use nao2gto_utils , only: deallocate_md_ftable

    implicit none

    ! Arguments
    type(hfx_system_type), intent(inout) :: hfx_sys

    ! -------------------------------------------------------------------------

    ! Get rid of internal tables
    call deallocate_md_ftable()
    call nao2gto_index_free()
    if ( associated(subshell) ) then
      deallocate(subshell)
      subshell => null()
    endif

    if ( associated(D2Sindx) ) then
      deallocate(D2Sindx)
      D2Sindx => null()
    endif

    if ( associated(eri_prescreen) ) then
      deallocate(eri_prescreen)
      eri_prescreen => null()
    endif

    ! Reset cell parameters
    hfx_sys%cell(:,:) = 0.0_dp
    hfx_sys%cell_r(:,:) = 0.0_dp

    ! Make system information unusable
    hfx_sys%maxnh => null()
    hfx_sys%na => null()
    hfx_sys%norb => null()
    hfx_sys%nspin => null()
    hfx_sys%nua => null()
    hfx_sys%nuo => null()
    hfx_sys%nuotot => null()
    hfx_sys%iaorb => null()
    hfx_sys%indxua => null()
    hfx_sys%iphorb => null()
    hfx_sys%isa => null()
    hfx_sys%listh => null()
    hfx_sys%listhptr => null()
    hfx_sys%nsc => null()
    hfx_sys%numh => null()
    hfx_sys%xa => null()

  end subroutine nao2gto_system_free

  ! ***************************************************************************
  !> \brief Updates the supercell used by NAO2GTO data structures
  !!
  !! \author Yann Pouillon
  !!
  !! \param[inout] hfx_sys: information about the Hartree-Fock system data
  !! \param[in] nsc: number of supercells in each direction (diagonal of mscell)
  !! \param[in] ucell: unit cell in real space
  ! ***************************************************************************
  subroutine nao2gto_system_update_cell(hfx_sys, nsc, ucell)

    use precision, only: dp
    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    type(hfx_system_type), intent(inout) :: hfx_sys
    integer, intent(in)  :: nsc(3)
    real(dp), intent(in) :: ucell(3,3)

    ! Local variables
    integer  :: ic, jc

    ! -------------------------------------------------------------------------

    ! Determine supercell both in real space and reciprocal space
    hfx_sys%cell(:,:) = 0.0_dp
    hfx_sys%cell_r(:,:) = 0.0_dp
    do ic=1,3
      do jc=1,3
        hfx_sys%cell(jc,ic) = ucell(jc,ic) * nsc(ic)
      enddo
    enddo
    call reclat(hfx_sys%cell, hfx_sys%cell_r, 0)

  end subroutine nao2gto_system_update_cell

end module nao2gto_wrappers
