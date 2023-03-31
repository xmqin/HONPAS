! *** Module: nao2gto_libint ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Fortran/C interfaces to the Libint and Libderiv libraries
!!
!! \author Manuel Guidon
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2007 Created [Manuel Guidon]
!!      - 10.2009 Refactored [Manuel Guidon]
!!      - 01.2016 Imported to SIESTA and edited [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_libint

  use, intrinsic :: iso_c_binding

  implicit none

  public

  ! ***************************************************************************
  ! *** Library: libint                                                     ***
  ! ***************************************************************************

  !> Maximum angular momentum considered by Libint
  !!
  !! \warning libint_max_am is equal to the value of Libint's configure script
  !!          plus one.
  integer, parameter :: libint_max_am = 5

  !> Dimension length of the 4D array of function pointers at the core of
  !! Libint (see \ref build_eri)
  integer, parameter :: build_eri_size = libint_max_am - 1

  !> Size of the VRR internal array of Libint
  integer, parameter :: libint_vrr_classes_size = 2*(libint_max_am - 1) + 1

  !> Size of the F internal array storing coefficients for Libint
  !! (see \ref prim_data)
  integer, parameter :: prim_data_f_size = 4*(libint_max_am - 1) + 1

  !> \brief Libint data structure
  type, bind(c) :: Libint_t
    type(c_ptr)    :: int_stack
    type(c_ptr)    :: PrimQuartet
    real(c_double) :: AB(3)
    real(c_double) :: CD(3)
    type(c_ptr)    :: &
      vrr_classes(libint_vrr_classes_size,libint_vrr_classes_size)
    type(c_ptr)    :: vrr_stack
  end type Libint_t

  !> \brief Internal coefficients for Libint and Libderiv
  type, bind(c) :: prim_data
    real(c_double) :: F(prim_data_f_size)
    real(c_double) :: U(3,6)
    real(c_double) :: twozeta_a
    real(c_double) :: twozeta_b
    real(c_double) :: twozeta_c
    real(c_double) :: twozeta_d
    real(c_double) :: oo2z
    real(c_double) :: oo2n
    real(c_double) :: oo2zn
    real(c_double) :: poz
    real(c_double) :: pon
    real(c_double) :: oo2p
    real(c_double) :: ss_r12_ss
  end type prim_data

  !> Libint is built around a matrix of function pointers
  type(c_funptr), dimension( &
&   0:build_eri_size, 0:build_eri_size, &
&   0:build_eri_size, 0:build_eri_size), &
&   bind(c) :: build_eri

  !> \brief Interfaces to the C routines of Libint
  interface

    subroutine init_libint_base() bind(c)
      import
    end subroutine init_libint_base

    function init_libint(libint_data, max_am, max_num_prim_comb) bind(c)
      import
      integer(kind=c_int)        :: init_libint
      type(Libint_t)             :: libint_data
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
    end function init_libint

    subroutine free_libint(libint_data) bind(c)
      import
      type(Libint_t) :: libint_data
    end subroutine free_libint

    function libint_storage_required(max_am, max_num_prim_comb) bind(c)
      import
      integer(kind=c_int)        :: libint_storage_required
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
    end function libint_storage_required

  end interface

  ! ***************************************************************************
  ! *** Library: libderiv                                                   ***
  ! ***************************************************************************

  !> Maximum angular momentum considered by Libderiv
  !!
  !! \warning libderiv_max_am1 is equal to the value of Libint's configure
  !!          script plus one.
  integer, parameter :: libderiv_max_am1 = 4
  integer, parameter :: libderiv_max_am12 = 3

  !> Dimension length of the 4D arrays of function pointers at the core of
  !! Libderiv (see \ref build_deriv1_eri)
  integer, parameter :: build_deriv1_eri_size = libderiv_max_am1 - 1
  integer, parameter :: build_deriv12_eri_size = libderiv_max_am12 - 1

  !> Size of the DVRR internal array of Libderiv
  integer, parameter :: libint_dvrr_classes_size = 2*(libderiv_max_am1 - 1) + 1

  !> \brief Libderiv data structure
  type, bind(c) :: Libderiv_t
    type(c_ptr)    :: int_stack
    type(c_ptr)    :: PrimQuartet
    type(c_ptr)    :: zero_stack
    type(c_ptr)    :: ABCD(156)
    real(c_double) :: AB(3)
    real(c_double) :: CD(3)
    type(c_ptr)    :: &
&     deriv_classes(12,libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: &
&     deriv2_classes(144,libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: &
&     dvrr_classes(libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: dvrr_stack
  end type Libderiv_t

  !> Libderiv(deriv1) is built around a matrix of function pointers
  type(c_funptr), dimension( &
&   0:build_deriv1_eri_size, 0:build_deriv1_eri_size, &
    0:build_deriv1_eri_size, 0:build_deriv1_eri_size), &
&   bind(c) :: build_deriv1_eri

  !> Libderiv(deriv12) is built around a matrix of function pointers
  type(c_funptr), dimension( &
&   0:build_deriv12_eri_size, 0:build_deriv12_eri_size, &
    0:build_deriv12_eri_size, 0:build_deriv12_eri_size), &
&   bind(c) :: build_deriv12_eri

  !> Interfaces to the C routines of Libderiv
  interface

    subroutine init_libderiv_base() bind(c)
      import
    end subroutine init_libderiv_base

    function init_libderiv1(libderiv_data, max_am, max_num_prim_quartets, &
               max_cart_class_size) bind(c)
      import
      integer(kind=c_int)        :: init_libderiv1
      type(Libderiv_t)           :: libderiv_data
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_quartets
      integer(kind=c_int), value :: max_cart_class_size
    end function init_libderiv1

    function init_libderiv12(libderiv_data, max_am, max_num_prim_quartets, &
               max_cart_class_size) bind(c)
      import
      integer(kind=c_int)        :: init_libderiv12
      type(Libderiv_t)           :: libderiv_data
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_quartets
      integer(kind=c_int), value :: max_cart_class_size
    end function init_libderiv12

    subroutine free_libderiv(libderiv_data) bind(c)
      import
      type(Libderiv_t) :: libderiv_data
    end subroutine free_libderiv

    function libderiv1_storage_required(max_am, max_num_prim_quartets, &
               max_cart_class_size) bind(c)
      import
      integer(kind=c_int)        :: libderiv1_storage_required
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_quartets
      integer(kind=c_int), value :: max_cart_class_size
    end function libderiv1_storage_required

    function libderiv12_storage_required(max_am, max_num_prim_quartets, &
               max_cart_class_size) bind(c)
      import
      integer(kind=c_int)        :: libderiv12_storage_required
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_quartets
      integer(kind=c_int), value :: max_cart_class_size
    end function libderiv12_storage_required

  end interface

end module nao2gto_libint
