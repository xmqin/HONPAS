! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_matel_registry
  !
  ! This module provides a general interface to the
  ! functions needed by (new_)matel
  ! 
  ! Sketch of usage (see examples in the code):
  !
  ! Each type of function is "registered", and a global
  ! index is obtained. The user is responsible for the correct
  ! bookeeping of the indexes (see routine register_rfs for an
  ! example appropriate for the orbitals, KB projectors, and Vna
  ! in Siesta, and routine "overlap").
  !
  ! Once registered, a function can be evaluated, its cutoff and 
  ! angular momentum obtained, etc.  So far only the functions
  ! needed by (new_)matel are implemented: rcut, lcut, evaluate (former phiatm).
  !
  ! The registry contains some meta-data for each function (its
  ! kind as a string, and the original indexes handled by the
  ! user at registration time).
  !
  use precision, only: dp
  use sys, only: die

  use radial, only: rad_func
  use trialorbitalclass, only: trialorbital

  implicit none

  private

  type id_t
     character(len=10) :: func_type = "null"
     integer           :: n_indexes = 0
     integer           :: indexes(4) = (/-1,-1,-1,-1/)
  end type id_t

  ! Two kinds of functions are supported:

  ! A) Finite-range functions with pure (l,m) angular part.
  ! Each function is represented by the radfunc record (for its radial data)
  ! and a meta-data section, which includes l, m, and id.
  !
  type ext_radfunc_t
     type(rad_func), pointer :: func => null()
     integer        :: l
     integer        :: m
     type(id_t)     :: id
  end type ext_radfunc_t

  !
  ! B) Trial "Wannier" orbitals, which are hybrids with a maximum l.
  ! They are defined in module trialorbitalclass, created by Richard Korytar.
  !

  ! This umbrella type can hold either function.
  ! The lcut and rcut fields are used to avoid the overhead of
  ! dispatch for queries

  type registered_function_t
     type(ext_radfunc_t), pointer :: rf => null()
     type(trialorbital), pointer  :: tf => null()
     real(dp)        :: rcut = 0.0_dp
     integer         :: lcut = -huge(1)
  end type registered_function_t
  !
  ! This is a dynamic array, allocated and extended as needed.
  type(registered_function_t), dimension(:), allocatable, save :: matel_pool

  integer, parameter :: initial_size = 100  ! initial size of pool
  integer, parameter :: delta_size = 50     ! chunk size for extensions
  integer, private   :: nmax_funcs = 0  ! current maximum size of the pool
  integer, private   :: nfuncs = 0      ! number of functions in pool

  public :: register_in_rf_pool, register_in_tf_pool
  public :: evaluate, rcut, lcut
  public :: evaluate_x, evaluate_y, evaluate_z
  public :: show_pool

CONTAINS

  function valid(gindex) result (ok)
    integer, intent(in)    :: gindex
    logical                :: ok

    ok = (gindex > 0 .AND. gindex <= nfuncs)
  end function valid

  !> Extends the size of the pool if needed
  subroutine check_size_of_pool(nfuncs)
    integer, intent(in) :: nfuncs

    type(registered_function_t), dimension(:), allocatable :: tmp_pool
    
    if (nfuncs > nmax_funcs) then
       if (.not. allocated(matel_pool)) then
          allocate(matel_pool(initial_size))
          nmax_funcs = initial_size
       else
          ! copy data to tmp array
          allocate(tmp_pool(nmax_funcs))

          ! Each element (derived type) involves just two scalars
          ! and two pointers, so the copy overhead is small.
          ! Note that we do not use the F2003 "move_alloc" idiom.
          
          tmp_pool(1:nmax_funcs) = matel_pool(1:nmax_funcs)
          deallocate(matel_pool)
          allocate(matel_pool(nmax_funcs + delta_size))
          matel_pool(1:nmax_funcs) = tmp_pool(1:nmax_funcs)
          nmax_funcs = nmax_funcs + delta_size

          deallocate(tmp_pool)
          
       endif
    endif
  end subroutine check_size_of_pool
               
  !
  !   This is the main entry to the registry for simple radial functions
  !
  subroutine register_in_rf_pool(func,l,m,func_type,indexes,gindex)
    type(rad_func), pointer :: func             ! function data
    integer, intent(in)     :: l, m           
    character(len=*), intent(in) :: func_type   ! mnemonic kind
    integer, intent(in)          :: indexes(:)  ! legacy indexes
    integer, intent(out)    :: gindex           ! global index

    integer :: n_indexes 
    type(ext_radfunc_t), pointer :: rf

    n_indexes = size(indexes)

    nfuncs = nfuncs + 1
    call check_size_of_pool(nfuncs)

    gindex = nfuncs

    allocate(matel_pool(gindex)%rf)
    rf => matel_pool(gindex)%rf
    rf%func => func
    rf%l = l
    rf%m = m

    rf%id%func_type = func_type
    rf%id%n_indexes = n_indexes
    rf%id%indexes(1:n_indexes) = indexes(:)

    ! To speed up queries
    matel_pool(gindex)%lcut = l
    matel_pool(gindex)%rcut = func%cutoff

  end subroutine register_in_rf_pool

  !
  !   This is the main entry to the registry for "Wannier" trial orbs
  !
  subroutine register_in_tf_pool(tf,gindex)
    use trialorbitalclass, only: print_trialorb
    use trialorbitalclass, only: gettrialrcut, gettriallmax

    type(trialorbital), intent(in) :: tf
    integer, intent(out)    :: gindex           ! global index

    nfuncs = nfuncs + 1
    call check_size_of_pool(nfuncs)

    gindex = nfuncs

    allocate(matel_pool(gindex)%tf)
    matel_pool(gindex)%tf = tf

    ! To speed up queries
    matel_pool(gindex)%lcut = gettriallmax(tf)
    matel_pool(gindex)%rcut = gettrialrcut(tf)

  end subroutine register_in_tf_pool

!--------------------------------------------------------------
  subroutine show_pool()
    use trialorbitalclass, only: print_trialorb

    integer :: gindex

    do gindex = 1, nfuncs
       if (associated(matel_pool(gindex)%rf)) then
          call print_rf(matel_pool(gindex)%rf)
       else if (associated(matel_pool(gindex)%tf)) then
          call print_trialorb(matel_pool(gindex)%tf)
       endif
    enddo
  end subroutine show_pool

!--------------------------------------------------------------

  function rcut(gindex) result(cutoff)
    integer, intent(in)    :: gindex
    real(dp)               :: cutoff

    if (valid(gindex)) then
       cutoff = matel_pool(gindex)%rcut
!       if (associated(matel_pool(gindex)%rf)) then
!          cutoff = matel_pool(gindex)%rf%func%cutoff
!       else if (associated(matel_pool(gindex)%tf)) then
!          cutoff = gettrialrcut(matel_pool(gindex)%tf)
!       endif
    else
       call die("Invalid gindex")
    endif
  end function rcut

!--------------------------------------------------------------
  function lcut(gindex) result(l)
    integer, intent(in)    :: gindex
    integer                :: l

    if (valid(gindex)) then
       l = matel_pool(gindex)%lcut
!       if (associated(matel_pool(gindex)%rf)) then
!          l = matel_pool(gindex)%rf%l
!       else if (associated(matel_pool(gindex)%tf)) then
!          l = gettriallmax(matel_pool(gindex)%tf)
!       endif
    else
       call die("Invalid gindex")
    endif
  end function lcut

!--------------------------------------------------------------
  subroutine evaluate(gindex,r,f,grad)
  use trialorbitalclass, only: gettrialwavefunction

    integer, intent(in)    :: gindex
    real(dp), intent(in)   :: r(3)
    real(dp), intent(out)  :: f
    real(dp), intent(out)  :: grad(3)

    if (valid(gindex)) then
       if (associated(matel_pool(gindex)%rf)) then
          call evaluate_ext_radfunc(matel_pool(gindex)%rf,r,f,grad)
       else if (associated(matel_pool(gindex)%tf)) then
          f = gettrialwavefunction(matel_pool(gindex)%tf,r)
          grad = 0.0_dp   ! do not know how to get this...
       endif
    else
       call die("Invalid gindex")
    endif
  end subroutine evaluate

  SUBROUTINE evaluate_x(gindex,r,xphi,grxphi)
    !      Calculates x*phiatm and its gradient

    integer, intent(in)   :: gindex
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: xphi, grxphi(3)

    real(dp) phi, grphi(3), x

    call evaluate(gindex,r,phi,grphi)
    x = r(1)
    xphi = x * phi
    grxphi(1) = x * grphi(1) + phi
    grxphi(2) = x * grphi(2)
    grxphi(3) = x * grphi(3)
  END SUBROUTINE evaluate_x

  SUBROUTINE evaluate_y(gindex,r,yphi,gryphi)
    !     Calculates y*phiatm and its gradient

    integer, intent(in)   :: gindex
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: yphi, gryphi(3)

    real(dp) phi, grphi(3), y

    call evaluate(gindex,r,phi,grphi)
    y = r(2)
    yphi = y * phi
    gryphi(1) = y * grphi(1)
    gryphi(2) = y * grphi(2) + phi
    gryphi(3) = y * grphi(3)
  END SUBROUTINE evaluate_y

  SUBROUTINE evaluate_z(gindex,r,zphi,grzphi)
    !     Calculates z*phiatm and its gradient

    integer, intent(in)   :: gindex
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: zphi, grzphi(3)

    real(dp) phi, grphi(3), z
    call evaluate(gindex,r,phi,grphi)
    z = r(3)
    zphi = z * phi
    grzphi(1) = z * grphi(1)
    grzphi(2) = z * grphi(2)
    grzphi(3) = z * grphi(3) + phi
  END SUBROUTINE evaluate_z

  ! helper function
  ! it is here because of the "vna" singularity. This should be
  ! removed.
  ! Also, there should be a module "radfuncClass" containing the
  ! extended type and the evaluator. The registration could be
  ! done on the basis of the extended type, or in the legacy form
  ! based on the radial part plus metadata.

  subroutine evaluate_ext_radfunc(rf,r,phi,grphi)

    use radial, only: rad_func, rad_get
    use spher_harm, only:  rlylm

    type(ext_radfunc_t), pointer :: rf
    real(dp), intent(in)   :: r(3)
    real(dp), intent(out)  :: phi
    real(dp), intent(out)  :: grphi(3)

    type(rad_func), pointer :: func

    integer, parameter   :: max_l = 5
    integer, parameter   :: max_ilm = (max_l+1)*(max_l+1)

    real(dp) rmod, phir, dphidr
    real(dp) rly(max_ilm), grly(3,max_ilm)
    integer i, l, m, ilm

    real(dp), parameter     :: tiny = 1.e-20_dp

    func => rf%func
    l = rf%l
    m = rf%m

    ! the addition of "tiny" avoids division by 0
    !
    rmod = sqrt(sum(r*r)) + tiny

    if(rmod > func%cutoff) then

       phi = 0.0_dp
       grphi(1:3) = 0.0_dp

    else

       call rad_get(func,rmod,phir,dphidr)

       if (rf%id%func_type(1:3) == "vna") then
          phi=phir
          grphi(1:3)=dphidr*r(1:3)/rmod
       else
          ilm = l*l + l + m + 1
          call rlylm( l, r, rly, grly )
          phi = phir * rly(ilm)
          do i = 1,3
             grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm)
          enddo
       endif
    endif
  end subroutine evaluate_ext_radfunc

  subroutine print_rf(rf)
    type(ext_radfunc_t), intent(in) :: rf
    print "(a,i2,i2,f8.4)", trim(rf%id%func_type), rf%l, rf%m, rf%func%cutoff
  end subroutine print_rf

end module m_matel_registry
  
