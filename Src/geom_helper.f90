!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.


! Many arrays in SIESTA are "not needed"
! In the sense that one can retrieve, by simple calculations
! the exact same information by other arrays.
! The "iaorb" can be substituted by checks in the lasto array
! The "indxuo" array can be substituted by a call to "x=MOD(i,no_u); if (x==0)x=no_u"
! etc.

! Here we provide a few of such routines for places
! in the code where such arrays just "clutter the namespace"

! It will be advised that these routines are only used
! in "one time calculations", i.e. for heavy duty calculations
! frequent calls to these routines may not be the best way.

! There are however a few routines which are encouraged
! to use through-out, and also in heavy duty calculations.
! These are:
!   - cell_abc which does not do any calculations that would not
!     be performed otherwise.
!   - ucorb, even though there are calculations in this routine
!     the lack of passing, and look-up in the indxuo array
!     seems like a good enough argument to fully refrain from
!     using that array.

module geom_helper

  use intrinsic_missing, only : UCORB  => MODP
  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)

  private

  ! We fetch this from the intrinsic_missing module
  ! It does the same thing
  public :: ucorb
  public :: iaorb
  interface cell_abc
     module procedure cell_abc_ij, cell_abc_c
  end interface cell_abc
  public :: cell_abc
  public :: cell_d, cell_a, cell_b, cell_c
  public :: xa_in_cell

contains

  function iaorb(iorb,lasto) result(ia)
! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: iorb
    integer, intent(in) :: lasto(:) ! actually lasto(0:na_u) !
! *********************
! * OUTPUT variables  *
! *********************
    integer :: ia, iiorb, na

    ! Could just do:
    ! ia = SFIND(lasto,iorb,NEAREST=+1) - 1

    ! retrieve number of atoms
    na = ubound(lasto,dim=1)-1

    ! We need to transfer to the unit-cell
    iiorb = ucorb(iorb,lasto(na+1))

    ! Start searching in the "best" guess of the placement
    ! Our best guess is the number of orbitals on the first atom
    ! 1+int(iorb / orb(1)) == atom (if all atoms have the same orb. count, this is exact)
    ! This is the most qualified quess we can make
    ! It makes no sense to take the mean of the orbital count...!
    ! If iiorb and lasto are both > 0 then this will always
    ! return in the range [1;na]
    ia = max(1,min(int(real(iiorb,dp) / lasto(2)),na))

    ia_loop: do
       if ( iiorb < lasto(ia) ) then
          ! 1. We have overestimated the orbital position
          ia = ia - 1
          cycle ia_loop
       else if ( iiorb > lasto(ia+1) ) then
          ! 2. We have overestimated the orbital position
          ia = ia + 1
          cycle ia_loop
       else if ( iiorb == lasto(ia) ) then
          ! 3. it is on the former atom
          ia = ia - 1
       end if
       ! We have found it!!!
       return
    end do ia_loop

    call die('SOMETHING WENT TERRIBLY WRONG IN IAORB')
  end function iaorb


  pure function cell_abc_ij(recell,xai,xaj,xij) result(nnn)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xai(3), xaj(3), xij(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer  :: nnn(3)
    real(dp) :: xd(3)
    ! The actual length between two atomic centered orbitals:
    xd (:) = xij(:) - (xaj(:)-xai(:))
    nnn(1) = nint(sum(xd*recell(:,1)))
    nnn(2) = nint(sum(xd*recell(:,2)))
    nnn(3) = nint(sum(xd*recell(:,3)))

  end function cell_abc_ij

  pure function cell_abc_c(recell,xR) result(nnn)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xR(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer :: nnn(3)
    ! The actual length between two atomic centered orbitals:
    nnn(1) = nint(sum(xR*recell(:,1)))
    nnn(2) = nint(sum(xR*recell(:,2)))
    nnn(3) = nint(sum(xR*recell(:,3)))

  end function cell_abc_c


  pure function cell_a(recell,xai,xaj,xij) result(n)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xai(3), xaj(3), xij(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer :: n
    n = cell_d(recell,xai,xaj,xij,1)
  end function cell_a

  pure function cell_b(recell,xai,xaj,xij) result(n)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xai(3), xaj(3), xij(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer :: n
    n = cell_d(recell,xai,xaj,xij,2)
  end function cell_b

  pure function cell_c(recell,xai,xaj,xij) result(n)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xai(3), xaj(3), xij(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer :: n
    n = cell_d(recell,xai,xaj,xij,3)
  end function cell_c

  pure function cell_d(recell,xai,xaj,xij,abc) result(n)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: recell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: xai(3), xaj(3), xij(3)
    integer,  intent(in) :: abc
! *********************
! * OUTPUT variables  *
! *********************
    integer  :: n
    real(dp) :: xd(3)
    ! The actual length between two atomic centered orbitals:
    xd (:) = xij(:) - (xaj(:)-xai(:))
    n = nint(sum(xd*recell(:,abc)))
  end function cell_d


  function cell_v_dir(cell,v) result(dir)
    use intrinsic_missing, only : VEC_PROJ, VNORM
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: cell(3,3) ! the reciprocal cell, WITHOUT 2Pi!!
    real(dp), intent(in) :: v(3)
! *********************
! * OUTPUT variables  *
! *********************
    integer  :: dir
    real(dp) :: c(3)
    c(1) = vnorm(VEC_PROJ(cell(:,1)/vnorm(cell(:,1)),v))
    c(2) = vnorm(VEC_PROJ(cell(:,2)/vnorm(cell(:,2)),v))
    c(3) = vnorm(VEC_PROJ(cell(:,3)/vnorm(cell(:,3)),v))
    if ( c(1) > c(2) .and. c(1) > c(3) ) dir = 1
    if ( c(2) > c(1) .and. c(2) > c(3) ) dir = 2
    if ( c(3) > c(1) .and. c(3) > c(2) ) dir = 3
  end function cell_v_dir

! shift an atomic coordinate to lie inside the cell
! This enables direct comparisons between different atomic 
! coordinates as they will be lying comparatively.
  function xa_in_cell(rcell,ucell, xa) result(uxa)
! *********************
! * INPUT variables   *
! *********************
    real(dp), intent(in) :: rcell(3,3), ucell(3,3), xa(3)
! *********************
! * OUTPUT variables  *
! *********************
    real(dp) :: uxa(3)

    ! Local variables
    integer :: n
    n = floor(sum(xa*rcell(:,1)))
    uxa(:) =  xa(:) - n * ucell(:,1)
    n = floor(sum(xa*rcell(:,2)))
    uxa(:) = uxa(:) - n * ucell(:,2)
    n = floor(sum(xa*rcell(:,3)))
    uxa(:) = uxa(:) - n * ucell(:,3)
  end function xa_in_cell

end module geom_helper
