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
!

! This module contains the geometrical object:
!  SQUARE

module m_geom_square

  use m_geom_aux
  use m_geom_plane, only : geo_plane_delta
  use m_geom_plane, only : geo_plane_gauss
  use m_geom_plane, only : geo_plane_exp
  use m_geom_plane, only : voxel_in_plane
  use m_geom_plane, only : voxel_val_plane
  use m_geom_plane, only : correct_plane

  implicit none

  private

  ! Types
  public :: geo_square_delta
  public :: geo_square_gauss
  public :: geo_square_exp
  ! Interfaces
  public :: fgeo_read_square
  public :: correct_square
  public :: voxel_in_square
  public :: voxel_val_square

  ! Explicit routines
  public :: fgeo_read_square_delta
  public :: fgeo_read_square_gauss
  public :: fgeo_read_square_exp
  public :: correct_square_delta
  public :: correct_square_gauss
  public :: correct_square_exp
  public :: voxel_in_square_delta
  public :: voxel_in_square_gauss
  public :: voxel_in_square_exp
  public :: voxel_val_square_delta
  public :: voxel_val_square_gauss
  public :: voxel_val_square_exp

  type geo_square_delta
     ! Describing a bounded plane which can be described by the infinite plane
     ! One can retrieve the normal vector of the plane spanned by: v1 and v2 by
     !    n = (v1 - o) \cross (v2 - o)
     ! 'o' describes the origo of the
     ! plane spanned by the vectors, v1 and v2. 
     ! With n we describe the infinite plane as in the previous case
     ! and subsequently bound it by:
     !    ---------------------     --------
     !    | v1_1 | v2_1 | n_1 |     | x1_1 |
     !    | v1_2 | v2_2 | n_2 | x = | x1_2 |
     !    | v1_3 | v2_3 | n_3 |     | x1_3 |
     !    ---------------------     --------
     ! and x = (a,b,c)^T.
     ! Now to be in the bounded plane, a and b must both satisfy:
     !    \{a,b\} \in [0;1], we check whether c is small enough via
     ! the inifinite plane method.
     sequence
     ! The coordinate from which the plane is spanned
     ! by the two vectors v1 and v2 is assigned to p%c
     ! The vectors which span the plane
     real(dp) :: v1(3), v2(3)
     ! Infinite plane
     type(geo_plane_delta) :: p
  end type geo_square_delta

  type geo_square_gauss
     sequence
     real(dp) :: v1(3), v2(3)
     type(geo_plane_gauss) :: p
  end type geo_square_gauss

  type geo_square_exp
     sequence
     real(dp) :: v1(3), v2(3)
     type(geo_plane_exp) :: p
  end type geo_square_exp

  interface fgeo_read_square
     module procedure fgeo_read_square_delta
     module procedure fgeo_read_square_gauss
     module procedure fgeo_read_square_exp
  end interface

  interface correct_square
     module procedure correct_square_delta
     module procedure correct_square_gauss
     module procedure correct_square_exp
  end interface

  interface voxel_in_square
     module procedure voxel_in_square_delta
     module procedure voxel_in_square_gauss
     module procedure voxel_in_square_exp
  end interface

  interface voxel_val_square
     module procedure voxel_val_square_delta
     module procedure voxel_val_square_gauss
     module procedure voxel_val_square_exp
  end interface

contains

  ! This subroutine returns .true. in has if the box spanned in
  !    lr,lr+dx,lr+dx+dy,..., has any points inside the plane.
  !  - lr is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_in_square_delta(plane,ll,d) result(has)
    type(geo_square_delta), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has

    ! We retrieve whether the box passes through the
    ! infinite plane spanned by the local plane
    has = voxel_in_plane(plane%p,ll,d)

    ! If it does not intersect, we can safely return
    if ( .not. has ) return
    
    has = P_voxel_square(ll,d,plane%v1,plane%v2,plane%p%n,plane%p%c)

  end function voxel_in_square_delta

  function voxel_in_square_gauss(plane,ll,d) result(has)
    type(geo_square_gauss), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has

    has = voxel_in_plane(plane%p,ll,d)

    if ( .not. has ) return
    
    has = P_voxel_square(ll,d,plane%v1,plane%v2,plane%p%n,plane%p%c)
    
  end function voxel_in_square_gauss

  function voxel_in_square_exp(plane,ll,d) result(has)
    type(geo_square_exp), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has

    has = voxel_in_plane(plane%p,ll,d)

    if ( .not. has ) return
    
    has = P_voxel_square(ll,d,plane%v1,plane%v2,plane%p%n,plane%p%c)
    
  end function voxel_in_square_exp


  function P_voxel_square(ll,d,v1,v2,n,o) result(has)
    real(dp), intent(in) :: ll(3), d(3), v1(3), v2(3), n(3), o(3)
    logical :: has
    integer :: ipiv(4)
    real(dp) :: sys(3,4)
    
    ! We need to check whether the bounds of the voxel
    ! lies within the plane...
    ! Thus we solve the linear equation:
    ! Ax=B
    ! where A = [v1,v2,n] and B = ll+.5*d
    sys(:,1) = v1
    sys(:,2) = v2
    ! We do need the normal vector, otherwise
    ! we would not be able to retrieve the
    ! solution to the system. (i.e. we do not
    ! know if the center cuts the plane)
    sys(:,3) = n
    ! Move the point down to origo
    sys(:,4) = ll + .5_dp * d - o
    call dgesv(3,1,sys(1,1),3,ipiv,sys(1,4),3,ipiv(4))
    if ( ipiv(4) /= 0 ) then
       write(*,*) ' LAPACK-error: Could not solve the linear &
            &equation for the voxel'
       ! This should never happen as any point in space
       ! should be reacheable by the vectors 
       ! (we however, need that v1 and v2 are linearly 
       ! independent)
    end if

    ! We already know that it lies in the plane
    ! so now we need to ensure that it lies within
    ! the two spanning vectors plane...
    has= 0._dp <= sys(1,4) .and. sys(1,4) <= 1.0_dp .and. &
         0._dp <= sys(2,4) .and. sys(2,4) <= 1.0_dp

  end function P_voxel_square

  function voxel_val_square_delta(plane,ll,d) result(vol)
    type(geo_square_delta), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: vol
    vol = 1._dp
  end function voxel_val_square_delta

  ! This subroutine returns the value of the charge
  ! that is attributed from the gaussian function about the
  ! finite plane
  !  - ll is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_val_square_gauss(plane,ll,d) result(vol)
    type(geo_square_gauss), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: vol

    ! It has already been checked that it exists
    vol = voxel_val_plane(plane%p,ll,d)
    
  end function voxel_val_square_gauss

  function voxel_val_square_exp(plane,ll,d) result(vol)
    type(geo_square_exp), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: vol

    vol = voxel_val_plane(plane%p,ll,d)
    
  end function voxel_val_square_exp

  
  ! Subroutine for initializing the planes by transforming them
  ! to the corresponding cell
  subroutine correct_square_delta(cell,plane)
    real(dp), intent(in) :: cell(3,3)
    type(geo_square_delta), intent(inout) :: plane
    real(dp), parameter :: EPS = 1e-5 
    real(dp) :: fac
    logical :: inde
    integer :: i

    ! The simple test is that if any-one coordinate
    ! is zero, and the other is not, then
    ! they are linearly independent
    inde = .false.
    do i = 1 , 3
       if ( abs(plane%v1(i)) < EPS .and. &
            abs(plane%v2(i)) > EPS ) then
          inde = .true.
       else if ( abs(plane%v2(i)) < EPS .and. &
            abs(plane%v1(i)) > EPS ) then
          inde = .true.
       else if ( abs(plane%v1(i)) > EPS .and. &
            abs(plane%v2(i)) > EPS ) then
          ! We create the factor here (as we must not divide by 0) 
          fac = plane%v1(i) / plane%v2(i)
       end if
    end do
       
    if ( .not. inde ) then
       inde = .not. ( &
            count( abs(fac * plane%v2 - plane%v1) <= EPS ) == 3 )
    end if


    ! Kill if the vectors are linearly dependent...
    if ( .not. inde ) then
       call die('The geometric finite plane has linearly &
            &dependent spanning vectors. This is not allowed.')
    end if
    
    ! Correct the infinite plane that the thing lies in
    call correct_plane(cell,plane%p)
    
  end subroutine correct_square_delta

  ! Subroutine for initializing the planes by transforming them
  ! to the corresponding cell
  subroutine correct_square_gauss(cell,plane)
    real(dp), intent(in) :: cell(3,3)
    type(geo_square_gauss), intent(inout) :: plane
    real(dp), parameter :: EPS = 1e-5 
    type(geo_square_delta) :: t_plane
    type(geo_plane_delta) :: ti_plane

    t_plane%v1 = plane%v1
    t_plane%v2 = plane%v2
    ti_plane%n = plane%p%n
    ti_plane%c = plane%p%c
    t_plane%p  = ti_plane
    call correct_square_delta(cell,t_plane)
    plane%p%n   = t_plane%p%n
    plane%p%d_A = t_plane%p%d + plane%p%cut_off
    plane%p%d_B = t_plane%p%d - plane%p%cut_off

  end subroutine correct_square_gauss

  ! Subroutine for initializing the planes by transforming them
  ! to the corresponding cell
  subroutine correct_square_exp(cell,plane)
    real(dp), intent(in) :: cell(3,3)
    type(geo_square_exp), intent(inout) :: plane
    real(dp), parameter :: EPS = 1e-5 
    type(geo_square_delta) :: t_plane
    type(geo_plane_delta) :: ti_plane

    t_plane%v1 = plane%v1
    t_plane%v2 = plane%v2
    ti_plane%n = plane%p%n
    ti_plane%c = plane%p%c
    t_plane%p  = ti_plane
    call correct_square_delta(cell,t_plane)
    plane%p%n   = t_plane%p%n
    plane%p%d_A = t_plane%p%d + plane%p%cut_off
    plane%p%d_B = t_plane%p%d - plane%p%cut_off

  end subroutine correct_square_exp

    
  ! Reading in information from block about the square-delta object
  subroutine fgeo_read_square_delta(bName, ngeom,geom,params,par_unit)

    use intrinsic_missing, only : VNORM
    use fdf

! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bName
    integer, intent(in) :: ngeom
    character(len=*), intent(in), optional :: par_unit

! ********************
! * OUTPUT variables *
! ********************
    type(geo_square_delta), intent(out) :: geom(ngeom)
    real(dp),  intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer :: t, ip
    real(dp) :: v

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_SQUARE_DELTA,t)

    if ( t /= ngeom ) &
         call die('Could not find any delta squares')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')
    
    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_SQUARE_DELTA ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step delta
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the delta mark')
       
      ! A plane is defined by the lower-left corner
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the lower-left corner &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%p%c,units=.true.)

       ! A plane is defined by two spanning vectors
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v1,units=.true.)

       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the second vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v2,units=.true.)

       ! We need to calculate the equivalent 
       ! infinite plane for this finite plane
       geom(ip)%p%n(1) = &
            geom(ip)%v1(2)*geom(ip)%v2(3) - &
            geom(ip)%v1(3)*geom(ip)%v2(2)
       geom(ip)%p%n(2) = &
            geom(ip)%v1(3)*geom(ip)%v2(1) - &
            geom(ip)%v1(1)*geom(ip)%v2(3)
       geom(ip)%p%n(3) = &
            geom(ip)%v1(1)*geom(ip)%v2(2) - &
            geom(ip)%v1(2)*geom(ip)%v2(1)
       ! Normalize the normal vector
       geom(ip)%p%n = geom(ip)%p%n / VNORM(geom(ip)%p%n)

    end do count_geom

  end subroutine fgeo_read_square_delta


  subroutine fgeo_read_square_gauss(bName, ngeom,geom,params,par_unit)

    use intrinsic_missing, only : VNORM
    use fdf

! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bName
    integer, intent(in) :: ngeom
    character(len=*), intent(in), optional :: par_unit

! ********************
! * OUTPUT variables *
! ********************
    type(geo_square_gauss), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer  :: t, ip
    real(dp) :: v

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_SQUARE_GAUSS,t)

    if ( t /= ngeom ) &
         call die('Could not find any Gaussian squares')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_SQUARE_GAUSS ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the gauss mark')

       ! We need to read in the \sigma and cut_off
       call fgeo_read_vals(pline,geom(ip)%p%c(1:2),units=.true.,unit_offset=1)
       v = geom(ip)%p%c(1)
       geom(ip)%p%inv_2var = 1._dp/(2._dp*v**2)
       geom(ip)%p%cut_off  = geom(ip)%p%c(2)
       
      ! A plane is defined by the lower-left corner
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the lower-left corner &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%p%c,units=.true.)

       ! A plane is defined by two spanning vectors
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v1,units=.true.)

       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the second vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v2,units=.true.)

       ! We need to calculate the equivalent 
       ! infinite plane for this finite plane
       geom(ip)%p%n(1) = &
            geom(ip)%v1(2)*geom(ip)%v2(3) - &
            geom(ip)%v1(3)*geom(ip)%v2(2)
       geom(ip)%p%n(2) = &
            geom(ip)%v1(3)*geom(ip)%v2(1) - &
            geom(ip)%v1(1)*geom(ip)%v2(3)
       geom(ip)%p%n(3) = &
            geom(ip)%v1(1)*geom(ip)%v2(2) - &
            geom(ip)%v1(2)*geom(ip)%v2(1)
       ! Normalize the normal vector
       geom(ip)%p%n = geom(ip)%p%n / VNORM(geom(ip)%p%n)

    end do count_geom

  end subroutine fgeo_read_square_gauss

  subroutine fgeo_read_square_exp(bName, ngeom,geom,params,par_unit)

    use intrinsic_missing, only : VNORM
    use fdf

! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bName
    integer, intent(in) :: ngeom
    character(len=*), intent(in), optional :: par_unit

! ********************
! * OUTPUT variables *
! ********************
    type(geo_square_exp), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer  :: t, ip
    real(dp) :: v

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_SQUARE_EXP,t)

    if ( t /= ngeom ) &
         call die('Could not find any Gaussian squares')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_SQUARE_EXP ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the gauss mark')

       ! We need to read in the half-length and cut_off
       call fgeo_read_vals(pline,geom(ip)%p%c(1:2),units=.true.,unit_offset=1)
       v = geom(ip)%p%c(1)
       geom(ip)%p%h        = log(2._dp) / v
       geom(ip)%p%cut_off  = geom(ip)%p%c(2)
       
      ! A plane is defined by the lower-left corner
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the lower-left corner &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%p%c,units=.true.)

       ! A plane is defined by two spanning vectors
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v1,units=.true.)

       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the second vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%v2,units=.true.)

       ! We need to calculate the equivalent 
       ! infinite plane for this finite plane
       geom(ip)%p%n(1) = &
            geom(ip)%v1(2)*geom(ip)%v2(3) - &
            geom(ip)%v1(3)*geom(ip)%v2(2)
       geom(ip)%p%n(2) = &
            geom(ip)%v1(3)*geom(ip)%v2(1) - &
            geom(ip)%v1(1)*geom(ip)%v2(3)
       geom(ip)%p%n(3) = &
            geom(ip)%v1(1)*geom(ip)%v2(2) - &
            geom(ip)%v1(2)*geom(ip)%v2(1)
       ! Normalize the normal vector
       geom(ip)%p%n = geom(ip)%p%n / VNORM(geom(ip)%p%n)

    end do count_geom

  end subroutine fgeo_read_square_exp

end module m_geom_square
