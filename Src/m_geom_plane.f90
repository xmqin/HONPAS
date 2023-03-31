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
!  PLANE

module m_geom_plane

  use m_geom_aux
  
  implicit none

  private

  ! Types
  public :: geo_plane_delta
  public :: geo_plane_gauss
  public :: geo_plane_exp
  ! Interfaces
  public :: fgeo_read_plane
  public :: correct_plane
  public :: voxel_in_plane
  public :: voxel_val_plane

  ! Explicit routines
  public :: fgeo_read_plane_delta
  public :: fgeo_read_plane_gauss
  public :: fgeo_read_plane_exp
  public :: correct_plane_delta
  public :: correct_plane_gauss
  public :: correct_plane_exp
  public :: voxel_in_plane_delta
  public :: voxel_in_plane_gauss
  public :: voxel_in_plane_exp
  public :: voxel_val_plane_delta
  public :: voxel_val_plane_gauss
  public :: voxel_val_plane_exp

  type geo_plane_delta
     ! Describing an infinite plane which will satisfy the following:
     !    n \dot ( x1 - c ) = 0
     ! where x1 is any coordinate (x,y,z)
     sequence
     ! The coordinate from which the normal vector
     ! is draw...
     real(dp) :: c(3)
     ! The NORMALIZED normal vector to the plane in point 'c'.
     real(dp) :: n(3)
     ! Auxillary distance length which is merely '- n dot c'
     ! This will ease the computations as it is used extensively
     real(dp) :: d = 0._dp
  end type geo_plane_delta

  type geo_plane_gauss
     sequence
     ! The gauss tail of the plane
     ! This is: 1/(2*\sigma^2)
     real(dp) :: inv_2var
     ! The cutoff-radii of the Gaussian function
     real(dp) :: cut_off

     real(dp) :: c(3)
     real(dp) :: n(3)
     real(dp) :: d_A = 0._dp ! d_Above
     real(dp) :: d_B = 0._dp ! d_Below
  end type geo_plane_gauss

  type geo_plane_exp
     sequence
     ! The gauss tail of the plane
     ! This is: h^{-1} =  + r_{1/2} / log(2.)
     real(dp) :: h
     ! The cutoff-radii of the exponential function
     real(dp) :: cut_off

     real(dp) :: c(3)
     real(dp) :: n(3)
     real(dp) :: d_A = 0._dp ! d_Above
     real(dp) :: d_B = 0._dp ! d_Below
  end type geo_plane_exp

  interface fgeo_read_plane
     module procedure fgeo_read_plane_delta
     module procedure fgeo_read_plane_gauss
     module procedure fgeo_read_plane_exp
  end interface

  interface correct_plane
     module procedure correct_plane_delta
     module procedure correct_plane_gauss
     module procedure correct_plane_exp
  end interface

  interface voxel_in_plane
     module procedure voxel_in_plane_delta
     module procedure voxel_in_plane_gauss
     module procedure voxel_in_plane_exp
  end interface

  interface voxel_val_plane
     module procedure voxel_val_plane_delta
     module procedure voxel_val_plane_gauss
     module procedure voxel_val_plane_exp
  end interface

contains

  ! This subroutine returns .true. in has if the box spanned in
  !    lr,lr+dx,lr+dx+dy,..., has any points inside the plane.
  !  - lr is the "lower-left" of the box.
  !  - d is the extension of the box

  ! TODO, for now this can only handle Cartesian boxes, it should be
  ! extended to be used in skewed axes....
  ! TODO this should be easy if we project the normal vector into
  ! the cell coordinates (i.e. do nb = (n dot bx, n dot by, n dot bz)


  function voxel_in_plane_delta(plane,ll,d) result(has)
    type(geo_plane_delta), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has
    ! The positive and negative vertices
    real(dp) :: p(3), n(3)

    ! We will use the (n,p)-vertex check
    ! We select 2 vertices which are those farthest (n=negative)
    ! and furthest (p=positiv) in the direction of the normal
    ! vector describing the plane.

    call voxel_pos_neg(ll,d,plane%n,p,n)

    has = .false.
    ! We can then consider whether the vertices lie both outside
    ! if not, we have an intersection.

    ! The positive vertex lies on the left side of the plane, 
    ! hence we know that it will not be crossing the plane
    if ( plane%n(1)*p(1)+plane%n(2)*p(2)+plane%n(3)*p(3) <  plane%d ) return

    if ( plane%n(1)*n(1)+plane%n(2)*n(2)+plane%n(3)*n(3) <= plane%d ) then
       ! The positive vertex lies on the right side of the plane, 
       ! This check ensures that the negative vertex lies on the left side
       ! of the plane. Hence we have an intersection.
       has = .true.
    end if

  end function voxel_in_plane_delta

  function voxel_in_plane_gauss(plane,ll,d) result(has)
    type(geo_plane_gauss), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has
    real(dp) :: p(3), n(3)

    call voxel_pos_neg(ll,d,plane%n,p,n)

    has = .false.
    if ( plane%n(1)*p(1)+plane%n(2)*p(2)+plane%n(3)*p(3) <  plane%d_B ) return
    if ( plane%n(1)*n(1)+plane%n(2)*n(2)+plane%n(3)*n(3) <= plane%d_A ) then
       has = .true.
    end if

  end function voxel_in_plane_gauss

  function voxel_in_plane_exp(plane,ll,d) result(has)
    type(geo_plane_exp), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has
    real(dp) :: p(3), n(3)

    call voxel_pos_neg(ll,d,plane%n,p,n)

    has = .false.
    if ( plane%n(1)*p(1)+plane%n(2)*p(2)+plane%n(3)*p(3) <  plane%d_B ) return 
    if ( plane%n(1)*n(1)+plane%n(2)*n(2)+plane%n(3)*n(3) <= plane%d_A ) then
       has = .true.
    end if

  end function voxel_in_plane_exp


  ! This subroutine returns the value of the charge
  ! that is attributed from the exponential function
  !  - ll is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_val_plane_gauss(plane,ll,d) result(vol)
    type(geo_plane_gauss), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: c(3)
    real(dp) :: vol
    real(dp) :: radii
    
    ! find the center of box
    c = ll + d * 0.5_dp

    radii = plane%n(1)*c(1)+plane%n(2)*c(2)+plane%n(3)*c(3) - &
         (plane%d_B + plane%cut_off)
    if ( abs(radii) <= plane%cut_off ) then
       vol = dexp(-radii**2*plane%inv_2var)
    else
       vol = 0._dp
    end if
    
  end function voxel_val_plane_gauss

  function voxel_val_plane_exp(plane,ll,d) result(vol)
    type(geo_plane_exp), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: c(3)
    real(dp) :: vol
    real(dp) :: radii
    
    ! find the center of box
    c = ll + d * 0.5_dp

    radii = abs(plane%n(1)*c(1)+plane%n(2)*c(2)+plane%n(3)*c(3) - &
         (plane%d_B + plane%cut_off))
    if ( radii <= plane%cut_off ) then
       vol = dexp(-radii*plane%h)
    else
       vol = 0._dp
    end if
    
  end function voxel_val_plane_exp

  function voxel_val_plane_delta(plane,ll,d) result(vol)
    type(geo_plane_delta), intent(in) :: plane
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: vol
    vol = 1._dp
  end function voxel_val_plane_delta

  ! Returns the positive and negative corner of the voxel
  pure subroutine voxel_pos_neg(ll,d,e,p,n)
    real(dp), intent(in) :: ll(3), d(3), e(3)
    real(dp), intent(out) :: p(3), n(3)
    integer :: i

    do i = 1 , 3
       if ( e(i) >= 0._dp ) then
          ! The normal-vector points in the i'th-direction
          p(i) = ll(i) + d(i)

          ! hence the negative vertex is the minimum coordinate
          n(i) = ll(i)
       else
          ! The normal-vector points in the negative i'th-direction
          p(i) = ll(i)
       
          ! hence the negative index lies in the positive i'th-direction
          n(i) = ll(i) + d(i)
       end if
    end do

  end subroutine voxel_pos_neg

  ! Subroutine for initializing the planes by transforming them
  ! to the corresponding cell
  subroutine correct_plane_delta(cell,plane)
    use intrinsic_missing, only : VNORM
    real(dp), intent(in) :: cell(3,3)
    type(geo_plane_delta), intent(inout) :: plane
    real(dp) :: tmp(3)
    real(dp) :: ncell(3)
    integer :: i

    ! Normalize the vector
    plane%n = plane%n / VNORM(plane%n)

    ! Correct the direction of the normal-vector to 
    ! be that of the cell-direction
    do i = 1 , 3 
       ncell = cell(:,i) / VNORM(cell(:,i))
       tmp(i) = &
            plane%n(1)*ncell(1) + &
            plane%n(2)*ncell(2) + &
            plane%n(3)*ncell(3)
    end do
    ! Copy the new normal vector in units of the cell-direction
    plane%n = tmp / VNORM(tmp)

    ! Calculate the auxillary distance
    plane%d = plane%n(1)*plane%c(1) &
         + plane%n(2)*plane%c(2) &
         + plane%n(3)*plane%c(3)

  end subroutine correct_plane_delta

  subroutine correct_plane_gauss(cell,plane)
    real(dp), intent(in) :: cell(3,3)
    type(geo_plane_gauss), intent(inout) :: plane
    type(geo_plane_delta) :: t_plane

    t_plane%n = plane%n
    t_plane%c = plane%c
    call correct_plane_delta(cell,t_plane)
    plane%n   = t_plane%n
    plane%d_A = t_plane%d + plane%cut_off
    plane%d_B = t_plane%d - plane%cut_off

  end subroutine correct_plane_gauss

  subroutine correct_plane_exp(cell,plane)
    real(dp), intent(in) :: cell(3,3)
    type(geo_plane_exp), intent(inout) :: plane
    type(geo_plane_delta) :: t_plane

    t_plane%n = plane%n
    t_plane%c = plane%c
    call correct_plane_delta(cell,t_plane)
    plane%n   = t_plane%n
    plane%d_A = t_plane%d + plane%cut_off
    plane%d_B = t_plane%d - plane%cut_off

  end subroutine correct_plane_exp


  subroutine fgeo_read_plane_delta(bName, ngeom,geom,params,par_unit)

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
    type(geo_plane_delta), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

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
    call fgeo_count(bName,GEOM_PLANE_DELTA,t)

    if ( t /= ngeom ) &
         call die('Could not find any delta planes')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')
    
    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_PLANE_DELTA ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step delta
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the delta mark')
       
       ! An infinite plane is defined by any point in the plane
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the point for the &
            &infinite plane of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%c,units=.true.)

       ! A plane is defined by the normal vector
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%n,units=.false.)

    end do count_geom

  end subroutine fgeo_read_plane_delta


  subroutine fgeo_read_plane_gauss(bName, ngeom,geom,params,par_unit)

    use intrinsic_missing, only : VNORM
    use fdf

! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bName
    character(len=*), intent(in), optional :: par_unit

! ********************
! * OUTPUT variables *
! ********************
    integer, intent(in) :: ngeom
    type(geo_plane_gauss), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

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
    call fgeo_count(bName,GEOM_PLANE_GAUSS,t)

    if ( t /= ngeom ) &
         call die('Could not find any Gaussian planes')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_PLANE_GAUSS ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the gauss mark')

       ! We need to read in the \sigma and cut_off
       call fgeo_read_vals(pline,geom(ip)%c(1:2),units=.true.,unit_offset=1)
       v = geom(ip)%c(1)
       geom(ip)%inv_2var = 1._dp/(2._dp*v**2)
       geom(ip)%cut_off  = geom(ip)%c(2)
       
       ! An infinite plane is defined by any point in the plane
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the point for the &
            &infinite plane of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%c,units=.true.)

       ! A plane is defined by the normal vector
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%n,units=.false.)

    end do count_geom

  end subroutine fgeo_read_plane_gauss

  subroutine fgeo_read_plane_exp(bName, ngeom,geom,params,par_unit)

    use intrinsic_missing, only : VNORM
    use fdf

! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bName
    character(len=*), intent(in), optional :: par_unit

! ********************
! * OUTPUT variables *
! ********************
    integer, intent(in) :: ngeom
    type(geo_plane_exp), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

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
    call fgeo_count(bName,GEOM_PLANE_EXP,t)

    if ( t /= ngeom ) &
         call die('Could not find any Exponential planes')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_PLANE_EXP ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the gauss mark')

       ! Calculate the exponential decay length
       ! this ensures that exp(-r/h) == .5 for r == half-length
       ! Note, we save the inverse as that is faster for computation
       call fgeo_read_vals(pline,geom(ip)%c(1:2),units=.true.,unit_offset=1)
       v = geom(ip)%c(1)
       geom(ip)%h       = log(2._dp) / v
       v = geom(ip)%c(2)
       geom(ip)%cut_off = v
       
       ! An infinite plane is defined by any point in the plane
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the point for the &
            &infinite plane of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%c,units=.true.)

       ! A plane is defined by the normal vector
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of a plane geometry object')
       call fgeo_read_vals(pline,geom(ip)%n,units=.false.)

    end do count_geom

  end subroutine fgeo_read_plane_exp

end module m_geom_plane
