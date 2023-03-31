!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.
!

! This module contains the geometrical object:
!  BOX

module m_geom_box

  use m_geom_aux

  implicit none

  private

  ! Types
  public :: geo_box_delta
  ! Interfaces
  public :: fgeo_read_box
  public :: voxel_in_box
  public :: voxel_val_box

  ! Explicit routines
  public :: fgeo_read_box_delta
  public :: voxel_in_box_delta
  public :: voxel_val_box_delta

  type geo_box_delta
     ! Describing a bounded box which can be described by 3 vectors.
     ! With the three vectors we bound it by:
     !    ----------------------     --------
     !    | v1_1 | v2_1 | v3_1 |     | x1_1 |
     !    | v1_2 | v2_2 | v3_2 | x = | x1_2 |
     !    | v1_3 | v2_3 | v3_3 |     | x1_3 |
     !    ----------------------     --------
     ! and x = (a,b,c)^T.
     ! Now to be in the bounded box, a, b and c must satisfy:
     !    \{a,b,c\} \in [0;1].
     sequence
     ! The coordinate from which the plane is spanned
     ! by the three vectors v1, v2 and v3.
     ! The vectors which span the box
     ! This is now the inverse of the box vectors,
     ! i.e. we have v^{-1} to easy calculate x
     !   v^{-1} x1 = x
     real(dp) :: v(3,3)
     ! The origin of the vectors spanning the box
     real(dp) :: c(3)
  end type geo_box_delta

  interface fgeo_read_box
     module procedure fgeo_read_box_delta
  end interface

  interface voxel_in_box
     module procedure voxel_in_box_delta
  end interface

  interface voxel_val_box
     module procedure voxel_val_box_delta
  end interface

contains

  ! This subroutine returns .true. in has if the box spanned in
  !    lr,lr+dx,lr+dx+dy,..., has any points inside the plane.
  !  - lr is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_in_box_delta(box,ll,d) result(has)
    type(geo_box_delta), intent(in) :: box
    real(dp), intent(in) :: ll(3), d(3)
    logical :: has

    has = P_voxel_box(ll,d,box%v,box%c)

  end function voxel_in_box_delta

  function P_voxel_box(ll,d,v,o) result(has)
    real(dp), intent(in) :: ll(3), d(3), v(3,3), o(3)
    logical :: has
    real(dp) :: p(3), s(3)

    ! We need to check whether the bounds of the voxel
    ! lies within the plane...
    ! Thus we solve the linear equation:
    ! Ax=B
    ! where A = [v1,v2,v3] and B = ll+.5*d
    ! Move the point down to origo
    p = ll + .5_dp * d - o
    s(1) = v(1,1)*p(1)+v(1,2)*p(2)+v(1,3)*p(3)
    s(2) = v(2,1)*p(1)+v(2,2)*p(2)+v(2,3)*p(3)
    s(3) = v(3,1)*p(1)+v(3,2)*p(2)+v(3,3)*p(3)

    ! We already know that it lies in the plane
    ! so now we need to ensure that it lies within
    ! the two spanning vectors plane...
    has= 0._dp <= s(1) .and. s(1) <= 1.0_dp .and. &
         0._dp <= s(2) .and. s(2) <= 1.0_dp .and. &
         0._dp <= s(3) .and. s(3) <= 1.0_dp

  end function P_voxel_box

  function voxel_val_box_delta(box,ll,d) result(vol)
    type(geo_box_delta), intent(in) :: box
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: vol
    vol = 1._dp
  end function voxel_val_box_delta

  ! Reading in information from block about the box-delta object
  subroutine fgeo_read_box_delta(bName, ngeom,geom,params,par_unit)

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
    type(geo_box_delta), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer :: t, ip
    integer :: ipiv(3), info
    real(dp) :: v, work(12)

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_BOX_DELTA,t)

    if ( t /= ngeom ) &
         call die('Could not find any delta boxes')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')
    
    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_BOX_DELTA ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step delta
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the delta mark')
       
       ! A box is defined by the lower-left corner
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the lower-left corner &
            &of the box geometry object')
       call fgeo_read_vals(pline,geom(ip)%c,units=.true.)

       ! A box is defined by three spanning vectors
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the first vector &
            &of the box geometry object')
       call fgeo_read_vals(pline,geom(ip)%v(:,1),units=.true.)

       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the second vector &
            &of the box geometry object')
       call fgeo_read_vals(pline,geom(ip)%v(:,2),units=.true.)

       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the third vector &
            &of the box geometry object')
       call fgeo_read_vals(pline,geom(ip)%v(:,3),units=.true.)

       ! Calculate the inverse of v to faster make a solution
       ! This will prevent lapack dgesv calls and rely on simple
       ! vector products.
       ! Although inverses are not as precise they should be
       ! more than adequate in these situations.
       call dgetrf(3,3,geom(ip)%v,3,ipiv,info)
       if ( info /= 0 ) then
          call die('Box geometry could not invert vectors')
       end if
       call dgetri(3,geom(ip)%v,3,ipiv,work,12,info)
       if ( info /= 0 ) then
          call die('Box geometry could not invert vectors')
       end if
       
    end do count_geom

  end subroutine fgeo_read_box_delta

end module m_geom_box
