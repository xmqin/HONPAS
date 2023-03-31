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
!  COORD

module m_geom_coord

  use m_geom_aux

  implicit none

  private

  ! Types
  public :: geo_coord_exp
  public :: geo_coord_gauss
  ! Interfaces
  public :: fgeo_read_coord
  public :: voxel_in_coord
  public :: voxel_val_coord

  ! Explicit routines
  public :: fgeo_read_coord_exp
  public :: voxel_in_coord_exp
  public :: voxel_val_coord_exp
  public :: fgeo_read_coord_gauss
  public :: voxel_in_coord_gauss
  public :: voxel_val_coord_gauss

  type geo_coord_exp
     ! The half-length in the exponential
     ! This is exactly where the charge will be halfed
     ! I.e. e^{-r/h} = .5.
     ! Hence the user specifies a length scale
     ! which will be calculated to be 
     ! h^{-1} =  - r_{1/2} / log(.5), log == natural logarithm.
     ! h^{-1} =  + r_{1/2} / log(2.), log == natural logarithm.
     ! NOTICE THAT WE SAVE THE INVERSE FOR EFFICIENCY
     real(dp) :: h
     ! The cutoff-radii of the exponential function
     real(dp) :: cut_off
     ! Number of atoms in this geometry
     integer :: na
     ! Atomic positions of the atoms
     real(dp), pointer :: xa(:,:) => null()
  end type geo_coord_exp

  type geo_coord_gauss
     ! The gauss tail of the plane
     ! This is: 1/(2*\sigma^2)
     real(dp) :: inv_2var
     ! The cutoff-radii of the exponential function
     real(dp) :: cut_off
     ! Number of atoms in this geometry
     integer :: na
     ! Atomic positions of the atoms
     real(dp), pointer :: xa(:,:) => null()
  end type geo_coord_gauss
  
  interface fgeo_read_coord
     module procedure fgeo_read_coord_exp
     module procedure fgeo_read_coord_gauss
  end interface

  interface voxel_in_coord
     module procedure voxel_in_coord_exp
     module procedure voxel_in_coord_gauss
  end interface 

  interface voxel_val_coord
     module procedure voxel_val_coord_exp
     module procedure voxel_val_coord_gauss
  end interface 

contains


  ! This subroutine returns .true. in has if the
  ! center of the box spanned in lr,lr+dx,lr+dx+dy,..., 
  ! is shorter than the cut-off radii of any of the
  ! atoms situated in the atomic configuration
  !  - lr is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_in_coord_exp(ce,ll,d) result(has)
    use intrinsic_missing, only : VNORM
    type(geo_coord_exp), intent(in) :: ce
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: rv(3)
    logical :: has
    integer :: ia
    real(dp) :: radii

    rv = ll + 0.5_dp * d
    ! Search through the atomic list
    do ia = 1 , ce%na
       ! Calculate the distance between the two points.
       radii = VNORM(ce%xa(:,ia) - rv)
       if ( radii <= ce%cut_off ) then
          has = .true.
          return
       end if
    end do
    has = .false.
    
  end function voxel_in_coord_exp

  function voxel_in_coord_gauss(cg,ll,d) result(has)
    use intrinsic_missing, only : VNORM
    type(geo_coord_gauss), intent(in) :: cg
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: rv(3)
    logical :: has
    integer :: ia
    real(dp) :: radii

    rv = ll + 0.5_dp * d
    ! Search through the atomic list
    do ia = 1 , cg%na
       ! Calculate the distance between the two points.
       radii = VNORM(cg%xa(:,ia) - rv)
       if ( radii <= cg%cut_off ) then
          has = .true.
          return
       end if
    end do
    has = .false.
    
  end function voxel_in_coord_gauss

  ! This subroutine returns the value of the charge
  ! that is attributed from the exponential function
  !  - ll is the "lower-left" of the box.
  !  - d is the extension of the box
  function voxel_val_coord_exp(ce,ll,d) result(vol)
    use intrinsic_missing, only : VNORM
    type(geo_coord_exp), intent(in) :: ce
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: rv(3)
    real(dp) :: vol
    integer :: ia
    real(dp) :: radii

    rv = ll + 0.5_dp * d
    vol = 0._dp
    ! Search through the atomic list
    do ia = 1 , ce%na
       ! Calculate the distance between the two points.
       radii = VNORM(ce%xa(:,ia) - rv)
       if ( radii <= ce%cut_off ) then
          ! Add the contribution from this geometrical object
          vol = vol + dexp(-radii*ce%h)
       end if
    end do
    
  end function voxel_val_coord_exp

  function voxel_val_coord_gauss(cg,ll,d) result(vol)
    use intrinsic_missing, only : VNORM
    type(geo_coord_gauss), intent(in) :: cg
    real(dp), intent(in) :: ll(3), d(3)
    real(dp) :: rv(3)
    real(dp) :: vol
    integer :: ia
    real(dp) :: radii

    rv = ll + 0.5_dp * d
    vol = 0._dp
    ! Search through the atomic list
    do ia = 1 , cg%na
       ! Calculate the distance between the two points.
       radii = VNORM(cg%xa(:,ia) - rv)
       if ( radii <= cg%cut_off ) then
          ! Add the contribution from this geometrical object
          vol = vol + dexp(-radii**2*cg%inv_2var)
       end if
    end do
    
  end function voxel_val_coord_gauss

    
  subroutine fgeo_read_coord_exp(bName, ngeom,geom,params,par_unit)
    
    use intrinsic_missing, only : VNORM
    use parallel, only : IONode
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
    type(geo_coord_exp), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer  :: t, ip, i
    real(dp) :: v, vv(2)

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_COORD_EXP,t)

    if ( t /= ngeom ) &
         call die('Could not find any coordinates')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_COORD_EXP ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the exp mark')
  
       ! The first value we encounter MUST be the
       ! half-length of the exponential function
       call fgeo_read_vals(pline,vv,units=.true.,unit_offset=1)
       if (any(vv <= 0._dp) ) call die("Exponential parameters &
            &are not positive and above 0. Please correct")
       ! Calculate the exponential decay length
       ! this ensures that exp(-r/h) == .5 for r == half-length
       ! Note, we save the inverse as that is faster for computation
       geom(ip)%h       = log(2._dp) / vv(1)
       geom(ip)%cut_off = vv(2)

       ! If the user requests an exponential decay
       ! where the cut-off radii is shorter than the 
       ! half-length radii, then the user should be
       ! notified. 
       ! Remark, that this would "simulate" spheres of constant charge
       ! rather than spheres of decaying charge, which could make sense...
       if ( IONode .and. vv(2) <= vv(1) ) then
          write(*,*)'The geometrical object has hard, (roughly) constant &
               &valued spheres. Your cut-off radii is shorter &
               &than your half-length constant'
       end if

       ! Now we need to read in the number of points
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the number of points &
            &in the geometrical object correctly. Please correct')
       
       ! Read in number of points
       geom(ip)%na = fdf_bintegers(pline,1)
       if ( geom(ip)%na < 1 ) &
            call die('You have defined a negative count (or 0) of points &
            &for your geometrical object?')

       ! Allocate to read in the stuff...
       nullify(geom(ip)%xa)
       allocate(geom(ip)%xa(3,geom(ip)%na))
       do i = 1 , geom(ip)%na
          if ( .not. fdf_bline(bfdf,pline) ) &
               call die('Could not read the correct number of coordinates.')
          call fgeo_read_vals(pline,geom(ip)%xa(:,i),units=.true.)
       end do

       ! End of reading the coord
       ! *******************

    end do count_geom

  end subroutine fgeo_read_coord_exp

  subroutine fgeo_read_coord_gauss(bName, ngeom,geom,params,par_unit)
    
    use intrinsic_missing, only : VNORM
    use parallel, only : IONode
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
    type(geo_coord_gauss), intent(out) :: geom(ngeom)
    real(dp), intent(out) :: params(ngeom)

! ********************
! * LOCAL variables  *
! ********************
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer  :: t, ip, i
    real(dp) :: v, vv(2)

    ! if none found, simply return
    if ( ngeom <= 0 ) return

    ! Count the number of planes
    call fgeo_count(bName,GEOM_COORD_GAUSS,t)

    if ( t /= ngeom ) &
         call die('Could not find any coordinates')

    ! Do the processing
    if ( .not. fdf_block(bName,bfdf) ) &
         call die('Could not find the block again...?')

    ip = 0
    count_geom: do while ( ip < ngeom )

       ! First loop until we have the geometry
       call fgeo_next(bfdf,t,v,unit=par_unit)

       ! If it is not the correct geometry, cycle...
       if ( t /= GEOM_COORD_GAUSS ) cycle
       
       ! Counter for geometry
       ip = ip + 1

       ! The parameter for the geometry
       params(ip) = v

       ! Step to gauss
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not step the gauss mark')
  
       ! The first value we encounter MUST be the
       ! std. dev gaussian
       call fgeo_read_vals(pline,vv,units=.true.,unit_offset=1)
       if (any(vv <= 0._dp) ) call die("Gaussian parameters &
            &are not positive and above 0. Please correct")
       ! Calculate the Gaussian variance
       ! Note, we save the inverse as that is faster for computation
       geom(ip)%inv_2var = 1._dp/(2._dp*vv(1)**2)
       geom(ip)%cut_off  = vv(2)

       ! If the user requests an std. dev
       ! where the cut-off radii is shorter than the 
       ! std. dev radii, then the user should be
       ! notified. 
       ! Remark, that this would "simulate" spheres of constant charge
       ! rather than spheres of decaying charge, which could make sense...
       if ( IONode .and. vv(2) <= vv(1) * 2._dp ) then
          write(*,*)'The geometrical object has hard, (roughly) constant &
               &valued spheres. Your cut-off radii is shorter &
               &than your std. deviation constant'
       end if

       ! Now we need to read in the number of points
       if ( .not. fdf_bline(bfdf,pline) ) &
            call die('Could not read the number of points &
            &in the geometrical object correctly. Please correct')
       
       ! Read in number of points
       geom(ip)%na = fdf_bintegers(pline,1)
       if ( geom(ip)%na < 1 ) &
            call die('You have defined a negative count (or 0) of points &
            &for your geometrical object?')

       ! Allocate to read in the stuff...
       nullify(geom(ip)%xa)
       allocate(geom(ip)%xa(3,geom(ip)%na))
       do i = 1 , geom(ip)%na
          if ( .not. fdf_bline(bfdf,pline) ) &
               call die('Could not read the correct number of coordinates.')
          call fgeo_read_vals(pline,geom(ip)%xa(:,i),units=.true.)
       end do

       ! End of reading the coord
       ! *******************

    end do count_geom

  end subroutine fgeo_read_coord_gauss

end module m_geom_coord
