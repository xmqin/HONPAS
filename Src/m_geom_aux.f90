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

module m_geom_aux

  implicit none

  integer, parameter :: dp = SELECTED_REAL_KIND(r=307,p=15)

  ! integer requests for types
  integer, parameter :: GEOM_NONE = 0

  ! Plane geometries
  integer, parameter :: GEOM_PLANE_DELTA = 100
  integer, parameter :: GEOM_PLANE_GAUSS = 101
  integer, parameter :: GEOM_PLANE_EXP   = 102

  ! Square geometries
  integer, parameter :: GEOM_SQUARE_DELTA = 200
  integer, parameter :: GEOM_SQUARE_GAUSS = 201
  integer, parameter :: GEOM_SQUARE_EXP   = 202

  ! Coordinate geometries
  integer, parameter :: GEOM_COORD_EXP = 300
  integer, parameter :: GEOM_COORD_GAUSS = 301

  ! 3D boxes
  integer, parameter :: GEOM_BOX_DELTA = 400

contains


  ! This is allowed for different block-names, so as to
  ! be used generally.
  ! This returns the next geometrical object in the FDF block
  subroutine fgeo_next(bfdf, t,val,unit)
    use fdf

! ********************
! * INPUT variables  *
! ********************
    type(block_fdf), intent(inout) :: bfdf
    character(len=*), intent(in), optional :: unit

! ********************
! * OUTPUT variables *
! ********************
    integer, intent(out) :: t
    real(dp), intent(out) :: val

! ********************
! * LOCAL variables  *
! ********************
    character(len=200) :: g, g_t
    type(parsed_line), pointer :: pline => null()

    t = GEOM_NONE

    ! First we need to determine the number of stuff in the block

    ! If the str is one of the following:
    val = 0._dp
    do while ( t == GEOM_NONE ) 
       ! Read next line, if not found, return
       if ( .not. fdf_bline(bfdf,pline) ) return

       if ( fdf_bnnames(pline) > 0 ) then

          ! Read in the line
          g = fdf_bnames(pline,1)

          if ( leqi(g,'plane') ) then ! We have an infinite plane
             ! pass
          else if ( leqi(g,'square') ) then ! We have a finite plane
             ! pass
          else if ( leqi(g,'coords') ) then ! We have coordinates
             ! pass
          else if ( leqi(g,'box') ) then ! We have a box
             ! pass
          else
             cycle
          end if

          ! Save the additional value...
          if ( fdf_bnvalues(pline) > 0 ) val = fdf_bvalues(pline,1)
          if ( present(unit) ) then
             if ( fdf_bnnames(pline) > 1 ) then
                val = val * fdf_convfac(fdf_bnames(pline,2),unit)
             else
                call die('Unit for parameter of geom-block is not provided, &
                     &please consult manual.')
             end if
          end if

          ! Figure out which kind of geometry we have
          if ( .not. fdf_bline(bfdf,pline) ) &
               call die('Structure of geom-block is errorneous, &
               &please correct.')

          g_t = fdf_bnames(pline,1)
          
          if ( leqi(g,'plane') ) then ! We have an infinite plane
             
             if ( leqi(g_t,'delta') ) then
                t = GEOM_PLANE_DELTA
             else if ( leqi(g_t,'gauss') ) then
                t = GEOM_PLANE_GAUSS
             else if ( leqi(g_t,'exp') ) then
                t = GEOM_PLANE_EXP
             end if

          else if ( leqi(g,'square') ) then ! We have a finite plane

             if ( leqi(g_t,'delta') ) then
                t = GEOM_SQUARE_DELTA
             else if ( leqi(g_t,'gauss') ) then
                t = GEOM_SQUARE_GAUSS
             else if ( leqi(g_t,'exp') ) then
                t = GEOM_SQUARE_EXP
             end if

          else if ( leqi(g,'box') ) then ! We have a box
          
             if ( leqi(g_t,'delta') ) then
                t = GEOM_BOX_DELTA
             end if

          else if ( leqi(g,'coords') ) then ! We have coordinates
          
             if ( leqi(g_t,'exp') ) then
                t = GEOM_COORD_EXP
             else if ( leqi(g_t,'gauss') ) then
                t = GEOM_COORD_GAUSS
             end if

          end if
          
       end if

    end do

    ! move back again so that the geometry routine will read in
    ! the line with type again (it may contain parsable values)
    if ( .not. fdf_bbackspace(bfdf) ) &
         call die('Could not move backwards in the FDF file')

  end subroutine fgeo_next

  subroutine fgeo_count(bName,t,count)
    use fdf
! ********************
! * INPUT variables  *
! ********************
    character(len=*), intent(in) :: bname
    integer, intent(in) :: t
    
! ********************
! * OUTPUT variables *
! ********************
    integer, intent(out) :: count

    type(block_fdf) :: bfdf
    integer :: nT
    real(dp) :: v

    count = 0
    if ( .not. fdf_block(bName,bfdf) ) return
    nT = 1
    do while ( nT /= GEOM_NONE ) 
       call fgeo_next(bfdf,nT,v)
       if ( nT == t ) then
          count = count + 1
       end if
    end do

  end subroutine fgeo_count


  ! Subroutine for doing generic read-ins of fdf values in a line.
  subroutine fgeo_read_vals(pline,v,units,val_offset,unit_offset)
    use fdf

! **********************
! * INPUT variables    *
! **********************
    type(parsed_line), pointer :: pline

! **********************
! * OUTPUT variables   *
! **********************
    real(dp), intent(out) :: v(:)

! **********************
! * OPTIONAL variables *
! **********************
    logical, intent(in), optional :: units
    integer, intent(in), optional :: val_offset, unit_offset

    logical :: lunits
    integer :: lval_off, lunit_off, i

    lunits = .true.
    if ( present(units) ) lunits = units
    lval_off = 0
    if ( present(val_offset) ) lval_off = val_offset
    lunit_off = 0
    if ( present(unit_offset) ) lunit_off = unit_offset
    lunit_off = lunit_off + lval_off
    do i = 1 , size(v)
       if ( lunits ) then
          v(i) = fdf_bvalues(pline,i+lval_off) &
               * fdf_convfac(fdf_bnames(pline,1+lunit_off),'Bohr') ! Correct for units
       else
          v(i) = fdf_bvalues(pline,i+lval_off)
       end if
    end do

  end subroutine fgeo_read_vals

end module m_geom_aux
