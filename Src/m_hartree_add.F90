!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
!

! This module adds the possibility of reading in a certain
! FDF-block which enables the addition of extra Hartree potential in the
! unit-cell in planes, boxes, or for creating corrugated Hartree potentials.

module m_hartree_add

  ! We add the reading of the geometrical objects.
  use m_geom_objects
  use precision, only : dp

  implicit none

  private
  save

  public :: init_hartree_add
  public :: read_hartree_add
  public :: hartree_add

  ! The different objects that this module is in hartree of handling

  ! We check for minus 1 to assert that we have not read the block
  integer :: N_geom = -1

  ! Size of the geometric object
  integer :: N_p_d = 0
  type h_plane_delta
     sequence
     ! These objects describe any hartree addition within
     ! any one full unit-cell plane points.
     type(geo_plane_delta) :: geo
     ! Each of these readings correspond to the additional hartree
     ! that will be added to the planes shown above.
     real(dp) :: hartree
  end type h_plane_delta
  type(h_plane_delta), allocatable :: p_d(:)

  integer :: N_p_g = 0
  type h_plane_gauss
     sequence
     type(geo_plane_gauss) :: geo
     real(dp) :: hartree
  end type h_plane_gauss
  type(h_plane_gauss), allocatable :: p_g(:)

  integer :: N_p_e = 0
  type h_plane_exp
     sequence
     type(geo_plane_exp) :: geo
     real(dp) :: hartree
  end type h_plane_exp
  type(h_plane_exp), allocatable :: p_e(:)

  integer :: N_s_d = 0
  type h_square_delta
     sequence
     type(geo_square_delta) :: geo
     real(dp) :: hartree
  end type h_square_delta
  type(h_square_delta), allocatable :: s_d(:)

  integer :: N_s_g = 0
  type h_square_gauss
     sequence
     type(geo_square_gauss) :: geo
     real(dp) :: hartree
  end type h_square_gauss
  type(h_square_gauss), allocatable :: s_g(:)

  integer :: N_s_e = 0
  type h_square_exp
     sequence
     type(geo_square_exp) :: geo
     real(dp) :: hartree
  end type h_square_exp
  type(h_square_exp), allocatable :: s_e(:)

  integer :: N_b_d = 0
  type h_box_delta
     sequence
     type(geo_box_delta) :: geo
     real(dp) :: hartree
  end type h_box_delta
  type(h_box_delta), allocatable :: b_d(:)

  integer :: N_c_e = 0
  type h_coord_exp
     type(geo_coord_exp) :: geo
     real(dp) :: hartree
  end type h_coord_exp
  type(h_coord_exp), allocatable :: c_e(:)

  integer :: N_c_g = 0
  type h_coord_gauss
     type(geo_coord_gauss) :: geo
     real(dp) :: hartree
  end type h_coord_gauss
  type(h_coord_gauss), allocatable :: c_g(:)

contains

  ! We need to initialize information about the grid which we distribute in
  subroutine init_hartree_add(ucell, mesh)
    use parallel, only : IONode
    use units, only : Ang, eV
    use m_mesh_node, only : offset_r, dL, dMesh, meshl
#ifdef MPI
    use mpi_siesta
#endif
    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of mesh divisions of each lattice vector
    integer, intent(in) :: mesh(3)
    ! Local counters of the grid-intersections
    integer, allocatable :: count_is(:)
    integer :: ix,iy,iz, iC, i
    real(dp) :: ll(3), llZ(3), llYZ(3), tot_hartree
    
#ifdef MPI
    integer :: MPIerror
#endif

    if ( N_geom < 1 ) return

    ! Start the printing of the message
    if ( IONode ) then
       write(*,'(/,a)') 'Geometric Hartree distributions:'
    end if

    ! Allocate to count the number of cross-sections (we make it
    ! twice as big so that we can re-use it for MPI
    allocate(count_is(2*N_geom))
    count_is(:) = 0

    ! We do a loop in the local grid
!$OMP parallel do default(shared), &
!$OMP firstprivate(offset_r,dMesh,dL), &
!$OMP private(iz,llZ,iy,llYZ,ix,iC,ll,i), &
!$OMP reduction(+:count_is)
    do iz = 0 , meshl(3) - 1
       llZ(:) = offset_r + iz*dL(:,3)
       do iy = 0 , meshl(2) - 1
          llYZ(:) = iy*dL(:,2) + llZ(:)
          do ix = 0 , meshl(1) - 1
             ! initialize the count-counter
             iC = 1
             ! The lower-left corner of the current box
             ll = ix*dL(:,1) + llYZ

             ! Count entries in each geometric object
             do i = 1 , N_p_d
                if ( voxel_in(p_d(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             do i = 1 , N_p_g
                if ( voxel_in(p_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             do i = 1 , N_p_e
                if ( voxel_in(p_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do

             do i = 1 , N_s_d
                if ( voxel_in(s_d(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             do i = 1 , N_s_g
                if ( voxel_in(s_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             do i = 1 , N_s_e
                if ( voxel_in(s_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do

             do i = 1 , N_b_d
                if ( voxel_in(b_d(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do

             do i = 1 , N_c_e
                if ( voxel_in(c_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             do i = 1 , N_c_g
                if ( voxel_in(c_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                end if
                iC = iC + 1
             end do
             
          end do
       end do
    end do
!$OMP end parallel do

#ifdef MPI
    ! Sum the counts
    call MPI_AllReduce(count_is(1),count_is(N_geom+1), &
         N_geom, MPI_Integer, MPI_Sum, MPI_Comm_World,MPIerror)
    count_is(1:N_geom) = count_is(N_geom+1:N_geom*2)
#endif

    ! We will calculate how much hartree that
    ! should be added in every voxel.
    iC = 1
    tot_hartree = 0._dp
    do i = 1 , N_p_d
       tot_hartree = tot_hartree + p_d(i)%hartree
       call hartree_print('Infinite delta-plane',p_d(i)%hartree,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_p_g
       tot_hartree = tot_hartree + p_g(i)%hartree
       call hartree_print('Infinite Gauss-plane',p_g(i)%hartree,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_p_e
       tot_hartree = tot_hartree + p_e(i)%hartree
       call hartree_print('Infinite Exp-plane',p_e(i)%hartree,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_s_d
       tot_hartree = tot_hartree + s_d(i)%hartree
       call hartree_print('Finite delta-plane',s_d(i)%hartree,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_s_g
       tot_hartree = tot_hartree + s_g(i)%hartree
       call hartree_print('Finite Gauss-plane',s_g(i)%hartree,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_s_e
       tot_hartree = tot_hartree + s_e(i)%hartree
       call hartree_print('Finite Exp-plane',s_e(i)%hartree,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_b_d
       tot_hartree = tot_hartree + b_d(i)%hartree
       call hartree_print('Finite delta-box',b_d(i)%hartree,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_c_e
       tot_hartree = tot_hartree + c_e(i)%hartree
       call hartree_print('Exp. spheres',c_e(i)%hartree,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_c_g
       tot_hartree = tot_hartree + c_g(i)%hartree
       call hartree_print('Gaussian spheres',c_g(i)%hartree,count_is(iC))
       iC = iC + 1
    end do

    if ( IONode ) then
       if ( N_geom > 1 ) then ! If there is only one geometry, it makes no sense to print
          write(*,'(a,e12.5,a)') 'Total geometric Hartree potential added: ', &
               tot_hartree/eV,' eV'
       end if
       write(*,*) ! new-line
    end if
    
    deallocate(count_is)

  contains 

    subroutine hartree_print(info,hartree,counts)
      use parallel, only : IONode
      character(len=*), intent(in) :: info
      real(dp), intent(in) :: hartree
      integer, intent(in) :: counts
      if ( counts <= 0 ) then
         call die('No elements are touching the '//trim(info)//', this &
              &is not allowed. Please correct your input.')
      end if

      if ( IONode ) then
         write(*,'(tr3,a,t26,a,e12.5,a,i0,a)') &
              trim(info),' adding ',hartree/eV, ' eV in ',counts,' grid points'
      end if
    end subroutine hartree_print

  end subroutine init_hartree_add


  ! Routine for adding the hartrees in the regions of space where any
  ! geometries are defined.
  subroutine hartree_add(ucell,npt_l,Vscf)

    use precision, only : grid_p
    use m_mesh_node, only : offset_r, dL, dMesh, meshl

! *********************
! * INPUT variables   *
! *********************
    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Total number of local points in Vscf, per spin
    integer, intent(in) :: npt_l
    ! The electronic Hartree potential local to this processor 
    real(grid_p), intent(inout) :: Vscf(npt_l)

! *********************
! * LOCAL variables   *
! *********************
    integer :: ix, iy, iz, imesh, i
    real(dp) :: ll(3), llZ(3), llYZ(3)

    ! If no geometries exists, return
    if ( N_geom < 1 ) return

    ! The mesh is running along:
    !  x-direction, then y, then z
!$OMP parallel do default(shared), &
!$OMP firstprivate(offset_r,dMesh,dL), &
!$OMP private(iz,llZ,iy,llYZ,ix,ll,i,imesh)
    do iz = 0 , meshl(3) - 1
       imesh = iz * meshl(2) * meshl(1)
       llZ(:) = offset_r(:) + iz*dL(:,3)
       do iy = 0 , meshl(2) - 1
          llYZ(:) = iy*dL(:,2) + llZ(:)
          do ix = 0 , meshl(1) - 1
             imesh = imesh + 1
             ll = ix*dL(:,1) + llYZ

             ! probably a spin loop would be appropriate

             ! Add the extra hartree
             do i = 1 , N_p_d
                if ( voxel_in(p_d(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(p_d(i)%geo,ll,dMesh) * p_d(i)%hartree
                end if
             end do
             do i = 1 , N_p_g
                if ( voxel_in(p_g(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(p_g(i)%geo,ll,dMesh) * p_g(i)%hartree
                end if
             end do
             do i = 1 , N_p_e
                if ( voxel_in(p_e(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(p_e(i)%geo,ll,dMesh) * p_e(i)%hartree
                end if
             end do

             do i = 1 , N_s_d
                if ( voxel_in(s_d(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(s_d(i)%geo,ll,dMesh) * s_d(i)%hartree
                end if
             end do
             do i = 1 , N_s_g
                if ( voxel_in(s_g(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(s_g(i)%geo,ll,dMesh) * s_g(i)%hartree
                end if
             end do
             do i = 1 , N_s_e
                if ( voxel_in(s_e(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(s_e(i)%geo,ll,dMesh) * s_e(i)%hartree
                end if
             end do

             do i = 1 , N_b_d
                if ( voxel_in(b_d(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(b_d(i)%geo,ll,dMesh) * b_d(i)%hartree
                end if
             end do

             do i = 1 , N_c_e
                if ( voxel_in(c_e(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(c_e(i)%geo,ll,dMesh) * c_e(i)%hartree
                end if
             end do
             do i = 1 , N_c_g
                if ( voxel_in(c_g(i)%geo,ll, dMesh) ) then
                   Vscf(imesh) = Vscf(imesh) + &
                        voxel_val(c_g(i)%geo,ll,dMesh) * c_g(i)%hartree
                end if
             end do

          end do
       end do
    end do 
!$OMP end parallel do

  end subroutine hartree_add

  ! Routine for initialization of the options
  ! provided in the FDF
  subroutine read_hartree_add( )

    use parallel, only: IONode
    use m_cite, only: add_citation
    use intrinsic_missing, only : EYE

    real(dp) :: cell(3,3)
    integer :: i

    ! If any of the hartrees are defined, it means we have
    ! already read the options
    if ( N_geom > -1 ) return
    N_geom = 0

    ! Initialize a cell
    ! As we run through the geometric stuff in
    ! Cartesian coordinates, we do not need to correct
    ! vectors etc., 
    call EYE(3,cell)

    ! We read in the options for the Geometry.Hartree block
    call fgeo_count('Geometry.Hartree', GEOM_PLANE_DELTA, N_p_d)
    if ( N_p_d > 0 ) then
       allocate(p_d(N_p_d))
       call fgeo_read('Geometry.Hartree', N_p_d, p_d(:)%geo, &
            p_d(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_p_d

       do i = 1 , N_p_d
          call correct_plane(cell,p_d(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_PLANE_GAUSS, N_p_g)
    if ( N_p_g > 0 ) then
       allocate(p_g(N_p_g))
       call fgeo_read('Geometry.Hartree', N_p_g, p_g(:)%geo, &
            p_g(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_p_g
       do i = 1 , N_p_g
          call correct_plane(cell,p_g(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_PLANE_EXP, N_p_e)
    if ( N_p_e > 0 ) then
       allocate(p_e(N_p_e))
       call fgeo_read('Geometry.Hartree', N_p_e, p_e(:)%geo, &
            p_e(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_p_e
       do i = 1 , N_p_e
          call correct_plane(cell,p_e(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_SQUARE_DELTA, N_s_d)
    if ( N_s_d > 0 ) then
       allocate(s_d(N_s_d))
       call fgeo_read('Geometry.Hartree', N_s_d, s_d(:)%geo, &
            s_d(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_s_d
       do i = 1 , N_s_d
          call correct_square(cell,s_d(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_SQUARE_GAUSS, N_s_g)
    if ( N_s_g > 0 ) then
       allocate(s_g(N_s_g))
       call fgeo_read('Geometry.Hartree', N_s_g, s_g(:)%geo, &
            s_g(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_s_g
       do i = 1 , N_s_g
          call correct_square(cell,s_g(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_SQUARE_EXP, N_s_e)
    if ( N_s_e > 0 ) then
       allocate(s_e(N_s_e))
       call fgeo_read('Geometry.Hartree', N_s_e, s_e(:)%geo, &
            s_e(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_s_e
       do i = 1 , N_s_e
          call correct_square(cell,s_e(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Hartree', GEOM_BOX_DELTA, N_b_d)
    if ( N_b_d > 0 ) then
       allocate(b_d(N_b_d))
       call fgeo_read('Geometry.Hartree', N_b_d, b_d(:)%geo, &
            b_d(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_b_d
    end if

    call fgeo_count('Geometry.Hartree', GEOM_COORD_EXP, N_c_e)
    if ( N_c_e > 0 ) then
       allocate(c_e(N_c_e))
       call fgeo_read('Geometry.Hartree', N_c_e, c_e(:)%geo, &
            c_e(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_c_e
    end if

    call fgeo_count('Geometry.Hartree', GEOM_COORD_GAUSS, N_c_g)
    if ( N_c_g > 0 ) then
       allocate(c_g(N_c_g))
       call fgeo_read('Geometry.Hartree', N_c_g, c_g(:)%geo, &
            c_g(:)%hartree , par_unit='Ry')
       N_geom = N_geom + N_c_g
    end if

    if ( N_geom > 0 .and. IONode ) then
       call add_citation("10.1039/C5CP04613K")
    end if

  end subroutine read_hartree_add

end module m_hartree_add
  
