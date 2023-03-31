!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
!

! This module adds the possibility of reading in a certain
! FDF-block which enables the addition of extra charge in the
! unit-cell in planes, or for creating corrugated charge densities.

module m_charge_add

  ! We add the reading of the geometrical objects.
  use m_geom_objects
  use precision, only : dp

  implicit none

  private
  save

  public :: init_charge_add
  public :: read_charge_add
  public :: charge_add

  ! The different objects that this module is in charge of handling

  ! We check for minus 1 to assert that we have not read the block
  integer :: N_geom = -1


  ! Size of the geometric object
  integer :: N_p_d = 0
  type c_plane_delta
     sequence
     ! These objects describe any charge addition within
     ! any one full unit-cell plane points.
     type(geo_plane_delta) :: geo
     ! Each of these readings correspond to the additional charge
     ! that will be added to the planes shown above.
     real(dp) :: charge
     ! This array will hold the division charges that will be adjusted
     ! on each grid-point, they can be viewed as
     !    dRho(i) = charge_inf_planes(i)/count(voxel_in(plane))
     real(dp) :: dRho
  end type c_plane_delta
  type(c_plane_delta), allocatable :: p_d(:)

  integer :: N_p_g = 0
  type c_plane_gauss
     sequence
     type(geo_plane_gauss) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_plane_gauss
  type(c_plane_gauss), allocatable :: p_g(:)

  integer :: N_p_e = 0
  type c_plane_exp
     sequence
     type(geo_plane_exp) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_plane_exp
  type(c_plane_exp), allocatable :: p_e(:)

  integer :: N_s_d = 0
  type c_square_delta
     sequence
     type(geo_square_delta) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_square_delta
  type(c_square_delta), allocatable :: s_d(:)

  integer :: N_s_g = 0
  type c_square_gauss
     sequence
     type(geo_square_gauss) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_square_gauss
  type(c_square_gauss), allocatable :: s_g(:)

  integer :: N_s_e = 0
  type c_square_exp
     sequence
     type(geo_square_exp) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_square_exp
  type(c_square_exp), allocatable :: s_e(:)

  integer :: N_b_d = 0
  type c_box_delta
     sequence
     type(geo_box_delta) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_box_delta
  type(c_box_delta), allocatable :: b_d(:)

  integer :: N_c_e = 0
  type c_coord_exp
     type(geo_coord_exp) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_coord_exp
  type(c_coord_exp), allocatable :: c_e(:)

  integer :: N_c_g = 0
  type c_coord_gauss
     type(geo_coord_gauss) :: geo
     real(dp) :: charge
     real(dp) :: dRho
  end type c_coord_gauss
  type(c_coord_gauss), allocatable :: c_g(:)

contains

  ! We need to initialize information about the grid which we distribute
  ! in
  subroutine init_charge_add(ucell, mesh)
    use parallel, only : IONode
    use units, only : Ang
    use m_mesh_node, only : offset_r, dL, dMesh, meshl
#ifdef MPI
    use mpi_siesta
#endif
    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of mesh distributevisions of each lattice vector
    integer, intent(in) :: mesh(3)
    ! Local counters of the grid-intersections
    integer, allocatable :: count_is(:)
    real(dp), allocatable :: val_sum(:)
    integer :: ix,iy,iz, iC, i
    real(dp) :: ll(3), llZ(3), llYZ(3), tot_char
    
    ! To calculate the volume of the mesh-elements
    real(dp) :: volcel, dVol
    external :: volcel
#ifdef MPI
    integer :: MPIerror
#endif

    if ( N_geom < 1 ) return

    ! Start the printing of the message
    if ( IONode ) then
       write(*,'(/,a)') 'Geometric charge distributions:'
    end if

    ! We calculate the volume of one mesh-element
    ! This is needed to correct the amount of charge put in
    ! every mesh element
    dVol = volcel(dL)

    ! We must reset the dRho that our planes, etc. adds
    ! to the grid.
    ! Hence, we will step through the grids and count

    ! Allocate to count the number of cross-sections (we make it
    ! twice as big so that we can re-use it for MPI
    allocate(count_is(2*N_geom))
    count_is(:) = 0
    allocate(val_sum(2*N_geom))
    val_sum(:) = 0._dp
    
    ! We do a loop in the local grid
!$OMP parallel do default(shared), &
!$OMP firstprivate(offset_r,dMesh,dL), &
!$OMP private(iz,llZ,iy,llYZ,ix,iC,ll,i), &
!$OMP reduction(+:count_is,val_sum)
    do iz = 0 , meshl(3) - 1
       llZ(:) = offset_r(:) + iz*dL(:,3)
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
                   val_sum(iC)  = val_sum(iC) + voxel_val(p_d(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do
             do i = 1 , N_p_g
                if ( voxel_in(p_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(p_g(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do
             do i = 1 , N_p_e
                if ( voxel_in(p_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(p_e(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do

             do i = 1 , N_s_d
                if ( voxel_in(s_d(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(s_d(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do
             do i = 1 , N_s_g
                if ( voxel_in(s_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(s_g(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do
             do i = 1 , N_s_e
                if ( voxel_in(s_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(s_e(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do

             do i = 1 , N_b_d
                if ( voxel_in(b_d(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(b_d(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do

             do i = 1 , N_c_e
                if ( voxel_in(c_e(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(c_e(i)%geo,ll,dMesh)
                end if
                iC = iC + 1
             end do
             do i = 1 , N_c_g
                if ( voxel_in(c_g(i)%geo,ll, dMesh) ) then
                   count_is(iC) = count_is(iC) + 1
                   val_sum(iC)  = val_sum(iC) + voxel_val(c_g(i)%geo,ll,dMesh)
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

    ! Sum the contributions from the methods
    call MPI_AllReduce(val_sum(1),val_sum(N_geom+1), &
         N_geom, MPI_Double_Precision, MPI_Sum, MPI_Comm_World,MPIerror)
    val_sum(1:N_geom) = val_sum(N_geom+1:N_geom*2)
#endif

    ! We will calculate how much charge that
    ! should be added in every voxel.
    iC = 1
    tot_char = 0._dp
    do i = 1 , N_p_d
       tot_char = tot_char + p_d(i)%charge
       call calc_dRho_print('Infinite delta-plane', &
            p_d(i)%charge,dVol,val_sum(iC),p_d(i)%dRho,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_p_g
       tot_char = tot_char + p_g(i)%charge
       call calc_dRho_print('Infinite Gauss-plane', &
            p_g(i)%charge,dVol,val_sum(iC),p_g(i)%dRho,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_p_e
       tot_char = tot_char + p_e(i)%charge
       call calc_dRho_print('Infinite Exp-plane', &
            p_e(i)%charge,dVol,val_sum(iC),p_e(i)%dRho,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_s_d
       tot_char = tot_char + s_d(i)%charge
       call calc_dRho_print('Finite delta-plane', &
            s_d(i)%charge,dVol,val_sum(iC),s_d(i)%dRho,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_s_g
       tot_char = tot_char + s_g(i)%charge
       call calc_dRho_print('Finite Gauss-plane', &
            s_g(i)%charge,dVol,val_sum(iC),s_g(i)%dRho,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_s_e
       tot_char = tot_char + s_e(i)%charge
       call calc_dRho_print('Finite Exp-plane', &
            s_e(i)%charge,dVol,val_sum(iC),s_e(i)%dRho,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_b_d
       tot_char = tot_char + b_d(i)%charge
       call calc_dRho_print('Finite delta-box', &
            b_d(i)%charge,dVol,val_sum(iC),b_d(i)%dRho,count_is(iC))
       iC = iC + 1
    end do

    do i = 1 , N_c_e
       tot_char = tot_char + c_e(i)%charge
       call calc_dRho_print('Exp. spheres', &
            c_e(i)%charge,dVol,val_sum(iC),c_e(i)%dRho,count_is(iC))
       iC = iC + 1
    end do
    do i = 1 , N_c_g
       tot_char = tot_char + c_g(i)%charge
       call calc_dRho_print('Gaussian spheres', &
            c_g(i)%charge,dVol,val_sum(iC),c_g(i)%dRho,count_is(iC))
       iC = iC + 1
    end do

    if ( IONode ) then
       if ( N_geom > 1 ) then ! If there is only one geometry, it makes no sense to print
          write(*,'(a,e12.5)') 'Total geometric charge added: ',tot_char
       end if
       write(*,*) ! new-line
    end if
    
    deallocate(count_is)
    deallocate(val_sum)

  contains 

    subroutine calc_dRho_print(info,charge,dVol,val_sum,dRho,counts)
      use parallel, only : IONode
      character(len=*), intent(in) :: info
      real(dp), intent(in) :: charge, dVol, val_sum
      integer, intent(in) :: counts
      real(dp), intent(out) :: dRho
      if ( counts > 0 ) then
         dRho = charge / dVol / val_sum
      else
         call die('No elements are touching the '//trim(info)//', this &
              &leads to incorrect doping. Please correct your input.')
      end if

      if ( IONode ) then
         write(*,'(tr3,a,t26,a,e12.5,a,i0,a)') &
              trim(info),' adding ',charge,' e divided in ',counts,' grid points'
      end if
    end subroutine calc_dRho_print

  end subroutine init_charge_add


  ! Routine for adding the charges in the regions of space where any
  ! geometries are defined.
  subroutine charge_add(sign,ucell,npt_l,DRho)

    use precision, only : grid_p
    use m_mesh_node, only : offset_r, dL, dMesh, meshl

! *********************
! * INPUT variables   *
! *********************
    ! The sign of the operation '+' for addition, '-' for subtraction
    character(len=1), intent(in) :: sign
    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Total number of local points in rho, per spin
    integer, intent(in) :: npt_l
    ! The electronic density local to this processor 
    ! (note that this is the total density, hence charge is total spin charge)
    real(grid_p), intent(inout) :: DRho(npt_l)

! *********************
! * LOCAL variables   *
! *********************
    integer :: ix, iy, iz, imesh, i
    integer :: lsign
    real(dp) :: ll(3), llZ(3), llYZ(3)

    ! If no geometries exists, return
    if ( N_geom < 1 ) return

    ! Update sign
    lsign = -1
    if ( sign == '+' ) lsign = 1

    ! The mesh is running along:
    !  x-direction, then y, then z
!$OMP parallel do default(shared), &
!$OMP firstprivate(offset_r,dMesh,dL,lsign), &
!$OMP private(iz,llZ,iy,llYZ,ix,ll,i,imesh)
    do iz = 0 , meshl(3) - 1
       imesh = iz * meshl(2) * meshl(1)
       llZ(:) = offset_r + iz*dL(:,3)
       do iy = 0 , meshl(2) - 1
          llYZ(:) = iy*dL(:,2) + llZ(:)
          do ix = 0 , meshl(1) - 1
             imesh = imesh + 1
             ll = ix*dL(:,1) + llYZ

             ! Add the extra charge
             do i = 1 , N_p_d
                if ( voxel_in(p_d(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(p_d(i)%geo,ll,dMesh) * p_d(i)%dRho
                end if
             end do
             do i = 1 , N_p_g
                if ( voxel_in(p_g(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(p_g(i)%geo,ll,dMesh) * p_g(i)%dRho
                end if
             end do
             do i = 1 , N_p_e
                if ( voxel_in(p_e(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(p_e(i)%geo,ll,dMesh) * p_e(i)%dRho
                end if
             end do

             do i = 1 , N_s_d
                if ( voxel_in(s_d(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(s_d(i)%geo,ll,dMesh) * s_d(i)%dRho
                end if
             end do
             do i = 1 , N_s_g
                if ( voxel_in(s_g(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(s_g(i)%geo,ll,dMesh) * s_g(i)%dRho
                end if
             end do
             do i = 1 , N_s_e
                if ( voxel_in(s_e(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(s_e(i)%geo,ll,dMesh) * s_e(i)%dRho
                end if
             end do

             do i = 1 , N_b_d
                if ( voxel_in(b_d(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(b_d(i)%geo,ll,dMesh) * b_d(i)%dRho
                end if
             end do

             do i = 1 , N_c_e
                if ( voxel_in(c_e(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(c_e(i)%geo,ll,dMesh) * c_e(i)%dRho
                end if
             end do
             do i = 1 , N_c_g
                if ( voxel_in(c_g(i)%geo,ll, dMesh) ) then
                   DRho(imesh) = DRho(imesh) + lsign * &
                        voxel_val(c_g(i)%geo,ll,dMesh) * c_g(i)%dRho
                end if
             end do

          end do
       end do
    end do 
!$OMP end parallel do
    
  end subroutine charge_add

  ! Routine for initialization of the options
  ! provided in the FDF
  subroutine read_charge_add(nspin,char_net)

    use parallel, only: IONode
    use fdf, only: fdf_deprecated, fdf_obsolete
    use m_cite, only: add_citation
    use intrinsic_missing, only : EYE

    integer, intent(in) :: nspin
    ! The net-charge of the system (i.e. a sum of all things)
    real(dp), intent(inout) :: char_net
    real(dp) :: cell(3,3)
    integer  :: i

    ! If any of the charges are defined, it means we have
    ! already read the options
    if ( N_geom > -1 ) return
    N_geom = 0

    ! Initialize a cell
    ! As we run through the geometric stuff in
    ! Cartesian coordinates, we do not need to correct
    ! vectors etc., 
    call EYE(3,cell)

    ! notify about deprecated routines
    call fdf_deprecated('ChargeGeometries','Geometry.Charge')
    call fdf_obsolete('ChargeGeometries')

    ! We read in the options for the Geometry.Charge block
    call fgeo_count('Geometry.Charge', GEOM_PLANE_DELTA, N_p_d)
    if ( N_p_d > 0 ) then
       allocate(p_d(N_p_d))
       call fgeo_read('Geometry.Charge', N_p_d, p_d(:)%geo, p_d(:)%charge)
       N_geom = N_geom + N_p_d

       ! Update the charge of the system
       ! As we put the same (opposite) charge in the system
       ! we will sum here
       char_net = char_net + sum(p_d(:)%charge)

       ! Correct the planes
       do i = 1 , N_p_d
          call correct_plane(cell,p_d(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_PLANE_GAUSS, N_p_g)
    if ( N_p_g > 0 ) then
       allocate(p_g(N_p_g))
       call fgeo_read('Geometry.Charge', N_p_g, p_g(:)%geo, p_g(:)%charge)
       N_geom = N_geom + N_p_g
       char_net = char_net + sum(p_g(:)%charge)

       do i = 1 , N_p_g
          call correct_plane(cell,p_g(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_PLANE_EXP, N_p_e)
    if ( N_p_e > 0 ) then
       allocate(p_e(N_p_e))
       call fgeo_read('Geometry.Charge', N_p_e, p_e(:)%geo, p_e(:)%charge)
       N_geom = N_geom + N_p_e
       char_net = char_net + sum(p_e(:)%charge)

       do i = 1 , N_p_e
          call correct_plane(cell,p_e(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_SQUARE_DELTA, N_s_d)
    if ( N_s_d > 0 ) then
       allocate(s_d(N_s_d))
       call fgeo_read('Geometry.Charge', N_s_d, s_d(:)%geo, s_d(:)%charge)
       N_geom = N_geom + N_s_d
       char_net = char_net + sum(s_d(:)%charge)

       do i = 1 , N_s_d
          call correct_square(cell,s_d(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_SQUARE_GAUSS, N_s_g)
    if ( N_s_g > 0 ) then
       allocate(s_g(N_s_g))
       call fgeo_read('Geometry.Charge', N_s_g, s_g(:)%geo, s_g(:)%charge)
       N_geom = N_geom + N_s_g
       char_net = char_net + sum(s_g(:)%charge)

       do i = 1 , N_s_g
          call correct_square(cell,s_g(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_SQUARE_EXP, N_s_e)
    if ( N_s_e > 0 ) then
       allocate(s_e(N_s_e))
       call fgeo_read('Geometry.Charge', N_s_e, s_e(:)%geo, s_e(:)%charge)
       N_geom = N_geom + N_s_e
       char_net = char_net + sum(s_e(:)%charge)

       do i = 1 , N_s_e
          call correct_square(cell,s_e(i)%geo)
       end do
    end if

    call fgeo_count('Geometry.Charge', GEOM_BOX_DELTA, N_b_d)
    if ( N_b_d > 0 ) then
       allocate(b_d(N_b_d))
       call fgeo_read('Geometry.Charge', N_b_d, b_d(:)%geo, b_d(:)%charge)
       N_geom = N_geom + N_b_d
       char_net = char_net + sum(b_d(:)%charge)
    end if

    call fgeo_count('Geometry.Charge', GEOM_COORD_EXP, N_c_e)
    if ( N_c_e > 0 ) then
       allocate(c_e(N_c_e))
       call fgeo_read('Geometry.Charge', N_c_e, c_e(:)%geo, c_e(:)%charge)
       N_geom = N_geom + N_c_e
       char_net = char_net + sum(c_e(:)%charge)
    end if

    call fgeo_count('Geometry.Charge', GEOM_COORD_GAUSS, N_c_g)
    if ( N_c_g > 0 ) then
       allocate(c_g(N_c_g))
       call fgeo_read('Geometry.Charge', N_c_g, c_g(:)%geo, c_g(:)%charge)
       N_geom = N_geom + N_c_g
       char_net = char_net + sum(c_g(:)%charge)
    end if

    if ( N_geom > 0 .and. IONode ) then
       call add_citation("10.1039/C5CP04613K")
    end if
    
  end subroutine read_charge_add

end module m_charge_add
  
