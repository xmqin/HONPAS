!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_ts_voltage

  ! Module for containing various ways of implementing the voltage in a
  ! TranSIESTA run.
  ! As a standard method TranSIESTA applies the voltage ramp across the
  ! entire unit cell.
  ! As a future progress the voltage ramp might be situated in the contact region
  ! only. This makes more physical sense as the voltage drop does not occur
  ! in the electrode regions.
  !
  ! Created and copyrighted by: Nick Papior Andersen, 2012
  ! The use of this program is allowed for not-for-profit research only.
  ! Copy or disemination of all or part of this package is not
  ! permitted without prior and explicit authorization by the author.

  use precision, only : dp
  use m_ts_tdir

  implicit none

  private

  ! The idea is to have sub routines in this module to do
  ! various voltage layouts
  public :: ts_voltage
  public :: ts_init_voltage

  ! The corresponding bias for either the left electrode
  real(dp), save :: V_low = 0._dp
  real(dp), save :: V_high = 0._dp

contains

  subroutine ts_init_voltage(cell, na_u, xa, nmesh)
    use parallel, only : IONode
    use m_ts_electype, only : TotUsedAtoms, Elec
    use m_ts_options, only : Elecs, N_Elec
    use m_ts_options, only : Volt, Hartree_fname
    use units, only : eV, ang

    use m_geom_box, only : voxel_in_box

    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: nmesh(3) ! total number of mesh-points.

    logical :: bool
    real(dp) :: tmp, ll(3), zero(3)
    integer :: iElL, iElR, iEl, ia, iia

    if ( IONode ) then
      write(*,*)
    end if

    if ( ts_tidx < 1 ) then

      if ( IONode .and. len_trim(Hartree_fname) == 0 ) then
        write(*,'(a)')'ts-voltage: Lifted locally on each electrode'
        write(*,'(a)')'ts-voltage: WARNING ****************************************'
        write(*,'(a)')'ts-voltage: WARNING Prefer to use a custom Poisson solution!'
        write(*,'(a)')'ts-voltage: WARNING ****************************************'

      else if ( IONode .and. len_trim(Hartree_fname) > 0 ) then
        write(*,'(3a)')'ts-voltage: ', &
            'User supplied Poisson solution in file ', &
            trim(Hartree_fname)
#ifdef NCDF_4
        call ts_ncdf_voltage_assert(Hartree_fname,cell,nmesh)
#else
        call die('NetCDF4 has not been compiled in, this functionality &
            &is not supported.')
#endif
      end if

      ! Find the lowest and highest chemical potential
      V_low = Elecs(1)%mu%mu
      V_high = Elecs(1)%mu%mu
      do iEl = 2 , size(Elecs)
        V_low = min(Elecs(iEl)%mu%mu,V_low)
        V_high = max(Elecs(iEl)%mu%mu,V_high)
      end do
      if ( Volt < 0._dp ) then
        ! with a negative bias, we have to reverse
        ! the high-low
        ! Remember that in this case V_high and V-low are
        ! only used for the external Hartree potential.
        tmp = V_high
        V_high = V_low
        V_low = tmp
      end if

      ! Check that all electrode atoms are residing in the boxes
      ! defined by the electrodes
      bool = .false.
      zero = 0._dp
      do iEl = 1 , N_Elec
        do ia = 1 , TotUsedAtoms(Elecs(iEl))
          iia = Elecs(iEl)%idx_a - 1 + ia
          ll = xa(:,iia)
          if (.not.voxel_in_box(Elecs(iEl)%box,ll,zero)) then
            if ( IONode ) & 
                write(*,'(3a,i0)')'Electrode ', &
                trim(Elecs(iEl)%name),&
                ' does not reside within the defined Hartree box &
                &for the potential. Please ensure unit-cell &
                &vectors have the _correct_ direction: ',iia
            bool = .true.
          end if
        end do
      end do

      if ( bool .and. len_trim(Hartree_fname) == 0 ) then
        call die('ts-voltage: Check output, an electrode cannot be &
            &correctly applied a bias.')
      end if

      return

    end if

    ! set the left chemical potential
    call get_elec_indices(Elecs, cell, na_u, xa, iElL, iElR)

    ! this will be the applied bias in the "lower"
    ! electrode in the unit-cell
    V_high = Elecs(iElL)%mu%mu
    V_low = Elecs(iElR)%mu%mu

    ! Print out the coordinates of the ramp placement
    call print_ts_voltage()

  contains

    subroutine get_elec_indices(Elecs, cell, na_u, xa, iElL, iElR)
      use m_ts_electype, only : Elec_frac

      ! ***********************
      ! * INPUT variables     *
      ! ***********************
      type(Elec), intent(in) :: Elecs(2)
      real(dp), intent(in) :: cell(3,3)
      integer,  intent(in) :: na_u
      real(dp), intent(in) :: xa(3,na_u)

      ! ***********************
      ! * OUTPUT variables    *
      ! ***********************
      integer, intent(out) :: iElL, iElR

      ! ***********************
      ! * LOCAL variables     *
      ! ***********************
      integer  :: i
      real(dp) :: r1, r2

      call Elec_frac(Elecs(1), cell, na_u, xa, ts_tidx, fmin = r1)
      call Elec_frac(Elecs(2), cell, na_u, xa, ts_tidx, fmin = r2)

      if ( r1 < r2 ) then
        iElL = 1
        iElR = 2
      else
        iElL = 2
        iElR = 1
      end if

    end subroutine get_elec_indices

  end subroutine ts_init_voltage

  subroutine ts_voltage(cell, nmesh, nmeshl, Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : Hartree_fname
    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: nmesh(3), nmeshl(3)
    ! ***********************
    ! * OUTPUT variables    *
    ! ***********************
    real(grid_p), intent(inout) :: Vscf(:)

    call timer('ts_volt',1)

    ! Voltage drop in between the electrodes
    ! The indices for the full cell is set
    ! correctly to not have two routines doing the
    ! same
    if ( ts_tidx > 0 ) then
      call ts_ramp(cell, nmesh, nmeshl, Vscf)
#ifdef NCDF_4
    else if ( len_trim(Hartree_fname) > 0 ) then
      call ts_ncdf_Voltage(cell, Hartree_fname, 'V', nmesh, nmeshl, Vscf)
#endif
    else
      call ts_elec_only(nmesh, nmeshl, Vscf)
    end if

    call timer('ts_volt',2)

  end subroutine ts_voltage

  ! Apply a ramp to the electrostatic potential.
  subroutine ts_ramp(cell, nmesh, nmeshl, Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : Volt
    use m_mesh_node,  only : offset_i
    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: nmesh(3)
    integer,       intent(in) :: nmeshl(3)
    ! ***********************
    ! * OUTPUT variables    *
    ! ***********************
    real(grid_p), intent(inout) :: Vscf(:)

    ! ***********************
    ! * LOCAL variables     *
    ! ***********************
    integer  :: i1, i2, i3, idT, imesh
    real(dp) :: dF, dV

    ! field in [0;end]: v = e*x = f*index
    dF = - (V_high - V_low) / real(nmesh(ts_tidx) - 1,dp)

    ! Find quantities in mesh coordinates
    if ( product(nmeshl) /= size(Vscf) ) &
        call die('ramp_elec: Vscf size not correct')

    ! Add the electric field potential to the input potential
    imesh = 0
    if ( ts_tidx == 1 ) then
      
      do i3 = 1 , nmeshl(3)
        do i2 = 1 , nmeshl(2)
          do i1 = offset_i(1), offset_i(1)+nmeshl(1)-1
            
            dV = V_high + dF*i1
            
            imesh = imesh + 1
            Vscf(imesh) = Vscf(imesh) + dV
            
          end do
        end do
      end do
      
    else if ( ts_tidx == 2 ) then
      
      do i3 = 1,nmeshl(3)
        do i2 = offset_i(2),offset_i(2)+nmeshl(2)-1
          
          dV = V_high + dF*i2
          
          do i1 = 1 , nmeshl(1)
            imesh = imesh + 1
            Vscf(imesh) = Vscf(imesh) + dV
          end do
          
        end do
      end do
      
    else
      
      do i3 = offset_i(3),offset_i(3)+nmeshl(3)-1

        dV = V_high + dF*i3

        do i2 = 1 , nmeshl(2)
          do i1 = 1 , nmeshl(1)
            imesh = imesh + 1
            Vscf(imesh) = Vscf(imesh) + dV
          end do
        end do

      end do

    end if

  end subroutine ts_ramp


  subroutine ts_elec_only(nmesh, nmeshl, Vscf)
    use precision,    only : grid_p
#ifdef MPI
    use mpi_siesta
#endif

    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs
    use m_mesh_node,  only : dL
    use m_mesh_node,  only : dMesh, offset_r, offset_i, mesh_correct_idx
    use m_geom_box, only : voxel_in_box

    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    integer,       intent(in) :: nmesh(3), nmeshl(3)
    ! ***********************
    ! * OUTPUT variables    *
    ! ***********************
    real(grid_p), intent(inout) :: Vscf(:)

    ! ***********************
    ! * LOCAL variables     *
    ! ***********************
    integer  :: i1, i2, i3, iEl, imesh
    integer :: imin(3), imax(3), idx(3)
    real(dp) :: ll(3)
#ifdef TRANSIESTA_BOX
    integer, allocatable :: n_V(:)

    allocate(n_V(N_Elec))
    n_V = 0
#endif

    ! We do a loop in the local grid
    do iEl = 1 , N_Elec

      call Elec_box2grididx(Elecs(iEl),nmesh,dL,imin,imax)

      ! Now we have the minimum index for the box encompassing
      ! the electrode

      ! Loop the indices, and figure out whether
      ! each of them lies in the local grid
!$OMP parallel do default(shared), private(i1,i2,i3,ll,idx,imesh)
      do i3 = imin(3) , imax(3)
        do i2 = imin(2) , imax(2)
          do i1 = imin(1) , imax(1)

            ! Transform the index to the unit-cell index
            idx(1) = i1
            idx(2) = i2
            idx(3) = i3
            call mesh_correct_idx(nmesh,idx)
            idx = idx - offset_i - 1

            if ( idx(1) >= 0 .and. idx(2) >= 0 .and. idx(3) >= 0 .and. &
                idx(1) < nmeshl(1) .and. idx(2) < nmeshl(2) .and. &
                idx(3) < nmeshl(3) ) then

              ll(:) = offset_r(:) + &
                  idx(1)*dL(:,1) + idx(2)*dL(:,2) + idx(3)*dL(:,3)

              if ( voxel_in_box(Elecs(iEl)%box, ll, dMesh)) then
                imesh = idx(1) + (idx(2) + idx(3)*nmeshl(2))*nmeshl(1) + 1 ! to correct for -1
#ifdef TRANSIESTA_BOX
                n_V(iEl) = n_V(iEl) + 1
#endif
                Vscf(imesh) = Vscf(imesh) + Elecs(iEl)%mu%mu
              end if
            end if

          end do
        end do
      end do
!$OMP end parallel do

    end do

#ifdef TRANSIESTA_BOX
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place,n_V,N_Elec,MPI_Integer,MPI_Sum, &
        MPI_Comm_World,i1)
#endif
    print '(10(tr1,i8))',n_V,size(vscf)

    if ( any(n_V == 0) ) then
      iEl = minloc(n_V,1)
      write(*,'(3a)') 'ts-voltage: Elec-Box ',trim(Elecs(iEl)%name), &
          ' is not within the device grid box.'
      !       call die('ts-voltage: Elec-Box has an electrode outside &
      !            &the grid. Check TS output')
    end if

    deallocate(n_V)
#endif

  end subroutine ts_elec_only


  ! Print out the voltage direction dependent on the cell parameters.
  subroutine print_ts_voltage()
    use parallel,     only : IONode
    use units,        only : eV

    if ( .not. IONode ) return

    ! Print the ramp coordinates
    write(*,'(a,2(f6.3,tr1,a),a)') &
        'ts-voltage: Ramp ', V_high/eV, 'eV to ', V_low/eV, 'eV ', &
        'placed in cell'

  end subroutine print_ts_voltage

#ifdef NCDF_4
  ! Read in a potential file from a NetCDF file.
  ! Thus the potential landscape can be fully customized by the user.
  ! We note that the potential landscape need only be calculated
  ! for one V, direct interpolation is possible as 
  ! the solution to the Poisson equation is linearly dependent on the BC
  subroutine ts_ncdf_voltage(cell, fname, V_name, nmesh, nmeshl, V)
    use precision, only: grid_p
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Bcast, MPI_Grid_Real
#endif
    use m_ncdf_io, only : cdf_r_grid
    use netcdf_ncdf

#ifdef TRANSIESTA_VOLTAGE_DEBUG
    use iogrid_netcdf, only: write_grid_netcdf
#endif

    real(dp), intent(in) :: cell(3,3)
    character(len=*), intent(in) :: fname, V_name
    ! global and local number of mesh-divisions
    integer, intent(in) :: nmesh(3), nmeshl(3)
    real(grid_p), intent(inout) :: V(:)

    type(hNCDF) :: ncdf
    real(grid_p) :: Vmm(2), fact
    real(grid_p), allocatable :: tmpV(:)
#ifdef MPI
    integer :: MPIerror
#endif

    ! All nodes should allocate an auxilliary grid
    allocate(tmpV(product(nmeshl)))
    if ( size(V) /= size(tmpV) ) &
        call die('ncdf_voltage: Vscf size not correct')

    ! Open the file
    call ncdf_open(ncdf,fname, mode=NF90_NOWRITE)

    ! Read in the grid (should be in Ry)
    call cdf_r_grid(ncdf,trim(V_name),nmeshl,tmpV)

    ! retrieve the max min from the file (should be in Ry)
    call ncdf_get_var(ncdf,trim(V_name)//'min',Vmm(1))
    call ncdf_get_var(ncdf,trim(V_name)//'max',Vmm(2))
#ifdef MPI
    call MPI_Bcast(Vmm,2,MPI_Grid_Real,0,MPI_Comm_World,MPIerror)
#endif

    call ncdf_close(ncdf)

    ! Correct the limits so that we align to the current potential
    fact = ( V_high - V_low ) / ( Vmm(2) - Vmm(1) ) 
    ! Align the bottom potentials so that the range becomes correct
    Vmm(1) = V_low - Vmm(1) * fact

#ifdef TRANSIESTA_VOLTAGE_DEBUG
    tmpV(:) = tmpV(:) * fact + Vmm(1)
    call write_grid_netcdf( cell, nmesh, 1, product(nmeshl), tmpV, &
        "TransiestaHartreePotential")
    call bye('transiesta debug for Hartree potential')
#endif

    V(:) = Vmm(1) + V(:) + tmpV(:) * fact

    deallocate(tmpV)

  end subroutine ts_ncdf_voltage


  subroutine ts_ncdf_voltage_assert(fname, cell, nmesh)

    use parallel, only : IONode

    use netcdf_ncdf
    use dictionary

    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: nmesh(3)

    type(hNCDF) :: ncdf
    type(dictionary_t) :: dic
    character(len=3), parameter :: XYZ = 'ABC'
    integer :: lnmesh(3), i
    logical :: found

    if ( .not. IONode ) return

    ! Ensure that it will fail if not found
    lnmesh = -1
    call ncdf_open(ncdf, trim(fname), mode=NF90_NOWRITE)

    ! Get the grid-size, here we just fetch all dimensions
    ! and find them afterwards
    call ncdf_inq(ncdf,dict_dim=dic)

    ! We allow the names to be n1/x/a
    if ( 'a' .in. dic ) then
      call assign(lnmesh(1),dic,'a')
    else if ( 'n1' .in. dic ) then
      call assign(lnmesh(1),dic,'n1')
    else if ( 'x' .in. dic ) then
      call assign(lnmesh(1),dic,'x')
    end if

    ! We allow the names to be n2/y/b
    if ( 'b' .in. dic ) then
      call assign(lnmesh(2),dic,'b')
    else if ( 'n2' .in. dic ) then
      call assign(lnmesh(2),dic,'n2')
    else if ( 'y' .in. dic ) then
      call assign(lnmesh(2),dic,'y')
    end if

    ! We allow the names to be n3/z/c
    if ( 'c' .in. dic ) then
      call assign(lnmesh(3),dic,'c')
    else if ( 'n3' .in. dic ) then
      call assign(lnmesh(3),dic,'n3')
    else if ( 'z' .in. dic ) then
      call assign(lnmesh(3),dic,'z')
    end if

    ! Clean up
    call ncdf_close(ncdf)
    call delete(dic)

    ! Check variables
    if ( any(lnmesh /= nmesh) ) then
      
      write(*,*)
      write(*,'(a)') 'TS.Poisson file cannot be used!'
      write(*,'(a)') 'Please carefully read the below error message:'

      write(*,'(/,a)') 'TranSiesta internal grid dimensions are:'
      do i = 1, 3
        write(*,'(a,i0)') '  '//trim(XYZ(i:i))//': ', nmesh(i)
      end do

      write(*,'(/,3a)') 'File: ', trim(fname), ' grid dimensions are:'
      found = .true.
      do i = 1, 3
        write(*,'(a,i0)') '  '//trim(XYZ(i:i))//': ', lnmesh(i)
        found = found .and. lnmesh(i) > 0
      end do

      write(*,*)
      if ( found ) then
        write(*,'(a)') 'The grid dimensions *MUST* be equivalent.'
      else
        write(*,'(3a)') 'The lattice dimensions in ', trim(fname), ' *MUST* be named:'
        write(*,'(a)') '  1st lattice vector: a|n1|x'
        write(*,'(a)') '  2nd lattice vector: b|n2|y'
        write(*,'(a)') '  3rd lattice vector: c|n3|z'
      end if
      write(*,*)
        
      call die('Incorrect grid in NetCDF file for user supplied &
          &Poisson solution.')
      
    else
      write(*,'(a)') 'ts-voltage: User Poisson solution asserted...'
    end if
    
  end subroutine ts_ncdf_voltage_assert
  
#endif
  
end module m_ts_voltage
