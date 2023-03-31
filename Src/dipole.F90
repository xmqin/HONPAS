!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

! Created by Nick Papior in June 2020.

!> Module intented for handling dipoles in the system
!!
!! This module does 2 things:
!! 1. Returns the origin of a dipole which defaults to the
!!    center of the system in case the user does not know where
!!    it is.
!! 2. Calculate the dipole based on a given origin.
!!
!! This module is highly connected with the E-field module
module dipole_m

  use precision, only: dp

  implicit none

  private

  integer, parameter, public :: SLABDIPOLECORR_NONE = 0
  integer, parameter, public :: SLABDIPOLECORR_CHARGE = 1
  integer, parameter, public :: SLABDIPOLECORR_VACUUM = 2
  integer, public, save :: SLABDIPOLECORR = SLABDIPOLECORR_NONE

  public :: init_dipole_correction
  public :: get_dipole_origin
  public :: dipole_charge
  public :: dipole_potential

  public :: get_field_from_dipole
  public :: get_dipole_from_field

  type :: dipole_vacuum_dir_t
    !< Whether `xyz` is a fixed user defined coordinate in the unit cell
    !!
    !! If this is `.false.` it means that the coordinate will be
    !! in the vacuum region.
    logical :: use_xyz = .false.
    !< User defined coordinate, only used if `%use_xyz` is true.
    real(dp) :: xyz(3) = 0._dp
    !< Dipole direction
    real(dp) :: dir(3) = 0._dp
    !< dE tolerance for determining a flat dE/dr = constant, default 1e-4 eV/Ang/e
    real(dp) :: dE_tol = 3.88937962930779e-6_dp
  end type dipole_vacuum_dir_t

  !< Determine position for vacuum dipole
  type(dipole_vacuum_dir_t), public, save :: dip_vacuum

contains

  !< Read options for the dipole correction
  subroutine init_dipole_correction(nbcell, bcell, acting_efield)
    use parallel,     only: ionode
    use siesta_cml
    use fdf
    use units,        only: Ang, eV
    use m_cite, only: add_citation
    use m_char, only: lcase
    use intrinsic_missing, only: vnorm, vec_proj_sca

    !< number of bulk cell-directions
    integer, intent(in) :: nbcell
    real(dp), intent(in) :: bcell(3,nbcell)
    !< Whether there is an intrinsic electric field
    !!
    !! In case the user request a dipole correction,
    !! this will *always* be true on exit.
    logical, intent(inout) :: acting_efield

    real(dp) :: proj
    integer :: i
    character(len=32) :: dipole_str
    logical :: old_slabdipole

    ! Only allow slab dipole corrections for slabs and chains and molecules
    if ( acting_efield ) then
      dipole_str = lcase(fdf_get("Slab.DipoleCorrection", "charge"))
    else
      dipole_str = lcase(fdf_get("Slab.DipoleCorrection", "none"))
    end if

    select case ( trim(dipole_str) )
    case ( "true", "t", ".true.", "yes", "charge", "" )
      ! All these options reflect a dipole correction
      ! based on the original implementation
      ! One *may* use the Slab.DipoleCorrection.Origin block to specify the
      ! origin of the dipole
      ! An empty string is equivalent to a "true"
      SLABDIPOLECORR = SLABDIPOLECORR_CHARGE
    case ( "none", "no", "n", "false", "f", ".false." )
      SLABDIPOLECORR = SLABDIPOLECORR_NONE
    case ( "vacuum" )
      SLABDIPOLECORR = SLABDIPOLECORR_VACUUM
    case default
      write(*,'(2a)') "Dipole option: ", trim(dipole_str)
      call die("init_dipole: Could not determine the dipole option &
          &[charge, vacuum, none]")
    end select

    ! Further, remove dipole corrections if this is not a slab calculation
    ! No dipole corrections for molecules, chains or bulk materials.
    if ( nbcell /= 2 ) then
      SLABDIPOLECORR = SLABDIPOLECORR_NONE
    end if

    select case ( SLABDIPOLECORR )
    case ( SLABDIPOLECORR_NONE )
      if ( acting_efield .and. IONode ) then
        write(6,'(/,(a))') &
            'efield: WARNING!', &
            'efield: There is no slab-dipole correction [Slab.DipoleCorrection none] &
            &although an external efield is present.', &
            'efield: For correct physics Slab.DipoleCorrection should be .true.', &
            'efield: This is only for backwards compatibility!'
        write(6,*) ! newline
      end if

    case ( SLABDIPOLECORR_CHARGE )
      if (cml_p) &
          call cmlAddParameter( xf=mainXML, name='Slab.DipoleCorrection', &
          value="charge", &
          dictRef="siesta:slab.dipole_correction")
      if ( IONode ) then
        call add_citation("10.1103/PhysRevB.59.12301")
        write(6,'(/,(a))') &
            'dipole: A dipole layer will be introduced in the vacuum', &
            'dipole: region to compensate the system dipole'
      end if
      acting_efield = .true.

    case ( SLABDIPOLECORR_VACUUM )
      if (cml_p) &
          call cmlAddParameter( xf=mainXML, name='Slab.DipoleCorrection', &
          value="vacuum", &
          dictRef="siesta:slab.dipole_correction")

      ! Retrieve vacuum information
      call get_vacuum(nbcell, bcell, dip_vacuum)

      if ( IONode ) then
        write(6,'(/,(a))') &
            'dipole: A dipole layer will be introduced in the vacuum', &
            'dipole: region to compensate the vacuum field'
        write(6,'(a,3(tr1,f12.6))') "dipole: Dipole direction", dip_vacuum%dir
        write(6,'(a,tr1,f12.6)') "dipole: Field tolerance [eV/Ang/e]", dip_vacuum%dE_tol/eV*Ang
        if ( dip_vacuum%use_xyz ) then
          write(6,'(a,3(tr1,f12.3))') "dipole: Vacuum plane point [Ang]", dip_vacuum%xyz/Ang
        else
          write(6,'(a)') "dipole: Vacuum plane point automatically determined"
        end if
        write(6,*) ! newline
      end if

      ! Check if the direction is orthogonal to the bulk lattice vectors
      do i = 1, nbcell
        proj = VEC_PROJ_SCA(bcell(:,i), dip_vacuum%dir)
        if ( abs(proj) > 1.e-7_dp ) then
          call die("dipole: Dipole direction (vacuum plane normal vector) &
              &must be orthogonal to the bulk lattice vectors.")
        end if
      end do

      acting_efield = .true.
    end select

  end subroutine init_dipole_correction

  !< Return the dipole origin
  !!
  !! If the user does not specify the dipole origin
  !! the dipole origin will be the center of the system (not center-of-mass).
  !!
  !! However, in some skewed molecules where one knows the dipole is
  !! not located at the center of the system one should specify the
  !! origin of the dipole through this block: Slab.DipoleCorrection.Origin
  !!
  !! NOTE: This routine does not take into account shifts
  !! in the atomic structure when doing Grid.CellSampling
  !! I.e. if users want a specific coordinate for the dipole placement
  !! one should not *trust* the average dipole.
  subroutine get_dipole_origin(ucell, na_u, xa, x0)
    use sys, only: die
    use fdf

    integer, intent(in) :: na_u
    real(dp), intent(in) :: ucell(3,3), xa(3,na_u)
    real(dp), intent(out) :: x0(3)

    ! Local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=132) :: length_unit
    real(dp) :: cfactor
    integer :: ia

    x0(:) = 0.0_dp

    if ( fdf_block("Slab.DipoleCorrection.Origin", bfdf) ) then

      ! TODO consider allowing atomic indices
      ! as input. Then the origin could be calculated
      ! with respect to these coordinates. Here
      ! multiple values of the same atom *could* be allowed.
      ! Say a single atom above a hexagon?

      if ( .not. fdf_bline(bfdf,pline) ) then
        call die("Slab.DipoleCorrection.Origin could not find origin")
      end if

      if ( .not. fdf_bmatch(pline,"vvvn") ) then
        call die("Slab.DipoleCorrection.Origin both a coordinate and unit *must* be present")
      end if

      length_unit = fdf_bnames(pline, 1)
      cfactor = fdf_convfac(length_unit, "Bohr")

      ! Retrieve value
      ! This will not work for MD simulations (unless
      ! the dipole origin is fixed!)
      do ia = 1,3
        x0(ia) = fdf_bvalues(pline,ia) * cfactor
      end do

      call fdf_bclose(bfdf)

    else

      ! Default to the average molecule position
      do ia = 1, na_u
        x0(1:3) = x0(1:3) + xa(1:3,ia)
      end do
      ! Average
      do ia = 1, 3
        x0(ia) = x0(ia) / na_u
      end do

    end if

  end subroutine get_dipole_origin

  !< Return direction and coordinate for vacuum region
  subroutine get_vacuum(nbcell, bcell, dip_vacuum)
    use siesta_cml
    use sys, only: die
    use fdf
    use m_char, only: lcase
    use intrinsic_missing, only: vnorm

    integer, intent(in) :: nbcell
    real(dp), intent(in) :: bcell(3,nbcell)
    type(dipole_vacuum_dir_t), intent(inout) :: dip_vacuum

    ! Local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=64) :: unit, name
    logical :: found_block
    real(dp) :: cfactor
    integer :: ix
    logical :: found_direction
    external :: cross

    found_block = fdf_block("Slab.DipoleCorrection.Vacuum", bfdf)

    dip_vacuum%use_xyz = .false.
    found_direction = .false.

    if ( found_block ) then

      do while ( fdf_bline(bfdf, pline) )

        ! First string is the value
        name = lcase(fdf_bnames(pline, 1))

        select case ( trim(name) )
        case ( "direction" )

          found_direction = .true.

          do ix = 1 , 3
            dip_vacuum%dir(ix) = fdf_bvalues(pline,ix)
          end do

        case ( "point" )

          dip_vacuum%use_xyz = .true.

          unit = fdf_bnames(pline, 2)
          cfactor = fdf_convfac(unit, "Bohr")
          do ix = 1 , 3
            dip_vacuum%xyz(ix) = fdf_bvalues(pline,ix) * cfactor
          end do

        case ( "tolerance", "tol" )

          unit = fdf_bnames(pline, 2)
          cfactor = fdf_convfac(unit, "Ry/Bohr/e")
          dip_vacuum%dE_tol = fdf_bvalues(pline, 1) * cfactor

        end select

      end do

      call fdf_bclose(bfdf)

    end if

    if ( .not. found_direction ) then
      if ( nbcell == 2 ) then
        ! we can actually figure out the direction our selves
        call cross(bcell(1,1), bcell(1,2), dip_vacuum%dir)
        found_direction = .true.
      else
        call die("dipole: Slab.DipoleCorrection.Vacuum not present but necessary, &
            &please fix input.")
      end if
    end if

    ! Normalize direction
    dip_vacuum%dir(:) = dip_vacuum%dir(:) / vnorm(dip_vacuum%dir)

    ! Check if the direction has been specified
    if ( .not. found_direction ) then
      call die("dipole: Slab.DipoleCorrection.Vacuum did not contain direction")
    end if

    if ( cml_p ) then
      call cmlStartPropertyList(xf=mainXML, &
          dictRef='siesta:slab.dipole_vacuum', &
          title='Potential field calculation in vacuum')
      if ( dip_vacuum%use_xyz ) then
        call cmlAddProperty(xf=mainXML, &
            dictref='siesta:slab.dipole_vacuum_point', &
            title='Vacuum point', value="auto")
      else
        call cmlAddProperty(xf=mainXML, &
            dictref='siesta:slab.dipole_vacuum_point', &
            units="siestaUnits:Bohr", &
            title='Vacuum point', value=dip_vacuum%xyz)
      end if
      call cmlAddProperty(xf=mainXML, &
          dictref='siesta:slab.dipole_vacuum_direction', &
          title='Dipole direction', value=dip_vacuum%dir)
      call cmlAddProperty(xf=mainXML, &
          dictref='siesta:slab.dipole_vacuum_tolerance', &
          units='siestaUnits:Ry_Bohr_e', &
          title='Field tolerance for vacuum', value=dip_vacuum%dE_tol)
      call cmlEndPropertyList(xf=mainXML)
    end if

  end subroutine get_vacuum


  !< Calculate the field resulting from a dipole
  function get_field_from_dipole(dipole, cell) result(efield)
    use units, only: pi
    real(dp), intent(in)  :: dipole(3)
    real(dp), intent(in)  :: cell(3,3)
    real(dp) :: efield(3)

    real(dp), external :: volcel

    ! 4 * 2 * pi * dipole / volcel(cell)
    efield(1:3) = -8.0_dp * pi * dipole(1:3) / volcel(cell)

  end function get_field_from_dipole

  !< Calculate the dipole that creates a given field
  function get_dipole_from_field(efield, cell) result(dipole)
    use units, only: pi
    real(dp), intent(in)  :: efield(3)
    real(dp), intent(in)  :: cell(3,3)
    real(dp) :: dipole(3)
    real(dp), external :: volcel
    dipole(:) = - efield(:) * volcel(cell) / (8._dp * Pi)

  end function get_dipole_from_field

  !< Calculate dipole for a given origin of the dipole
  subroutine dipole_charge( cell, dvol, ntm, ntml, nsm, drho, X0, dipole)
! ********************************************************************
! Finds the electric dipole
! Written by J.M.Soler. July 1997.
! Modified for distributed drho matrix using a 2D grid of processors.
! Routine now is based on intrinsic structure of grid distribution
! in order to calculate position of grid points in space for local
! matrix. J.D.Gale March 1999.
! *********** INPUT **************************************************
! real*8  cell(3,3)     : Unit cell vectors
! real*8  dvol          : Voxel volume
! integer ntm(3)        : Global number of divisions of each lattice vector
! integer ntml(3)       : Local number of divisions of each lattice vector
! integer nsm           : Number of sub-points for each mesh point
! real*?  drho(*ntml)   : Minus neutral charge density at mesh points
! real*8  X0(3)         : Origin in cartesian coordinates of the dipole
! *********** OUTPUT *************************************************
! real*8 dipole(3)   : Electric dipole
! *********** UNITS **************************************************
! cell   in atomic units (Bohr)
! dvol   in atomic units (Bohr**3)
! drho   in atomic units (electrons/Bohr**3)
! X0     in atomic units (Bohr)
! dipole in atomic units (electrons*Bohr)
! ********************************************************************

    use precision,  only : grid_p
    use sys,        only : die
    use mesh,       only : meshLim
#ifdef MPI
    use mpi_siesta
#endif

    integer, intent(in) :: ntm(3), ntml(3), nsm
    real(grid_p), intent(in) :: drho(ntml(1),ntml(2),ntml(3))
    real(dp), intent(in) :: cell(3,3), dvol, X0(3)
    real(dp), intent(out) :: dipole(3)

    external :: reclat

    ! Internal variables and arrays
    integer :: I, I1, I2, I3, I10, I20, I30
    integer :: MG1, MG2, MG3, Xoffset, Yoffset, Zoffset
    real(dp) :: D(3), DD(3), DX(3,2:3), Rcell(3,3), X0L(3)
#ifdef MPI
    integer :: MPIerror
#endif

    ! Assign local variables
    MG1 = ntm(1)
    MG2 = ntm(2)
    MG3 = ntm(3)

    do I = 1 , 3
      DD(I) = 1._dp / ntm(I)
    end do

    ! Find reciprocal cell vectors (without the factor 2*pi)
    call reclat(cell, Rcell, 0 )

    ! Find origin in lattice coordinates
    do I = 1,3
      X0L(I) = X0(1)*Rcell(1,I) + X0(2)*Rcell(2,I) + X0(3)*Rcell(3,I)
    end do

    ! Initialize dipole
    dipole(1:3) = 0.0_dp

    ! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

    ! Find dipole by direct integration allowing for block distributed
    ! structure of drho
    I30 = Zoffset - 1
    do I3 = 1 , ntml(3)
      I30 = I30 + 1
      D(3) = DD(3) * I30 - X0L(3)
      if ( D(3) < - 0.5_dp ) then
        D(3) = D(3) + 1.0_dp
      else if ( D(3) > + 0.5_dp ) then
        D(3) = D(3) - 1.0_dp
      end if
      DX(1,3) = cell(1,3)*D(3)
      DX(2,3) = cell(2,3)*D(3)
      DX(3,3) = cell(3,3)*D(3)
      I20 = Yoffset - 1
      do I2 = 1 , ntml(2)
        I20 = I20 + 1
        D(2) = DD(2) * I20 - X0L(2)
        if ( D(2) < - 0.5_dp ) then
          D(2) = D(2) + 1.0_dp
        else if ( D(2) > + 0.5_dp ) then
          D(2) = D(2) - 1.0_dp
        end if
        DX(1,2) = cell(1,2)*D(2) + DX(1,3)
        DX(2,2) = cell(2,2)*D(2) + DX(2,3)
        DX(3,2) = cell(3,2)*D(2) + DX(3,3)
        I10 = Xoffset - 1
        do I1 = 1 , ntml(1)
          I10 = I10 + 1
          D(1) = DD(1) * I10 - X0L(1)
          if ( D(1) < - 0.5_dp ) then
            D(1) = D(1) + 1.0_dp
          else if ( D(1) > + 0.5_dp ) then
            D(1) = D(1) - 1.0_dp
          end if
          dipole(1) = dipole(1) - (cell(1,1)*D(1) + DX(1,2)) * drho(I1,I2,I3)
          dipole(2) = dipole(2) - (cell(2,1)*D(1) + DX(2,2)) * drho(I1,I2,I3)
          dipole(3) = dipole(3) - (cell(3,1)*D(1) + DX(3,2)) * drho(I1,I2,I3)
        end do
      end do
    end do

    ! Correct for volume
    dipole(1) = dipole(1) * dvol
    dipole(2) = dipole(2) * dvol
    dipole(3) = dipole(3) * dvol

#ifdef MPI
    call MPI_AllReduce(dipole, d, 3, MPI_double_precision, MPI_sum, &
        MPI_Comm_World, MPIerror)
    dipole(1:3) = d(1:3)
#endif

  end subroutine dipole_charge


  subroutine dipole_potential(cell, na_u, isa, xa, ntm, ntml, nsm, dip_vacuum, V, dipole)
! ********************************************************************
! Calculates the dipole field in the Vacuum region
! This routine only works for FFT calculated Hartree potentials
! since any dipoles will induce a compensating field in the solution.
! This method will calculate the resulting dipole by averaging in a vacuum
! region where the double derivative is below a given tolerance
! Written by N.Papior. June 2020.
! Requires distributed V matrix using a 2D grid of processors.
! *********** INPUT **************************************************
! real*8  cell(3,3)     : Unit cell vectors
! integer na_u          : Atoms in unit cell
! integer isa(na_u)     : Species index for atom
! real*8  xa(3,na_u)    : Atomic coordinates
! integer ntm(3)        : Global number of divisions of each lattice vector
! integer ntml(3)       : Local number of divisions of each lattice vector
! integer nsm           : Number of sub-points for each mesh point
! dipole_vacuum_dir_t
!     dip_vacuum        : type containing vacuum information (dipole direction)
!                         and the tolerance for the electric field
! real*?  V(*ntml)      : Hartree potential as calculated from FFT
! *********** OUTPUT *************************************************
! real*8 dipole(3)      : A global dipole that creates a potential ramp
! *********** UNITS **************************************************
! cell, xa  in atomic units (Bohr)
! dE_tol  in atomic units (Ry/Bohr/e) (electric field change)
! V     in atomic units (Ry)
! dipole in atomic units (Ry*Bohr)
! ********************************************************************
!
!  Modules
!
    use precision,  only : grid_p
    use parallel,   only : Node
    use sys,        only : die
    use mesh,       only : meshLim
    use atmfuncs,     only : rcut
#ifdef MPI
    use mpi_siesta
#endif
    use units, only: Pi, eV, Ang
    use intrinsic_missing, only: MODULOP, vnorm

    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u, isa(na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: ntm(3), ntml(3), nsm
    type(dipole_vacuum_dir_t), intent(in) :: dip_vacuum
    real(grid_p), intent(in) :: V(ntml(1),ntml(2),ntml(3))
    real(dp), intent(out) :: dipole(3)

    real(dp), external :: volcel
    external :: reclat

    ! Internal variables and arrays
    real(dp) :: rc, rc_dplane, xfrac, field, rcell(3,3)
    real(dp) :: xmin, xmax, dplane
    logical :: valid_x0
    integer :: ix, i0, is, ia, idir
    integer :: i1, i2, i3, i, offset
    real(dp), allocatable :: Vavg(:), E(:), dE(:)

#ifdef MPI
    integer :: MPIerror
#endif

    ! Figure out which lattice vector dir corresponds to
    idir = 0
    do ix = 1, 3
      xfrac = dot_product(cell(:,ix), dip_vacuum%dir) / vnorm(cell(:,ix))
      if ( abs(xfrac) > 0.001_dp ) then
        ! We are having at least a projection of 0.1% onto this
        ! lattice vector
        if ( idir == 0 ) then
          idir = ix
        else
          print *, "Direction 1: ", idir
          print *, "Direction 2: ", ix
          call die("dipole_from_potential: found multiple directions of &
              &parallel lattice vectors to the dipole direction, this is not implemented.")
        end if
      end if
    end do
    if ( idir == 0 ) then
      call die("dipole_from_potential: no direction was found parallel to the dipole &
          &direction, how is this possible?")
    end if

    ! Find the origin of a shifted cell, with the system centered in it
    ! This is done at every call, because of possible atomic movements

    ! Find reciprocal unit cell and distance between lattice planes
    call reclat(cell, rcell, 0)

    dplane = 1.0_dp / sqrt(dot_product(rcell(:,idir), rcell(:,idir)))
    xmin = 1.0e30_dp
    xmax = -1.0e30_dp
    do ia = 1 , na_u
      is = isa(ia)
      rc = rcut(is,0)
      rc_dplane = rc / dplane
      ! Only take into account the dipole direction
      xfrac = dot_product(xa(:,ia), rcell(:,idir))
      xmin = min( xmin, xfrac - rc_dplane )
      xmax = max( xmax, xfrac + rc_dplane )
    end do
    ! Now we have min and max along the dipole direction

    if ( dip_vacuum%use_xyz ) then
      ! Find vacuum position
      rc = dot_product(dip_vacuum%xyz, rcell(:,idir))
      rc = modulo(rc, 1._dp)
      if ( xmin <= rc .and. rc <= xmax ) then
        call die("dipole_potential: Vacuum point lies in atomic range.")
      end if
    else
      rc = (xmin + xmax) / 2 - 0.5_dp
    end if
    i0 = MODULOP(nint(rc*ntm(idir)), ntm(idir))

    ! Initialize calculated field
    ! We will finally convert to dipole
    field = 0.0_dp
    dipole(:) = 0.0_dp

    ! Some processor *must* have the contribution
    ! distribute ghost-layers

    allocate(Vavg(ntm(idir)))
    Vavg(:) = 0._dp

    offset = (meshLim(1,idir)-1)*nsm

    ! Calculate average for all planes with normal along dir
    select case ( idir )
    case ( 1 )
      ! Calculate plane average constant
      rc = 1._dp / (ntm(2) * ntm(3))
      do i3 = 1 , ntml(3)
        do i2 = 1 , ntml(2)
          do i1 = 1 , ntml(1)
            i = offset + i1
            Vavg(i) = Vavg(i) + V(i1,i2,i3)
          end do
        end do
      end do
    case ( 2 )
      rc = 1._dp / (ntm(1) * ntm(3))
      do i3 = 1 , ntml(3)
        do i2 = 1 , ntml(2)
          i = offset + i2
          do i1 = 1 , ntml(1)
            Vavg(i) = Vavg(i) + V(i1,i2,i3)
          end do
        end do
      end do
    case ( 3 )
      rc = 1._dp / (ntm(1) * ntm(2))
      do i3 = 1 , ntml(3)
        i = offset + i3
        do i2 = 1 , ntml(2)
          do i1 = 1 , ntml(1)
            Vavg(i) = Vavg(i) + V(i1,i2,i3)
          end do
        end do
      end do
    end select

    ! We will do two averages
    ! 1. number of elements in plane (current rc, calculated above)
    ! 2. Reduce calculations in next step by adding voxel length along
    !    dipole direction (corrected for here)
    !    This step removes the rc factor in E(:) calculation
    rc = ntm(idir) * rc / VNORM(cell(:,idir))
    do i = offset + 1 , offset + ntml(idir)
      Vavg(i) = Vavg(i) * rc
    end do

    ! Now perform a reduction
    ! There is no need to do this on any other than 1 core
#ifdef MPI
    allocate(E(ntm(idir)))
    E(:) = Vavg(:)
    call MPI_Reduce(E, Vavg, ntm(idir), MPI_double_precision, MPI_sum, 0, &
        MPI_Comm_World, MPIerror)
    deallocate(E)
#endif

    ! Now perform average calculation of field along idir
    ! This is a 3-point stencil
    if ( Node == 0 ) then

      allocate(E(ntm(idir)), dE(ntm(idir)))
      E(1) = 0.5_dp * (Vavg(2) - Vavg(ntm(idir)))
      dE(1) = Vavg(2) + Vavg(ntm(idir)) - 2 * Vavg(1)
      do i = 2, ntm(idir) - 1
        E(i) = 0.5_dp * (Vavg(i+1) - Vavg(i-1))
        dE(i) = Vavg(i+1) + Vavg(i-1) - 2 * Vavg(i)
      end do
      E(ntm(idir)) = 0.5_dp * (Vavg(1) - Vavg(ntm(idir)-1))
      dE(ntm(idir)) = Vavg(1) + Vavg(ntm(idir)-1) - 2 * Vavg(ntm(idir))

      !! DEBUG
      ! write out table with coordinate E and dE
      !rc = vnorm(cell(:,idir)) / ntm(idir)
      !do i = 1 , ntm(idir)
      !  write(8639,'(3(tr1,f14.6))') i * rc, E(i), dE(i)
      !end do

      ! Now we have everything
      ! Initialize field sign at i0(idir)
      i = i0
      rc = E(i)
      field = abs(E(i))

      ! step back until we hit a non-valid point
      do while ( abs(dE(i)) <= dip_vacuum%dE_tol )
        i = i - 1
        if ( i < 1 ) i = ntm(idir)
        if ( i == i0 ) then
          call die("dipole_from_potential: dE tolerance too high, all points valid")
        end if
      end do
      ! step one line back
      ! this is ok even if dE(i0) > dip_vacuum%dE_tol
      ! since we check at the final stage whether it is ok
      i = i + 1

      ! Find minimum field value in the vacuum
      ! The reason for the minimum value is that the actual
      ! vacuum is *far* from the system and in this case
      ! for a FFT solved Hartree potential there will be
      ! a local minima in the vacuum point.
      ! The tolerance is just to cutoff regions of space
      ! where we have charge.
      i1 = 0
      do while ( abs(dE(i)) <= dip_vacuum%dE_tol )
        field = min(field, abs(E(i)))
        i = i + 1
        i1 = i1 + 1
        ! do it round-robin
        if ( i > ntm(idir) ) i = 1
      end do

      ! Get correct sign of field
      field = sign(field, rc)

      ! We shouldn't use the negative sign here
      ! since the dipole will be reverted in dhscf
      dipole(idir) = field * volcel(cell) / (Pi * 2 * 4)

      deallocate(E, dE)
    end if

    deallocate(Vavg)

#ifdef MPI
    ! Bcast result from root node
    call MPI_Bcast(dipole(idir),1,Mpi_Double_Precision, 0, &
        MPI_Comm_World, MPIerror)
#endif

  end subroutine dipole_potential

end module dipole_m
