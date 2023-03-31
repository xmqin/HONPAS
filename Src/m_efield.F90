! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_efield

  ! This module implements routines to deal with an external
  ! electric field.
  !
  ! The field required can be specified in two possible ways:
  !
  ! 1. Through and FDF block 'ExternalElectricField'. For example:
  !
  !   %block ExternalElectricField
  !     0.000  0.000  3.000  V/Ang
  !   %endblock ExternalElectricField
  
  ! 2. Through the FDF Option 'SlabDipoleCorrection'. 
  !    If 'true', the program will then calculate at every SCF step
  !    the electric field required to compensate the dipole of the
  !    system. The potential added to the grid corresponds to that
  !    of a dipole layer in the middle of the vacuum layer. For slabs,
  !    this exactly compensates the electric field at the vacuum created
  !    by the dipole moment of the system, thus allowing to deal with
  !    asymmetric slabs (and compute properties such as the work funcion
  !    of each of the surfaces).
  !    See L. Bengtsson, PRB 59, 12301 (1999), DOI: 10.1103/PhysRevB.59.12301
  
  use precision, only: dp
  use sys,       only: die

  implicit none

  logical, public  :: acting_efield = .false.
  real(dp), public :: user_specified_field(3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

  private
  public :: initialize_efield
  public :: add_potential_from_field

contains

  subroutine get_user_specified_field(input_field)

    use fdf

    real(dp), intent(out) :: input_field(3)

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    character(len=132) :: eunits

    real(dp) :: cfactor
    integer :: ix

    input_field(1:3) = 0.0_dp

    if (.not. fdf_block('ExternalElectricField',bfdf) ) RETURN
    do while ( fdf_bline(bfdf,pline) )
      if (.not. fdf_bmatch(pline,"vvvn")) &
          call die("Wrong format in ElectricField block")
      eunits = fdf_bnames(pline,1)
      cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
      do ix = 1,3
        input_field(ix) = fdf_bvalues(pline,ix) * cfactor
      end do
    end do

    call fdf_bclose(bfdf)

  end subroutine get_user_specified_field

  !< Orthogonalize a given field for the periodic directions
  subroutine orthogonalize_efield(input_field,orthog_field,orthog,nbcell)

    use siesta_geom, only : ucell, xa, na_u, isa

    real(dp), intent(in)  :: input_field(3)
    real(dp), intent(out) :: orthog_field(3)
    logical,  intent(out) :: orthog ! originally orthogonal?
    integer,  intent(out) :: nbcell ! shape of cell

    ! tolerance for bulk components of the electric field
    real(dp), parameter :: tol = 1.0d-12

    real(dp) :: eb1, eb2, eb3, b1xb2(3), bcell(3,3)
    integer :: ix
    character(len=8) :: shape

    real(dp), external :: ddot

    call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
    orthog = .true.

    select case ( nbcell )
    case ( 1 )
      eb1 = ddot(3,input_field,1,bcell,1) / ddot(3,bcell,1,bcell,1)
      if (abs(eb1) .gt. tol) then
        orthog = .false.
        do ix = 1,3
          orthog_field(ix) = input_field(ix) - eb1 * bcell(ix,1)
        enddo
      endif

    case ( 2 )
      eb1 = ddot(3,input_field,1,bcell(1,1),1)/ &
          ddot(3,bcell(1,1),1,bcell(1,1),1)
      eb2 = ddot(3,input_field,1,bcell(1,2),1)/ &
          ddot(3,bcell(1,2),1,bcell(1,2),1)
      if ( abs(eb1) > tol .or. abs(eb2) > tol ) then
        orthog = .false.
        call cross( bcell(1,1), bcell(1,2), b1xb2 )
        eb3 = ddot(3,input_field,1,b1xb2,1)/ddot(3,b1xb2,1,b1xb2,1)
        do ix = 1,3
          orthog_field(ix) = eb3 * b1xb2(ix)
        end do
      end if

    case ( 3 )
      orthog = .false.
      do ix = 1,3
        orthog_field(ix) = 0.0_dp
      end do

    end select

  end subroutine orthogonalize_efield

  !< Initialize electric-field options
  !!
  !! It sets the module variables
  !! - `user_specified_field`
  !! - `acting_efield`
  !! - `dipole_correction`
  subroutine initialize_efield()

    use parallel,     only: ionode
    use siesta_cml,   only: cml_p, cmlAddProperty, mainXML
    use fdf
    use units,        only: Ang, eV
    use intrinsic_missing, only: vnorm

    real(dp) :: input_field(3), orthog_field(3)
    logical  :: orthog
    integer :: nbcell

    call get_user_specified_field(input_field)
    acting_efield = vnorm(input_field) /= 0.0_dp

    ! always orthogonalize retrieve the shape of the cell
    call orthogonalize_efield(input_field,orthog_field,orthog,nbcell)
      
    if ( acting_efield ) then
      if ( orthog ) then
        if ( ionode ) then
          write(6,'(/,a,3f12.6,a)') &
              'efield: External electric field =', &
              input_field/eV*Ang, ' eV/Ang/e'
        end if
        user_specified_field(1:3) = input_field(1:3)
      else
        if ( ionode ) then
          write(6,'(a,(/,a,3f12.6))') &
              'efield: WARNING: Non zero bulk electric field.', &
              'efield: Input field (eV/Ang/e) =', input_field/eV*Ang, &
              'efield: Orthogonalized field   =', orthog_field/eV*Ang
        end if
        user_specified_field(1:3) = orthog_field(1:3)
      end if
      if (cml_p) &
          call cmlAddProperty(xf=mainXML, &
          value=user_specified_field, &
          dictref='siesta:elfield', &
          units='siestaUnits:Ry_Bohr_e')
    end if

    ! finally determine whether we do have a field
    acting_efield = vnorm(user_specified_field) /= 0.0_dp

  end subroutine initialize_efield

  subroutine add_potential_from_field(efield, cell, na, isa, xa, ntm, ntml, nsm, V)

! Adds the potential created by an external electric field.
! Written by J.M.Soler. Feb. 1998.
! Modified to operate only on the sub-matrix of the potential stored
! locally. J.D.Gale March 1999.
! Modularized by A. Garcia, Nov. 2009
! ********* Input ******************************************************
! real*8  cell(3,3) : Unit cell vectors
! integer na        : Number of atoms
! integer isa(na)   : Atomic species indexes
! real*8  xa(3,na)  : Atomic positions (cartesian coordinates)
! integer ntm(3)    : Number of mesh divisions in each cell direction
! integer ntml(3)   : Number of local mesh divisions in each cell direction
! integer nsm       : Number of sub-mesh points along each axis
! real*8  field(3)  : Electric field
! ********* Input and output *******************************************
! real*?    v(*ntml) : Electron potential, to which that created by the
! electric field is added
! ********* Units ******************************************************
! Distances in Bohr
! Energies in Rydbergs
! Electric field in Ry/Bohr
! Dipoles in electrons*Bohr
! ********* Behaviour **************************************************
! The sign of the potential is that for electrons (v=+E*x), i.e. 
! opposite to that of the conventional electrostatic potential.
! Notice that the potential is not initialized.
! Bulk electric fields are not allowed. If the specified electric field
! is not orthogonal to all bulk directions, it is orthogonalized, and
! a warning message is printed.
! The electric field produces a discontinuity of the potential in the
! periodic cell, which is automatically placed in the middle of the
! vacuum region.

    use precision, only: grid_p
    use atmfuncs,     only : rcut
    use mesh,         only : meshLim

    integer, intent(in)          ::  na, isa(na), ntm(3), ntml(3), nsm
    real(dp), intent(in)         ::  cell(3,3), efield(3), xa(3,na)
    real(grid_p), intent(inout)  ::  v(ntml(1),ntml(2),ntml(3))

    integer :: i0(3), i1, i2, i3, ia
    integer :: is, ix
    integer :: j1, j2, j3
    integer :: Xoffset, Yoffset, Zoffset, i30, i20, i10

    real(dp) :: dplane(3), f(3), rc, rcell(3,3), v0
    real(dp) :: xfrac, xmax(3), xmean, xmin(3)

    real(dp), external :: ddot
    external :: reclat

    ! Find the origin of a shifted cell, with the system centered in it
    ! This is done at every call, because of possible atomic movements

    ! Find reciprocal unit cell and distance between lattice planes
    call reclat( cell, rcell, 0 )

    ! Find the geometric center of the system
    do ix = 1,3
      dplane(ix) = 1._dp / sqrt(ddot(3,rcell(1,ix),1,rcell(1,ix),1))
      xmin(ix) =  1.0e30_dp
      xmax(ix) = -1.0e30_dp
    end do
    do ia = 1,na
      is = isa(ia)
      rc = rcut(is,0)
      do ix = 1,3
        xfrac = ddot(3,xa(1,ia),1,rcell(1,ix),1)
        xmin(ix) = min( xmin(ix), xfrac-rc/dplane(ix) )
        xmax(ix) = max( xmax(ix), xfrac+rc/dplane(ix) )
      end do
    end do

    ! Find the mesh index of the origin of the shifted cell
    ! Ensure that all indices become positive in the below loop
    do ix = 1,3
      xmean = (xmin(ix) + xmax(ix)) / 2
      i0(ix) = nint( (xmean-0.5_dp) * ntm(ix) ) - 10 * ntm(ix)

      ! Find the electric field in mesh coordinates, so that
      ! v = efield*x = f*index
      f(ix) = ddot(3,efield,1,cell(1,ix),1) / max( ntm(ix), 1 )
    end do

    ! Find the potential at the origin of the shifted cell, so that
    ! the potential is zero at the center of the cell
    v0 = - 0.5_dp*(f(1)*ntm(1) + f(2)*ntm(2) + f(3)*ntm(3))

    ! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

    ! Add the electric field potential to the input potential
    i30 = Zoffset - 1
    do i3 = 1 , ntml(3)
      i30 = i30 + 1
      j3 = mod( i30-i0(3), ntm(3) )
      i20 = Yoffset - 1
      do i2 = 1 , ntml(2)
        i20 = i20 + 1
        j2 = mod( i20-i0(2), ntm(2) )
        i10 = Xoffset - 1
        do i1 = 1 , ntml(1)
          i10 = i10 + 1
          j1 = mod( i10-i0(1), ntm(1) )
          v(i1,i2,i3) = v(i1,i2,i3) + v0 + f(1)*j1 + f(2)*j2 + f(3)*j3
        enddo
      enddo
    enddo

  end subroutine add_potential_from_field

end module m_efield
