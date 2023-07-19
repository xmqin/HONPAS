! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_efield

!     This module implements routines to deal with an external
!     electric field.
!
!     The field required can be specified in two possible ways:
!
!     1. Through and FDF block 'ExternalElectricField'. For example:
!
!       %block ExternalElectricField
!           0.000  0.000  3.000  V/Ang
!       %endblock ExternalElectricField

!    2. Through the FDF Option 'SlabDipoleCorrection'. 
!    If 'true', the program will then calculate at every SCF step
!    the electric field required to compensate the dipole of the
!    system. The potential added to the grid corresponds to that
!    of a dipole layer in the middle of the vacuum layer. For slabs,
!    this exactly compensates the electric field at the vacuum created
!    by the dipole moment of the system, thus allowing to deal with
!    asymmetric slabs (and compute properties such as the work funcion
!    of each of the surfaces).

      use precision, only: dp
      use sys,       only: die

      implicit none
      logical, public  :: acting_efield = .false.
      logical, public  :: dipole_correction = .false.
      real(dp), public :: user_specified_field(3) =
     $                                (/ 0.0_dp, 0.0_dp, 0.0_dp /)

      public   :: initialize_efield
      public   :: get_field_from_dipole, add_potential_from_field

      private

      CONTAINS

      subroutine get_user_specified_field(input_field)

      use fdf,       only: fdf_convfac, fdf_block, fdf_bline, block
      use parse,     only: parsed_line, digest, match
      use parse,     only: names, values

      use parallel,  only: ionode
      use m_mpi_utils, only: broadcast

      real(dp),   intent(out)     :: input_field(3)

      logical                     :: found
      type(block), pointer        :: bp => null()
      type(parsed_line), pointer  :: p  => null()
      character(len=132)          :: line, eunits

      real(dp)                    :: cfactor
      integer                     :: ix

      input_field(1:3) = 0.0_dp

      if (ionode) then
       found = fdf_block('ExternalElectricField',bp)
       if (found) then
         loop: DO
           if (.not. fdf_bline(bp,line)) exit loop
           p => digest(line)
           if (.not. match(p,"vvvn"))
     $       call die("Wrong format in ElectricField block")
           eunits = names(p,1)
           cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
           do ix = 1,3
              input_field(ix) = values(p,ix) * cfactor
           enddo
         enddo loop
       endif ! (found)
      endif ! (ionode)
      call broadcast(input_field)

      end subroutine get_user_specified_field
!---------------------------------------------------------------
       
      subroutine orthogonalize_efield(input_field,orthog_field,orthog)

!     Orthogonalizes the electric field to the bulk directions

      use siesta_geom, only : ucell, xa, na_u, isa

      real(dp), intent(in)  :: input_field(3)
      real(dp), intent(out) :: orthog_field(3)
      logical,  intent(out) :: orthog           ! originally orthogonal?

!     tolerance for bulk components of the electric field
      real(dp), parameter    :: tol = 1.0e-12_dp

      real(dp)  :: eb1, eb2, eb3, b1xb2(3), bcell(3,3)
      integer   :: nbcell, ix
      character(len=8) :: shape

      real(dp), external :: ddot

      call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
      orthog = .true.
      if (nbcell .eq. 1) then
         eb1 = ddot(3,input_field,1,bcell,1) / ddot(3,bcell,1,bcell,1)
         if (abs(eb1) .gt. tol) then
            orthog = .false.
            do ix = 1,3
               orthog_field(ix) = input_field(ix) - eb1 * bcell(ix,1)
            enddo
         endif
      elseif (nbcell .eq. 2) then
         eb1 = ddot(3,input_field,1,bcell(1,1),1)/
     .        ddot(3,bcell(1,1),1,bcell(1,1),1)
         eb2 = ddot(3,input_field,1,bcell(1,2),1)/
     .        ddot(3,bcell(1,2),1,bcell(1,2),1)
         if (abs(eb1).gt.tol .or. abs(eb2).gt.tol) then
            orthog = .false.
            call cross( bcell(1,1), bcell(1,2), b1xb2 )
            eb3 = ddot(3,input_field,1,b1xb2,1)/ddot(3,b1xb2,1,b1xb2,1)
            do ix = 1,3
               orthog_field(ix) = eb3 * b1xb2(ix)
            enddo
         endif
      elseif (nbcell .eq. 3) then
         orthog = .false.
         do ix = 1,3
            orthog_field(ix) = 0.0_dp
         enddo
      endif
      end subroutine orthogonalize_efield
!---------------------------------------------------------------

      subroutine initialize_efield()
!
!     Initializations that can be done once and for all
!     It sets the module variables:
!          user_specified_field
!          acting_efield
!          dipole_correction
!
      use parallel,     only: ionode
      use siesta_cml,   only: cml_p, cmlAddProperty, mainXML
      use m_fdf_global, only: fdf_global_get

      real(dp) :: input_field(3), orthog_field(3)
      logical  :: orthog

      real(dp), external :: ddot

      call get_user_specified_field(input_field)
      acting_efield = (ddot(3,input_field,1,input_field,1) /= 0.0_dp)
      if (acting_efield) then
         call orthogonalize_efield(input_field,orthog_field,orthog)

         if (orthog) then
            if (ionode) then
               write(6,'(/,a,3f12.6,a)')
     $              'efield: External electric field =',
     $              input_field, ' Ry/Bohr/e'
            endif
            user_specified_field(1:3) = input_field(1:3)
         else
            if (ionode) then
               write(6,'(a,(/,a,3f12.6))')
     .              'efield: WARNING: Non zero bulk electric field.',
     .              'efield: Input field (Ry/Bohr/e) =', input_field,
     .              'efield: Orthogonalized field    =', orthog_field
            endif
            user_specified_field(1:3) = orthog_field(1:3)
         endif
         if (cml_p) call cmlAddProperty(xf=mainXML,
     $           value=user_specified_field, 
     .           dictref='siesta:elfield', 
     .           units='siestaUnits:Ry_Bohr_e')
      endif
!
      acting_efield =
     $ (ddot(3,user_specified_field,1,user_specified_field,1) /= 0.0_dp)

      call fdf_global_get(dipole_correction,
     $                    'SlabDipoleCorrection',.false.)

      if (dipole_correction) then
         if (ionode) write(6,'(/,(a))')
     .        'efield: SlabDipoleCorrection = .true. in input file',
     .        'efield: A dipole layer will be introduced in the vacuum',
     .        'efield: region to compensate the system dipole'
         acting_efield = .true.
      endif

      end subroutine initialize_efield
!---------------------------------------------------------------

      function get_field_from_dipole(dipole, cell) result(efield)

      use units,       only: pi

      real(dp), intent(in)  :: dipole(3)
      real(dp), intent(in)  :: cell(3,3)
      real(dp)              :: efield(3)

      real(dp), external  :: volcel

      efield(1:3) = -4.0_dp * pi * dipole(1:3) * 2.0_dp / volcel(cell)

      end function get_field_from_dipole
!---------------------------------------------------------------

      subroutine add_potential_from_field(efield, cell, na, isa, xa,
     $                                    mesh, nsm, V)

!     Adds the potential created by an external electric field.
!     Written by J.M.Soler. Feb. 1998.
!     Modified to operate only on the sub-matrix of the potential stored
!     locally. J.D.Gale March 1999.
!     Modularized by A. Garcia, Nov 2009.
!
! Distances in Bohr
! Energies in Rydbergs
! Electric field in Ry/Bohr
! Dipoles in electrons*Bohr
!
! The sign of the potential is that for electrons (v=+E*x), i.e. 
! opposite to that of the conventional electrostatic potential.
! Notice that the potential is not initialized.
! Bulk electric fields are not allowed. If the specified electric field
! is not orthogonal to all bulk directions, it is orthogonalized, and
! a warning message is printed.
! The electric field produces a discontinuity of the potential in the
! periodic cell, which is automatically placed in the middle of the
! vacuum region.

      use precision, only: dp, grid_p
      use parallel,     only : Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use atmfuncs,     only : rcut
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

!     Electric field
      real(dp), intent(in)        :: efield(3)
!     Number of atoms, species indexes
      integer, intent(in)          ::  na, isa(na)
!     Number of mesh divisions in each cell direction
!     Number of sub-mesh points along each axis
      integer, intent(in)          ::  mesh(3), nsm
!     Unit cell and atomic positions
      real(dp), intent(in)         ::  cell(3,3), xa(3,na)
!     Electron potential, to which that created by the
!     electric field is added
      real(grid_p), intent(inout)  ::  v(*)

      external          reclat

      integer           i0(3), i1, i2, i3, ia, imesh,
     $                  is, iu, ix, j1, j2, j3, last, 
     .                  nbcell, ni, nn, nr, nv, meshl(3),
     .                  ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset,
     .                  Zoffset, Py, Pz, i30, i20, meshnsm(3),
     .                  NRemY, NRemZ
#ifdef MPI
      integer           MPIerror, npl
#endif
      real(dp)          dplane(3), f(3), pi, rc, rcell(3,3), v0,
     .                  xfrac, xmax(3), xmean, xmin(3)

      real(dp), external :: ddot

! Find the origin of a shifted cell, with the system centered in it
! This is done at every call, because of possible atomic movements

! Find reciprocal unit cell and distance between lattice planes
      call reclat( cell, rcell, 0 )
      do ix = 1,3
         dplane(ix) = 1.0_dp /
     $        sqrt(ddot(3,rcell(1,ix),1,rcell(1,ix),1))
      enddo

! Find the geometric center of the system
      do ix = 1,3
         xmin(ix) =  1.0e30_dp
         xmax(ix) = -1.0e30_dp
      enddo
      do ia = 1,na
         is = isa(ia)
         rc = rcut(is,0)
         do ix = 1,3
            xfrac = ddot(3,xa(1,ia),1,rcell(1,ix),1)
            xmin(ix) = min( xmin(ix), xfrac-rc/dplane(ix) )
            xmax(ix) = max( xmax(ix), xfrac+rc/dplane(ix) )
         enddo
      enddo

! Find the mesh index of the origin of the shifted cell
      do ix = 1,3
         xmean = (xmin(ix) + xmax(ix)) / 2
         i0(ix) = nint( (xmean-0.5_dp) * mesh(ix) )
      enddo

! Find the electric field in mesh coordinates, so that
! v = efield*x = f*index
      do ix = 1,3
         f(ix) = ddot(3,efield,1,cell(1,ix),1) / max( mesh(ix), 1 )
      enddo

! Find the potential at the origin of the shifted cell, so that
! the potential is zero at the center of the cell
      v0 = (- 0.5_dp) * (f(1)*mesh(1) + f(2)*mesh(2) + f(3)*mesh(3))

! Find local number of mesh points
      meshnsm(1:3) = mesh(1:3)/nsm
#ifdef MPI
      call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
      meshl(1:3) = meshl(1:3)*nsm
#else
      meshl(1:3) = mesh(1:3)
#endif

! Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      ProcessorZ = Nodes/ProcessorY

! Calculate blocking sizes
      BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
      NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
      BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
      NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

! Calculate coordinates of current node in processor grid
      Py = (Node/ProcessorZ)+1
      Pz = Node - (Py - 1)*ProcessorZ + 1

! Calculate starting point for grid
      Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
      Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

! Add the electric field potential to the input potential
      imesh = 0
      i30 = Zoffset - 1
      do i3 = 0,meshl(3)-1
         i30 = i30 + 1
         i20 = Yoffset - 1
         do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            do i1 = 0,meshl(1)-1
               imesh = imesh + 1
               j1 = mod( i1-i0(1)+10*mesh(1), mesh(1) )
               j2 = mod( i20-i0(2)+10*mesh(2), mesh(2) )
               j3 = mod( i30-i0(3)+10*mesh(3), mesh(3) )
               v(imesh) = v(imesh) + v0 + f(1)*j1 + f(2)*j2 + f(3)*j3
            enddo
         enddo
      enddo

      end subroutine add_potential_from_field
!---------------------------------------------------------------

      end module m_efield
