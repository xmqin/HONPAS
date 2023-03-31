! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine write_orb_indx( na_u, na_s, no_u, no_s, isa, xa, &
                           iaorb, iphorb, indxuo, nsc, ucell )

! Writes a list of basis orbital indexes and labels in file 
! (SistemLabel).ORB_INDX J.M.Soler. Dec.2009

  ! Used module routines and objects
  use atmfuncs,  only: cnfigfio ! Returns principal quantum number
  use precision, only: dp       ! Double precision real kind
  use atmfuncs,  only: izofis   ! Returns atomic number
  use atmfuncs,  only: labelfis ! Returns atom label
  use atmfuncs,  only: lofio    ! Returns angular mumentum number
  use atmfuncs,  only: mofio    ! Returns magnetic quantum number
  use atmfuncs,  only: pol      ! Returns whether an orbital is polarized
  use atmfuncs,  only: rcut     ! Returns orbital cutoff radius
  use cellsubs,  only: reclat   ! Finds reciprocal unit cell vectors
  use files,     only: slabel   ! Label of simulated system
  use atmfuncs,  only: symfio   ! Returns angular symmetry of orbital
  use atmfuncs,  only: zetafio  ! Returns zeta number of orbital

  ! Input routine arguments
  implicit none
  integer, intent(in):: na_u         ! Total number of atoms in unit cell
  integer, intent(in):: na_s         ! Total number of atoms in super cell
  integer, intent(in):: no_u         ! Total number of orbitals in unit cell
  integer, intent(in):: no_s         ! Total number of orbitals in super cell
  integer, intent(in):: isa(na_s)    ! Species index of each atom
  real(dp),intent(in):: xa(3,na_s)   ! Atomic positions (in a.u.)
  integer, intent(in):: iaorb(no_s)  ! Atom to which each orbital belongs
  integer, intent(in):: iphorb(no_s) ! Orbital index within atom
  integer, intent(in):: indxuo(no_s) ! Index of equiv. atom in unit cell
  integer, intent(in):: nsc(3)       ! Supercell size along each cell vector
  real(dp),intent(in):: ucell(3,3)   ! Unit cell vectors: ucell(ixyz,ivec)

  ! Used external routines
  external:: io_assign   ! Finds and reserves a I/O file unit
  external:: io_close    ! Closes a reserved I/O file unit

  ! Internal variables and arrays
  character(len=32):: atom_label, orb_sym
  logical :: polarized
  integer :: atomic_number, i, ia, iao, io, is, isc(3), isum, iu, iua, iuo, &
             l, m, n, z
  real(dp):: dxa(3), rc, rcell(3,3)

  ! Find reciprocal unit cell vectors (without 2*pi factor)
  call reclat( ucell, rcell, 0 )

  ! Find an available I/O unit and open file for output
  call io_assign(iu)
  open( iu, file=trim(slabel)//'.ORB_INDX' )

  ! Loop on atoms in supercell
  do io = 1,no_s

    ! Find orbital indexes
    ia = iaorb(io)                ! Atom to which orbital belongs
    is = isa(ia)                  ! Atomic species index
    iuo = indxuo(io)              ! Equivalent orbital in first unit cell
    iua = iaorb(iuo)              ! Equivalent atom in first unit cell
    atomic_number = izofis( is )  ! Atomic number
    atom_label = labelfis( is )   ! Atomic label
    iao = iphorb( io )            ! Orbital index within atom
    n = cnfigfio( is, iao)        ! Orbital's principal quantum number
    l = lofio( is, iao )          ! Orbital's angular mumentum number
    m = mofio( is, iao )          ! (Real) orbital's magnetic quantum number
    z = zetafio( is, iao )        ! 'Zeta' index of orbital
    polarized = pol( is, iao )    ! Is this a polarization orbital?
    orb_sym = symfio( is, iao )   ! Name of orbital's angular symmetry
    rc = rcut( is, iao )          ! Orbital's cutoff radius
    dxa(:) = xa(:,ia) - xa(:,iua) ! Cell vector of atom ia
    isc(:) = nint( matmul(dxa,rcell) )  ! Cell index of atom ia
    do i = 1,3
      if (isc(i)>nsc(i)/2) isc(i) = isc(i) - nsc(i) ! Same centered in isc=0
    end do

    ! Write orbital indexes
    if (io==1) write(iu,'(i7,i8,a,/)') no_u, no_s, &
      ' = orbitals in unit cell and supercell. See end of file.'
    if (io==1) & ! Write header
    write(iu,'(a6,a6,a3,1x,a6,a4,4a3,a3,a13,a8, a9,a6)') &
      'io', 'ia', 'is', 'spec', 'iao', &
      'n', 'l', 'm', 'z', 'p', 'sym', 'rc', 'isc  ', 'iuo'
    write(iu,'(i6,i6,i3,1x,a6,i4,4i3,l3,a13,f8.3,3i3,i6)') &
      io, ia, is, trim(atom_label), iao, &
      n, l, m, z, polarized, trim(orb_sym), rc, isc(:), iuo

  end do ! io

  ! Write explanations
  write(iu,'(/,(a))') &
    'Column codes:', &
    '  io = Orbital index in supercell', &
    '  ia = Atom to which orbital belongs', &
    '  is = Atomic species index', &
    'spec = Atomic species label', &
    ' iao = Orbital index within atom', &
    '   n = Principal quantum number', &
    '   l = Angular mumentum quantum number', &
    '   m = Magnetic quantum number of (real) orbital:', &
    '       m<0 => sin(m*phi), m>=0 => cos(m*phi)', &
    '   z = Zeta index of orbital', &
    '   p = Is this a polarization orbital? (False|True)', &
    ' sym = Symmetry name of real orbital', &
    '  rc = Cutoff radius of orbital (Bohr)', &
    ' isc = Unit cell indexes to which orbital belongs:', &
    '       center(io) = center(iuo) + sum_(i=1:3) cell_vec(i) * isc(i)', &
    ' iuo = Equivalent orbital in first unit cell', ' '

  ! Close output file
  call io_close(iu)

end subroutine write_orb_indx

