! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine outcoor(cell, xa, na, cohead, writec) 
c *******************************************************************
c Writes atomic coordinates in format given by fdf:AtomCoorFormatOut
c Default: same format as was used for the reading.
c
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996. 
C Modified by DSP, Aug. 1998.
c ********* INPUT ***************************************************
c integer na                : number of atoms
c double precision cell(3,3): Lattice (supercell) vectors
c double precision alat     : Lattice constant (in Bohr) 
c double precision xa(3,na) : atomic coordinates in Bohr cartesian
c integer isa(na)           : atomic species of different atoms
c character*(*) cohead      : special phrase to include in heading
c logical writec            : writing coor only if true
c *******************************************************************
c inversion of lattice vectors for fractional coordinates is done
c through subroutine reclat.
c *******************************************************************

      use atmfuncs,  only : labelfis
      use siesta_geom,  only : isa
      use fdf, only : fdf_physical, fdf_string, leqi
      use precision, only : dp
      use units, only : Ang
      use sys,   only: die
      use alloc, only: re_alloc, de_alloc

      implicit          none

      integer,          intent(in) :: na
      real(dp),         intent(in) :: cell(3,3)
      real(dp),         intent(in) :: xa(3,na)
      character(len=*), intent(in) :: cohead
      logical,          intent(in) :: writec

c Internal variables and arrays

      character         acf*22, acf_default*22, acfout*22, 
     .                  pieceh*20, titl*60
      logical           frstme
      integer           ia, ix
      real(dp)          recell(3,3), alat

      real(dp), dimension(:,:), pointer :: xap

      save              frstme, acfout, alat

      data              frstme /.true./
c--------------------------------------------------------------------

      if (frstme) then

C read format for output of atomic coordinates

        acf_default = 'NotScaledCartesianBohr'
        acf = fdf_string('AtomicCoordinatesFormat',acf_default)
        acfout = fdf_string('AtomCoorFormatOut',acf)

c Scale atomic coordinates
c   Coord. option Bohr       => Do nothing
c   Coord. option Ang        => Multiply by 0.529177 (Bohr --> Ang)
c   Coord. option Scaled     => Divide by lattice constant
c   Coord. option Fractional => Multiply by inverse of lattice vectors

        alat = fdf_physical('LatticeConstant',0._dp,'Bohr')
        if (alat.eq.0._dp .and. leqi(acfout,'ScaledCartesian')) then
           write(6,"(/,2a)") 'outcoor: WARNING: Explicit lattice ',
     .         'constant is needed for ScaledCartesian output.'
           write(6,"(2a)")   'outcoor:          NotScaledCartesianAng',
     .         'format being used instead.'
           acfout = 'NotScaledCartesianAng'
        endif

        frstme = .false.
      endif

c Write coordinates at every time or relaxation step?

      if ( (cohead .eq. ' ') .and. ( .not. writec) ) return

C Allocate local memory
      nullify( xap )
      call re_alloc( xap, 1, 3, 1, na, name='xap', routine='outcoor' )

c write coordinates according to format 

      if (leqi(acfout,'NotScaledCartesianBohr') .or. 
     .    leqi(acfout,'Bohr') ) then
        pieceh = '(Bohr):'
        xap = xa
      else if (leqi(acfout,'NotScaledCartesianAng') .or.
     .         leqi(acfout,'Ang') ) then
        pieceh = '(Ang):'
        xap = xa / Ang
      else if (leqi(acfout,'ScaledCartesian')) then
        pieceh = '(scaled):'
        xap = xa / alat
      else if (leqi(acfout,'ScaledByLatticeVectors') .or. 
     .         leqi(acfout,'Fractional') ) then
        pieceh = '(fractional):'
        call reclat(cell, recell, 0)
        xap = matmul(transpose(recell),xa)

      else
        write(6,"(/,'outcoor: ',a)") repeat('*',72)
        write(6,"('outcoor:                  INPUT ERROR')")
        write(6,'(a)') 'outcoor: '
        write(6,'(2a)') 'outcoor: You must use one of the following',
     .                            ' coordinate output options:'
        write(6,'(a)') 'outcoor:     - NotScaledCartesianBohr (or Bohr)'
        write(6,'(a)') 'outcoor:     - NotScaledCartesianAng (or Ang) '
        write(6,'(a)') 'outcoor:     - ScaledCartesian                '
        write(6,'(2a)') 'outcoor:     - ScaledByLatticeVectors ',
     .                                               '(or Fractional)'
        write(6,"('outcoor: ',a)") repeat('*',72)
        call
     $   die('outcoor: ERROR: Wrong atomic-coordinate output format')
      endif

c writing a heading for the coordinates

      if (cohead .eq. ' ') then
         titl = 'outcoor: Atomic coordinates ' // pieceh 
      else
         titl = 'outcoor: ' // trim(cohead) // 
     .          ' atomic coordinates ' // pieceh
      endif

      write(6,'(/,a)') titl

c writing the coordinates

      write(6,'(3f14.8,i4,2x,i6,2x,a)')
     .  ((xap(ix,ia),ix=1,3),isa(ia),ia,trim(labelfis(isa(ia))),ia=1,na)

C Deallocate local memory
      call de_alloc( xap, name='xap', routine='outcoor')

      end subroutine outcoor
