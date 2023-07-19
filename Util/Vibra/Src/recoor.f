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
      subroutine recoor(overflow, cell, alat, xa, isa, xmass, na) 

c *******************************************************************
c Reads atomic coordinates and format in which they are given, to be
c transformed into Bohr cartesian for internal handling. 
c It also shifts them all according to AtomicCoordinatesOrigin.
c
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996.
c ********* INPUT ***************************************************
c logical overflow          : true if some dimension is too small
c integer na                : number of atoms
c double precision cell(3,3): Lattice (supercell) vectors
c double precision alat     : Lattice constant (in Bohr) 
c ********* OUTPUT **************************************************
c double precision xa(3,na) : atomic coordinates in Bohr cartesian
c integer isa(na)           : atomic species of different atoms
c integer xmass(na)         : atomic masses of different atoms
c *******************************************************************

      use fdf

      implicit          none
      logical           overflow
      integer           na
      integer           isa(na)
      double precision  xa(3,na), cell(3,3), alat, xmass(na)

c Internal variables and arrays

      character         acf*22, acf_defect*22
      logical           leqi
      integer           iscale, ia, i, ix, iunit
      double precision  origin(3), xac(3)



c enable FDF input/output

      data origin /3*0.d0/

C format of atomic coordinates

      acf_defect = 'NotScaledCartesianBohr'
      acf = fdf_string('AtomicCoordinatesFormat',acf_defect)
      if (leqi(acf,'NotScaledCartesianBohr')) then
        iscale = 0
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   'Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '(in Bohr units)'
      else if (leqi(acf,'NotScaledCartesianAng')) then
        iscale = 1
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   'Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '(in Angstroms)'
      else if (leqi(acf,'ScaledCartesian')) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'recoor: ERROR: Explicit lattice ',
     .       'constant is needed for ScaledCartesian format'
           stop 'recoor: ERROR: Explicit lattice constant needed'
        endif
        iscale = 2
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   'Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '(in units of alat)'
      else if (leqi(acf,'ScaledByLatticeVectors') .or. 
     .         leqi(acf,'Fractional') ) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'recoor: ERROR: Explicit lattice ',
     .       'constant is needed for Fractional format'
           stop 'recoor: ERROR: Explicit lattice constant needed'
        endif
        iscale = 3
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   'Ref. to lattice vectors'
      else
        write(6,"(/,'recoor: ',72(1h*))")
        write(6,"('recoor:                  INPUT ERROR')")
        write(6,'(a)') 'recoor: '
        write(6,'(2a)') 'recoor: You must use one of the following',
     .                            ' coordinate scaling options:'
        write(6,'(a)') 'recoor:     - NotScaledCartesianBohr         '
        write(6,'(a)') 'recoor:     - NotScaledCartesianAng          '
        write(6,'(a)') 'recoor:     - ScaledCartesian                '
        write(6,'(2a)') 'recoor:     - ScaledByLatticeVectors ',
     .                                               '(or Fractional)'
        write(6,"('recoor: ',72(1h*))")
        stop 'recoor: ERROR: Wrong atomic-coordinate input format'
      endif


c read atomic coordinates and species

      if (.not. overflow) then

        if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .    then
          do ia = 1,na
            read(iunit,*) (xa(i,ia), i=1,3), isa(ia), xmass(ia)
          enddo
        else
          write(6,"(/,'recoor: ',72(1h*))")
          write(6,"('recoor:                  INPUT ERROR')")
          write(6,'(a)')
     .    'recoor:   You must specify the atomic coordinates'
          write(6,"('recoor: ',72(1h*))")
          stop 'recoor: ERROR: Atomic coordinates missing'
        endif

        if ( fdf_block('AtomicCoordinatesOrigin',iunit) ) then
          read(iunit,*) (origin(i),i=1,3)
          do ia = 1,na
            do i = 1,3
              xa(i,ia) = xa(i,ia) + origin(i)
            enddo
          enddo
        endif


c Scale atomic coordinates
c   Coord. option = 0 => Do nothing
c   Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)
c   Coord. option = 2 => Multiply by lattice constant
c   Coord. option = 3 => Multiply by lattice vectors

        if (iscale .eq. 1) then
          do ia = 1,na
            do ix = 1,3
              xa(ix,ia) = 1.d0 / 0.529177d0 * xa(ix,ia)
            enddo
          enddo
        elseif (iscale .eq. 2) then
          do ia = 1,na
            do ix = 1,3
              xa(ix,ia) = alat * xa(ix,ia)
            enddo
          enddo
        elseif (iscale .eq. 3) then
          do ia = 1,na
            do ix = 1,3
              xac(ix) = xa(ix,ia)
            enddo
            do ix = 1,3
              xa(ix,ia) = cell(ix,1) * xac(1) +
     .                    cell(ix,2) * xac(2) +
     .                    cell(ix,3) * xac(3)
            enddo
          enddo
        endif

        write(6,'(a)') 'recoor: Atomic coordinates (Bohr) and species'
        do ia = 1,na
          write(6,"('recoor: ',i4,2x,3f10.5,i3)")
     .                    ia,(xa(ix,ia),ix=1,3),isa(ia)
        enddo

      endif

      return
      end
