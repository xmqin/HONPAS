! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine coor(na,cell)

C Reads atomic coordinates and format in which they are given, to be
C transformed into Bohr cartesian for internal handling. 
C It also shifts them all according to AtomicCoordinatesOrigin.
C Written by E. Artacho, December 1997, on the original piece of the
C redata subroutine written by P. Ordejon in December 1996.
C Modified by J.M.Soler. August 1998.

!! Modified by Alberto Garcia, May 16, 2000

C integer na  : Number of atoms in unit cell
!               Note: the user might specify a supercell in the
!                     input file, and the coordinates are generated
!                     accordingly in this routine, but this does
!                     not have anything to do with the 'virtual'
!                     supercell used to deal with k-point sampling.

      use precision, only : dp
      use sys, only: die
      use fdf
      use alloc
      use parallel,  only : Node
      use units, only: Ang
      use zmatrix
      use siesta_geom, only : xa, isa, cisa
      implicit none

      integer,  intent(out) :: na
      real(dp), intent(out) :: cell(3,3)

C Internal variables and arrays

      character             :: acf*22, acf_default*22
      logical               :: isdiag
      integer               :: ia, ic, iscale, iua, ix, i, j,
     .                         coor_mscell(3,3), ncells, 
     &                         coor_nsc(3), nua
      real(dp)              :: alat, dcell(3,3), dscell(3,3), 
     .                         origin(3), unit_cell(3,3),
     .                         volcel, xac(3), super_cell(3,3)
      external          digcel, redcel, volcel

      logical               :: lOrigin

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline


C Read unit cell and supercell
      call redcel( alat, unit_cell, super_cell, coor_mscell )

C Format of atomic coordinates
      acf_default = 'NotScaledCartesianBohr'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)
      if (leqi(acf,'NotScaledCartesianBohr') .or.
     .    leqi(acf,'Bohr') ) then
        iscale = 0
        if (Node.eq.0) then
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                      ',
     .     '      (in Bohr units)'
        endif
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .         leqi(acf,'Ang') ) then
        iscale = 1
        if (Node.eq.0) then
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                    ',
     .     '      (in Angstroms)'
        endif
      else if (leqi(acf,'ScaledCartesian')) then
        if (alat.eq.0.0_dp) then
          if (Node.eq.0) then
            write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .        'constant is needed for ScaledCartesian format'
          endif
          call die()
        endif
        iscale = 2
        if (Node.eq.0) then
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                    ',
     .     '      (in units of alat)'
        endif
      else if (leqi(acf,'ScaledByLatticeVectors') .or. 
     .         leqi(acf,'Fractional') ) then
        if (alat.eq.0.0_dp) then
          if (Node.eq.0) then
            write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .        'constant is needed for Fractional format'
          endif
          call die()
        endif
        iscale = 3
        if (Node.eq.0) then
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Fractional'
        endif
      else
        if (Node.eq.0) then
          write(6,"(/,'coor: ',72('*'))")
          write(6,"('coor:                  INPUT ERROR')")
          write(6,'(a)') 'coor: '
          write(6,'(2a)') 'coor: You must use one of the following',
     .                              ' coordinate scaling options:'
          write(6,'(a)') 'coor:     - NotScaledCartesianBohr (or Bohr)'
          write(6,'(a)') 'coor:     - NotScaledCartesianAng (or Ang) '
          write(6,'(a)') 'coor:     - ScaledCartesian                '
          write(6,'(2a)') 'coor:     - ScaledByLatticeVectors ',
     .                                                 '(or Fractional)'
          write(6,"('coor: ',72('*'))")
        endif

        call die()

      endif

C Read atomic coordinates and species
      na = fdf_integer('NumberOfAtoms',0)
      if ( na == 0 ) then
       ! estimate the number of atoms by reading the block
       ! Currently we do not try to estimate from a Z-matrix
        na = fdf_block_linecount('AtomicCoordinatesAndAtomicSpecies',
     &      'vvvi')
      end if
      if ( na == 0 ) then
         call die("Must specify number of atoms AND coordinates!")
      end if
!
!     Check if we need more space to accommodate supercell
!     (still a "real" supercell, specified for convenience!)
!
      if (volcel(unit_cell) .lt. 1.d-8) then
        ncells = 1
      else
          ncells = nint( volcel(super_cell) / volcel(unit_cell) )
      endif

      nua = na
      na  = na * ncells

C Find origin with which to translate all coordinates
      lOrigin = fdf_block('AtomicCoordinatesOrigin',bfdf)
      if (lOrigin) then 
         if (.not. fdf_bline(bfdf,pline))
     .        call die('coor: ERROR in AtomicCoordinatesOrigin block')
         origin(1) = fdf_bvalues(pline,1)
         origin(2) = fdf_bvalues(pline,2)
         origin(3) = fdf_bvalues(pline,3)
         call fdf_bclose(bfdf)
      else
         origin = 0._dp
      endif

C Set array dimensions
      nullify( isa, cisa, xa )
      call re_alloc( isa, 1, na, 'isa', 'coor' )
      call re_alloc( xa, 1, 3, 1, na, 'xa', 'coor' )
      allocate(cisa(na))
!      call re_alloc( cisa, 1, na, 'cisa', 'coor' )

C Attempt to read a Z matrix
      call read_Zmatrix( nua, isa, alat, unit_cell, lOrigin, origin )

      if (lUseZmatrix) then
        
        cell = unit_cell
        if (ncells > 1 ) call die("Cannot use SuperCell with Zmatrix")

C       Generate Cartesian coordinates from Z-matrix
        call Zmat_to_Cartesian(xa) 

      else

C If Z matrix hasn't been found, read regular atomic coordinates
        if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',bfdf) ) then
          do ia= 1, nua
            if (.not. fdf_bline(bfdf,pline))
     .        call die('coor: ERROR in ' //
     .                 'AtomicCoordinatesAndAtomicSpecies block')
            xa(1,ia) = fdf_bvalues(pline,1)
            xa(2,ia) = fdf_bvalues(pline,2)
            xa(3,ia) = fdf_bvalues(pline,3)
            isa(ia)  = fdf_bintegers(pline,1)
            if (lOrigin) then
              xa(1:3,ia) = xa(1:3,ia) + origin(1:3)
            endif

          enddo
          call fdf_bclose(bfdf)
        else
          call die("coor: You must specify the atomic coordinates")
        endif

C Scale atomic coordinates
C   Coord. option = 0 => Do nothing
C   Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)
C   Coord. option = 2 => Multiply by lattice constant
C   Coord. option = 3 => Multiply by lattice vectors

        if (iscale .eq. 1) then
          xa = xa * Ang
        elseif (iscale .eq. 2) then
          xa = xa * alat
        elseif (iscale .eq. 3) then
          do ia = 1,nua
            do ix = 1,3
              xac(ix) = xa(ix,ia)
            enddo
            do ix = 1,3
                xa(ix,ia) = unit_cell(ix,1) * xac(1) +
     .                      unit_cell(ix,2) * xac(2) +
     .                      unit_cell(ix,3) * xac(3)
            enddo
          enddo
        endif
   

C Expand the coordinates to whole supercell

        if (ncells.gt.1) then

C Find equivalent diagonal combination of cell/supercell
          call digcel( unit_cell, coor_mscell, dcell, dscell,  
     &                 coor_nsc, isdiag )
          if ( .not. isdiag) then
             !
             ! Inform the user but terminate the program
             !
             if (Node == 0) then
                write(6,*) "Supercell is not diagonal."
                write(6,*) "Equiv. diagonal unit cell:"
                do i=1,3
                   write(6,"(3f12.6)") (dcell(j,i),j=1,3)
                enddo
                write(6,*) "Original unit cell:"
                do i=1,3
                   write(6,"(3f12.6)") (unit_cell(j,i),j=1,3)
                enddo
                write(6,*) "Supercell multipliers:", coor_nsc(:)
!!                coor_nsc(:) = abs(coor_nsc(:))
             endif
             call die("Not safe to use non-diagonal supercells")
          endif

C Expand coordinates
          call superx_coor( dcell, coor_nsc, nua, na, xa, dscell )
          cell = dscell

C Expand index isa
          ia = 0
          do ic = 1,ncells
            do iua = 1,nua
              ia = ia + 1
              isa(ia) = isa(iua)
            enddo
          enddo

        else

           cell = unit_cell

        endif  ! (ncells  > 1 ?)

      endif    ! (use Zmatrix?)

      ! Construct references: cannot use isa alone, since 
      ! refs may not start with a number.

      do ia=1, na
        write(cisa(ia), '("siesta:e",i3.3)') isa(ia)
      enddo

      end subroutine coor

      subroutine superx_coor( unit_cell,NSC,nua,maxa,XA,super_cell)
C **********************************************************************
C Generates supercell vectors and atomic positions.
C Written by J.M.Soler, August 1998
C *************** Input ************************************************
C Real*8  unit_cell(3,3)  : Unit cell vectors unit_cell(Ixyz,Ivector)
C Integer NSC(3)      : Number of cells in each supercell direction:
C                         super_cell(ix,i) = unit_cell(ix,i) * NSC(i)
C Integer NUA          : Number of atoms in unit cell
C Integer MAXA        : Second dimension of XA
C *************** Input and output *************************************
C Real*8  XA(3,MAXA)  : Atomic positions in unit cell (input) and
C                       in supercell (output), in cartesian coord.
C Real*8  super_cell(3,3)  : Supercell vectors
C *********** Units ****************************************************
C Units of CELL and XA are arbitrary but must be the same
C *********** Behavior *************************************************
C - If NUA*NCELLS > MAXA (where NCELLS is the total number of cells),
C   the supercell atomic coordinates are not generated.
C - The first supercell atoms are those of the initial unit cell, i.e.
C   the positions XA(i,ia) for (ia.le.NUA) are not modified.
C - The remaining atoms are ordered by unit cells, i.e. the atom ia
C   is equivalent to the unit-cell atom ja=MOD(ia-1,NUA)+1
C **********************************************************************

      use precision, only : dp
      use sys, only: die

      implicit          none

      integer           maxa, nua, NSC(3)
      real(dp)          super_cell(3,3), unit_cell(3,3), XA(3,maxa)

C Internal variables
      integer           I, I1, I2, I3, IA, IX, JA, ncells
      real(dp)          XC(3)

C Find supercell vectors
      do I = 1,3
        do IX = 1,3
          super_cell(IX,I) = unit_cell(IX,I) * NSC(I)
        enddo
      enddo

C Expand atomic positions to supercell
      ncells = NSC(1) * NSC(2) * NSC(3)
      if (nua*ncells .le. maxa) then
        IA = 0
        do I3 = 0,NSC(3)-1
        do I2 = 0,NSC(2)-1
        do I1 = 0,NSC(1)-1
          do IX = 1,3
            XC(IX) = unit_cell(IX,1)*I1 + unit_cell(IX,2)*I2 +
     $                                    unit_cell(IX,3)*I3
          enddo
          do JA = 1,nua
            IA = IA + 1
            do IX = 1,3
              XA(IX,IA) = XA(IX,JA) + XC(IX)
            enddo
          enddo
        enddo
        enddo
        enddo
      else
         call die("Not enough space for supercell extension")
      endif

      end subroutine superx_coor
