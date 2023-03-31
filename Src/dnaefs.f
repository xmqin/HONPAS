! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_dnaefs
      public :: dnaefs
      CONTAINS
      subroutine dnaefs( nua, na, scell, xa, indxua, rmaxv,
     .                   isa, DEna, fa, stress, forces_and_stress)
C *********************************************************************
C Correction of Neutral Atom energies, forces and stress due to the
C overlap between ionic (bare pseudopotential) charges.
C Energies in Ry. Lengths in Bohr.
C Written by J.Soler and P.Ordejon, August-October'96, June'98 
C Modified by DSP Aug., 1998
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxv            : Maximum cutoff for NA potential
C integer isa(na)          : Species index of each atom
C logical forces_and_stress   Determines whether fa and stress are touched
C **************************** OUTPUT *********************************
C real*8 DEna              : NA energy correction
C ********************** INPUT and OUTPUT *****************************
C real*8 fa(3,na)          : Atomic forces (NA part added to input)
C real*8 stress(3,3)       : Stress tensor (NA part added to input)
C *********************************************************************
C The following function must exits:
C
C    INTEGER  FUNCTION  IZOFIS(IS) : Returns the atomic number
C Input:
C      INTEGER IS : Specie index
C
C *********************************************************************

      use precision
      use atmfuncs,  only: izofis, psover
      use neighbour, only: mneighb, jna=>jan, xij, r2ij,
     &                     reset_neighbour_arrays

      implicit none

      integer  na, nua

      integer  indxua(na), isa(na)

      real(dp)
     .  scell(3,3), DEna, fa(3,nua), rmaxv, stress(3,3), xa(3,na)

      logical, intent(in)  :: forces_and_stress

C Internal variables ......................................................
      integer
     .  ia, is, ix, ja, jn, js, jx, jua, nnia

      real(dp)
     .  dvdr, fij(3), rij, r2min, vij, volcel, volume
       
      parameter ( r2min = 1.d-15 )
C ......................

      call mneighb( scell, 2.0d0*rmaxv, na, xa, 0, 0, nnia )

      volume = nua * volcel(scell) / na
      DEna = 0.0d0
 
      do ia = 1,nua

C Find neighbour atoms
        call mneighb( scell, 2.0d0*rmaxv, na, xa, ia, 0, nnia )
        do jn = 1,nnia
          ja = jna(jn)
          jua = indxua(ja)
          is = isa(ia)
          js = isa(ja) 

          if (r2ij(jn).gt.r2min .and.
     .        izofis(is).gt.0 .and. izofis(js).gt.0) then
            rij = sqrt( r2ij(jn) )
            call psover( is, js, rij, vij, dvdr )
            DEna = DEna + vij / 2.0d0
            if (forces_and_stress) then
               do ix = 1,3
                  fij(ix) = dvdr * xij(ix,jn) / rij / 2.0d0
                  fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                  fa(ix,jua) = fa(ix,jua) - fij(ix)
                  do jx = 1,3
                     stress(jx,ix) = stress(jx,ix) +
     .                    xij(jx,jn) * fij(ix) / volume
                  enddo
               enddo
            endif
          endif
        enddo
      enddo
      call reset_neighbour_arrays( )
      end subroutine dnaefs
      end module m_dnaefs
