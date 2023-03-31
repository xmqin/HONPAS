! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_naefs
      public :: naefs
      CONTAINS
      subroutine naefs(nua, na, scell, xa, indxua, rmaxv,
     .                 isa, Ena, fa, stress, forces_and_stress)
C *********************************************************************
C Neutral Atom (NA) energy, forces and stress.
C This is the self energy of rho_na=-Laplacian(v_na(Ry))/(8*pi)
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, August-October'96, June'98 
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
C real*8 Ena               : NA energy
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
      use atmfuncs,  only: izofis, vna_gindex
      use neighbour, only: jna=>jan, xij, mneighb,
     &                     reset_neighbour_arrays
      use m_new_matel,   only : new_matel

      implicit none

      integer na, nua

      integer
     .  indxua(na), isa(na)

      real(dp)
     .  scell(3,3), Ena, fa(3,nua), rmaxv, stress(3,3), xa(3,na)

      logical, intent(in)  :: forces_and_stress

C Internal variables ......................................................
      integer  ia, is, ix, ja, jn, js, jx, jua, nnia, ig, jg

      real(dp)  fij(3), pi, vij, volcel, volume 
      
C ......................
      call timer( 'naefs', 1 )

C Initialize neighb subroutine
      call mneighb( scell, 2.d0*rmaxv, na, xa, 0, 0, nnia )

      pi = 4.0d0 * atan(1.0d0)
      volume = nua * volcel(scell) / na
      Ena = 0.0d0
 
      do ia = 1,nua

C Find neighbour atoms
        call mneighb( scell, 2.0d0*rmaxv, na, xa, ia, 0, nnia )
        do jn = 1,nnia
          ja = jna(jn)
          jua = indxua(ja)
          is = isa(ia)
          js = isa(ja)
          if (izofis(is).gt.0 .and. izofis(js).gt.0) then
            ig = vna_gindex(is)
            jg = vna_gindex(js)
            call new_MATEL( 'T', ig, jg, xij(1:3,jn), vij, fij )
            Ena = Ena + vij / (16.0d0*pi)
            if (forces_and_stress) then
               do ix = 1,3
                  fij(ix) = fij(ix) / (16.0d0*pi)
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

C     Free local memory
!      call new_MATEL( 'T', 0, 0, 0, 0, xij, vij, fij )
      call reset_neighbour_arrays( )
      call timer( 'naefs', 2 )
      end subroutine naefs
      end module m_naefs


