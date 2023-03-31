! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_intramol_pressure
CONTAINS
subroutine remove_intramol_pressure(ucell,stress,na_u,xa,fa,mstress)
 use precision, only: dp
 use zmatrix, only: nZmol, nZmolStartAtom
 use zmatrix, only: lUseZmatrix, nZmolAtoms

 implicit none

 real(dp), external :: volcel

 integer, intent(in)   :: na_u
 real(dp), intent(in)  :: ucell(3,3), xa(3,na_u), fa(3,na_u)
 real(dp), intent(in)  :: stress(3,3)
 real(dp), intent(out) :: mstress(3,3)

 real(dp) :: virial, fmean, volume
 integer  :: ifirst, ilast, natoms_mol, imol, ix, ia

! Find intramolecular contributions to the pressure 
        volume = volcel(ucell)
        virial = 0.0_dp
        if (lUseZmatrix) then
           ! Molecule by molecule
           do imol = 1, nZmol
              ifirst = nZmolStartAtom(imol)
              natoms_mol = nZmolAtoms(imol)
              ilast  = ifirst + natoms_mol - 1
              do ix = 1,3
               fmean = sum(fa(ix,ifirst:ilast)) / natoms_mol
               do ia = ifirst, ilast
                 virial = virial + xa(ix,ia) * (fa(ix,ia) - fmean)
               enddo
              enddo
           enddo
        else  
           ! No Zmatrix
           ! Consider the whole system as a molecule, even if it is meaningless
           do ix = 1,3
              fmean = sum(fa(ix,1:na_u)) / na_u
              do ia = 1,na_u
                 virial = virial + xa(ix,ia) * (fa(ix,ia) - fmean)
              enddo
           enddo
        endif
!
        mstress(1:3,1:3) = stress(1:3,1:3)
        do ix = 1, 3
           mstress(ix,ix) = mstress(ix,ix) + virial / volume / 3.0_dp
        enddo

end subroutine remove_intramol_pressure
end module m_intramol_pressure
