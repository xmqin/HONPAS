! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module basis_enthalpy
  implicit none
  private
  public :: write_basis_enthalpy

CONTAINS

  subroutine write_basis_enthalpy(FE,FEharris)
    use precision, only: dp
    use fdf,       only: fdf_physical
    use units,     only: GPa, eV

    real(dp), intent(in) :: FE  ! Electronic (Free)Energy
    real(dp), intent(in) :: FEharris  ! Electronic (Free)HarrisEnergy

    real(dp)  orb_vol, basis_enthalpy, basis_harris_enthalpy, basis_pressure
    integer   iu

    call orb_volumes(orb_vol)
    basis_pressure = GPa * fdf_physical("BasisPressure",0.2_dp,"GPa")
    basis_enthalpy = FE + basis_pressure * orb_vol
    basis_harris_enthalpy = FEHarris + basis_pressure * orb_vol

    write(6,"(a,f18.6)") "(Free)E+ p_basis*V_orbitals  = ", basis_enthalpy/eV
    write(6,"(a,f18.6)") "(Free)Eharris+ p_basis*V_orbitals  = ", basis_harris_enthalpy/eV

    call io_assign(iu)
    open(iu,file="BASIS_ENTHALPY",form="formatted",status="replace")
    write(iu,*) basis_enthalpy/eV
    write(iu,*) "The above number is the electronic (free)energy:", FE/eV
    write(iu,*) "Plus the pressure : ", basis_pressure,  &
                                    " ( ", basis_pressure/GPa, " GPa)"
    write(iu,*) "     times the orbital volume (in Bohr**3): ", orb_vol
    call io_close(iu)

    call io_assign(iu)
    open(iu,file="BASIS_HARRIS_ENTHALPY",form="formatted",status="replace")
    write(iu,*) basis_harris_enthalpy/eV
    write(iu,*) "The above number is the electronic (free) harris energy:", FEHarris/eV
    write(iu,*) "Plus the pressure : ", basis_pressure,  &
                                    " ( ", basis_pressure/GPa, " GPa)"
    write(iu,*) "     times the orbital volume (in Bohr**3): ", orb_vol
    call io_close(iu)

  end subroutine write_basis_enthalpy
!--------------------------------------------------------------------
  subroutine orb_volumes(orb_vol)
    !
    ! Computes the total volume of the basis orbitals
    ! in the unit cell, for the purposes of calculating
    ! the "basis enthalpy" used in the optimization
    ! procedure. See:
    ! E. Anglada et al, PRB ...
    !
    ! Notes: 
    ! 1. Only a single representative of each "nlz" shell is
    !    used (e.g., just one and not 5 for a 'd' shell)
    ! 2. All the species are treated equally. Presumably
    !    one might want to give different weights to different
    !    species (for example, to avoid the case in which
    !    an impurity atom's basis details are drowned by the
    !    much larger influence of the host orbitals).
    !    This will be implemented via an fdf block:
    !
    !    %block PAO.BasisEnthalpyWeights
    !      Species_number_i  Weigth_i
    !      ...
    !    %endblock PAO.BasisEnthalpyWeights.

    use precision, only: dp
    use atmfuncs,  only: lofio, rcut
    use atomlist,  only: lasto, iphorb
    use siesta_geom, only: isa, na_u
    use units,     only: pi

    real(dp), intent(out) :: orb_vol

    integer ::  ia, ioa, is, l, io
    real(dp) :: r

    orb_vol = 0.0_dp

    do ia = 1,na_u
       is = isa(ia)
       do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          l = lofio(is,ioa)
          r = rcut(is,ioa)
          orb_vol = orb_vol + ((4.0_dp/3.0_dp)*pi*r**3)/(2*l+1) 
          !Just one per nl shell
       enddo
    enddo
  end subroutine orb_volumes

end module basis_enthalpy
