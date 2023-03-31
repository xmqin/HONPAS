! Auxilary functions for transiesta

module m_ts_aux
      
  use precision, only : dp

  implicit none

  private

  interface nf
     module procedure nf_z
     module procedure nf_d
     module procedure nf2
     module procedure nf22
  end interface nf
  public :: nf

contains
  
  ! Fermi Function
  elemental function nf_d(d) result(nf)
    real(dp), intent(in) :: d
    real(dp) :: nf
    nf = 1._dp/(1._dp + exp(d))
  end function nf_d
  elemental function nf_z(z) result(nf)
    complex(dp), intent(in) :: z
    complex(dp) :: nf
    nf = dcmplx(1._dp,0._dp)/(dcmplx(1._dp,0._dp) + exp(z))
  end function nf_z

  ! Double fermi function
  ! calculates
  !   nF(E-E1) - nF(E-E2) at the same temperature [kT]
  elemental function nf2(E,E1,E2,kT) result(nf)
    real(dp), intent(in) :: E,E1,E2,kT
    real(dp) :: nf
    nf = nf22(E,E1,kT,E2,kT)
  end function nf2

  ! Double fermi function
  ! calculates
  !   nF(E-E1,kT1) - nF(E-E2,kT2) at separate temperatures
  elemental function nf22(E,E1,kT1,E2,kT2) result(nf)
    real(dp), intent(in) :: E,E1,kT1,E2,kT2
    real(dp) :: nf
    nf = nf_d((E-E1)/kT1) - nf_d((E-E2)/kT2)
  end function nf22

end module m_ts_aux
