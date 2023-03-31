!!@LICENSE
!
!******************************************************************************
! MODULE xcmod
! Stores the information about the XC functional to use,
! and provides routines to set and to get it.
!
!******************************************************************************
! subroutine setXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Sets the xc functional(s) to be used by atomxc and/or cellxc
! ------------------------- INPUT ---------------------------------------------
!     integer,         :: n       ! Number of functionals
!     character(len=*),:: func(n) ! Functional name labels
!     character(len=*),:: auth(n) ! Functional author labels
!     real(dp),        :: wx(n)   ! Functional weights for exchange
!     real(dp),        :: wc(n)   ! Functional weights for correlation
!                                   It should be sum(wx)=sum(wc)=1
!
! Allowed functional/author values:
! XCfunc: 
!   'LDA' or 'LSD' => Local density approximation
!            'GGA' => Generalized gradients approx.
!            'VDW' => Van der Waals functional
! XCauth:
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
!           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
!                     the local density limit of the next:
!            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
!           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
!         'revPBE' => GGA Zhang & Yang, PRL 80,890(1998)
!            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
!             'WC' => GGA Wu-Cohen (see subroutine wcxc)
!         'PBESOL' => GGA Perdew et al, PRL, 100, 136406 (2008)
!           'AM05' => GGA Mattsson & Armiento, PRB, 79, 155101 (2009)
!      'PBEJsJrLO' => GGA Reparametrizations of the PBE functional by
!     'PBEJsJrHEG' => GGA   L.S.Pedroza et al, PRB 79, 201106 (2009) and
!      'PBEGcGxLO' => GGA   M.M.Odashima et al, JCTC 5, 798 (2009)
!     'PBEGcGxHEG' => GGA using 4 different combinations of criteria
!          'DRSLL' => VDW Dion et al, PRL 92, 246401 (2004)
!          'LMKLL' => VDW K.Lee et al, PRB 82, 081101 (2010)
!            'KBM' => VDW optB88-vdW of J.Klimes et al, JPCM 22, 022201 (2010)
!            'C09' => VDW V.R. Cooper, PRB 81, 161104 (2010)
!             'BH' => VDW K. Berland and Per Hyldgaard, PRB 89, 035412 (2014)
!             'VV' => VDW Vydrov-VanVoorhis, JCP 133, 244103 (2010)
!
! ------------------------ USAGE ----------------------------------------------
!   use siestaXC, only: setXC
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!
! --------------------- BEHAVIOUR ---------------------------------------------
! - Stops with an error message if n is larger than internal parameter maxFunc
! - Prints a warning message if sum(wx)/=1 or sum(wc)/=1
!
!******************************************************************************
! subroutine getXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Returns the xc functional(s) that has been previously set
! --------------------- OPTIONAL OUTPUT ---------------------------------------
!     integer         :: n       ! Number of functionals
!     character(len=*):: func(n) ! Functional name labels
!     character(len=*):: auth(n) ! Functional author labels
!     real(dp)        :: wx(n)   ! Functional weights for exchange
!     real(dp)        :: wc(n)   ! Functional weights for correlation
!
! ------------------------ USAGE ----------------------------------------------
!   use precision, only: dp
!   use siestaXC,  only: getXC
!   integer,parameter:: maxFunc = 10
!   character(len=20):: func(maxFunc), auth(maxFunc)
!   real(dp):: wx(maxFunc), wc(maxFunc)
!   call getXC( n, func, auth, wx, wc )
!
! --------------------- BEHAVIOUR ---------------------------------------------
! - Does not change any output array whose size is smaller than nFunc
!
!******************************************************************************

module xcmod

  use precision, only: dp              ! Double precision real kind
  use sys,       only: die             ! Termination routine
  use m_vdwxc,   only: vdw_set_author  ! Sets vdW functional flavour

  implicit none

public:: &
  setXC, &! Sets XC functional(s) to be used
  getXC   ! Returns the XC functional(s) being used

private ! Nothing is declared public beyond this point

  integer, parameter :: maxFunc = 20
  integer,           save :: nXCfunc=0
  character(len=20), save :: XCauth(MaxFunc), XCfunc(MaxFunc)
  real(dp),          save :: XCweightX(MaxFunc), XCweightC(MaxFunc)

contains

  subroutine setXC( n, func, auth, wx, wc )
    implicit none
    integer,         intent(in):: n       ! Number of functionals
    character(len=*),intent(in):: func(n) ! Functional name labels
    character(len=*),intent(in):: auth(n) ! Functional author labels
    real(dp),        intent(in):: wx(n)   ! Functl weights for exchng
    real(dp),        intent(in):: wc(n)   ! Functl weights for correl
    integer:: i, j
    if (n>maxFunc) call die('setXC: ERROR: parameter maxFunc too small')
    nXCfunc = n
    XCfunc(1:n) = func(1:n)
    XCauth(1:n) = auth(1:n)
    XCweightX(1:n) = wx(1:n)
    XCweightC(1:n) = wc(1:n)
    do i = 1,n
      if (XCfunc(i)=='VDW' .or. XCfunc(i)=='vdw' .or. XCfunc(i)=='vdW') then
        XCfunc(i) = 'VDW'
        do j = 1,i-1
          if (XCfunc(j)=='VDW' .and. XCauth(j)/=XCauth(i)) &
            call die('setXC ERROR: mixing different VDW authors not allowed')
        end do ! j
        call vdw_set_author( XCauth(i) )
      end if ! (XCfunc(i)=='VDW')
    end do ! i
  end subroutine setXC

  subroutine getXC( n, func, auth, wx, wc )
    implicit none
    integer,         optional,intent(out):: n       ! Number of functionals
    character(len=*),optional,intent(out):: func(:) ! Functional name labels
    character(len=*),optional,intent(out):: auth(:) ! Functional author labels
    real(dp),        optional,intent(out):: wx(:)   ! Functl weights for exchng
    real(dp),        optional,intent(out):: wc(:)   ! Functl weights for correl
    integer:: nf
    nf = nXCfunc
    if (present(n)) n = nf
    if (present(func)) then
       if (size(func)>=nf) func(1:nf) = XCfunc(1:nf)
    end if
    if (present(auth)) then
       if (size(auth)>=nf) auth(1:nf) = XCauth(1:nf)
    end if
    if (present(wx)) then
       if (size(wx)  >=nf) wx(1:nf)   = XCweightX(1:nf)
    end if
    if (present(wc)) then
       if (size(wc)  >=nf) wc(1:nf)   = XCweightC(1:nf)
    end if
  end subroutine getXC

end module xcmod
