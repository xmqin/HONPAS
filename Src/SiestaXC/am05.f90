module am05
! $Header:$
!**********************************************************************
! Armiento Mattsson am05 functional for exchange and correlation,
!  
! Spin polarization type: exchange index scales with the separate 
! spin densities; correlation index and prefactor both scales with
! the separate spin densities. (xscss)
! 
! Version: 7 (xscss, web)
!
! The latest version of this subroutine file is available at
!   http://dft.sandia.gov/functionals/AM05.html
!
!
! Usage:
!
! Below follows first a number of usage examples on how to calculate energy
! and potential in codes that use different schemes for calculating
! the potential, and therefore uses different input quantities. After the
! usage examples follows the main am05 routine.
!
! Examples for the following schemes are included:
!
!   am05trads : The traditional scheme used e.g. by PRB 33, 8800 (1986).
!   am05wbs : The White and Bird scheme [PRB 50, 4954 (1994)].
!   am05pgjs : The Pople, Gill, and Johnson scheme [CPL 199, 557 (1992)].
!   am05wbnums : The White and Bird scheme with numerical derivatives.
!
!
! Citation request: 
!
! When using this functional, please cite:
! "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
!  R. Armiento and A. E. Mattsson (unpublished)."
! (The first paper for the AM05 functional, the second for the
! spin-polarized version)
!
!
! License and copyright:
!
! (c) Rickard Armiento 2005-2008
!
! Permission is hereby granted, free of charge, to any person obtaining 
! a copy of this software and associated documentation files (the 
! "Software"), to deal in the Software without restriction, including 
! without limitation the rights to use, copy, modify, merge, publish, 
! distribute, sublicense, and/or sell copies of the Software, and to 
! permit persons to whom the Software is furnished to do so, subject 
! to the following conditions:
!
! The above copyright notice and this permission notice shall be 
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
! OTHER DEALINGS IN THE SOFTWARE.
!
!**********************************************************************
  use precision, only : dp
  use sys,       only : die

  implicit none

  public :: am05wbs

CONTAINS

!**********************************************************************
! saferecp
! Helper function for making divisions with built in cutoff 
! for very small denominators
!
! output
!   set   variable to assign the value of the reciprocal
!
! input
!   nom   nominator
!   denom denominator
!**********************************************************************

  subroutine saferecp(set, nom, denom)

!     ** Input parameters
  real(dp) nom, denom

!     ** Output parameters
  real(dp) set

  if (denom .ge. 1e-30) then
    set = nom/denom
  elseif (nom .le. 1e-30) then
    set = 0.0_dp
  elseif (nom .ge. 1e-30) then
    set = 1d30
  endif

  return
  end subroutine saferecp

!**********************************************************************
! am05wbs
! Usage example for White and Bird scheme [PRB 50, 4954 (1994)].
!
! input
!   nup    electron density [bohr**(-3)]
!   ndn    electron density [bohr**(-3)]
!   gup    abs of gradient of upspin density, |grad(nup)| [bohr**(-4)]
!   gdn    abs of gradient of downspin density, |grad(ndn)| [bohr**(-4)]
!
! output
!   fxc       exchange-correlation energy density [hartree]
!   dfxcup    d(fxc)/dnup [hartree]
!   dfxcdn    d(fxc)/dndn [hartree]
!   dfxcdgup  d(fxc)/d(|grad(nup)|) [hartree * bohr] 
!   dfxcdgdn  d(fxc)/d(|grad(ndn)|) [hartree * bohr] 
!
! Note that dfxcdgtot = d(fxc)/d(|grad(nup+ndn)|) = 0      
!
!**********************************************************************

  subroutine am05wbs(nup, ndn, gup, gdn, fx, fc, &
        dfxup, dfxdn, dfcup, dfcdn, dfxdgup, dfxdgdn, dfcdgup, dfcdgdn)

  implicit none

! ** Input parameters
  real(dp) nup, gup
  real(dp) ndn, gdn

! ** Output parameters
  real(dp) fx, fc, dfxup, dfxdn, dfcup, dfcdn
  real(dp) dfxdgup,dfxdgdn,dfcdgup,dfcdgdn

! ** Internal parameters
  real(dp) kFup, sup
  real(dp) kFdn, sdn
  real(dp) ex,ec 

  integer pot
  real(dp) pi
  parameter (pi = 3.141592653589793238462643383279502884197_dp)

! ** Dummy parameters (not used)
  real(dp) uup, udn, tup, tdn, vxup, vxdn, vcup, vcdn, stot

  kFup = (3.0_dp*pi**2*2.0_dp*nup)**(1.0_dp/3.0_dp)
  kFdn = (3.0_dp*pi**2*2.0_dp*ndn)**(1.0_dp/3.0_dp)
  call saferecp(sup,gup,(2.0_dp*kFup*nup))
  call saferecp(sdn,gdn,(2.0_dp*kFdn*ndn))

! ** Not needed input for pot=1
  pot = 1 
  stot = 0.0_dp
  uup = 0.0_dp
  udn = 0.0_dp
  tup = 0.0_dp
  tdn = 0.0_dp

  call am05_xscss(nup,ndn,sup,sdn,stot,uup,udn,tup,tdn, &
                  fx,fc,vxup,vxdn,vcup,vcdn, &
                  dfxup,dfxdgup,dfcup,dfcdgup, &
                  dfxdn,dfxdgdn,dfcdn,dfcdgdn,pot)

!     ** Form the output 
!     For SIESTA don't multiply by gradient since we need to divide
!     by this quantity on return to get components.
!      dfxdgup = gup*dfxdgup
!      dfxdgdn = gdn*dfxdgdn
!      dfcdgup = gup*dfcdgup
!      dfcdgdn = gdn*dfcdgdn

  return
  end subroutine am05wbs

!**********************************************************************
! am05
! Calculate the Armiento Mattsson AM05 exchange-correlation energy 
! functional and various functional derivatives for different
! potential schemes.
!
! Spin polarization type: exchange index scales with the separate 
! spin densities; correlation index and prefactor both scales with
! the separate spin densities. (xscss)
!
! Input:
!   nup     electron upspin density [bohr**(-3)]
!   ndn     electron downspin density [bohr**(-3)]
!
!   sup      scaled gradient of upspin-density
!   sdn      scaled gradient of downspin-density
!   stot     scaled gradient of total density
!
!   uup      scaled grad nup * grad | grad nup |
!   udn      scaled grad ndn * grad | grad ndn |
!
!   tup      scaled laplacian of upspin-density
!   tdn      scaled laplacian of downspin-density
!
!   pot      integer: 
!              2 = calculate potential in the traditional scheme
!                  (all input needed and all output well defined)
!              1 = calculate quantities for the White and Bird,
!                  and the Pople, Gill, and Johnson schemes for
!                  potentials (u and t and stot are never touched
!                  and vx and vc give undefined output)
!              0 = don't calculate potential (u and t and stot are never 
!                  touched and only ex and ec output are well defined)
!
! Output:
!
!   fx          exchange energy density [hartree], 
!                 Total exchange Ex = Integrate[fx]
!   fc          correlation energy density [hartree]
!                 Total correlation Ec = Integrate[fc]
!    
! if pot = 1 or 2:
!
!   dfxup     d(fx)/d(nup)
!   dfxdn     d(fx)/d(ndn)
!   dfcup     d(fc)/d(nup)
!   dfcdn     d(fc)/d(ndn)
!   dfxdgup   d(fx)/d(|grad(nup)|) * 1/|grad(nup)|
!   dfxdgdn   d(fx)/d(|grad(ndn)|) * 1/|grad(ndn)|
!   dfcdgup   d(fc)/d(|grad(nup)|) * 1/|grad(nup)|
!   dfcdgdn   d(fc)/d(|grad(ndn)|) * 1/|grad(ndn)|
!
! Note that dfxcdgtot = d(fxc)/d(|grad(nup+ndn)|) * 1/|grad(nup+ndn)| = 0
!
! if pot = 2:
!
!   vxup      upspin exchange potential
!   vxdn      downspin exchange potential
!   vcup      upspin correlation potential
!   vcdn      downspin correlation potential
!
! Citation request: when using this functional, please cite:
! "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
!  R. Armiento and A. E. Mattsson (unpublished)."
!
! (The first paper for the AM05 functional, the second for the
! spin-polarized version)
!
!**********************************************************************
  subroutine am05_xscss(nup,ndn,sup,sdn,stot,uup,udn,tup,tdn, &
                        fx,fc,vxup,vxdn,vcup,vcdn, &
                        dfxup,dfxdgup,dfcup,dfcdgup, &
                        dfxdn,dfxdgdn,dfcdn,dfcdgdn,pot)

! ** Input parameters
  real(dp) nup,ndn,sup,sdn,stot,uup,udn,tup,tdn
  integer pot

! ** Output parameters
  real(dp) fx,fc,vxup,vxdn,vcup,vcdn
  real(dp) dfxup,dfxdgup,dfcup,dfcdgup
  real(dp) dfxdn,dfxdgdn,dfcdn,dfcdgdn

! ** Constants
  real(dp) pi, g, a, c
  parameter (pi = 3.141592653589793238462643383279502884197_dp)
  parameter (g = 0.8098_dp, a = 2.804_dp)
  parameter (c = 0.7168_dp)

! ** Local variables
  real(dp) s2, kF, nn, ntot
  real(dp) n(2), s(2), t(2), u(2)
  real(dp) dfx(2), dfxdg(2), dfc(2), dfcdg(2), vx(2), vc(2)
  real(dp) exlda, vxlda, eclda, vclda(2), X(2), Xsos, Xsossos
  real(dp) Hx, Hxsos, Hxsossos, Hc(2), Hcsos, Hcsossos
  real(dp) F, Fsos, s2Fsossos 
  real(dp) szsoz, mixder
  real(dp) denom, denomsos, sdenomsoss
  real(dp) zfac, zosn, w

  integer i

! ** Initialization
  fx = 0.0_dp
  fc = 0.0_dp

  n(1) = nup
  n(2) = ndn
  s(1) = sup
  s(2) = sdn
  t(1) = tup
  t(2) = tdn
  u(1) = uup
  u(2) = udn

  ntot = DMAX1(n(1)+n(2),1e-16+n(2),n(1)+1e-16)

  do i = 1,2

    dfx(i) = 0.0_dp
    dfxdg(i) = 0.0_dp
    dfc(i) = 0.0_dp 
    dfcdg(i) = 0.0_dp 
    vx(i) = 0.0_dp
    vc(i) = 0.0_dp

!   ** Avoid floating point exceptions
    if(s(i) .ge. 1.0d12) then
      X(i) = 0.0_dp
      Hc(i) = g
    elseif(s(i) .le. 1.0d-30) then
      X(i) = 1.0_dp
      Hc(i) = 1.0_dp
    else
      X(i) = 1.0_dp/(1.0_dp + a*s(i)**2)
      Hc(i) = X(i) + g*(1.0_dp - X(i))
    endif
  enddo

! *******************
!   LDA correlation
! *******************
  call am05_xscss_ldapwc(n,eclda,vclda)

! ** Loop over spin up and down
  do 20 i = 1,2

! *************
!    Cutoffs
! *************
    if (n(i) .le. 1.0d-16) then
!
!     ** specifically handle a small density in the LDA limit (s=0) n->0
!     ** A small s in AM05 is always indication that we are in the
!     ** interior LDA limit, regardless of what density we have
      if (s(i) .le. 1.0d-30) then
! The following line is the original code
!       fc = fc + n(i)*eclda
! The following line is in the form needed for SIESTA
        fc = fc + n(i)*eclda/ntot
        vc(i) = vclda(i)*Hc(3-i) + eclda*(1-Hc(3-i))
        dfc(i) = vclda(i)*Hc(3-i) + eclda*(1-Hc(3-i))
        goto 20
      endif

!     ** otherwise, go on as usual but with a lowest cutoff density
      n(i) = 1.0d-16
    endif

!   ** Scaling density and gradient
!   ** Ex = 0.5 * (Ex[2*spin up n] + Ex[2*spin down n] ) 
    nn = 2.0_dp*n(i)
    s2 = s(i)**2
    kF = (3.0_dp*pi**2*2.0_dp*n(i))**(1.0_dp/3.0_dp)
     
!   *******************
!     LDA exchange
!   *******************
    call am05_xscss_ldax(nn,exlda,vxlda)

!   *****************************
!     Exchange energy density
!   *****************************

!   ** Airy LAA refinement function
    call am05_xscss_lambertw(s(i)**(3.0_dp/2.0_dp)/sqrt(24.0_dp),w)

!   ** am05_lambertw give back argument if it is < 1.0e-20
!   ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
!   ** zosn = normalized z/s
    if (s(i) < 1.0e-14) then
      zosn = 1.0_dp
    else
      zosn = 24.0_dp**(1.0_dp/3.0_dp)*w**(2.0_dp/3.0_dp)/s(i)
    end if
    zfac = s2*(zosn*27.0_dp/32.0_dp/pi**2)**2

!   ** denom = denominator of Airy LAA refinement function
    denom = 1.0_dp + c*s2*zosn*(1.0_dp + zfac)**(1.0_dp/4.0_dp)
    F = (c*s2 + 1.0_dp)/denom
      
!   ** Exchange refinement function
    Hx = X(i) + (1.0_dp - X(i))*F

!   ** Exchange energy density, Ex = Integrate[fx]
    fx = fx + 0.5_dp*nn*exlda*Hx

!   ********************
!      Correlation
!   ********************

!   ** Correlation energy density, Ec = Integrate[fc]
! The following line is the original code
!   fc = fc + 0.5_dp*nn*eclda*Hc(i)
! The following line is in the form needed for SIESTA
   fc = fc + 0.5_dp*nn*eclda*Hc(i)/ntot

!   ** goto next spin if we are only calculating energies
    if (pot .eq. 0) goto 20

!   ***************************
!     Exchange derivatives for White and Bird and 
!     Pople, Gill, and Johnson schemes
!   ***************************

!   ** Interpolation index derivatives: 1/s dX/ds
    Xsos = -2.0_dp*a*X(i)**2

!   ** Airy LAA refinement function derivatives, 1/s dF/ds 
!   ** szsoz = s*(dz/ds)/z
    szsoz = 1.0_dp/(1.0_dp + w)
         
    Fsos = c/denom**2*(2.0_dp - zosn* &
              ((1.0_dp - c*s2)*(1.0_dp + zfac)**(1.0_dp/4.0_dp) + &
              (1.0_dp + c*s2)*(1.0_dp + 3.0_dp/2.0_dp*zfac)/ &
              (1.0_dp + zfac)**(3.0_dp/4.0_dp)*szsoz))

!   ** Refinement function derivatives, 1/s dHx/ds
!   ** We use that (1 - X) = a*X*s2
    Hxsos = (1.0_dp - X(i))*Fsos - (F - 1.0_dp)*Xsos

    dfx(i) = vxlda*Hx - 4.0_dp/3.0_dp*exlda*s2*Hxsos
    dfxdg(i) = exlda*Hxsos/((2.0_dp*kF)**2*n(i))

!   *****************************
!     Correlation derivatives for White and Bird and 
!     Pople, Gill, and Johnson schemes
!   *****************************

!   ** Correlation refinement function derivatives, 1/s dF/ds 
    Hcsos = Xsos*(1.0_dp - g)

!   ** Note: n(3-i) gives the density of the *other* spin, etc.
    dfc(i) = vclda(i)/ntot*(n(1)*Hc(1)+n(2)*Hc(2)) + &
              eclda*n(3-i)/ntot*(Hc(i) - Hc(3-i)) -  &
              4.0_dp/3.0_dp*eclda*s2*Hcsos

    dfcdg(i) = eclda*Hcsos/((2.0_dp*kF)**2*n(i))

!   ** goto next spin if only doing W&B/Pople et. al.
    if (pot .eq. 1) goto 20

!   ***************************
!     Exchange potential in traditional scheme
!   ***************************
!
!   ** Interpolation index derivatives: 1/s d/ds(1/s dX/ds)                   
    Xsossos = 8.0_dp*a**2*X(i)**3

!   ** Airy LAA refinement function derivatives s^2 1/s d/ds (1/s dF/ds) 
!   **  mixder = szsoz + s^2*(d^2z/ds^2)/z
    mixder = (2.0_dp - w)/(2.0_dp*(1.0_dp + w)**3)

!   ** denomsos = 1/s d(denom)/ds,  sdenomsoss = s*d/ds(1/s d(denom)/ds))
    denomsos = c*zosn/(1.0_dp + zfac)**(3.0_dp/4.0_dp)* &
              (1.0_dp + zfac + (1.0_dp + 3.0_dp/2.0_dp*zfac)*szsoz)

    sdenomsoss = c*zosn/(1.0_dp + zfac)**(7.0_dp/4.0_dp)* &
              (-1.0_dp - zfac*(2.0_dp + zfac) &
              + (1.0_dp + zfac/2.0_dp*(5.0_dp + 3.0_dp*zfac))*mixder &
              + 3.0_dp/2.0_dp*zfac*(1.0_dp + zfac/2.0_dp)*szsoz**2)

    s2Fsossos = (-4.0_dp*c*s2*denom*denomsos + (c*s2 + 1.0_dp)* &
              (2.0_dp*s2*denomsos**2 - denom*sdenomsoss))/denom**3

!   ** Refinement function derivatives 1/s d/ds (1/s dHx/ds) 
!   ** We use that (1 - X) = a*X*s2
    Hxsossos = - 2.0_dp*Fsos*Xsos + a*X(i)*s2Fsossos - (F - 1.0_dp)*Xsossos

!   ** vx formula for gradient dependent functional,
    vx(i) = vxlda*(Hx - s2*Hxsos) + &
              exlda*((4.0_dp/3.0_dp*s2-t(i))*Hxsos +  &
              (4.0_dp/3.0_dp*s(i)**3-u(i))*s(i)*Hxsossos)

!   *****************************
!     Correlation potential in traditional scheme
!   *****************************
         
!   ** Correlation refinement function derivatives, 1/s d/ds (1/s dHc/ds)
    Hcsossos = Xsossos*(1.0_dp - g)

!   ** vc formula for gradient dependent functional,
    vc(i) = vclda(i)*(Hc(i) - s2*Hcsos) + &
              eclda*((4.0_dp/3.0_dp*s2 - t(i))*Hcsos +  &
              (4.0_dp/3.0_dp*s(i)**3-u(i))*s(i)*Hcsossos) + &
           (eclda-vclda(i))*n(3-i)/ntot* &
              (Hc(i)-Hc(3-i)-s(i)**2*Hcsos)+ &
           (eclda-vclda(3-i))* &
              0.5_dp* &
             (stot**2*ntot**(8.0_dp/3.0_dp)*0.5_dp**(2.0_dp/3.0_dp) -  &
                s(1)**2*n(1)**(8.0_dp/3.0_dp) -  &
                s(2)**2*n(2)**(8.0_dp/3.0_dp))/ &
             (n(i)**(5.0_dp/3.0_dp))/ntot*Hcsos

 20 continue

  dfxup = dfx(1)
  dfxdgup = dfxdg(1)
  dfcup = dfc(1)
  dfcdgup = dfcdg(1)
  dfxdn = dfx(2)
  dfxdgdn = dfxdg(2)
  dfcdn = dfc(2)
  dfcdgdn = dfcdg(2)

  vxup = vx(1) 
  vxdn = vx(2)
  vcup = vc(1) 
  vcdn = vc(2)

  return
  end subroutine am05_xscss

! ******************************************
!   Local density approximation exchange
!
!   input
!   n        electron density [bohr**(-3)]
!
!   output
!   ex       exchange energy per electron [hartree]
!   vx       exchange potential [hartree]
!
!   Copyright (c) 2005, Rickard Armiento
! ******************************************

  subroutine am05_xscss_ldax(n,ex,vx)
! ** Input parameters
  real(dp) n

! ** Output parameters
  real(dp) ex, vx

! ** Constants
  real(dp) pi
  parameter (pi = 3.141592653589793238462643383279502884197_dp)

  vx = -(3.0_dp*n/pi)**(1.0_dp/3.0_dp)
  ex = (3.0_dp/4.0_dp*vx)

  return
  end subroutine am05_xscss_ldax

! ***********************************************
! Local density approximation correlation
!
! input
! n(2)     electron upspin/downspin density [bohr**(-3)]
!
! output
! ec       correlation energy per electron [hartree]
! vc(2)    correlation upspin/downspin potential [hartree]
!
! As parameterized by Perdew Wang,
!   Phys. Rev. B 45, 13244 (1992) 
! Based on Monte Carlo data by Ceperley Alder, 
!   Phys. Rev. Lett. 45, 566 (1980)
!
! (Clean room implementation from paper)
!
! Copyright (c) 2005, Rickard Armiento
! ***********************************************
  subroutine am05_xscss_ldapwc(n,ec,vc)

! ** Input parameters
  real(dp) n(2)

! ** Output parameters
  real(dp) ec, vc(2)

! ** Constants
  real(dp) pi
  real(dp) A0,a01,b01,b02,b03,b04
  real(dp) A1,a11,b11,b12,b13,b14
  real(dp) Aa,aa1,ba1,ba2,ba3,ba4
  parameter (pi = 3.141592653589793238462643383279502884197_dp)
  parameter (a01 = 0.21370_dp)
  parameter (b01 = 7.5957_dp)
  parameter (b02 = 3.5876_dp)
  parameter (b03 = 1.6382_dp)
  parameter (b04 = 0.49294_dp)
  parameter (a11 = 0.20548_dp)
  parameter (b11 = 14.1189_dp)
  parameter (b12 = 6.1977_dp)
  parameter (b13 = 3.3662_dp)
  parameter (b14 = 0.62517_dp)
  parameter (aa1 = 0.11125_dp)
  parameter (ba1 = 10.357_dp)
  parameter (ba2 = 3.6231_dp)
  parameter (ba3 = 0.88026_dp)
  parameter (ba4 = 0.49671_dp)
! ** Paper actually use this:
! parameter (A0 = 0.031091_dp)
! parameter (A1 = 0.015545_dp)
! parameter (Aa = 0.016887_dp)
! ** But routines now "defacto standard" was distributed using:
  parameter (A0 = 0.0310907_dp)
  parameter (A1 = 0.01554535_dp)
  parameter (Aa = 0.0168869_dp)

! ** Local variables
  real(dp) xi, rsq, f, fp, fb0, mac, ec0, ec1, ecrs, ecxi
  real(dp) Q0, Q1, Q1p, ec0rs, ec1rs, acrs, fdenom

! ** Actual values from paper
! fdenom = (2.0_dp**(4.0_dp/3.0_dp)-2.0_dp)
! fb0 = 4.0_dp/(9.0_dp*(2.0_dp**(1.0_dp/3.0_dp)-1.0_dp))
! ** Replaced with "defacto standard" approximations
  fdenom = 0.5198421_dp
  fb0 = 1.709921_dp

! ** Cutoff
  if((n(1)+n(2)) .le. 1e-30) then
    ec = 0.0_dp
    vc(1) = 0.0_dp
    vc(2) = 0.0_dp
    return
  endif

  xi = (n(1)-n(2))/(n(1)+n(2))
  rsq = (3.0_dp/(4.0_dp*pi*(n(1)+n(2))))**(1.0_dp/6.0_dp)
  f = ((1.0_dp+xi)**(4.0_dp/3.0_dp)+(1.0_dp-xi)**(4.0_dp/3.0_dp)-2.0_dp)/fdenom

  mac = -2.0_dp*Aa*(1.0_dp + aa1*rsq**2)*log(1.0_dp + 1.0_dp/ &
           (2.0_dp*Aa*rsq*(ba1 + rsq*(ba2 + rsq*(ba3 + ba4*rsq)))))

  ec0 = -2.0_dp*A0*(1.0_dp + a01*rsq**2)*log(1.0_dp + 1.0_dp/ &
           (2.0_dp*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))))

  ec1 = -2.0_dp*A1*(1.0_dp + a11*rsq**2)*log(1.0_dp + 1.0_dp/ &
           (2.0_dp*A1*rsq*(b11 + rsq*(b12 + rsq*(b13 + b14*rsq)))))

  ec = ec0 - mac*f/fb0*(1.0_dp-xi**4) + (ec1-ec0)*f*xi**4

  Q0 = -2.0_dp*A0*(1.0_dp + a01*rsq**2)
  Q1 = 2.0_dp*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))
  Q1p = A0*(b01/rsq+2.0_dp*b02+3.0_dp*b03*rsq+4.0_dp*b04*rsq**2)
  ec0rs = -2.0_dp*A0*a01*log(1.0_dp + 1.0_dp/Q1)-Q0*Q1p/(Q1**2+Q1)

  Q0 = -2.0_dp*A1*(1.0_dp + a11*rsq**2)
  Q1 = 2.0_dp*A1*rsq*(b11 + rsq*(b12 + rsq*(b13 + b14*rsq)))
  Q1p = A1*(b11/rsq+2.0_dp*b12+3.0_dp*b13*rsq+4.0_dp*b14*rsq**2)
  ec1rs = -2.0_dp*A1*a11*log(1.0_dp + 1.0_dp/Q1)-Q0*Q1p/(Q1**2+Q1)

  Q0 = -2.0_dp*Aa*(1.0_dp + aa1*rsq**2)
  Q1 = 2.0_dp*Aa*rsq*(ba1 + rsq*(ba2 + rsq*(ba3 + ba4*rsq)))
  Q1p = Aa*(ba1/rsq+2.0_dp*ba2+3.0_dp*ba3*rsq+4.0_dp*ba4*rsq**2)
  acrs = -2.0_dp*Aa*aa1*log(1.0_dp + 1.0_dp/Q1)-Q0*Q1p/(Q1**2+Q1)

  ecrs = ec0rs*(1.0_dp - f*xi**4) + ec1rs*f*xi**4 - acrs*f/fb0*(1.0_dp - xi**4)
  fp = 4.0_dp/3.0_dp*((1.0_dp+xi)**(1.0_dp/3.0_dp) - (1.0_dp-xi)**(1.0_dp/3.0_dp))/fdenom
  ecxi = 4.0_dp*xi**3*f*(ec1-ec0+mac/fb0)+fp*(xi**4*ec1-xi**4*ec0+ (1.0_dp-xi**4)*(-mac)/fb0)

  vc(1)=ec - rsq**2/3.0_dp*ecrs - (xi-1.0_dp)*ecxi
  vc(2)=ec - rsq**2/3.0_dp*ecrs - (xi+1.0_dp)*ecxi
  
  end subroutine am05_xscss_ldapwc


! ***********************************************
! LambertW function. 
!
! Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), 
!   Adv. in Comp. Math. 5(4):329-359. 
! Implementation approach loosely inspired by the 
! GNU Octave version by N. N. Schraudolph, but this 
! implementation is only for real values and 
! principal branch.
!
! Copyright (c) 2005, Rickard Armiento
! ***********************************************

  subroutine am05_xscss_lambertw(z,result)
! input
  real(dp) z
! output
  real(dp) result
! local variables
  real(dp) e,t,p    
  integer i

! ** If z too low, go with the first term of the power expansion, z
  if( z .lt. 1.0d-20) then
    result = z
    return
  endif

  e = exp(1.0_dp)

! ** Inital guess
  if( abs(z + 1.0_dp/e) .gt. 1.45_dp ) then
!   ** Asymptotic expansion at 0 and Inf
    result = log(z)
    result = result - log(result)
  else
!   ** Series expansion about -1/e to first order
    result = 1.0_dp*sqrt(2.0_dp*e*z + 2.0_dp) - 1.0_dp
  endif

! ** Find result through iteration
  do i = 1,10
    p = exp(result)
    t = result*p - z
    if ( result .ne. -1.0_dp ) then
      t = t/(p*(result + 1.0_dp) - 0.5_dp*(result + 2.0_dp)*t/(result + 1.0_dp))
    else
      t = 0.0_dp
    endif
    result = result - t
    if (abs(t) < (2.48_dp*1.0d-14)*(1.0_dp + abs(result))) then
      return
    endif
  enddo

! ** This should never happen!;
  call die('am05_xscss_lambertw: iteration limit reached.')

  return
  end subroutine am05_xscss_lambertw

end module am05
