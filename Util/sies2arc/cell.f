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
      subroutine cell(rv,a,b,c,alpha,beta,gamma)
      implicit real*8(a-h,o-z)
      include 'constants'
      dimension rv(3,3)
C
C  Convert cell parameters into cartesian frame
C
C  Julian Gale, Sept 1991
C
      if (alpha.eq.90.0) then
        cosa=0.0d0
      else
        alp=alpha*degtorad
        cosa=cos(alp)
      endif
      if (beta.eq.90.0) then
        cosb=0.0d0
      else
        bet=beta*degtorad
        cosb=cos(bet)
      endif
      if (gamma.eq.90.0) then
        sing=1.0d0
        cosg=0.0d0
      else
        gam=gamma*degtorad
        sing=sin(gam)
        cosg=cos(gam)
      endif
      rv(2,1)=0.0d0
      rv(3,1)=0.0d0
      rv(3,2)=0.0d0
      rv(1,1)=a
      rv(1,2)=b*cosg
      rv(2,2)=b*sing
      rv(1,3)=c*cosb
      rv(2,3)=c*(cosa-cosg*cosb)/sing
      trm1=rv(2,3)/c
      rv(3,3)=c*sqrt(1.0d0-cosb**2-trm1**2)
      return
      end
