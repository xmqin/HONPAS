! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine uncell(rv,a,b,c,alpha,beta,gamma)
      implicit real*8(a-h,o-z)
      include 'constants'
      dimension rv(3,3),temp(6)
C
C  Convert cell vectors to parameters
C
C  Julian Gale, April 1992
C
      do i=1,3
        temp(i)=0.0d0
        do j=1,3
          temp(i)=temp(i)+rv(j,i)**2
        enddo
        temp(i)=sqrt(temp(i))
      enddo
      a=temp(1)
      b=temp(2)
      c=temp(3)
      do i=1,3
        temp(3+i)=0.0d0
      enddo
      do j=1,3
        temp(4)=temp(4)+rv(j,2)*rv(j,3)
        temp(5)=temp(5)+rv(j,1)*rv(j,3)
        temp(6)=temp(6)+rv(j,1)*rv(j,2)
      enddo
      temp(4)=temp(4)/(temp(2)*temp(3))
      temp(5)=temp(5)/(temp(1)*temp(3))
      temp(6)=temp(6)/(temp(1)*temp(2))
      alpha=radtodeg*acos(temp(4))
      beta=radtodeg*acos(temp(5))
      gamma=radtodeg*acos(temp(6))
C
C  Avoid round off errors for 90.0 and 120.0 degrees
C
      if (abs(alpha-90.0).lt.0.00001) alpha=90.0
      if (abs(alpha-120.0).lt.0.00001) alpha=120.0
      if (abs(beta-90.0).lt.0.00001) beta=90.0
      if (abs(beta-120.0).lt.0.00001) beta=120.0
      if (abs(gamma-90.0).lt.0.00001) gamma=90.0
      if (abs(gamma-120.0).lt.0.00001) gamma=120.0
      return
      end
