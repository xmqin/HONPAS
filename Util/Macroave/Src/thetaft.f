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
	subroutine thetaft(n,L,lav,ft)
C
C Calculates the Fourier coefficients of the convoluting function
C    w(x) = (1/lav) theta(lav/2 - abs(x))
C asuming it is periodic with period L:
C
C       |                                          |
C  1/lav|__________                      __________|
C       |          |                     |         |
C       |__________|_____________________|_________|___
C       0         lav/2               L-lav/2     L
C

        implicit none

	integer n,j
        real*8 L,lav
	real*8 ft(2*n),x,phi
        real*8 pi

	pi=4.0d0*datan(1.0d0)

        ft(1)=n/L
	ft(2)=0.0

        do j=2,n/2+1
          ft(2*(j-1)+1) = (n/(pi*lav*(j-1)/2.))*
     .                    sin(pi*lav*(j-1)/2./L)*
     .                    cos(pi*(L-lav/2.)*(j-1)/L)*
     .                    cos(pi*L*(j-1)/L)
          ft(2*(j-1)+2) = (n/(pi*lav*(j-1)/2.))*
     .                    sin(pi*lav*(j-1)/2./L)*
     .                    cos(pi*(L-lav/2.)*(j-1)/L)*
     .                    sin(pi*L*(j-1)/L)
        enddo

        do j=n/2+2,n
          ft(2*(j-1)+1) = (n/(pi*lav*(-n+j-1)/2.))*
     .                    sin(pi*lav*(-n+j-1)/2./L)*
     .                    cos(pi*(L-lav/2.)*(-n+j-1)/L)*
     .                    cos(pi*L*(-n+j-1)/L)
          ft(2*(j-1)+2) = (n/(pi*lav*(-n+j-1)/2.))*
     .                    sin(pi*lav*(-n+j-1)/2./L)*
     .                    cos(pi*(L-lav/2.)*(-n+j-1)/L)*
     .                    sin(pi*L*(-n+j-1)/L)
        enddo


	return
	end
