! ---
! Copyright (C) 1996-2014 The SIESTA group
! This file is distributed under the terms of the
! GNU General Public License: see COPYING in the top directory
! or http:
! See Docs/Contributors.txt for a list of contributors.
! ---
!
!
      subroutine fft3d(f, Mesh, isn)
C
C 3d FFT based on FFT routine of Temperton
C
C On input :
C
C real*8 f() : contains data to be Fourier transformed
C integer Mesh(3) : contains dimensions of grid
C integer isn : indicates the direction of the transform
C
C On exit :
C
C real*8 f() : contains data Fourier transformed
C
C Julian Gale, July 1999
C Serial version simplified by A. Garcia, May 2015
C
C
!     Use Temperton's GPFA package
!
      use m_fft_gpfa, only: fft_gpfa

      implicit none

      integer, parameter :: dp = selected_real_kind(10,100)
C
C Passed arguments
C
      real(dp) ::  f(*)
      integer  ::  Mesh(3), isn
C
C Local variables
C
      integer  :: n1, n2, n3, i, n, IOffSet
      real(dp) ::  scale
C
C Set mesh size variables
C
      n1 = Mesh(1)
      n2 = Mesh(2)
      n3 = Mesh(3)
      n = n1*n2*n3

C
C FFT in X direction
C
      call fft_gpfa(f,f(2),2,2*n1,n1,n2*n3,isn)
C
C FFT in Y direction
C
      do i=0,n3-1
         IOffSet=2*n1*n2*i
         call fft_gpfa(f(IOffSet+1),f(IOffSet+2),
     . 2*n1,2,n2,n1,isn)
      enddo
C
C FFT in Z direction
C
      call fft_gpfa(f,f(2),2*n1*n2,2,n3,n1*n2,isn)
C
C Scale values
C
      if (isn.gt.0) then
        scale=1.0_dp/dble(n)
        do i=1,2*n
          f(i)=f(i)*scale
        enddo
      endif
      end
