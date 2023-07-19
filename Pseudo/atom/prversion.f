      subroutine prversion
c
c     Simple routine to print the version string. Could be extended to
c     provide more information, if needed.
c
c     Alberto Garcia, Feb. 23, 1998
c
      implicit none
      include 'version.h'

      write(6,'(a)') version

      end
