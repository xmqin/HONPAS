c
      subroutine potrw(vd,r,nr,k,kj,ist,rc)
c
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonalble plot.
c
C     .. Parameters ..
      double precision zero, pzf
      parameter (zero=0.D0,pzf=0.01D0)
C     ..
C     .. Scalar Arguments ..
      integer ist, k, kj, nr
      double precision rc
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*7
C     ..
c
      if (kj .eq. 0) then
         write(filename,9900) k
 9900 format('PSWFNR',i1)
      else if (kj .eq. -1) then
         write(filename,9930) k
 9930 format('PTWFNR',i1)
      else
         write(filename,9910) k
 9910 format('AEWFNR',i1)
      endif
c
      open(unit=3,file=filename,form='formatted',status='unknown')
c
c     Write out r, the wavefunction, and rc (kludge to pass it to
c     the plotting program)
c
      step = zero
      do 10 j = 2, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j)*ist, rc
            step = step + pzf
         end if
   10 continue
 9000 format(1x,f7.4,3x,f12.8,2x,f8.4)
c
      close(unit=3)
cag
      return
c
      end


