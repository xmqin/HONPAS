c
      subroutine potrv(vd,r,nr,k,zion)
c
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonable plot.
c
C     .. Scalar Arguments ..
      integer k, nr
      double precision zion
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*7
C     ..
cag
      write(filename,9900) k
 9900 format('PSPOTR',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
cag
c     Write out r, V(r) and the Coulomb potential (for reference)
c
      step = 0.0D0
      do 10 j = 5, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j), -2*zion/r(j)
            step = step + 0.01D0
         end if
   10 continue
 9000 format(1x,f7.4,3x,f10.5,g20.10)
c
      close(unit=3)
cag
      return
c
      end
