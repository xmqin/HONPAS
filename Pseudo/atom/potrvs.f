c
      subroutine potrvs(vd,r,nr,k)
c
c  Generates file for plotting of Screened ionic pseudopotentials
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonable plot.
c
C     .. Scalar Arguments ..
      integer k, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*10
C     ..
      write(filename,9900) k
 9900 format('SCRPSPOTR',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')

c     Write out r and the screened Vps(r) 
c
      step = 0.0D0
      do 10 j = 5, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j)
            step = step + 0.01D0
         end if
   10 continue
 9000 format(1x,f7.4,3x,f10.5)
c
      close(unit=3)
      return
c
      end


