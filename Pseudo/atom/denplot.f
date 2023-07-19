c
      subroutine denplot
c
c  Prints the charge density
c     
      include 'radial.h'
      include 'charge.h'
      include 'param.h'
c
      double precision pi
      parameter (pi=3.141592653589d0)
c
      double precision fx, step, delta
      integer j
c
c     Minimum step-size for plotting. This cuts down on file size.
c     Set it to zero to recover old behavior.
c
      parameter (delta = 0.000d0)
ccc      parameter (delta = 0.005d0)
c
c     Specify name according to type of job.
c
      if (job .eq. 4) then
         open(unit=3,file='PTCHARGE',form='formatted',status='unknown')
      else
         open(unit=3,file='AECHARGE',form='formatted',status='unknown')
      endif
c
c     Keep for backwards compatibility
c
      open(unit=66,file='CHARGE',form='formatted',status='unknown')
c
      open(unit=4,file='RHO',form='formatted',status='unknown')
c
c     Write out r, cdu, cdd and cdc
c
      step = 0.0d0
      do 10 j = 2, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), cdu(j), cdd(j), cdc(j)
            write(66,9000) r(j), cdu(j), cdd(j), cdc(j)
            fx = 1.d0 / ( 4.d0 * pi * r(j)**2)
            write(4,9000) r(j), fx*cdu(j), fx*cdd(j), fx*cdc(j)
            step = step + delta
         endif
   10 continue
 9000 format(1x,f15.10,3x,3f18.5)
c
      close(unit=3)
      close(unit=66)
      close(unit=4)
c
      return
c
      end

