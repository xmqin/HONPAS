c $Id: auxibm.f,v 1.5 2002/07/11 21:47:00 wdpgaara Exp $
c
c     Auxiliary file for the plane-wave program. IBM RS/6000 Version
c     It is now very close to "BSD" flavors, thanks to the introduction
c     by IBM of idate_, itime_, flush_, etime_, and friends.
C
          SUBROUTINE cal_date(BDATE)
C
C    GETS THE DATE (DAY-MONTH-YEAR)
C    IBM RS/6000   VERSION
c
          CHARACTER*10 BDATE
          CHARACTER*3 MONTH(12)
          CHARACTER*1 DASH
c
          integer iarray(3)
          integer day, mon, year
c
          DATA DASH/'-'/
      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
C
      call idate_(iarray)

      day = iarray(1)
      mon = iarray(2)
      year = mod(iarray(3),100)

          WRITE(BDATE,102) day,DASH,MONTH(mon),DASH,year
 102  FORMAT(I2,A1,A3,A1,I2.2,' ')
          RETURN
          END
c
      subroutine time_of_day(a)
c
c     IBM RS/6000 version
c
      character a*10
      integer iarray(3)
c
          integer hour, min, sec
c
          call itime_(iarray)

          hour = iarray(1)
          min = iarray(2)
          sec = iarray(3)
c
      write(a,10) hour, min, sec 
10    format(1x,i2,':',i2.2,':',i2.2,1x)
c
      return
c
      end
c
      double precision function ran(i)
c
c     Emulates DEC FORTRAN ran(i) on the IBM RS/6000
c
C     .. Scalar Arguments ..
      integer i
c
      real rand
C     ..
      ran = rand() 
c
      return
c
      end
c
      subroutine sort(n,arrin,indx)
c     sorts an array by the heapsort method
c     w. h. press et al. numerical recipes
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision arrin(n)
      integer indx(n)
C     ..
C     .. Local Scalars ..
      double precision q
      integer i, indxt, ir, j, l
C     ..
      do 10 j = 1, n
         indx(j) = j
 10   continue
      l = n/2 + 1
      ir = n
 20   continue
      if (l .gt. 1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir .eq. 1) then
            indx(1) = indxt
c
            return
c
         end if
      end if
      i = l
      j = l + l
 30   continue
      if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
c
         go to 30
c
      end if
      indx(i) = indxt
c
      go to 20
c
      end
c
      logical function env_var(name,value)
c
c     IBM RS/6000 Version . Alberto Garcia, June 1991
c
      character*(*) name, value
c
      integer valuelg
c
      value = ' '
      call getenv(name,value)
c
c     Test the length of the returned string
c
      call ChrLen(value,0,valuelg)
      env_var = (valuelg .gt. 0)
c
      return
c
      end
c
      subroutine flush(n)
      integer n
C_Uncomment when it works properly...      call flush_(n)
      return
      end

      double precision function second() 
c
       real tarray(2)
c
       external etime_
       real etime_
c
       second = dble(etime_(tarray))
       return
c
       end


