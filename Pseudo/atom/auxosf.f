c
          subroutine cal_date(bdate)
c
c    gets the date (day-MONTH-year)
c    sun   version   (native fuction idate used)
c
          character*10 bdate
          character*3 month(12)
          character*1 dash
          
          external idate

          integer iarray(3)
c
      data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
      data dash/'-'/
c
          call idate(iarray)
          write(bdate,102) iarray(1),dash,month(iarray(2)),dash,
     &                     mod(iarray(3),100)
 102  format(i2,a1,a3,a1,i2.2,' ')
          return
          end
c
      double precision function second() 
c
      real tarray(2)
c
      external etime
      real etime
c
      second = dble(etime(tarray))
      return
c
      end
c

