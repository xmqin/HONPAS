c
          subroutine cal_date(bdate)
c
c    gets the date (day-MONTH-year)
c
          character*10 bdate
          character*3 month(12)
          character*1 dash
          
          integer values(8)
c
      data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
      data dash/'-'/
c
      call date_and_time(values=values)

          write(bdate,102) values(3),dash,month(values(2)),dash,
     &                     mod(values(1),100)
 102  format(i2,a1,a3,a1,i2.2,' ')

          end
c
      double precision function second() 
c
      real tt
      call cpu_time(tt)
      second = dble(tt)
      end
c

