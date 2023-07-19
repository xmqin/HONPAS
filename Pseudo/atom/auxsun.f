c
c $Id: auxsun.f,v 1.4 2002/07/11 21:47:00 wdpgaara Exp $
c
c $Log: auxsun.f,v $
c Revision 1.4  2002/07/11 21:47:00  wdpgaara
c Fix year problem
c
c Revision 1.3  1997/05/22 17:32:01  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1992/02/28  00:34:47  alberto
c New date routine.
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
          subroutine cal_date(bdate)
c
c    gets the date (day-MONTH-year)
c    sun   version   (native fuction idate used)
c
          character*10 bdate
          character*3 month(12)
          character*1 dash
          data dash/'-'/
c
          integer iarray(3)
c
      data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
c
          call idate(iarray)
          write(bdate,102) iarray(1),dash,month(iarray(2)),dash,
     &                     iarray(3)
 102  format(i2,a1,a3,a1,i2.2,' ')
          return
          end

