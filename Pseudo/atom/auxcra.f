c
c $Id: auxcra.f,v 1.3 1997/05/22 17:32:00 wdpgaara Exp $
c
c $Log: auxcra.f,v $
c Revision 1.3  1997/05/22 17:32:00  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1992/02/28  00:33:57  alberto
c New date routine.
c .
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
c  ********************************************************
c  *                                                      *
c  *   njtj                                               *
c  *     These are machine dependent routines.            *
c  *   Included are routine for Apollo, Sun,               *
c  *   Vax, and Cray systems.  The user must              *
c  *   1)compile with their systems lines uncommented     *
c  *   or 2)supply their own                              *
c  *   or 3)remove-comment out all references to          *
c  *   these calls in the program.                        *
c  *                                                      *
c  ********************************************************
c
C
cag      subroutine zedate(bdate)
c
c   Apollo version, this should
c   work on all apollo machines.
c
c   Gets the data (DAY-MONTH-YEAR)
c
cag% include '/sys/ins/base.ins.ftn'
cag% include '/sys/ins/cal.ins.ftn'
cag      character*10 bdate
cag      integer*2 loctim(6),lyear
cag      character*3 month(12)
cag      character*1 dash
cag      data dash/'-'/
cag      data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
cag     1 'OCT','NOV','DEC'/
cag      call cal_$decode_local_time(loctim)
cag      lyear = loctim(1) - (loctim(1)/100)*100
cag      write(bdate,100) loctim(3),dash,month(loctim(2)),dash,lyear
cag 100  format(i2,a1,a3,a1,i2,' ')
cag      return
cag      end
C
C **************Apollo end*************************
C
          SUBROUTINE CAL_DATE(BDATE)
C
C    GETS THE DATE (DAY-MONTH-YEAR)
C    CRAY   VERSION
C
          CHARACTER*10 BDATE
          CHARACTER*8 ADATE
          CHARACTER*3 MONTH(12)
          CHARACTER*1 DASH,DUM1,DUM2
          DATA DASH/'-'/
      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
C
          WRITE(ADATE,100) DATE()
          READ(ADATE,101) LMONTH,DUM1,LDAY,DUM2,LYEAR
          WRITE(BDATE,102) LDAY,DASH,MONTH(LMONTH),DASH,LYEAR
 100  FORMAT(A8)
 101  FORMAT(I2,A1,I2,A1,I2)
 102  FORMAT(I2,A1,A3,A1,I2,' ')
          RETURN
          END
C
C  *****************Cray end***********************
C
C   Gets the data (DAY-MONTH-YEAR)
C   VAX version
C
C     .. Scalar Arguments ..
c     character bdate*10
C     ..
C     .. Local Scalars ..
c     character adate*9
C     ..
c     call date(adate)
c     write(bdate,FMT=9000) adate
c9000 format(a9,' ')
c
c     return
c
c     end



