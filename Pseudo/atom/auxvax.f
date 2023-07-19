c
c $Id: auxvax.f,v 1.3 1997/05/22 17:32:02 wdpgaara Exp $
c
c $Log: auxvax.f,v $
c Revision 1.3  1997/05/22 17:32:02  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1992/02/28  00:36:18  alberto
c New date routine.
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
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
c  ****************Apollo start***********************
c
cag      subroutine zesec(t)
c
c   Apollo version, this should
c   work on all apollo machines.
c
c   Gets the cpu time is seconds
c
cag% INCLUDE '/sys/ins/base.ins.ftn'
cag% INCLUDE '/sys/ins/cal.ins.ftn'
cag% INCLUDE '/sys/ins/proc1.ins.ftn'
cag       REAL*8 T
cag       CHARACTER*6 CLOCK
cag       CALL PROC1_$GET_CPUT(CLOCK)
cag       CALL CAL_$FLOAT_CLOCK(CLOCK,T)
cag       RETURN
cag       END
C
C **************Apollo end*************************
C
C **************Cray start***********************
C
Cray      SUBROUTINE ZESEC(T)
C
C   GETS CPU TIME IN SECONDS
C   CRAY-2 VERSION
C
Cray      T = SECOND()
Cray      RETURN
Cray      END
C
C  *****************Vax start**********************
C
      subroutine zesec(t)
C
C   CALCULATES THE ELAPSED CPU TIME SINCE
C   THE FIRST CALL IN A VAX/VMS SYSTEM
C
C     .. Scalar Arguments ..
      double precision t
C     ..
C     .. Scalars in Common ..
      integer iflag
C     ..
C     .. Local Scalars ..
      integer its
C     ..
C     .. External Subroutines ..
      external lib$init_timer, lib$stat_timer
C     ..
C     .. Intrinsic Functions ..
      intrinsic dble
C     ..
C     .. Common blocks ..
      common /zesec2/iflag
C     ..
C     .. Data statements ..
c     The value 0 is stored only at compile time.
      data iflag/0/
C     ..
      if (iflag .eq. 0) then
         call lib$init_timer
         iflag = 1
         t = 0.0D0
      else
         call lib$stat_timer(2,its)
         t = 0.01D0*dble(its)
      end if
c
      return
c
      end
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
Cray      SUBROUTINE ZEDATE(BDATE)
C
C    GETS THE DATE (DAY-MONTH-YEAR)
C    CRAY-2 VERSION
C
Cray      CHARACTER*10 BDATE
Cray      CHARACTER*8 ADATE
Cray      CHARACTER*3 MONTH(12)
Cray      CHARACTER*1 DASH,DUM1,DUM2
Cray      DATA DASH/'-'/
Cray      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
Cray     2  'OCT','NOV','DEC'/
C
Cray      WRITE(ADATE,100) DATE()
Cray      READ(ADATE,101) LMONTH,DUM1,LDAY,DUM2,LYEAR
Cray      WRITE(BDATE,102) LDAY,DASH,MONTH(LMONTH),DASH,LYEAR
Cray 100  FORMAT(A8)
Cray 101  FORMAT(I2,A1,I2,A1,I2)
Cray 102  FORMAT(I2,A1,A3,A1,I2,' ')
Cray      RETURN
Cray      END
C
C  *****************Cray end***********************
      subroutine cal_date(bdate)
C
C   Gets the data (DAY-MONTH-YEAR)
C   VAX version
C
C     .. Scalar Arguments ..
      character bdate*10
C     ..
C     .. Local Scalars ..
      character adate*9
C     ..
      call date(adate)
      write(bdate,FMT=9000) adate
 9000 format(a9,' ')
c
      return
c
      end


