C
c $Id: ext.f,v 1.2 1997/05/22 17:32:13 wdpgaara Exp $
c
c $Log: ext.f,v $
c Revision 1.2  1997/05/22 17:32:13  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
C
C
      subroutine ext(i)
c
c  Stops program in case of errors or completion.
c
c  i is a stop parameter
c   000-099 main (0 is normal exit)
c   100-199 input
c   200-299 charge
c   300-399 vionic
c   400-499 velect
c   500-599 dsolv1
c   600-699 dsolv2 (including difnrl and difrel)
c   700-799 etotal
c   800-899 pseudo, pseudk, pseudt and pseudv
c
C     .. Scalar Arguments ..
      integer i
C     ..
      if (i .ne. 0) write(6,FMT=9000) i
 9000 format('stop parameter =',i3)
      close(unit=1)
      close(unit=3)
      close(unit=5)
      close(unit=6)
c
      stop
c
      end
