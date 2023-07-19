c
c $Id: string.f,v 1.3 1999/02/26 14:26:47 wdpgaara Exp $
c
c $Log: string.f,v $
c Revision 1.3  1999/02/26 14:26:47  wdpgaara
c Cosmetic changes.
c
c Revision 1.2  1997/05/22  17:32:33  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      SUBROUTINE CHRLEN(STRING,NCHAR,LCHAR)
C
C***********************************************************************
C
C  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
C  the length of the string up to the last nonblank, nonnull.
C
      CHARACTER STRING*(*)
      Integer nchar, lchar
c
      Integer i, ncopy 
C
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
C
      DO 10 I=1,NCOPY
        LCHAR=NCOPY+1-I
        IF(STRING(LCHAR:LCHAR).NE.' '.AND.
     *     STRING(LCHAR:LCHAR).NE.CHAR(0))RETURN
 10     CONTINUE
      LCHAR=0
      RETURN
      END
c
      subroutine chrcap(string,nchar)
C
C***********************************************************************
C
C  CHRCAP accepts a STRING of NCHAR characters and replaces
C  any lowercase letters by uppercase ones.
C
C
C     .. Scalar Arguments ..
      integer nchar
      character string*(*)
C     ..
C     .. Local Scalars ..
      integer i, itemp, ncopy
C     ..
C     .. Intrinsic Functions ..
      intrinsic char, ichar, len, lge, lle
C     ..
      ncopy = nchar
      if (ncopy .le. 0) ncopy = len(string)
      do 10 i = 1, ncopy
C
         if (lge(string(i:i),'a') .and. lle(string(i:i),'z')) then
            itemp = ichar(string(i:i)) + ichar('A') - ichar('a')
            string(i:i) = char(itemp)
         end if
   10 continue
c
      return
c
      end
c
      logical function leqi(strng1,strng2)
C
C***********************************************************************
C
C  Case-insensitive lexical equal-to comparison
C
C
C     .. Scalar Arguments ..
      character strng1*(*), strng2*(*)
C     ..
C     .. Local Scalars ..
      integer i, len1, len2, lenc
      character s1*1, s2*1
C     ..
C     .. External Subroutines ..
      external chrcap
C     ..
C     .. Intrinsic Functions ..
      intrinsic len, min
C     ..
      len1 = len(strng1)
      len2 = len(strng2)
      lenc = min(len1,len2)
C
      leqi = .FALSE.
      do 10 i = 1, lenc
         s1 = strng1(i:i)
         s2 = strng2(i:i)
         call chrcap(s1,1)
         call chrcap(s2,1)
         if (s1 .ne. s2) return
   10 continue
C
      if (len1 .gt. lenc .and. strng1(lenc+1:len1) .ne. ' ') return
      if (len2 .gt. lenc .and. strng2(lenc+1:len2) .ne. ' ') return
      leqi = .TRUE.
c
      return
c
      end
c
      subroutine loc_des(message)
c
c     Processes message to locate the delimiters % and $ used
c     in the warnp routine.
c
c     Alberto Garcia, Feb 1, 1991
c
C     .. Scalar Arguments ..
      character message*(*)
C     ..
C     .. Scalars in Common ..
      integer form_length
      character form_spec*200
C     ..
C     .. Local Scalars ..
      integer dol_pos, pct_pos
      character work*200
C     ..
C     .. Intrinsic Functions ..
      intrinsic index
C     ..
C     .. Common blocks ..
      common /fordes/form_spec, form_length
C     ..
      pct_pos = index(message,'%')
c
      form_spec = form_spec(1:form_length)//message(1:(pct_pos-1))//
     &            ''','
      form_length = form_length + (pct_pos-1) + 2
      work = message(pct_pos+1:)
c
      dol_pos = index(work,'$')
      form_spec = form_spec(1:form_length)//work(1:dol_pos-1)//','''
      form_length = form_length + (dol_pos-1) + 2
c
c        Return the rest of message
c
      message = ' '
      message = work(dol_pos+1:)
c
      return
c
      end
