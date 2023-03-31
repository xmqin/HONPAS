! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine wtof(string,rnum,nlen,ndp)
      implicit real*8(a-h,o-z)
C
C  Convert string to floating point number
C
C  nlen = length of string
C  ndp  = position of decimal point
C
C   9/92 Created
C
C  Julian Gale, Imperial College, March 1997
C
      character*1 c
      character*40 string
      logical lneg,lnegp
      npower=0
      i0=ichar('0')
      ineg=ichar('-')
      ipos=ichar('+')
C
C  Look for exponentiation
C
      lnegp=.false.
      nexp=index(string,'e')
      if (nexp.eq.0) nexp=index(string,'E')
      if (nexp.eq.0) nexp=index(string,'d')
      if (nexp.eq.0) nexp=index(string,'D')
      if (nexp.gt.0.and.nexp.lt.nlen) then
C
C  Exponential found
C
        nse=nexp+1
        if (index(string(nse:nse),'+').eq.1) then
          nse=nse+1
        elseif (index(string(nse:nse),'-').eq.1) then
          nse=nse+1
          lnegp=.true.
        endif
        nfct=1
        do i=nlen,nse,-1
          c=string(i:i)
          ic=ichar(c)
          n=ic-i0
          npower=npower+nfct*n
          nfct=nfct*10
        enddo
        nlen=nexp-1
        if (lnegp) npower=-npower
      endif
      rnum=0.0d0
      if (nlen.eq.0) return
      factor=0.1d0
      if (ndp.eq.0) ndp=nlen+1
      nend=1
      c=string(1:1)
      ic=ichar(c)
      if (ic.eq.ineg.or.ic.eq.ipos) nend=2
      lneg=(ic.eq.ineg)
      if (lneg) nend=2
C
C  In front of decimal point
C
      do i=ndp-1,nend,-1
        c=string(i:i)
        ic=ichar(c)
        n=ic-i0
        factor=factor*10.0d0
        rnum=rnum+factor*n
      enddo
C
C  After decimal point
C
      factor=1.0d0
      do i=ndp+1,nlen
        c=string(i:i)
        ic=ichar(c)
        n=ic-i0
C
C  The following condition is to avoid problems with tab characters
C
        if (n.ge.0.and.n.le.9) then
          factor=factor*0.1d0
          rnum=rnum+factor*n
        endif
      enddo
      if (lneg) rnum=-rnum
      rnum=rnum*(10.0d0**npower)
      return
      end
