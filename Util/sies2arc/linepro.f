! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine linepro(line_in,nword,words,nfloat,floats,norder,
     *  iline,maxword)
      implicit real*8(a-h,o-z)
C
C  Subroutine searches a line for words and numbers
C  Now amended to handle full 80 character string.
C
C   2/95 Ability to handle continuation lines added
C   4/98 Now reads from channel 4 as input is preprocessed
C
C  Julian Gale, Imperial College, February 1995
C
      common /warning/ nwarn
      dimension words(maxword),floats(maxword),norder(2*maxword)
      character*81 line
      character*80 line_in
      character*30 words
      character*40 string,blank
      character*1 c1,c2,space,dp,hash,slash
      logical lspace,lfirst,lread,lword
      data blank/'                                        '/
      dp='.'
      hash='#'
      space=' '
      slash='/'
      i0=ichar('0')
      i9=ichar('9')
      ineg=ichar('-')
      ipos=ichar('+')
      ia=ichar('a')
      iz=ichar('z')
      ica=ichar('A')
      icz=ichar('Z')
      nword=0
      nfloat=0
      iorder=0
      lread=.false.
   10 if (lread) then
        read(4,'(a)',end=100)line_in
        iline=iline+1
      endif
      line=line_in
      line(81:81)=' '
      lfirst=(line(1:1).ne.space)
      lread=.true.
C*******************
C  Locate strings  *
C*******************
      ilower=1
      iupper=index(line,'#')
      if (iupper.eq.0) iupper=80
      do while (ilower.lt.iupper)
        lspace=.false.
        i=ilower-1
C
C  Find space before alphanumeric character
C
        if (lfirst) then
          i=0
          nstart=1
          lfirst=.false.
          lspace=.true.
        else
          do while (.not.lspace.and.i.lt.iupper)
            i=i+1
            c1=line(i:i)
            c2=line(i+1:i+1)
            if (c1.eq.space.and.c2.ne.space) then
              lspace=.true.
            endif
          enddo
          nstart=i+1
        endif
C
C  Find last alphanumeric character
C
        do while (lspace.and.i.lt.iupper)
          i=i+1
          c1=line(i:i)
          c2=line(i+1:i+1)
          if (c1.ne.space.and.(c2.eq.space.or.c2.eq.hash)) then
            lspace=.false.
          endif
        enddo
        nend=i
        nlen=nend-nstart+1
        if (nstart.lt.iupper) then
C**********************************************
C  Decide what type of string has been found  *
C**********************************************
          string=blank
          string=line(nstart:nend)
          c1=string(1:1)
          if (c1.eq.dp) then
            c1=string(2:2)
          endif
          ic=ichar(c1)
          lword=((ic.lt.i0.or.ic.gt.i9).and.ic.ne.ineg.and.ic.ne.ipos)
C
C  Are there more than two letters in the string - if so this is a word
C
          nsw=0
          do i=1,nlen
            ii=ichar(string(i:i))
            if (ii.ge.ia.and.ii.le.iz) nsw=nsw+1
            if (ii.ge.ica.and.ii.le.icz) nsw=nsw+1
          enddo
          if (nsw.ge.2) lword = .true.
C
          if (lword) then
C*********
C  Word  *
C*********
            nword=nword+1
            if (nword.gt.maxword) then
              write(6,
     *'(/,''  **** Error - too many words - increase maxword ****'',/)')
              stop
            endif
            words(nword)=string
            iorder=iorder+1
            norder(iorder)=1
          else
C***********
C  Number  *
C***********
            nfloat=nfloat+1
            if (nfloat.gt.maxword) then
              write(6,
     *'(/,''  **** Error - too many words - increase maxword ****'',/)')
              stop
            endif
            islash=(index(string,slash))
            if (islash.ne.0) then
C
C  Fractions
C
              ndp=islash
            else
C
C  Floating point
C
              ndp=index(string,'.')
            endif
            call wtof(string,rnum,nlen,ndp)
            floats(nfloat)=rnum
            iorder=iorder+1
            norder(iorder)=2
          endif
        endif
C
C  Set values for next search
C
        ilower=nend+1
      enddo
      return
  100 write(6,'(/,a,a,/)')"  **** Warning - continuation character",
     *  " on last line of input ****"
      nwarn=nwarn+1
      return
      end
      subroutine linepronc(line_in,nword,words,nfloat,floats,
     *  norder,iline,maxword)
      implicit real*8(a-h,o-z)
C
C  Subroutine searches a line for words and numbers
C  Now amended to handle full 80 character string.
C
C  Version of linepro which ignores continuation characters for
C  use by firstpass.f
C
C   2/95 Ability to handle continuation lines added
C
C  Conditions of use:
C
C  GULP is available free of charge to academic institutions
C  and non-commerical establishments only. Copies should be
C  obtained from the author only and should not be distributed
C  in any form by the user to a third party without the express
C  permission of the author. This notice applies to all parts
C  of the program, except any library routines which are
C  distributed with the code for completeness. All rights for
C  such routines remain with the original distributor.
C
C  No claim is made that this program is free from errors and
C  no liability will be accepted for any loss or damage that
C  may result. The user is responsible for checking the validity
C  of their results.
C
C  Copyright Imperial College 1997
C
C  Julian Gale, Imperial College, February 1995
C
      dimension words(maxword),floats(maxword),norder(2*maxword)
      character*81 line
      character*80 line_in
      character*30 words
      character*40 string,blank
      character*1 c1,c2,space,dp,hash,slash
      logical lspace,lfirst,lword,lread
      data blank/'                                        '/
      dp='.'
      hash='#'
      space=' '
      slash='/'
      i0=ichar('0')
      i9=ichar('9')
      ineg=ichar('-')
      ipos=ichar('+')
      ia=ichar('a')
      iz=ichar('z')
      ica=ichar('A')
      icz=ichar('Z')
      nword=0
      nfloat=0
      iorder=0
      line=line_in
      line(81:81)=' '
      lfirst=(line(1:1).ne.space)
C*******************
C  Locate strings  *
C*******************
      ilower=1
      iupper=index(line,'#')
      if (iupper.eq.0) iupper=80
      do while (ilower.lt.iupper)
        lspace=.false.
        i=ilower-1
C
C  Find space before alphanumeric character
C
        if (lfirst) then
          i=0
          nstart=1
          lfirst=.false.
          lspace=.true.
        else
          do while (.not.lspace.and.i.lt.iupper)
            i=i+1
            c1=line(i:i)
            c2=line(i+1:i+1)
            if (c1.eq.space.and.c2.ne.space) then
              lspace=.true.
            endif
          enddo
          nstart=i+1
        endif
C
C  Find last alphanumeric character
C
        do while (lspace.and.i.lt.iupper)
          i=i+1
          c1=line(i:i)
          c2=line(i+1:i+1)
          if (c1.ne.space.and.(c2.eq.space.or.c2.eq.hash)) then
            lspace=.false.
          endif
        enddo
        nend=i
        nlen=nend-nstart+1
        if (nstart.lt.iupper) then
C**********************************************
C  Decide what type of string has been found  *
C**********************************************
          string=blank
          string=line(nstart:nend)
          c1=string(1:1)
          if (c1.eq.dp) then
            c1=string(2:2)
          endif
          ic=ichar(c1)
          lword=((ic.lt.i0.or.ic.gt.i9).and.ic.ne.ineg.and.ic.ne.ipos)
C
C  Are there more than two letters in the string - if so this is a word
C
          nsw=0
          do i=1,nlen
            ii=ichar(string(i:i))
            if (ii.ge.ia.and.ii.le.iz) nsw=nsw+1
            if (ii.ge.ica.and.ii.le.icz) nsw=nsw+1
          enddo
          if (nsw.ge.2) lword = .true.
C
          if (lword) then
C*********
C  Word  *
C*********
            nword=nword+1
            if (nword.gt.maxword) then
              write(6,'(/,2a,/)')
     *          "  **** Error - too many words - incr",
     *          "ease maxword ****"
              stop
            endif
            words(nword)=string
            iorder=iorder+1
            norder(iorder)=1
          else
C***********
C  Number  *
C***********
            nfloat=nfloat+1
            if (nfloat.gt.maxword) then
              write(6,'(/,2a,/)')
     *          "  **** Error - too many words - incr",
     *          "ease maxword ****"
              stop
            endif
            islash=(index(string,slash))
            if (islash.ne.0) then
C
C  Fractions
C
              ndp=islash
            else
C
C  Floating point
C
              ndp=index(string,'.')
            endif
            call wtof(string,rnum,nlen,ndp)
            floats(nfloat)=rnum
            iorder=iorder+1
            norder(iorder)=2
          endif
        endif
C
C  Set values for next search
C
        ilower=nend+1
      enddo
      return
      end
