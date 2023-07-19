c
c Symbols ------
c
c    Alberto Garcia, February 1993--April 1995--October 1995
c
c     This package filters the input stream, detecting and acting upon
c     directives of the form
c
c     %define NAME
c     %delete NAME
c     %NAME = value
c     %show
c
c     and ignoring comment lines (with a '#' in column one).
c
c     The only user-callable routines are:
c
c        check_directives ( call check_directives(unit_no: integer) )
c        getline    ( call getline(unit_no: integer, line: string) )
c        set_value  ( new_value = set_value(NAME: string, default: string))
c        defined    ( if (defined('NAME')) ...   )
c        insert     [experimental...]
c
c
c     Check_directives reads comment and directive lines, if present,
c     leaving the file ready for further processing.
c
c     Getline will provide the caller with the next input line that is
c     not a comment or a directive.
c
c     Set_value will return the value associated with the name NAME, or
c     a default value (which MUST BE PASSED AS A STRING) if NAME is not
c     in the symbol table.
c
c     %define and %delete are implemented as special cases of
c     the more general form
c
c     %NAME = value
c
c     by simply associating a value of '1' to the symbol to be %define(d)
c     and a value of '0' to the symbol to be %delete(d), and the routine
c     Defined simply calls Set_value.
c
c     
      subroutine check_directives(unit_no)
      
      implicit none

c     Reads any lines containing comments or directives, and readies
c     the file connected to unit_no for reading the next block of
c     information

      integer unit_no
      character*10 dummy_string

      call getline(unit_no,dummy_string)
      backspace(unit_no)

      return
      end
c
      subroutine getline(unit_no,string)

      implicit none

c     Provides the caller with a line that does not contain comments
c     or directives.

      integer unit_no
      character*(*) string

      character*132 line
      external directive

c 'line' is used as an internal file (it should be big enough to
c        accomodate requests of different sizes).

c Special lines:
c
c ---> Comment lines:
c      The character '#' must appear in column 1
c ---> A line containing '%' in column 1 is treated as a 
c      directive.
c
c      Blank lines are skipped
c

  10  continue

      read(unit_no,'(a132)') line

c  ...skip comment lines ( those include blank lines )

      if (line .eq. ' ' .or. line(1:1) .eq. '#') then
        write(6,'(a78)') line
c
        go to 10
c
      endif
c
c  ...process directives
c
      if (line(1:1) .eq. '%') then
c
	 write(6,'(a78)') line
         call directive(line(2:))
c
         go to 10
c
      endif
c
      string = line

      return
c
      end   
  
c
      subroutine directive(str)

      implicit none

      character*(*) str
c
      character*30 name, val_str
      integer eqloc
      
      logical success
      logical leqi, put_pair
      external leqi, put_pair
c
      success = .true.

      eqloc = index(str,'=')
      
      if (eqloc .ne. 0) then
         name = str(1:eqloc-1)
         val_str = str(eqloc+1:)
         success = put_pair(name,val_str)
      else 
         name = str(8:)
         if (leqi(str(1:6),'define')) then
            success = put_pair(name,'1.')
         else if (leqi(str(1:6),'delete')) then
            success = put_pair(name,'0.')
         else if (leqi(str(1:4),'show')) then
            call print_pairs
         else
            write(6,'(/,a,1x,a,/)') ' Unrecognized directive:', str
         endif

      endif
      
      if (.not. success) write(6,'(/,a,2x,a,/)') 
     &     'Set full... Could not process ',name
      
      return
      end
c
      logical function defined(name)
c
      implicit none
c
      character*(*) name
      double precision set_value
      external set_value
c
      defined = (nint(set_value(name,'0')) .eq. 1)
c
      return
c
      end
c
      subroutine insert(module,name)
c
      implicit none
c
      character*(*) name, module
      logical success
      logical put_pair
c
      success = put_pair(name,'1.')
      write(6,*) ' <---- defining ', name, ' in ',module, success
c
      return
c
      end
c
      logical function put_pair(name,num_str)
	
      implicit none

      character*(*) name, num_str

      include 'set2.h'
      
      double precision value
      integer i

      logical leqi
      external leqi, get_real

      call get_real(num_str,value)

      do 10 i = 1, nels

         if (leqi(name,el_name(i))) then
            put_pair = .true.
            el_val(i) = value
cag            write(6,*) name, ' reset. New value: ',num_str
            return
         endif
 10      continue

         if (nels .ne. nmax) then
            nels = nels + 1
            el_name(nels) = name
            el_val(i) = value
            put_pair = .true.
cag            write(6,*) name, ' set to: ', num_str
         endif

         return

         end

         double precision function set_value(name,default_str)

         implicit none
         
         character*(*) name, default_str
         double precision default_value

         include 'set2.h'

         logical leqi
         external leqi
	 external get_real

         integer i
c
         call get_real(default_str,default_value)
         set_value = default_value

         do 10 i=1, nels
            if (leqi(el_name(i),name)) set_value = el_val(i)
 10      continue

         return

         end
         
c     
      subroutine get_real(string,value)
         
c     Alberto Garcia, April 2, 1995
c     
c     It turns out that compiler behavior is not uniform when it comes
c     to interpreting real input format statements. This routine takes
c     the string "string" and extracts the value it represents as a real
c     number. String can contain the letters e,d,E,D signifying an
c     exponent, and embedded blanks everywhere. 
c     The bn descriptor is used in case some compiler's default behavior is
c     to pad the mantissa or exponent fields with significant zeroes.
c     
c     I have tested the routine against all kinds of *reasonable* strings.
c     Please let me know if you find any bugs.
c     
         implicit none
         
         character*(*) string
         double precision value
c     
c     80 should be enough...
c     
         character*80 mantissa, exponent
c     
         integer exp_loc, exp_val
c     
c     Paranoid section.
c     Make sure that the routine still works if the case of the source
c     is changed! Sorry, ASCII only...
c  
      character*1 low_e, cap_e, low_d, cap_d
      
      cap_d = char(68)
      cap_e = char(69)
      low_d = char(100)
      low_e = char(101)

cdebug      write(6,*) cap_d, cap_e, low_d, low_e
c
c     We are counting on finding only one exponent... (reasonable?)
c

      exp_loc = index(string,cap_e) + index(string,cap_d) +
     +          index(string,low_e) + index(string,low_d)

      if (exp_loc .ne. 0) then
	mantissa = string(1:exp_loc-1)
	exponent = string(exp_loc+1:)
      else
	mantissa = string
      endif
c
c    BN means "Treat blanks as no-significant"
c
	read(mantissa,'(bn,f80.0)') value
	
	if (exp_loc .ne. 0) then
		read(exponent,'(bn,i80)') exp_val
		value = value * 10.d0 ** exp_val
	endif

	return

	End
c
      subroutine print_pairs
c
      implicit none
c
      include 'set2.h'
c
      integer i
c
      write(6,9000) nels, nmax
 9000 format(/,1x,'Symbol set contains ',i2,' pairs,',
     +         ' out of a maximum of ',i3,':')
c 
      write(6,9010) (el_name(i),el_val(i),i=1,nels)
 9010 format(1x,a30,2x,g25.15)
      write(6,'(/)')
c
      return
c
      end    
c
c $Id: symbols.f,v 1.1 1997/05/22 18:12:35 wdpgaara Exp $
c
c $Log: symbols.f,v $
c Revision 1.1  1997/05/22 18:12:35  wdpgaara
c *** empty log message ***
c
c Revision 1.1.1.1  1997/01/07 08:38:56  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1995/07/01  17:10:43  wdpgaara
c Initial revision
c
c Revision 1.2  1993/10/07  17:37:43  garcia
c Check-in 10/7/93
c
c Revision 1.1  91/12/16  00:07:11  alberto
c Initial revision
c 
      SUBROUTINE CHRLEN(STRING,NCHAR,LCHAR)
C
C***********************************************************************
C
C  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
C  the length of the string up to the last nonblank, nonnull.
C
      CHARACTER CHAR*1
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
