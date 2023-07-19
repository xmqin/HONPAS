c
c     Copyright (c) 1993, 1995, 1997 
c     Alberto Garcia, wdpgaara@lg.ehu.es
c
c     Redistribution and use, with or without modification, are 
c     permitted provided that the above copyright notice is retained.
c
c Symbols ------
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
      if (dummy_string .ne. '#') backspace(unit_no)

      return
      end
c
      subroutine getline(unit_no,string)

      implicit none

c     Provides the caller with a line that does not contain comments
c     or directives.
c
c     If it encounters the end of file, returns '#'
c
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

      read(unit_no,'(a132)',end=999) line

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
      RETURN
 999  continue
      string = '#'
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








