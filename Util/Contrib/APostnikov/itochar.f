!

      character*10 function itochar(i)

C Returns a character string with a copy of input integer 
C with spaces filling on the right, if the integer
C has less than 10 digits.
C
C P. Ordejon, July 2003

      implicit none

      integer i, io, idiv, itest, iascii, ipos, iscan, nvac
      logical null
      character char*10


      null=.true.
      nvac=0

      io=i
      do iscan = 9, 0, -1
        idiv = 10**iscan
        itest = io/idiv
        iascii = itest+48
        ipos = 10-iscan
        if (null) then 
          if (achar(iascii) .eq. '0') then 
            nvac = nvac + 1
            goto 10
          else 
            null = .false.
          endif
        endif
        char(ipos:ipos) = achar(iascii)
10      continue
        io = io - itest*idiv 
      enddo

      if (nvac .eq. 0)  then
        itochar = char
      else
        itochar(10-nvac+1:10) = ' '
        itochar(1:10-nvac) = char(nvac+1:10)
      endif

      end

