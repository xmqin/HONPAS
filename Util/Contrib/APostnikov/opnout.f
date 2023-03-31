      subroutine opnout(io1,outfil)
C
C     Opens output file, preventing its accidental overrid
C
      implicit none
      integer io1
      character outfil*60,owrite*1
      logical filexist

      inquire (file=outfil, exist=filexist)
      if (filexist) then
        write (6,*) ' File ',trim(outfil),' exists.',
     .              ' Overwrite? (Y/N)'
        read (5,*) owrite
        if (owrite.eq.'Y'.or.owrite.eq.'y') then
          open (io1,file=outfil,form='formatted',
     .              status='REPLACE',err=801)
        else
          write (6,*) ' Then rename is first.  Bye...'
          stop
        endif
      else
        open (io1,file=outfil,form='formatted',status='NEW',err=802)
      endif
      return

  801 write (6,*) ' Error replacing formatted file ',trim(outfil)
      stop
  802 write (6,*) ' Error opening file ',trim(outfil),
     .            ' as new formatted'
      stop

      end

