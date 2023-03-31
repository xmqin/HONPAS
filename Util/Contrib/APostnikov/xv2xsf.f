C
C    xf2xsf,  a quick script to reformat the XV file 
C             with lattice/coordinates into a XCrysden xsf file.
C
C             Written by Andrei Postnikov, Mar 2006   Vers_0.3
C             apostnik@uos.de
C
      program xv2xsf
      implicit none
      integer ii1,io1
      parameter (ii1=11,io1=14)
      integer ii,jj,iat,nat,ityp,nz
      character inpfil*60,outfil*60,syslab*30,suffix*6
      double precision b2ang,cc_bohr(3,3),cc_ang(3,3),coord(3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      external opnout
C
C     string manipulation functions in Fortran used below:
C     len_trim(string): returns the length of string 
C                       without trailing blank characters,
C     trim(string)    : returns the string with railing blanks removed
      
      write (6,701,advance="no")
  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
C --- open output file:
C     outfil = syslab(1:len_trim(syslab))//'.XSF'
      outfil = trim(syslab)//'.XSF'
      call opnout(io1,outfil)
C --  read in translation vectors, convert into Ang:
      do ii=1,3
        read  (ii1,*,end=803,err=803)  (cc_bohr(jj,ii),jj=1,3)
      enddo
      cc_ang = cc_bohr*b2ang

C     write (io1,*) '# crystal structure from ',
C    .              syslab(1:len_trim(syslab)),'.XV'
      write (io1,*) '# crystal structure from ',trim(syslab),'.XV'
      write (io1,*) 'CRYSTAL'
      write (io1,*) '# Cell vectors in Angstroem:'
      write (io1,*) 'PRIMVEC'
      do ii=1,3
        write  (io1,'(3f12.7)')  (cc_ang(jj,ii),jj=1,3)
      enddo
      write (io1,*) '# Atom coordinates in Angstroem:'
      write (io1,*) 'PRIMCOORD'
      read  (ii1,*,end=804,err=804)  nat
      write (io1,'(2i5)') nat,1
      do iat=1,nat
        read  (ii1,*,end=805,err=805) ityp, nz, (coord(ii),ii=1,3)
        write (io1,'(i4,3f20.8)') nz, (coord(ii)*b2ang,ii=1,3)
      enddo

      close (ii1)
      close (io1)
      stop

  801 write (6,*) ' Error opening file ',
     .              trim(inpfil),' as old formatted'
      stop
  802 write (6,*) ' Error opening file ',
     .              trim(outfil),' as new formatted'
      stop
  803 write (6,*) ' End/Error reading XV for cell vector ',ii
      stop
  804 write (6,*) ' End/Error reading XV for number of atoms line'
      stop
  805 write (6,*) ' End/Error reading XV for atom number ',iat
      stop
      end

