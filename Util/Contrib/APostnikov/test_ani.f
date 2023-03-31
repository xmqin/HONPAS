C...............................................................
C
      subroutine test_ani(ii2,nat,nstep)        
C
C     reads from ii2 (ANI file), makes some consistency tests,
C     returns the number of MD steps.
C
      implicit none
      integer ii2,nat,nstep,na,iat,ii
      logical varcel
      double precision dummy
C     double precision coord(3,nat)
      character symbol*2

      rewind (ii2)
C --- read through input file, to get the total number of MD steps.
C     It must be specified at the beginning of AXSF file.
C       read from unformatted MD, as written in subr. pixmol (coord. in Ang):
      nstep = 1
  101 continue
      read (ii2,301,err=106,end=107) na
      if (na.ne.nat) then
        write (6,*) ' Error reading ANI file, step=',nstep,
     .              ' : na=',na,' not matching Nr of atoms in XV =',nat
        stop
      endif
      do iat=1,na
C       read (ii2,302,err=108,end=109) symbol,(coord(ii,iat),ii=1,3)
        read (ii2,302,err=108,end=109) symbol,(dummy,ii=1,3)
      enddo
      nstep = nstep + 1
      goto 101
  106 continue
C --- Error reading ANI file: 
      write (6,*) ' Error reading header of record ',nstep,
     .            ' in the ANI file.'
      nstep = nstep - 1  !  No. of full records
      if (nstep.gt.0) then
        write (6,*) ' Keep ',nstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif
  107 continue
C --- regular end of records ANI file: 
      nstep = nstep - 1  !  Attempt to read record Nr. nstep failed
      write (6,*) ' Cleanly read in ',nstep,'  MD steps'
      return
  108 continue
C --- Error reading ANI file: 
      write (6,*) ' Error reading coordinates block in record ',
     .              nstep,', for atom ',iat,' in the ANI file.'
      nstep = nstep - 1  !  No. of full records
      if (nstep.gt.0) then
        write (6,*) ' Keep ',nstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif
  109 continue
C --- irregular end of records in ANI file: 
      write (6,*) ' Uncomplete record in MD step No.',
     .              nstep,', for atom ',iat,' in the ANI file.'
      nstep = nstep - 1  !  No. of full records
      if (nstep.gt.0) then
        write (6,*) ' Keep ',nstep,' records.'
        return
      else
        write (6,*) ' Check the ANI file; bye'
        stop
      endif

  301 format(i5,/)
  302 format(a2,2x,3f12.6)

      end 
