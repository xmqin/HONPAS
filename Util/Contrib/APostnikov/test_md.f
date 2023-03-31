C...............................................................
C
      subroutine test_md(ii2,nat,varcel,nstep)        
C
C     reads from ii2 (MD file), makes some consistency tests,
C     finds out whether the records are for fixed cell or variable cell,
C     returns the number of MD steps.
C
C     MD file is written in iomd.f as follows (for formt = F ):
C       write(iupos) istep, xa(1..3,1..nat), va(1..3,1..nat)
C       if ( varcel ) write(iupos) cell(1..3,1..3), vcell(1..3,1..3)
C
      implicit none
      integer ii2,nat,nstep,istep,istep1,istep2,ii
      logical varcel
      double precision dummy

      rewind (ii2)
C --- assume MD was written with fixed cell...
      read (ii2,err=801,end=107) istep1,(dummy,ii=1,6*nat)
      read (ii2,err=102,end=107) istep2,(dummy,ii=1,6*nat)
      if (istep2.eq.istep1+1) then
        varcel = .false.
        write (6,*) ' .MD seems to be written without variable cell'
        goto 101
      endif
  102 rewind (ii2)
C --- if it fails, assume MD was written with variable cell...
      read (ii2,err=801,end=107) istep1,(dummy,ii=1,6*nat)
      read (ii2,err=803,end=107) (dummy,ii=1,18)
      read (ii2,err=802,end=107) istep2,(dummy,ii=1,6*nat)
      if (istep2.eq.istep1+1) then
        write (6,*) ' .MD seems to be written with variable cell'
        varcel = .true.
        goto 101
      endif
      write (6,*) ' Fail to identify variable cell in the .MD '
      stop
  101 continue
      nstep = 0
      rewind (ii2)
C --- read through input file, to get the total number of MD steps.
  104 continue
      nstep = nstep + 1
C ... variable format ( formt=F in subr. iomd):
      read (ii2,err=804,end=107) istep,(dummy,ii=1,6*nat)
      if (varcel) read (ii2,err=805,end=108) (dummy,ii=1,18)
      write (6,*) ' MD record No. ',nstep,',  istep=',istep
      goto 104
  108 continue
C --- irregular end of records in MD file: 
      write (6,*) ' Uncomplete record in MD step No.',nstep
      nstep = nstep - 1  !  No. of full records
      if (nstep.gt.0) then
        write (6,*) ' Keep ',nstep,' records.'
        return
      else
        write (6,*) ' Check the MD file; bye'
        stop
      endif
  107 continue
C --- regular end of records in MD or ANI file: 
      nstep = nstep - 1  !  Attempt to read record Nr. nstep failed
      if (nstep.gt.0) then
        write (6,*) '  Cleanly read in ',nstep,'  MD steps'
        return
      else
        write (6,*) ' End of record in the firt step: MD file empty?'
        stop
      endif
  801 continue
      write (6,*) ' Error reading coordinates ',
     .            ' for MD step No. 1,  istep1=',istep1
      stop
  802 continue
      write (6,*) ' Error reading coordinates ',
     .            ' for MD step No. 1,  istep2=',istep2
      stop
  803 continue
      write (6,*) ' Error reading variable cell ',
     .            ' for MD step No. 1'
      stop
  804 continue
      write (6,*) ' Error reading coordinates ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
  805 continue
      write (6,*) ' Error reading variable cell ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
      
      end 
