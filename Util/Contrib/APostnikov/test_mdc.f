C...............................................................
C
      subroutine test_mdc(ii2,nat,varcel,nstep)        
C
C     reads from ii2 (MD_CAR file), makes some consistency tests,
C     returns the number of MD steps.
C
C     MD_CAR file is written in md_out.F90 as follows:
C  write(iomd,"(a)") "---" // trim(slabel) //"---"
C  write(iomd,"(f10.1)") 1.0
C  do i=1, 3
C     write(iomd,"(3f16.9)") cell(:,i)/Ang
C  enddo
C  write(iomd,"(30i6)") natoms(:)
C  write(iomd,"(a)") "Direct"
C  call reclat(cell,celli,0)
C    do i=1, na
C      frac(:) =  matmul(transpose(celli),xa(:,i))
C      write(iomd,"(3f16.9)") frac(:)
C    enddo
C  call pxfflush(iomd)
C
      implicit none
      integer ii2,nat,nstep,istep,istep1,istep2,ii,jj
      logical varcel
      double precision dummy
      character lett*1

      nstep = 0
      rewind (ii2)
C --- read through input file, to get the total number of MD steps.
  104 continue
      nstep = nstep + 1
      read (ii2,'(a1)',err=801,end=107) lett !  "---" // trim(slabel) //"---"
      read (ii2,'(a1)',err=802,end=108) lett !  "       1.0"
      read (ii2,'(3f16.9)',err=803,end=108) (dummy,ii=1,9)
      read (ii2,'(a1)',err=804,end=108) lett !  List of atoms: modify!!!
      read (ii2,'(a1)',err=804,end=108) lett !  "Direct"
      do ii = 1,nat
        read (ii2,'(3f16.9)',err=805,end=108) (dummy,jj=1,3)
      enddo
      write (6,*) ' MD record No. ',nstep,',  istep=',istep
      goto 104
  108 continue
C --- irregular end of records in MD_CAR file: 
      write (6,*) ' Uncomplete record in MD_CAR step No.',nstep
      nstep = nstep - 1  !  No. of full records
      if (nstep.gt.0) then
        write (6,*) ' Keep ',nstep,' records.'
        return
      else
        write (6,*) ' Check the MD_CAR file; bye'
        stop
      endif
  107 continue
C --- regular end of records in MD_CAR file: 
      nstep = nstep - 1  !  Attempt to read record Nr. nstep failed
      if (nstep.gt.0) then
        write (6,*) '  Cleanly read in ',nstep,'  MD steps'
        return
      else
        write (6,*) ' End of record in the 1st step: MD_CAR file empty?'
        stop
      endif
  801 continue
      write (6,*) ' Error reading header ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
  802 continue
      write (6,*) ' Error reading 1.0 (the 2d) line ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
  803 continue
      write (6,*) ' Error reading lattice vectors ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
  804 continue
      write (6,*) ' Error reading coordinates ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
  805 continue
      write (6,*) ' Error reading lattice vectors ',
     .            ' for MD step No. ',nstep,'  istep=',istep
      stop
      
      end 
