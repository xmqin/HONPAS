! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      program siesta2arc
      implicit real*8(a-h,o-z)
      parameter (maxat=1000)
      parameter (maxsp=10)
      parameter (maxword=100)
      include 'constants'
C
C  Generates an archive file for visualising SIESTA output
C
C  Excute program as :
C
C  siesta2arc < rootname.out > rootname.arc
C
C  species symbols have to be truncated to 4 characters due
C  to fixed format of arc files
C
C  Julian Gale, Imperial College, March 1999
C
      dimension rvl(3,3)
      dimension x(maxat),y(maxat),z(maxat),n(maxat)
      dimension nat(maxsp),label(maxsp)
      dimension words(maxword),floats(maxword),nlorder(2*maxword)
      character*1 blank
      character*2 asym
      character*4 wtype1,lab2
      character*5 lab,label
      character*30 words
      character*80 line
      logical eof,fractional,first
      iout = 6
      eof=.false.
      etot=0.0d0
      first=.true.
      blank=' '
      wtype1='CORE'
      q=0.0d0
      ndimen=0
      iline=0
      fractional=.false.
      scale = 1.0d0
C***************************
C  Initialisation of file  *
C***************************
C
C  Write out header
C
      write(iout,'(''!BIOSYM archive 2'')')
C*********************************************
C  Loop over file looking for relevant data  *
C*********************************************
      do while (.not.eof)
        do i = 1,maxword
          words(i) = ' '
          floats(i) = 0.0d0
        enddo
C#############################
C  Read line and process it  #
C#############################
        iline = iline + 1
        read(5,'(a)',err=10,end=10) line
        call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *    maxword)
C######################
C  Coordinate format  #
C######################
        if (index(words(1),'AtomicCoord').ne.0) then
          if (index(words(2),'Frac').ne.0) then
            fractional = .true.
            ndimen = 3
          endif
          if (index(words(2),'Cart').ne.0) fractional = .false.
          if (index(words(2),'NotSc').ne.0) then
            fractional = .false.
            scale = 1.0d0
          endif
          if (index(words(2),'Bohr').ne.0) then
            scale = autoangs*scale
          endif
        endif
C#########################
C  Lattice Scale Factor  #
C#########################
        if (index(words(1),'LatticeCon').ne.0) then
          scale = floats(1)
          if (index(words(2),'Bohr').ne.0) scale=scale*autoangs
        endif
C#######################
C  Lattice Parameters  #
C#######################
        if (index(words(2),'LatticePar').ne.0) then
          ndimen = 3
          read(5,'(a)',err=10,end=10) line
          call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *      maxword)
          a = floats(1)*scale
          b = floats(2)*scale
          c = floats(3)*scale
          alpha = floats(4)*scale
          beta  = floats(5)*scale
          gamma = floats(6)*scale
C
C  Generate cell in standard orientation for coordinate output
C
          call cell(rvl,a,b,c,alpha,beta,gamma)
          ndimen = 3
          read(5,'(a)',err=10,end=10) line
        endif
C####################
C  Lattice Vectors  #
C####################
        if (index(words(2),'LatticeVec').ne.0) then
          ndimen = 3
          read(5,'(a)',err=10,end=10) line
          call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *      maxword)
          rvl(1,1)=floats(1)*scale
          rvl(2,1)=floats(2)*scale
          rvl(3,1)=floats(3)*scale
          read(5,'(a)',err=10,end=10) line
          call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *      maxword)
          rvl(1,2)=floats(1)*scale
          rvl(2,2)=floats(2)*scale
          rvl(3,2)=floats(3)*scale
          read(5,'(a)',err=10,end=10) line
          call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *      maxword)
          rvl(1,3)=floats(1)*scale
          rvl(2,3)=floats(2)*scale
          rvl(3,3)=floats(3)*scale
          print *,' Rv : ',rvl(1,1),rvl(2,1),rvl(3,1)
          print *,' Rv : ',rvl(1,2),rvl(2,2),rvl(3,2)
          print *,' Rv : ',rvl(1,3),rvl(2,3),rvl(3,3)
          call uncell(rvl,a,b,c,alpha,beta,gamma)
          read(5,'(a)',err=10,end=10) line
        endif
C####################
C  Number of atoms  #
C####################
        if (index(words(1),'NumberOfAt').ne.0) then
          numat = nint(floats(1))
          if (numat.gt.maxat) then
            write(6,'(''ERROR : Too many atoms'')')
            stop
          endif
        endif
C######################
C  Number of species  #
C######################
        if (index(words(1),'NumberOfSp').ne.0) then
          nsp = nint(floats(1))
          if (nsp.gt.maxsp) then
            write(6,'(''ERROR : Too many species'')')
            stop
          endif
        endif
C###################
C  Species labels  #
C###################
        if (index(words(2),'ChemicalSp').ne.0) then
          do i = 1,nsp
            read(5,'(a)',err=10,end=10) line
            call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *        maxword)
            nat(i) = nint(floats(2))
            label(i) = words(1)(1:5)
          enddo
          read(5,'(a)',err=10,end=10) line
        endif
C***********************
C  Configuration dump  *
C***********************
        if (index(words(1),'outcoor:').ne.0) then
C
C  Read coordinates
C
          do i = 1,numat
            read(5,'(a)',err=10,end=10) line
            call linepro(line,nword,words,nfloat,floats,nlorder,iline,
     *        maxword)
            x(i) = floats(1)
            y(i) = floats(2)
            z(i) = floats(3)
            n(i) = nint(floats(4))
          enddo
          if (first) then
            if (ndimen.eq.3) then
              write(iout,'(''PBC=ON'')')
            else
              write(iout,'(''PBC=OFF'')')
            endif
            first=.false.
          endif
C
C  Write out first line
C
          write(iout,'(64a1,f16.6)')(blank,i=1,64),etot*evtokcal
          write(iout,'(''!DATE'')')
          if (ndimen.eq.3) then
C****************
C  Bulk output  *
C****************
            write(iout,'(''PBC'',6f10.4)')a,b,c,alpha,beta,gamma
            do i=1,numat
              if (fractional) then
                xci=x(i)*rvl(1,1)+y(i)*rvl(1,2)+z(i)*rvl(1,3)
                yci=x(i)*rvl(2,1)+y(i)*rvl(2,2)+z(i)*rvl(2,3)
                zci=x(i)*rvl(3,1)+y(i)*rvl(3,2)+z(i)*rvl(3,3)
              else
                xci=x(i)*scale
                yci=y(i)*scale
                zci=z(i)*scale
              endif
              lab2=label(n(i))(1:4)
              asym=lab2(1:2)
              write(iout,
     *          '(a4,1x,3f15.9,1x,a4,1x,i4,2(1x,a2),1x,f8.4,1x,i4)')
     *          lab2,xci,yci,zci,wtype1,i,asym,asym,q,i
            enddo
            write(iout,'(''end'')')
          else
C*******************
C  Cluster output  *
C*******************
            do i=1,numat
              lab2=label(n(i))(1:4)
              asym=lab2(1:2)
              write(iout,
     *          '(a4,1x,3f15.9,1x,a4,1x,i4,2(1x,a2),1x,f8.4,1x,i4)')
     *          lab2,x(i),y(i),z(i),wtype1,i,asym,asym,q,i
            enddo
            write(iout,'(''end'')')
          endif
          write(iout,'(''end'')')
        endif
      enddo
   10 continue
      end
