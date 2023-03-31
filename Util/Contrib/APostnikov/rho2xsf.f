C
C    rho2xsf,  a script to transform 3-dim grid function
C             (i.e. LDOS, RHO, DRHO, etc.) written in SIESTA by subr. iorho
C             into an arbitrarily chosen grid for XCrysden,
C             using linear 4-point interpolation 
C
C   !!! --------------------  IMPORTANT  ---------------------- !!!
C   compile this code with the same compiler switches as Siesta,
C         in what regards using single/double precision,
C       otherwise reading the data from unformatted files 
C              written by iorho.F  might be spoiled.  
C                Don't say you haven't been warned.
C                                !!!
C
C             Written by Andrei Postnikov, Mar 2006   Vers_0.3
C             apostnik@uos.de
C
      program rho2xsf
      implicit none
      integer ii1,ii2,io1,is1
      parameter (ii1=11,ii2=12,io1=14,is1=13)
      integer mesh0(3),mesh1(3),ip,nspin,is,ii,jj,n1,n2,n3,
     .        iat,nat,nz,iix,iiy,iiz,ind,mn,ibox,nbox,
     .        ixmax,iymax,izmax,ixmin,iymin,izmin,ialloc,
     .        ishift,jshift,kshift,n1div,n2div,n3div,i1,i2,i3
      integer limit,maxa,maxo,maxuo,maxnh,maxna,il,ia,nrela(3),idim
      parameter (limit=5)  !  tried translations along each lattice vector
      character inpfil*60,outfil*60,syslab*30,suffix*6,
     .          unitlab*1,labunit*9,labbox*1,owrite*1
      logical unitb,charge,waves,filexist
      double precision b2ang,cc_bohr(3,3),cc_ang(3,3),cc_inv(3,3),
     .                 coort(3),obox(3),rbox(3,3),rinv(3,3),
     .                 cell(3,3),dum,rmaxo,rela,modu,rmesh(3),drela(3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      integer, allocatable :: ityp(:),iz(:)
      double precision, allocatable :: mass(:),coor_ang(:,:)
      character(len=2), allocatable :: label(:)
      real, allocatable :: func(:) !  NB! single precision, as in iorho.F
      real   fintp,fmax,fmin       !  NB! single precision, as in iorho.F
      external test_xv,read_xv,fillbox,inver3,intpl04
C
C     string manipulation functions in Fortran used below:
C     len_trim(string): returns the length of string 
C                       without trailing blank characters,
C     trim(string)    : returns the string with railing blanks removed
C     char(integer)   : returns the character in the specified position
C                       of computer's ASCII table, i.e. char(49)=1
      
      write (6,701,advance="no")
  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
      call test_xv(ii1,nat)
      allocate (ityp(nat))
      allocate (iz(nat))
      allocate (mass(nat))
      allocate (label(nat))
      allocate (coor_ang(1:3,1:nat),STAT=ialloc)
      if (ialloc.ne.0) then
         write (6,*) ' Fails to allocate space for ',nat,' atoms.'
         stop
      endif
      call read_xv(ii1,nat,ityp,iz,cc_ang,mass,label,coor_ang)
      call inver3(cc_ang,cc_inv)
      close (ii1)
C --- set up and fill output box:     
      call makebox(obox,rbox)
C --- invert the box vectors; will need it in a minute...
      call inver3(rbox,rinv)
C --- write selected atoms first into a scratch file (is1), for the case
C     there are zero. Then the label 'ATOMS' with no  atoms following
C     will crash XCrySDen.
C     open (is1, file='tmpfil',form='formatted',status='scratch')
      open (is1,               form='formatted',status='scratch')
      call fillbox(is1,obox,rbox,rinv,cc_ang,nat,coor_ang,nbox)

C --- open output file:
C     outfil = syslab(1:len_trim(syslab))//'.XSF'
      outfil = trim(syslab)//'.XSF'
      inquire (file=outfil, exist=filexist)
      if (filexist) then
        write (6,*) ' File ',trim(outfil),' exists. Overwrite? (Y/N)'
        read (5,*) owrite
        if (owrite.eq.'Y'.or.owrite.eq.'y') then
          open (io1,file=outfil,form='formatted',status='REPLACE')
        else
          write (6,*) ' Then rename is first. Bye...'
          stop
        endif
      else
        open (io1,file=outfil,form='formatted',status='NEW')
      endif
C
      write (6,*) ' The box contains ',nbox,' atoms.'
      if (nbox.gt.0) then
C --- add atoms to XSF; either as periodic cell or as for molecule
C       write (6,702)
C 101   write (6,703,advance="no")
C       read (5,*) labbox
C       if (labbox.eq.'Y'.or.labbox.eq.'y') then        
C         write (io1,'(a7)') 'CRYSTAL'
C         write (io1,'(a7)') 'PRIMVEC'
C         write (io1,201) (rbox(ii,1),ii=1,3)
C         write (io1,201) (rbox(ii,2),ii=1,3)
C         write (io1,201) (rbox(ii,3),ii=1,3)
C         write (io1,'(a9)') 'PRIMCOORD'
C         write (io1,'(i5,i2)') nbox,1
C       else if (labbox.eq.'N'.or.labbox.eq.'n') then
          write (io1,'(a5)') 'ATOMS'
C       else
C         write (6,*) "I do not understand:"
C         goto 101
C       endif
        rewind (is1)
        do ibox = 1,nbox
          read  (is1,202) iat,     (coort(jj),jj=1,3)
          write (io1,202) iz(iat), (coort(jj),jj=1,3)
        enddo
      else
        write (6,*) ' This is OK, just be warned that your XSF file',
     .              ' will have no ATOMS secton!'
      endif 
      close (is1)

      write (6,704) 
  102 write (6,705,advance="no") 
      read (5,*) n1,n2,n3
      if (n1.le.0.or.n2.le.0.or.n3.le.0) then
        write (6,*) ' Numbers must be positive, try again.'
        goto 102
      endif
      idim = 3    !   dimensionalty of output grid
      if (n1.eq.1) idim = idim - 1
      if (n2.eq.1) idim = idim - 1
      if (n3.eq.1) idim = idim - 1

C --- Look for grid data files to include:
  103 write (6,706,advance="no")
      read (5,*) suffix
C     inpfil = syslab(1:len_trim(syslab))//
C    .      '.'//suffix(1:len_trim(suffix))
      inpfil = trim(syslab)//'.'//trim(suffix)
      open (ii2,file=inpfil,form='unformatted',status='old',err=806)
      write (6,*) 'Found and opened: ',trim(inpfil)
      read (ii2,err=807) cell
      read (ii2,err=808) mesh0, nspin
      allocate (func(1:mesh0(1)*mesh0(2)*mesh0(3)),STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for ',
     .                mesh0(1)*mesh0(2)*mesh0(3),' grid points.'
        stop
      endif 
      write (6,*) 'mesh0 = (',mesh0,'),   nspin=',nspin
      write (io1,"('BEGIN_BLOCK_DATAGRID_',i1,'D')") idim
      write (io1,*) 'DATA_from:'//inpfil(1:len_trim(inpfil))
      do is=1,nspin
        fmax= -9.999E+10
        fmin=  9.999E+10
C       write (io1,*) 'BEGIN_DATAGRID_'//char(48+idim)//'D_'//
C    .        suffix(1:len_trim(suffix))//':spin_'//char(48+is)
        write (io1,*) 'BEGIN_DATAGRID_'//char(48+idim)//'D_'//
     .        trim(suffix)//':spin_'//char(48+is)
        if (n1.ne.1) write (io1,'(i6)',advance="no") n1
        if (n2.ne.1) write (io1,'(i6)',advance="no") n2
        if (n3.ne.1) write (io1,'(i6)',advance="no") n3
        write (io1,'()') 
        write (io1,'(1p,3e15.7)') (obox(ii),ii=1,3)
        if (n1.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,1),ii=1,3)
        if (n2.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,2),ii=1,3)
        if (n3.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,3),ii=1,3)
        
        ind=0                        
        do iiz=1,mesh0(3)
          do iiy=1,mesh0(2)
            read (ii2,err=809,end=810) (func(ind+iix),iix=1,mesh0(1))
            ind = ind + mesh0(1)
          enddo
        enddo

C ---   loop over mesh points
C       avoid division by zero if only 1 point is selected: 
        n1div=max(n1-1,1)  !  if (n1.eq.1) n1div=1 else n1div=n1-1
        n2div=max(n2-1,1)
        n3div=max(n3-1,1)
        ind = 0
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              ind = ind+1
              do ii=1,3
                rmesh(ii) = obox(ii) + 
     +                      rbox(ii,1)*(i1-1)/n1div +
     +                      rbox(ii,2)*(i2-1)/n2div +
     +                      rbox(ii,3)*(i3-1)/n3div
              enddo
C ---         rmesh(1..3) are absolute Cartesian coordinates
C             of the mesh point (i1,i2,i3) in Angstroem
C             Find its relative coordinates on the unit cell grid:
              do ii=1,3
                rela = 0.0
                do jj=1,3
                  rela = rela + cc_inv(ii,jj)*rmesh(jj)
                enddo
C               take modulo to make sure that it falls within [0,1]:
                modu = modulo( rela*mesh0(ii), dble(mesh0(ii)) )
                nrela(ii) = floor(modu) + 1
                drela(ii) = modu - nrela(ii) + 1
              enddo
C ---         mesh point rmesh(1..3) falls within the grid microcell
C             originating at the grid point nrela(1..3),
C             its relative coordinates within this microcell are drela(1..3)
C             Select neighboring grid points and make the interpolation:
              call intpl04 (func(1),fintp,
     .                      mesh0(1),mesh0(2),mesh0(3),nrela,drela)
              if (fintp.gt.fmax) then
                fmax = fintp
                ixmax = i1
                iymax = i2
                izmax = i3
              endif
              if (fintp.lt.fmin) then
                fmin = fintp
                ixmin = i1
                iymin = i2
                izmin = i3
              endif
              write (io1,'(1p,e13.5)',advance="no")
     $                          fintp        ! write without linebreak
              if (mod(ind,6).eq.0) write (io1,'()') ! linebreak after 6 numbers
            enddo    !  do i1 =
          enddo    !  do i2 =
        enddo    !  do i3 =
        write (io1,'()')  ! linebreak after all numbers
        write (io1,"(' END_DATAGRID_',i1,'D')") idim
        write (6,206) is,'max',fmax,ixmax,iymax,izmax
        write (6,206) is,'min',fmin,ixmin,iymin,izmin
      enddo  ! do is=1,nspin
      deallocate (func)
      write (io1,"('END_BLOCK_DATAGRID_',i1,'D')") idim
      close (ii2)
      goto 103                                                                  

  201 format (3f20.8)
  202 format (i4,3f20.8)
  203 format (3i6)
  204 format (3f12.7)
  205 format (1p,6e13.6)
C 205 format (1p,8e10.3)
  206 format (' For is=',i1,': ',a3,'. grid value =',1p,e12.5,
     .        ' at iix,iiy,iiz=',3i4)

  702 format(" The atom section may appear in the XSF file",
     .       " as for periodic structure, that will allow you",
     .       " to replicate the units within the XCrySDen.",
     .       " This only makes sense if your selected box",
     .       " coincides with the true periodic cell. Otherwise",
     .       " the atoms will be written non-periodically, as for",
     .       " molecule. ")
  703 format(" Do you want atom section as for periodic ",
     .       " structure, Y or N ? ")
  704 format (" Now define the grid. If you want it two-dimensional,",/
     .   " give 1 as number of grid points along one spanning vector.") 
  705 format (" Enter number of grid points along three vectors: ")
  706 format (' Add grid property (LDOS, RHO, ...;',
     .        ' or BYE if none): ')
  707 format ("A pity, indeed. In 99% of cases this happens because of",
     .   " incompatibility",/"of binary file format. Make sure",
     .   " that you run rho2xsf on THE SAME platform",/"where you have",
     .   " run Siesta and use the same compiler settings.")

  801 write (6,*) ' Error opening file ',
     .              trim(inpfil),' as old formatted'
      stop
  806 write (6,*) ' A wild guess! There is no file ',
     .              trim(inpfil),'; close XSF and quit.'
      close (io1)
      stop
  807 write (6,*) ' Error reading cell vectors from binary file !'
      write (6,707)
      stop
  808 write (6,*) ' Error reading n1,n2,n3,nspin from binary file !'
      write (6,707)
      stop
  809 write (6,*) ' Error reading function values on the grid'
      write (6,707)
      stop
  810 write (6,*) ' Unexpected end of data for iiy=',iiy,
     .            ' iiz=',iiz,'  is=',is
      stop
  811 write (6,*) ' Error opening file ',
     .              trim(outfil),' as new unformatted'
      stop
      end
