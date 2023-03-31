C
C    md2axsf,  a script to transform molecular dynamics files 
C             (MD or ANI) written in SIESTA by subr. iomd
C             into an XCrysden animation file (axsf format).
C    Three options are possible:
C    1. Use .MD or .MD_CAR written with MD.VariableCell = T,
C       then the output is variable-cell animated xsf
C    2. Use .MD or .MD_CAR written with MD.VariableCell = F,
C       then the output is fixed-cell (as in XV file) animated xsf
C    3. Use .ANI which contains no information about the cell;
C       then the output is fixed-cell (as in XV file) animated xsf
C
C       Written by Andrei Postnikov, Jul 2007   Vers_0.4
C       apostnik@uos.de
C
      program md2axsf
      implicit none
      integer ii1,ii2,io1,is1,iat,nat,na,ialloc,ii,jj,mdmod,mdstep,
     .        istep,nstep,idum,mdfirst,mdlast
      parameter (ii1=11,ii2=12,io1=14,is1=13)

      character inpfil*60,outfil*60,syslab*30,suffix*6,
     .          unitlab*1,labunit*9,labbox*1
      logical varcel
      double precision b2ang,cc_bohr(3,3),cc_ang(3,3),cc_velo(3,3),
     .                 obox(3),rbox(3,3),rinv(3,3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      double precision,  allocatable :: coord(:,:),veloc(:,:)
      integer,           allocatable :: nz(:), ityp(:)
      external test_md,test_mdc,test_ani,makebox,inver3,
     .         write_axsf1,write_axsf2,opnout
C
C     string manipulation functions in Fortran used below:
C     len_trimd(string): returns the length of string 
C                       without trailing blank characters,
C     trim(string)    : returns the string with railing blanks removed
      
      write (6,701,advance="no")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
C --  read in translation vectors, convert into Ang:
      do ii=1,3
        read  (ii1,702,end=803,err=803)  (cc_bohr(jj,ii),jj=1,3)
      enddo
      cc_ang = cc_bohr*b2ang
      read  (ii1,*,end=804,err=804)  nat
      allocate (coord(3,nat),STAT=ialloc)
      allocate (veloc(3,nat),STAT=ialloc)
      allocate (nz(nat),STAT=ialloc)
      allocate (ityp(nat),STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for ',nat,' atoms.'
        stop
      endif
      do iat=1,nat
        read (ii1,704,end=805,err=805) ityp(iat), nz(iat), 
     .                               (coord(ii,iat),ii=1,3)
      enddo
      close (ii1)
      write (6,*) ' with ',nat,'  atoms '
C --- finished with .XV
  103 write (6,705,advance="no")
      read (5,*) suffix
      inpfil = trim(syslab)//'.'//trim(suffix)
      if (trim(suffix).eq.'MD') then
        mdmod = 1
        open (ii2,file=inpfil,form='unformatted',status='old',err=806)
        write (6,*) 'Found and opened: ',inpfil
        call test_md(ii2,nat,varcel,nstep)
      elseif (trim(suffix).eq.'MD_CAR') then
        mdmod = 2
        open (ii2,file=inpfil,form='formatted',status='old',err=801)
        write (6,*) 'Found and opened: ',inpfil
        call test_mdc(ii2,nat,varcel,nstep)
      elseif (trim(suffix).eq.'ANI') then
        mdmod = 3
        open (ii2,file=inpfil,form='formatted',status='old',err=801)
        write (6,*) 'Found and opened: ',inpfil
        call test_ani(ii2,nat,nstep)
      else
        write (6,*) ' Only MD, MD_CAR or ANI allowed, try again:'
        goto 103
      endif
      write (6,710,advance="no")
      mdfirst = 0
      mdlast  = 0
      mdstep  = 0
      read (5,*) mdfirst, mdlast, mdstep     
      if (mdfirst.le.0) mdfirst = 1
      if (mdlast.le.0.or.mdlast.gt.nstep)  mdlast  = nstep
      if (mdstep.le.0)  mdstep  = 1
      if (mdlast.lt.mdfirst) mdlast = mdfirst
      write (6,711) mdfirst,mdfirst+mdstep,
     .              mdfirst+((mdlast-mdfirst)/mdstep)*mdstep,
     .              (mdlast-mdfirst)/mdstep+1
C --- open output file:
C     outfil = syslab(1:len_trim(syslab))//'.AXSF'
      outfil = trim(syslab)//'.AXSF'
      call opnout(io1,outfil)
C --- provide an option for defining output box. 
C     Without the box, AXSF generated from MD will contain unit cell vectors
C     (variable or not) and atoms in he nit cell; AXSF generated from ANI
C     will contain just the list of atoms (as for molecule) without cell.
C     With he box, the AXSF file will contain only the list of atoms
C     (no periodic cell), extended or reduced to the actual size
C     of the output box.
  104 write (6,712,advance="no") 
      read (5,*) labbox
      if (labbox.eq.'Y'.or.labbox.eq.'y') then
        call makebox(obox,rbox)
        call inver3(rbox,rinv)
        call write_axsf1(ii2,io1,nstep,mdfirst,mdlast,mdstep,
     .                   nat,nz,mdmod,varcel,cc_ang,coord,veloc,
     .                   obox,rbox,rinv)
      elseif (labbox.eq.'N'.or.labbox.eq.'n') then
        call write_axsf2(ii2,io1,nstep,mdfirst,mdlast,mdstep,
     .                   nat,nz,mdmod,varcel,cc_ang,coord,veloc)
      else
        goto 104
      endif

      close (io1) 
      write (6,*) 'Written to and closed ',outfil
      deallocate (coord,veloc,nz,ityp)

      stop

  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
  702 format(3x,3f18.9)
  703 format(i13)
  704 format(i3,i6,3f18.9)
  705 format(' Suffix of molecular dynamics file',
     .       ' (MD, MD_CAR, or ANI): ')
  706 format(' Was this MD run done with MD.VariableCell, T or F : ')
  707 format(' MD record No. ',i6,'  istep =',i6)
  708 format(i5,/)
  709 format(a2,2x,3f12.6)    
  710 format(' You may wish to keep only some of these steps,',/,
     .       ' MDfirst, MDfirst+MDstep, MDfirst+2*MDstep etc.', 
     .       ' till (not exceeding) MDlast',/,' from the list above.',
     .       ' Specify three numbers MDfirst, MDlast, MDstep -',/,
     .       ' or 0 for any of them as default : ')
  711 format(" OK, I'll keep steps Nr.",i6,','i6,', etc. till',i6,
     .       ' ( total',i6,' )')   
  712 format(" Do you want to define output box (Y or N): ")

  801 write (6,*) ' Error opening file ',
     .            trim(inpfil),' as old formatted'
      stop
  803 write (6,*) ' End/Error reading XV for cell vector ',ii
      stop
  804 write (6,*) ' End/Error reading XV for number of atoms line'
      stop
  805 write (6,*) ' End/Error reading XV for atom number ',iat
      stop
  806 write (6,*) ' Error opening file ',
     .            trim(outfil),' as old unformatted'
      stop
  807 write (6,*) ' Error reading file ',
     .            trim(outfil),' for nstep=',nstep
      write (6,*) ' This can be (but not necessarily) due to',
     .            ' your wrong guess about VariableCell.'
      stop
  808 write (6,*) ' Error reading file ',
     .            trim(outfil),' for nstep=',nstep
      stop

      end
