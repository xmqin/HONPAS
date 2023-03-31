C
C   ene2bxsf, a script to make the band-XSF file (for plotting
C             the Fermi surface) from siesta KP and EIG files.
C             The k-mesh must be generated without shift
C             (requirement of the BXSF format).
C
C             Written by Andrei Postnikov, Jun 2007, modif Jul 2010   Vers_0.4
C             postnikov@univ-metz.fr
C
      program ene2bxsf
      implicit none
      integer ii1,ii2,ii3,io1,il1
      parameter (ii1=11,ii2=12,ii3=13,io1=14,il1=15)
      integer ii,jj,nbands,nbmin,nbmax,nkp,ikp,ndum,iis,mdiv(3),ind,
     .        nspin,is,iband,iik,idiv1,idiv2,idiv3,itry1,itry2,itry3,
     .        ib,mkp,homo(2),lumo(2),ik,in2,in3,ialloc
      character inpfil*60,outfil*60,syslab*30,suffix*6
      double precision cell(3,3),efermi,twopi,rcell(3,3),small,
     .       relmin(3),kkp(3),dum
      logical unshift
      double precision,  allocatable :: relk(:,:), eneb(:,:), enek(:)
      integer,           allocatable :: ndiv(:,:),irrek(:)
      parameter (small=1.0d-04)
      external inver3,opnout
C
C     string manipulation functions in Fortran used below:
C     len_trim(string): returns the length of string 
C                       without trailing blank characters,
C     trim(string)    : returns the string with railing blanks removed
C     char(integer)   : returns the character in the specified position
C                       of computer's ASCII table, i.e. char(49)=1
C
      twopi = 8.d0*atan(1.d0)
      write (6,701,advance="no")
  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
      read (5,*) syslab
C --- open .XV and .EIG : 
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
C     inpfil = syslab(1:len_trim(syslab))//'.EIG'
      inpfil = trim(syslab)//'.EIG'
      open (ii2,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
C --- open LOG file:
      open (il1,file='eig2bxsf.LOG',form='formatted',status='unknown')
C --- read Fermi energy and total number of bands from .EIG :
      read  (ii2,*) efermi
      read  (ii2,*) nbands, nspin, nkp
      if (nspin.ne.1.and.nspin.ne.2) then
        write (6,*) 'A problem encountered: nspin=',nspin
        stop
      endif

C --- finds bands crossing the Fermi energy:
      do is=1,nspin
        homo(is)=0           !  higest (partially?) occupied band
        lumo(is)=nbands+1    !  lowest (partially?) unoccupied band
      enddo
      write (6,*) '  nkp=',nkp,'  nbands=',nbands
      allocate (eneb(nbands,nspin),STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for nbands=',
     .                nbands,', nspin=',nspin  
        stop
      endif
      do ik=1,nkp
        read (ii2,"(i5,10f12.5,/,(5x,10f12.5))") iik,   !  written in
     .  ((eneb(ib,is),ib=1,nbands),is=1,min(nspin,2))   !  ioeif.f
        if (iik.ne.ik) then
          write (6,*) ' iik=',iik,'.ne. ik=',ik,' for spin ',is
          stop
        endif
        do is=1,nspin
          do ib=1,nbands
            if (eneb(ib,is).lt.efermi.and.ib.gt.homo(is)) homo(is) = ib
            if (eneb(ib,is).gt.efermi.and.ib.lt.lumo(is)) lumo(is) = ib
          enddo
        enddo
      enddo
      deallocate (eneb)
      do is=1,nspin
        write (6,*) ' is=',is,'  homo, lumo=',homo(is),lumo(is)
        if (homo(is).lt.lumo(is)) write (6,201) is, homo(is),lumo(is)
      enddo

C --- read cell vectors from .XV, find reciprocal:
      do ii=1,3
        read  (ii1,*,end=803,err=803)   (cell(ii,jj),jj=1,3)
      enddo
      close (ii1)
      call inver3(cell,rcell)
C +++++++++++++++++++++++++
      write(il1,*) 'Direct cell vectors:'
      do ii=1,3
        write(il1, 502) (cell(ii,jj),jj=1,3)
  502   format(3f14.8)
      enddo
      write(il1,*) 'Reciprocal cell vectors without 2pi:'
      do ii=1,3
        write(il1, 502) (rcell(jj,ii),jj=1,3)
      enddo
      write(il1,*) 'Reciprocal cell vectors times 2pi:'
      do ii=1,3
        write(il1, 502) (rcell(jj,ii)*twopi,jj=1,3)
      enddo
C +++++++++++++++++++++++++

C --- open .KP as old:
C     inpfil = syslab(1:len_trim(syslab))//'.KP'
      inpfil = trim(syslab)//'.KP'
      open (ii3,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
C --- read k-points from the .KP file and recover their fractional
C     coordinates in terms of reciprocal vectors:
      read (ii3,*) nkp
      relmin(:)=1.d0/small
      allocate (relk(3,nkp))
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for relk(3,',nkp,')'
        stop
      endif
      do ikp=1,nkp
        read (ii3,*) iik,(kkp(jj),jj=1,3)
        if (iik.ne.ikp) then
          write (6,*) ' a mess in KP list: read in iik=',iik,
     .                ', but expected ikp=',ikp
          stop
        endif
C ---   find relative coordinates of k-point:
        do ii=1,3
          relk(ii,ikp)=( cell(ii,1)*kkp(1) +
     +                   cell(ii,2)*kkp(2) +
     +                   cell(ii,3)*kkp(3) )/twopi 
          if (abs(relk(ii,ikp)).lt.relmin(ii).and.
     .        abs(relk(ii,ikp)).gt.small) 
     .        relmin(ii)=abs(relk(ii,ikp))
        enddo
      enddo
      close (ii3)
C --- relmin(ii) is now the smallest relative coordinate of k points 
C     along the reciprocal vector (ii). Good chance that its inverse 
C     is number of divisions. Round up to the closest integer:
      mdiv(:)=floor(1.d0/relmin(:)+ 0.5)
      write (6,*) ' No. of divisions seems to be ',mdiv
      mkp = (mdiv(1)+1)*(mdiv(2)+1)*(mdiv(3)+1)
      write (6,*) ' Full No. of k-points on general grid:', mkp
      allocate (ndiv(3,mkp),STAT=ialloc)
      allocate (irrek(mkp), STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for mkp=',mkp
        stop
      endif
C --- decifer all k-point coordinates as integer fractions, ndiv=relk/mdiv, 
C     with  0 <= relk  < mdiv :
      unshift = .false.
      do ikp=1,nkp
        ndiv(:,ikp)=modulo(floor(relk(:,ikp)/relmin(:)+0.5),mdiv(:)) 
C          floor(relk/relmin)      : max. integer not exceeding relk/relmin
C          floor(relk/relmin)+0.5  : the integer closest to relk/relmin
C          modulo(..)              : an always non-negative coord. relk/mdiv
C++++++++++++++++++
        write (il1,501) ikp, (relk(ii,ikp),ii=1,3),
     .                     (ndiv(ii,ikp),mdiv(ii),ii=1,3)
  501   format (i4,'  relk=',3f10.6,'   ndiv=',3(i4,'/',i2))
C++++++++++++++++++
        if (.not.unshift.and.
     .       (ndiv(1,ikp)**2+ndiv(2,ikp)**2+ndiv(3,ikp)**2).eq.0) then
          unshift = .true.
          write (6,*) ' k-mesh contains Gamma-point; this is good...'
        endif
      enddo
      if (.not.unshift) then
        write (6,207) 
        stop
      endif
C --- attribute irreducible k-points 
C     to k-points over the full reciprocal cell:
      ind=0
C -  The right sequnce ("column-major") of writing energies into .BXSF file:
      do idiv1 = 0,mdiv(1)
      do idiv2 = 0,mdiv(2)
      div3: do idiv3 = 0,mdiv(3) 
        ind = ind + 1      !  global address on the k-mesh
        itry1=mod(idiv1,mdiv(1))
        itry2=mod(idiv2,mdiv(2))
        itry3=mod(idiv3,mdiv(3))
        do ikp=1,nkp
          if( (ndiv(1,ikp).eq.itry1).and.
     .        (ndiv(2,ikp).eq.itry2).and.
     .        (ndiv(3,ikp).eq.itry3) ) then
! ---       k-point is identified!
            irrek(ind) = ikp
            write (il1,503) idiv1,idiv2,idiv3,ikp   
            cycle div3 ! - skip the rest of ikp loop,
                       !   try the next idiv3. This is the regular termination
                       !   of the ikp search if the k-point is found.
          endif
        enddo 
! ---   haven't find any matching k-point; take a second chance, try inversion:
        itry1=mod(mdiv(1)-idiv1,mdiv(1))
        itry2=mod(mdiv(2)-idiv2,mdiv(2))
        itry3=mod(mdiv(3)-idiv3,mdiv(3))
        do ikp=1,nkp
          if( (ndiv(1,ikp).eq.itry1).and.
     .        (ndiv(2,ikp).eq.itry2).and.
     .        (ndiv(3,ikp).eq.itry3) ) then
C ---        k-point is identified!
            irrek(ind) = ikp
            write (il1,503) idiv1,idiv2,idiv3,ikp   
            cycle div3 ! - skip the rest of ikp loop,
                       !   try the next idiv3. This is the regular termination
                       !   of the ikp search if the k-point is found.
          endif
        enddo
! ---   Loops over k-points with and without inversion terminated
!       yet identification failed. Stop and complain:
        write (6,504) idiv1,mdiv(1),idiv2,mdiv(2),idiv3,mdiv(3)
        stop
      enddo div3
      enddo
      enddo
      allocate (enek(nkp),  STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for enek(1:',nkp,')'
        stop
      endif

C --- The writing sequence (into two BXSF files if two spins)
      do is=1,nspin
        if (nspin.eq.1) then
C         outfil = syslab(1:len_trim(syslab))//'.BXSF'
          outfil = trim(syslab)//'.BXSF'
        else
C         outfil = 
C    .    syslab(1:len_trim(syslab))//'_spin_'//char(48+is)//'.BXSF'
          outfil = trim(syslab)//'_spin_'//char(48+is)//'.BXSF'
        endif
        call opnout(io1,outfil)
        write (io1, "(a10)") 'BEGIN_INFO'
        write (io1, "(a5)")  '  #  '
        write (io1, "(a33)") '  #  Band-XCRYSDEN-Structure-File'
        write (io1, "(a5)")  '  #  '
        write (io1, "(a16,f10.4)")  '  Fermi energy: ',efermi
        write (io1, "(a8,/)")  'END_INFO'
        write (io1, "(a23)") 'BEGIN_BLOCK_BANDGRID_3D'
C       write (io1, *)  ' ',syslab(1:len_trim(syslab))
        write (io1, *)  ' ',trim(syslab)
        write (io1, "(a19,a1)") '  BANDGRID_3D_spin_',char(48+is)
        write (io1, "(i5)")     homo(is)-lumo(is)+1  !  No. of bands 
        write (io1, "(3i5)")    mdiv(1)+1, mdiv(2)+1, mdiv(3)+1
        write (io1, "(3f16.8)") 0.0, 0.0, 0.0
        write (io1, "(3f16.8)") (rcell(jj,1)*twopi,jj=1,3)
        write (io1, "(3f16.8)") (rcell(jj,2)*twopi,jj=1,3)
        write (io1, "(3f16.8)") (rcell(jj,3)*twopi,jj=1,3)
        do iband=lumo(is),homo(is)
          write (io1, '(a7,i5)') '  BAND:',iband
C --- again read band energies of the needed band and spin from .ENE 
C     and write them in the correct order into .BXSF
          rewind (ii2)
          read  (ii2,*) dum
          read  (ii2,*) ndum, ndum, ndum
C --- read in all energy values over all k points for the given spin band:
          do ik=1,nkp
            read (ii2,"(i5,10f12.5,/,(5x,10f12.5))") iik, 
     .      ((dum,ib=1,nbands),iis=1,is-1),    ! dummy read prev. spin, if any
     .       (dum,ib=1,iband-1),enek(ik),(dum,ib=iband+1,nbands),
     .      ((dum,ib=1,nbands),iis=is+1,nspin) ! dummy read next spin, if any
          enddo
C --- write into .BXSF file:
            write (io1, "(7f11.5)") 
     .            (enek(irrek(ii)),ii=1,mkp)
        enddo   !  do iband=lumo(is),homo(is)
        write (io1, '(a17)') '  END_BANDGRID_3D'
        write (io1, '(a21)') 'END_BLOCK_BANDGRID_3D'
        close (io1)
      enddo   !  do is=1,ispin
      deallocate (enek,relk,ndiv,irrek)
      close (ii2)
C     close (il1)
      stop

  201 format (' spin',i2,': band gap between bands ',i5,'  and ',i5)
  204 format (3f12.7)
  205 format (1p,6e13.6)
C 205 format (1p,8e10.3)
  206 format (' For is=',i1,': ',a3,'. grid value =',1p,e12.5,
     .        ' at ix,iy,iz=',3i4)
  207 format (' Gamma point not found in the k-mesh.', /
     .        ' Fermi surface calculation needs energy bands',
     .        ' calculated with unshifted mesh!')

  503 format('  idiv1,2,3 = (',3i5,' )   is  ikp Nr.',i6)
  504 format(' No irreducible k-point found for k-grid point (',
     .       i4,'/',i4,',   ',i4,'/',i4,',   ',i4,'/',i4,' )')

  801 write (6,*) ' Error opening file ',
C    .            inpfil(1:len_trim(inpfil)),' as old formatted'
     .            trim(inpfil),' as old formatted'
      stop
  803 write (6,*) ' End/Error reading XV for cell vector ',ii
      stop

      end
