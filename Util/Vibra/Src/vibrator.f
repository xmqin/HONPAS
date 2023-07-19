! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!

      program vibrator

c *********************************************************************
c Calculation of vibrational modes for clusters, linear chains, slabs 
c and 3D xtals.
c 
c Uses the FDF (Flexible Data Format) package (version 0.66.1.5)
c of J.M.Soler and A. Garcia, 
c
c Written by P.Ordejon, August'98
c
c Modified to add Infra-red intensity calculation by J.D. Gale, July'06
c
c **********************************************************************

      use fdf

      implicit none

      include 'vibra.h'

c Internal variables ...

      logical overflow, eigen, intensity

      character
     .  filein*20, fileout*20, fname*33

      character
     .  slabel*20, sname*150, slabel_defect*150, sname_defect*20, 
     .  paste*33

      integer 
     .  i, i1, i2, iatom, imass(maxa), imasssc(maxasc), iunit, 
     .  ik, iunit2, iunit3, ind, ix, j,
     .  ij, ii, iq, jj, jx, icall, 
     .  lx, ly, lz, lxmax, lymax, lzmax, 
     .  llx, lly, llz,
     .  lx_defect, ly_defect, lz_defect,
     .  na_defect, natoms, ncells, nnat

      integer maxlin, maxk, maxnq
      parameter (maxlin =  1000)
      parameter (maxk   =  5000)
      parameter (maxnq  =  1000)
     
      integer 
     .  in, neq, nk, nlines, lastk(maxlin)

      character label(maxlin)*8

      real*8 kpoint(3,maxk), ek(maxd,maxk)

      real*8 
     .  dx, alat, alp, blp, clp, alplp, betlp, gamlp, pi, xxx

      real*8
     .  zpe, planck

      real*8 phi0(3,maxa,3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz),
     .        phi(3,maxa,3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz)
      real*8 pp(3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz),
     .       pn(3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz)
      real*8 IRinten(3*maxa), BornQ(3,3,maxa)

      real*8
     .  b(3,maxa), cell(3,3), r(3), scell(3,3), xa(3,maxasc),xmass(maxa)

      real*8
     .  correct, dmin, q(3), qr(maxnq), r2, rl(3), rmass

c Correction terms to satisfy translational modes
      real*8 zero(3,3,maxa), zeroo(3,3)

c Work space for diagonalization.
      complex dc(maxd,maxd),phase,IRtrm,vecmw(3)
      real*8 work(maxd),work2(2,maxd)
      real*8 dd(maxd,maxd),zr(maxd,maxd),zi(maxd,maxd),omega(maxd)

c Conversion factor from dsqrt(K/M) in eV and Ang to cm**-1 is 519.6
      real*8 xmagic
      parameter (xmagic=519.6d0)

      external
     .  paste
c     .  io_assign, io_close, paste

      data pi / 3.1415926d0 /
      data overflow /.false./
      data nk / 0 /

c ZPE contains the zero point energy of the system
      zpe = 0.0d0

c Planck contains planck's constant in eV.cm
      planck = 6.62618d-5*2.997925d0/1.602189d0
c ...
     

C ****************** READ DATA FROM FDF FILE *********************

c Set up fdf ...
      filein = 'stdin'
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)
c ...

c Defile Name of the system ...
      sname_defect = ' '
      sname = fdf_string('SystemName',sname_defect)
      write(6,'(a,a)')
     . 'redata: System Name                      = ',sname
c ...

c Defile System Label (short name to label files) ...
      slabel_defect  = 'vibra'
      slabel = fdf_string('SystemLabel',slabel_defect)
      write(6,'(a,a)')
     . 'redata: System Label                     = ',slabel
c ...

c Read Number of Atoms in Unit cell ...
      na_defect = 0
      natoms = fdf_integer('NumberOfAtoms',na_defect)
      if (natoms .le. 0) then
        write(6,'(a)') 
     . 'ERROR: Number of atoms must be larger than zero.'
        write(6,'(a)') 
     . '       You MUST specify NumberOfatoms in fdf file.'
        stop
      endif
      write(6,'(a,i5)') 
     . 'Number of Atoms                  = ',natoms
c     check if dimension of number of atomos is large enough
      call chkdim('fcbuild', 'maxa', maxa, natoms, 1)
c ...

c Lattice constant of unit cell...
      alat = fdf_physical('LatticeConstant',0.d0,'Bohr')
      if (alat .eq. 0.d0) then
        write(6,'(a)') 
     . 'ERROR: No valid lattice constant specified.'
        write(6,'(a)') 
     . '       You MUST specify LatticeConstant in fdf file.'
        stop
      endif
      write(6,'(a,f10.5,a)') 'Lattice Constant    = ',alat,'  Bohr'
c ...

c Lattice vectors of unit cell...
      if ( fdf_block('LatticeParameters',iunit) .and.
     .     fdf_block('LatticeVectors',iunit) ) then
         write(6,'(2a)')'ERROR: Lattice vectors doubly ',
     .     'specified: by LatticeVectors and by LatticeParameters.' 
         stop 
      endif

      if ( fdf_block('LatticeParameters',iunit) ) then
         read(iunit,*) alp, blp, clp, alplp, betlp, gamlp
         write(6,'(a)')
     .    'Lattice Parameters (units of Lattice Constant) ='
         write(6,'(a,3f10.5,3f9.3)')
     .    '    ',alp,blp,clp,alplp,betlp,gamlp
         alplp = alplp * pi/180.d0
         betlp = betlp * pi/180.d0
         gamlp = gamlp * pi/180.d0
         cell(1,1) = alp
         cell(2,1) = 0.d0
         cell(3,1) = 0.d0
         cell(1,2) = blp * cos(gamlp)
         cell(2,2) = blp * sin(gamlp)
         cell(3,2) = 0.d0
         cell(1,3) = clp * cos(betlp)
         xxx = (cos(alplp) - cos(betlp)*cos(gamlp))/sin(gamlp)
         cell(2,3) = clp * xxx
         cell(3,3) = clp * sqrt(sin(betlp)*sin(betlp) - xxx*xxx)
      elseif ( fdf_block('LatticeVectors',iunit) ) then
        do i = 1,3
          read(iunit,*) (cell(j,i), j=1,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.d0
          enddo
          cell(i,i) = 1.d0
        enddo
      endif
      write(6,'(a)') 
     .   'Lattice vectors (in units of Lattice Constant) ='
      do i = 1,3
        write(6,'(a,3f10.5)')
     .   '        ',(cell(j,i), j=1,3)
      enddo
c ...

c Multiply cell vectors by lattice constant ...
      do i = 1,3
        do ix = 1,3
          cell(ix,i) = alat * cell(ix,i)
        enddo
      enddo
      write(6,'(a)') 
     .   'Lattice vectors (in Bohr) ='
      do i = 1,3
        write(6,'(a,3f10.5)') '        ',(cell(j,i), j=1,3)
      enddo
c ...

c Read atomic coordinates and species of unit cell...
      call recoor(overflow,cell,alat,b,imass,xmass,natoms)
c ...

c Define number of unit cells in the supercell ...
      lx_defect = 0
      ly_defect = 0
      lz_defect = 0
      lxmax = fdf_integer('SuperCell_1',lx_defect)
      lymax = fdf_integer('SuperCell_2',ly_defect)
      lzmax = fdf_integer('SuperCell_3',lz_defect)
      call chkdim('fcbuild', 'maxx', maxx, lxmax, 1)
      call chkdim('fcbuild', 'maxy', maxy, lymax, 1)
      call chkdim('fcbuild', 'maxz', maxz, lzmax, 1)
      ncells = (2*lxmax+1)*(2*lymax+1)*(2*lzmax+1)
      write(6,'(a,i5)') 'lxmax    = ',lxmax
      write(6,'(a,i5)') 'lymax    = ',lymax
      write(6,'(a,i5)') 'lzmax    = ',lzmax
      write(6,'(a,i5)') 'Number of unit cells in Supercell  = ',ncells
c ...

c Determine q points to compute the phonon dispersion relations ...
      call klines(maxk, nk, nlines, lastk, label, kpoint)
c ...

c Determine whether eigenvectors are computed, or only eigenvalues ...
      eigen = fdf_boolean('Eigenvectors',.false.)
      if (eigen) then
         icall = 3
         write(6,'(a,i5)') 'Eigenvectors =   True'
         write(6,'(a,i5)') 'Computing Eigenvalues and Eigenvectors'
      else
         icall = 1
         write(6,'(a,i5)') 'Eigenvectors =   False'
         write(6,'(a,i5)') 'Computing Eigenvalues only'
      endif

c Determine whether IR intensities are to be computed
      intensity = fdf_boolean('Intensities',.false.)

c If intensities are requested then eigenvectors must be computed
      if (intensity) then
        eigen = .true.
        icall = 3
      endif

C *************** END READ DATA FROM FDF FILE ********************

c Build lattice vector of the supercell ...
      do ix=1,3
        scell(ix,1) = (2*lxmax+1)*cell(ix,1)
        scell(ix,2) = (2*lymax+1)*cell(ix,2)
        scell(ix,3) = (2*lzmax+1)*cell(ix,3)
      enddo
c ...

C Build atomic coordinates in the supercell ...
c loop over unit cells within supercell
      iatom=0
      do lx=-lxmax,lxmax
      do ly=-lymax,lymax
      do lz=-lzmax,lzmax
        r(1) = lx*cell(1,1) + ly*cell(1,2) + lz*cell(1,3)
        r(2) = lx*cell(2,1) + ly*cell(2,2) + lz*cell(2,3)
        r(3) = lx*cell(3,1) + ly*cell(3,2) + lz*cell(3,3)
c loop over atoms in unit cell
        do i=1,natoms
          iatom=iatom+1
          imasssc(iatom)=imass(i)
          do ix=1,3
            xa(ix,iatom) = b(ix,i) + r(ix)
          enddo
        enddo
      enddo
      enddo
      enddo

      nnat = natoms * (2*lxmax+1) * (2*lymax+1) * (2*lzmax+1)
      if (iatom .ne. nnat) stop 'Error computing number of atoms'
c ...

c  Determine the indices of the atoms in the central (lx=ly=lz=0) cell
c  (those that need to be displaced to calculate the force constants)  ...
      i1 = natoms * (4*lxmax*lymax*lzmax + 
     .               2*(lxmax*lymax + lxmax*lzmax + lymax*lzmax) +
     .               lxmax + lymax + lzmax) + 1
      i2 = i1 + natoms - 1
c ...

c Read Force Constants Matrix ...
      call io_assign(iunit2)
      fname = paste(slabel,'.FC')
      open(iunit2,file=fname,status='old')
      read(iunit2,*)
c Negative displacements
      do j=1,natoms
        do ij=1,3
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do i=1,natoms
              read(iunit2,*) (pp(ii,i,lx,ly,lz),ii=1,3)
            enddo
          enddo
          enddo
          enddo
c Positive displacements
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do i=1,natoms
              read(iunit2,*) (pn(ii,i,lx,ly,lz),ii=1,3)
            enddo
          enddo
          enddo
          enddo
c Average  displacements
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do i=1,natoms
              do ii=1,3
                phi0(ii,i,ij,j,lx,ly,lz)= 0.5*
     .                (pp(ii,i,lx,ly,lz)+pn(ii,i,lx,ly,lz))
              enddo
            enddo
          enddo
          enddo
          enddo
        enddo
      enddo
      call io_close(iunit2)
c ...

c If intensities are required then read Born effective charges
      if (intensity) then
        call io_assign(iunit2)
        fname = paste(slabel,'.BC')
        open(iunit2,file=fname,status='old')
        read(iunit2,*)
        do j = 1,natoms
          do ii = 1,3
            read(iunit2,*) (BornQ(ij,ii,j),ij=1,3)
          enddo
        enddo
        call io_close(iunit2)
      endif

c =================================================================
c Now form phibar(1). This force constant matrix is hermitian.

      do i=1,natoms
        do j=1,natoms
          do ii=1,3
            do ij=1,3
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
c                phi(ii,i,ij,j,lx,ly,lz) = phi0(ii,i,ij,j,lx,ly,lz)
                phi(ii,i,ij,j,lx,ly,lz) = 
     .            0.5d0 * (phi0(ii,i,ij,j,lx,ly,lz) 
     .            + phi0(ij,j,ii,i,-lx,-ly,-lz))
              enddo
              enddo
              enddo
            enddo
          enddo
        enddo
      enddo

c      goto 100
c =================================================================
c Now form zero(ix,jx,j), defined as zijbeta in the notes.

      do j=1,natoms
        do ii=1,3
          do ij=1,3
            zero(ii,ij,j) = 0.0d0
            do i=1,natoms
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                zero(ii,ij,j) = zero(ii,ij,j) + phi(ii,i,ij,j,lx,ly,lz)
              enddo
              enddo
              enddo
            enddo
            zero(ii,ij,j) = zero(ii,ij,j) / (natoms*ncells)
          enddo
        enddo
      enddo
         

c =================================================================
c Now form zeroo(ix,jx), the sum ofer beta of z(i,j,beta)

      do ii=1,3
        do ij=1,3
          zeroo(ii,ij) = 0.0d0
          do j=1,natoms
            zeroo(ii,ij) = zeroo (ii,ij) + zero(ii,ij,j)
          enddo
          zeroo(ii,ij) = zeroo(ii,ij) / natoms
        enddo
      enddo

c =================================================================
c Now form phibar(2)

      do i=1,natoms
        do j=1,natoms
          do ii=1,3
            do ij=1,3
              correct = (zeroo(ii,ij)+zeroo(ij,ii))/2.0d0 -
     .                  (zero(ii,ij,j)+zero(ij,ii,i))
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                phi(ii,i,ij,j,lx,ly,lz) = phi(ii,i,ij,j,lx,ly,lz) +
     .                                    correct
              enddo
              enddo
              enddo
            enddo
          enddo
        enddo
      enddo

100   continue

c =================================================================
c Now finally compute the dynamical matrix. We store in real-hermitian 
c form because hermdp.f wants input this way.

c Loop over k points 

      do ik=1,nk

        q(1) = kpoint(1,ik)
        q(2) = kpoint(2,ik)
        q(3) = kpoint(3,ik)

        do ix=1,3*natoms
        do jx=1,3*natoms
            dc(ix,jx)=(0.0,0.0)
        enddo
        enddo

        do i=1,natoms
        do j=1,natoms
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do jj=1,3
              r(jj) = b(jj,i)-b(jj,j)+
     .                 lx*cell(jj,1) + ly*cell(jj,2) + lz*cell(jj,3)
            enddo
            dmin=1.e8
            do llx=-1,1
            do lly=-1,1
            do llz=-1,1
              R2=0.
              do jj=1,3
                rl(jj)=llx*scell(jj,1)+lly*scell(jj,2)+llz*scell(jj,3)
                R2=R2+(rl(jj)+r(jj))**2
              enddo
              if (abs(R2-dmin) .gt. 0.0001) then
                if (R2 .lt. dmin) then
                  neq = 1
                  dmin = R2
                  qr(1)=q(1)*(rl(1)+r(1))+
     .                  q(2)*(rl(2)+r(2))+
     .                  q(3)*(rl(3)+r(3))
                  goto 10
                endif
              endif
              if (abs(R2-dmin) .lt. 0.0001) then
                neq = neq+1
                qr(neq)=q(1)*(rl(1)+r(1))+
     .                  q(2)*(rl(2)+r(2))+
     .                  q(3)*(rl(3)+r(3))
              endif
10            continue
            enddo
            enddo
            enddo
c            qr = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
            do in = 1,neq
              phase = cos(qr(in))*(1.0,0.0) + sin(qr(in))*(0.0,1.0)
              do ii=1,3
              do ij=1,3
                ix = (i-1)*3+ii
                jx = (j-1)*3+ij
                dc(ix,jx)=dc(ix,jx)+phi(ii,i,ij,j,lx,ly,lz)*phase/neq
              enddo
              enddo
            enddo
          enddo
          enddo
          enddo
        enddo
        enddo

        do ix=1,3*natoms
        do jx=1,3*natoms
          i = (ix-mod(ix-1,3))/3+1
          j = (jx-mod(jx-1,3))/3+1
          dc(ix,jx) = dc(ix,jx) / dsqrt(xmass(i)*xmass(j))
        enddo
        enddo

        do ix=1,3*natoms
          do jx=ix,3*natoms
            dd(ix,jx)=imag(dc(jx,ix))
            dd(jx,ix)=real(dc(jx,ix))
          enddo
        enddo

Cc =================================================================

        call  hermdp(dd,zr,zi,omega,maxd,3*natoms,work,work2,icall)

c the eigenvalues omega are actually omega**2.
c zr(i,j)=i'th component of j'th eigenvector. zr is the real part, and
c zi is the imaginary part. Our matrix is real, so we should get zero inmaginary
c part.!!
c =================================================================
c Compute the IR intensities if requested
      if (intensity) then

c Loop over modes
        do i = 1,3*natoms
          IRinten(i) = 0.0d0
          do ij = 1,3
            IRtrm = (0.0d0,0.0d0)
            ind = - 3
            do j = 1,natoms
              ind = ind + 3

c Mass weight eigenvectors
              rmass = 1.0d0/sqrt(xmass(j))
              do ii = 1,3
                vecmw(ii) = cmplx(zr(ind+ii,i),zi(ind+ii,i))
                vecmw(ii) = vecmw(ii)*rmass
              enddo

c Multiply eigenvectors by Born effective charges
              do ii = 1,3
                IRtrm = IRtrm + BornQ(ii,ij,j)*vecmw(ii)
              enddo

            enddo
            IRinten(i) = IRinten(i) + IRtrm*IRtrm
          enddo
        enddo

      endif

c =================================================================
c write out eigenvalues.
        do 5 i=1,3*natoms
c convert to cm**-1. the conversion factor is xmagic (511**2)
          omega(i)=xmagic*xmagic*omega(i)
          if(omega(i).gt.0.0d0)then
            omega(i)=dsqrt(omega(i))
          else
            omega(i)=-omega(i)
            omega(i)=-dsqrt(omega(i))
c        write(*,*)' Caution: omega**2 .lt.0.0 .....',
c     1	' sqrt(abs(omega2))=',omega(i)
          end if
          write(*,*)' eigenvalue #',i,' omega=',omega(i)

c Add zero point energy contribution for positive modes
          if (omega(i).gt.0.0d0) then
            zpe = zpe + omega(i)
          endif

          ek(i,ik) = omega(i)
5       continue

Cc =================================================================
c write the eigenvalues and eigenvectors to output data file.
        if (icall .eq. 3) then
          call io_assign(iunit3)
          fname = paste(slabel,'.vectors')
          open(iunit3,file=fname,status='unknown')
          if (ik.eq.1) then
	    write(6,'(/,a)')' Writing eigenvalues and eigenvectors'
	    write(6,'(2a,/)')' to output file ', fname
          else
c         go to end of file
18          read(iunit3,*,end=28)
            goto 18
          endif
28        continue
          write(iunit3,'(/,a,3f12.6)') 'k            = ',
     .                                  q(1),q(2),q(3)
	  do 20 i=1,3*natoms
	    write(iunit3,'(a,i5)')     'Eigenvector  = ',i
	    write(iunit3,'(a,f12.6)')  'Frequency    = ',omega(i)
            if (intensity) then
	      write(iunit3,'(a,f12.6)')  'IR Intensity = ',IRinten(i)
            endif
	    write(iunit3,'(a)')        'Eigenmode (real part)'
	    write(iunit3,'(3e12.4)') (zr(j,i),j=1,3*natoms)
	    write(iunit3,'(a)')        'Eigenmode (imaginary part)'
	    write(iunit3,'(3e12.4)') (zi(j,i),j=1,3*natoms)
20	  continue
	  call io_close(iunit3)
        endif
Cc =================================================================

      enddo

c Finish computation of zero point energy (averaged across k points) and output
      if (nk.gt.0) then
        zpe = 0.5d0*planck*zpe/dble(nk)
        write(*,'('' Zero point energy = '',f20.6,'' eV'')') zpe
      endif

c Write eigenvalues ...
      if (nk .gt. 0) then
        call outbands(0, 1, maxd, 3*natoms, maxk, nk, nlines, lastk,
     .                label, kpoint, ek, 0.0d0)
      endif
C ...


   33 format(7(f10.4,2x))
Cc
      stop
      end
