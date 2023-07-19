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

      program fcbuild

c *********************************************************************
c Build supercell coordinates for Force Constant Matrix calculation,
c for clusters, linear chains, slabs and 3D xtals.
c Compatible with vibrator.f
c
c Uses the FDF (Flexible Data Format) package (version 0.66.1.5)
c of J.M.Soler and A. Garcia, 
c for input and output
c
c Produces an output file in FDF format, to be used by SIESTA
c
c Written by P.Ordejon, August'98
c
c **********************************************************************

      use fdf
      implicit none

      include 'vibra.h'

c Internal variables ...

      logical overflow

      character
     .  filein*20, fileout*20

      character
     .  slabel*20, sname*150, slabel_defect*150, sname_defect*20


      integer 
     .  i, i1, i2, iatom, imass(maxa), imasssc(maxasc), iunit, ix, j,
     .  lx, ly, lz, lxmax, lymax, lzmax, 
     .  lx_defect, ly_defect, lz_defect,
     .  na_defect, natoms, ncells, nnat

      real*8 
     .  dx, alat, alp, blp, clp, alplp, betlp, gamlp, pi, xxx

      real*8
     .  b(3,maxa), cell(3,3), r(3), scell(3,3), xa(3,maxasc), 
     .  xmass(maxa)

      data pi / 3.1415926d0 /
      data overflow /.false./

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

c Multiply cell vectors by by lattice constant ...
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

c Read the atomic displacement for Force Constant calculation ...
      dx = fdf_physical('AtomicDispl',0.04d0,'Bohr')
      write(6,'(a,f10.5,a)') 'Atomic displacements = ',dx,'  Bohr'
c ...

C *************** END READ DATA FROM FDF FILE ********************

c Build lattice vector of the supercell ...
      do ix=1,3
        scell(ix,1) = (2*lxmax+1)*cell(ix,1)/alat
        scell(ix,2) = (2*lymax+1)*cell(ix,2)/alat
        scell(ix,3) = (2*lzmax+1)*cell(ix,3)/alat
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

C ************** WRITE FDF FORMAT FILE WITH SUPERCELL DATA ***********

c  Open file to write ...

      call io_assign(iunit)
      open(unit=iunit,file='FC.fdf',status='new',err=100)

c  Write Number of atoms in Supercell ...
      write(iunit,'(a,i5)') 'NumberOfAtoms       ',nnat
      write(iunit,*) 
c ...

c  Write Lattice Constant...
      write(iunit,'(a,f15.10,a)') 'LatticeConstant ',alat,'  Bohr'
      write(iunit,*) 
c ...

c  Write Lattice vectors of supercell...
      write(iunit,'(a)') '%block LatticeVectors'
      do i=1,3
        write(iunit,10) (scell(ix,i),ix=1,3)
      enddo
      write(iunit,'(a)') '%endblock LatticeVectors'
      write(iunit,*)
c ...

c  Write atomic coordinates format ...
      write(iunit,'(a)') 
     .  'AtomicCoordinatesFormat  NotScaledCartesianBohr'
      write(iunit,*)
c ...

c  Write atomic coordinates and species in supercell ...
      write(iunit,'(a)') 
     .   '%block AtomicCoordinatesAndAtomicSpecies'
      do i=1,nnat
        write(iunit,11) (xa(ix,i),ix=1,3),imasssc(i)
      enddo
      write(iunit,'(a)') 
     .   '%endblock AtomicCoordinatesAndAtomicSpecies'
      write(iunit,*)
c ...

c  Write MD options: Force Constant Calculation ...
      write(iunit,'(2a)')    'MD.TypeOfRun       ','    FC'
      write(iunit,'(a,i5)') 'MD.FCfirst          ',i1
      write(iunit,'(a,i5)') 'MD.FClast           ',i2
      write(iunit,'(a,f15.10,a)') 'MD.FCdispl          ',dx,'  Bohr'
c ...
      
      write(6,'(/,3(/,a,/),/)') 
     . '************************************************************',
     . '    Output (in FDF format) has been written in FC.fdf',
     . '************************************************************'

10    format(3(f15.10,2x))
11    format(3(f15.10,2x),i5)

      stop

100   continue
      write(6,'(/,3(/,a,/),/)') 
     . '********************************************',
     . '    ERROR: File FC.fdf already exists!',
     . '********************************************'
      stop
      end

