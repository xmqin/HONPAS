! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine setspatial(na,xa,vecs,rcut,lspatialok)
C
C  Sets up the spatial decomposition of the atoms into boxes
C
C  Based on the equivalent routine from GULP, except that it 
C  now uses fractional coordinates to handle low symmetry cells.
C
C  On entry :
C
C  na            = number of atoms in cell
C  xa(3,na)      = Cartesian coordinates of atoms
C  vecs(3,3)     = cell vectors
C  rcut          = maximum potential cut-off radius
C  Node          = local Node number
C  Nodes         = total number of processors
C
C  On exit :
C
C  lspatialok    = if .true. then the cell can be decomposed
C                  into boxes
C
C  Set in modules:
C
C  nspcellat(ncell)        = number of atoms per spatial cell
C  nspcellatptr(maxspcellat, ncell)  = pointer to global atom number from atom number in cell
C  nspcellatptrcell =  pointer to image (1-27) of atom in cell
C  spmin(3)            = min spatial cell length in x, y, z
C  spcell           =  
C  nspcell(3)          = no. spatial cells in x, y, z
C  maxspcell        = 
C  maxspcellat      = maxval(nspcellat) + a bit extra
C  
C  nspcellat(ncell)              = number of atoms per cell
C  nspcellatptr(natom,ncell)     = pointer to atom from atom number in cell
C  nspcellatptrcell(natom,ncell) = pointer to image (1-27) of atom in cell

C  
C
C  Julian Gale, NRI, Curtin University, March 2004
C
      use parallel,   only : Node
      use precision,  only : dp
      use spatial,    only : nspcellat, nspcellatptr, nspcellatptrcell
     .                     , nbufferx, nbuffery, nbufferz
     .                     , spmin, spcell, nspcell
     .                     , maxspcell, maxspcellat
      use alloc,      only : re_alloc, de_alloc
      implicit none
C
C  Passed variables
C
      integer,     intent(in)     :: na
      logical,     intent(out)    :: lspatialok
      real(dp),    intent(in)     :: xa(3,na)
      real(dp),    intent(in)     :: vecs(3,3)
      real(dp),    intent(in)     :: rcut
C
C  Local variables
C
      integer                     :: i
      integer                     :: ii
      integer                     :: ind
      integer                     :: ixvec1cell(27)
      integer                     :: iyvec1cell(27)
      integer                     :: izvec1cell(27)
      integer                     :: maxxy
      integer                     :: maxx
      integer                     :: n
      integer                     :: n2bufferplus1x
      integer                     :: n2bufferplus1y
      integer                     :: n2bufferplus1z
      integer                     :: np
      integer                     :: npc
      integer                     :: nx
      integer                     :: ny
      integer                     :: nz
      integer                     :: nspcelltot
      logical                     :: ldebug
      logical,               save :: firstcall = .true.
      real(dp)                    :: a
      real(dp)                    :: b
      real(dp)                    :: c
      real(dp)                    :: alpha
      real(dp)                    :: beta
      real(dp)                    :: gamma
      real(dp)                    :: degtorad
      real(dp)                    :: small
      real(dp)                    :: xdivision
      real(dp)                    :: ydivision
      real(dp)                    :: zdivision
      real(dp), pointer           :: xyzfrac(:,:)
      real(dp)                    :: xi
      real(dp)                    :: yi
      real(dp)                    :: zi
      real(dp)                    :: xvec1cell(27)
      real(dp)                    :: yvec1cell(27)
      real(dp)                    :: zvec1cell(27)
      real(dp), dimension(:),   pointer       :: xinbox
      real(dp), dimension(:),   pointer       :: yinbox
      real(dp), dimension(:),   pointer       :: zinbox
C
C  Set the following according to whether debugging is required
C
c      ldebug = .false.
      ldebug = (Node.eq.0)
C
C  On first call initialise spatial decomposition arrays
C
      if (firstcall) then
        nullify(nspcellat)
        nullify(nspcellatptr)
        nullify(nspcellatptrcell)
        nullify(xinbox)
        nullify(yinbox)
        nullify(zinbox)
        firstcall = .false.
      endif
C
C  Allocate global memory
C
      call re_alloc(xinbox,1,na,name='xinbox')
      call re_alloc(yinbox,1,na,name='yinbox')
      call re_alloc(zinbox,1,na,name='zinbox')
C
C  Allocate local workspace
C
      nullify(xyzfrac)
      call re_alloc(xyzfrac,1,3,1,na,name="xyzfrac")
C
      small = 0.000001d0
      degtorad = 4.0d0*atan(1.0d0)/180.0d0
C
C  Find cell parameters
C
      call uncell(vecs,a,b,c,alpha,beta,gamma,degtorad)
C
C  Check that cut-off is greater than zero
C
      lspatialok = .true.
      if (rcut.lt.1.0d-10) lspatialok = .false.
C
C  If cell is not compatible return
C
      if (.not.lspatialok) return
C
C  Copy current atomic coordinates into box coordinate arrays
C
      xinbox(1:na) = xa(1,1:na)
      yinbox(1:na) = xa(2,1:na)
      zinbox(1:na) = xa(3,1:na)
C
C  Find cell lengths in each direction 
C
      spmin(1:3) = - small
      spcell(1) = a
      spcell(2) = b
      spcell(3) = c
C
C  Find numbers of cells in each direction - subtract a small amount for safety
C  Add 2*nbuffer to number of cells in each direction to allow for buffer on each side
C
      nspcell(1) = 2*nbufferx + (spcell(1) - 1.0d-6)/rcut
      nspcell(2) = 2*nbuffery + (spcell(2) - 1.0d-6)/rcut
      nspcell(3) = 2*nbufferz + (spcell(3) - 1.0d-6)/rcut
C
C  Check that minimum number of cells in any direction is at least 2*nbuffer + 1
C
      n2bufferplus1x = 2*nbufferx + 1
      n2bufferplus1y = 2*nbuffery + 1
      n2bufferplus1z = 2*nbufferz + 1
      nspcell(1) = max(nspcell(1),n2bufferplus1x)
      nspcell(2) = max(nspcell(2),n2bufferplus1y)
      nspcell(3) = max(nspcell(3),n2bufferplus1z)
C
C  Set up cell images
C
      call cellimagelist(vecs,xvec1cell,yvec1cell,zvec1cell,
     .                   ixvec1cell,iyvec1cell,izvec1cell)
C
C  Create an array of fractional coordinates to decide which box each atom belongs to
C  and place Cartesian coordinates within central cell
C
      call cart2frac(na,xinbox,yinbox,zinbox,vecs,xyzfrac)
C
C  Calculate the total number of cells and check that arrays are allocated to this size
C
      nspcelltot = nspcell(1)*nspcell(2)*nspcell(3)
      if (nspcelltot.gt.maxspcell) then
        maxspcell = nspcelltot
        call re_alloc(nspcellat,1,maxspcell,name='nspcellat')
      endif
C
C  Initialise atom counters
C
      nspcellat(1:nspcelltot) = 0
C
C  Find which box each atom belongs in. This is described by the arrays:
C
C  nspcellat(ncell)              = number of atoms per cell
C  nspcellatptr(natom,ncell)     = pointer to atom from atom number in cell
C  nspcellatptrcell(natom,ncell) = pointer to image (1-27) of atom in cell
C
C  Have to loop over cell images in order to fill buffer regions too
C
      xdivision = nspcell(1) - 2*nbufferx
      ydivision = nspcell(2) - 2*nbuffery
      zdivision = nspcell(3) - 2*nbufferz
      maxxy = nspcell(1)*nspcell(2)
      maxx  = nspcell(1)
      do i = 1,na
        do ii = 1,27
          nx = (xyzfrac(1,i) + ixvec1cell(ii) - spmin(1))*xdivision + 
     .         2*nbufferx
          if (nx.ge.1.and.nx.le.nspcell(1)) then
            ny = (xyzfrac(2,i) + iyvec1cell(ii) - spmin(2))*ydivision + 
     .           2*nbuffery
            if (ny.ge.1.and.ny.le.nspcell(2)) then
              nz = (xyzfrac(3,i) + izvec1cell(ii) - spmin(3))*zdivision+
     .             2*nbufferz
              if (nz.ge.1.and.nz.le.nspcell(3)) then
                ind = (nz-1)*maxxy + (ny-1)*maxx + nx
                nspcellat(ind) = nspcellat(ind) + 1
                if (nspcellat(ind).gt.maxspcellat) then
                  maxspcellat = nspcellat(ind) + 50
                  call re_alloc(nspcellatptr,1,maxspcellat,
     .              1,maxspcell,name='nspcellatptr')
                  call re_alloc(nspcellatptrcell,1,maxspcellat,
     .              1,maxspcell,name='nspcellatptrcell')
                endif
                nspcellatptr(nspcellat(ind),ind) = i
                nspcellatptrcell(nspcellat(ind),ind) = ii
              endif
            endif
          endif
        enddo
      enddo
C
C  Debugging information
C
      if (ldebug) then
        write(6,'(/,''  Spatial decomposition data : '',/)')
        write(6,'(''  Cutoff = '',f8.4,'' Angstroms'',/)') 
     .    rcut*0.529177d0
        write(6,'(''  Spatial decomposition cells (containing atoms) :''
     .    ,/)')
        write(6,'(6x,''Ncells'',4x,''Cell min'',4x,''Cell length'')')
        write(6,'(''  x : '',i6,2(3x,f9.4))') nspcell(1),spmin(1),
     .    spcell(1)
        write(6,'(''  y : '',i6,2(3x,f9.4))') nspcell(2),spmin(2),
     .    spcell(2)
        write(6,'(''  z : '',i6,2(3x,f9.4))') nspcell(3),spmin(3),
     .    spcell(3)
        write(6,'(/)')
        write(6,'(''----------------------------------------------'',
     .    ''----------------------------------'')')
        write(6,'('' Cell   No. of atoms    Atom No.      x       '',
     .    ''     y            z'')')
        write(6,'(''----------------------------------------------'',
     .    ''----------------------------------'')')
        do i = 1,nspcelltot
          if (nspcellat(i).gt.0) then
            np = nspcellatptr(1,i)
            npc = nspcellatptrcell(1,i)
            xi = xinbox(np) + xvec1cell(npc)
            yi = yinbox(np) + yvec1cell(npc)
            zi = zinbox(np) + zvec1cell(npc)
            write(6,'(i6,4x,i8,4x,i8,3(1x,f12.4))') i,nspcellat(i),
     .        np,xi,yi,zi
            do n = 2,nspcellat(i)
              np = nspcellatptr(n,i)
              npc = nspcellatptrcell(n,i)
              xi = xinbox(np) + xvec1cell(npc)
              yi = yinbox(np) + yvec1cell(npc)
              zi = zinbox(np) + zvec1cell(npc)
              write(6,'(22x,i8,3(1x,f12.4))') np,xi,yi,zi
            enddo
            write(6,'(''--------------------------------------------'',
     .        ''------------------------------------'')')
          endif
        enddo
        write(6,'(/)')
      endif
C
C  Deallocate local workspace
C
      call de_alloc(xyzfrac,name='xyzfrac')
      call de_alloc(xinbox,name='xinbox')
      call de_alloc(yinbox,name='yinbox')
      call de_alloc(zinbox,name='zinbox')
C
      return
      end
