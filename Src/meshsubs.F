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
C Includes the following subroutines connected to the mesh :
C
C InitMesh          = initialises the mesh
C SameMeshAndAtoms  = checks whether the mesh or atoms have changed
C InitAtomMesh      = initialises quantities relating to the atoms
C                     and the mesh
C PartialCoreOnMesh = calculates the partial core density on mesh
C NeutralAtomOnMesh = calculates the neutral atom potential on mesh
C PhiOnMesh         = calculates the orbital values on the mesh
C

      module mesh
C
C Stores quantities that are connected with the mesh
C
      use precision, only : dp

      implicit none

C ----------------------------------------------------------------------
C Mesh variables
C ----------------------------------------------------------------------
C real*8  cmesh(3,3)    : Mesh-cell vectors
C real*8  dxa(3,na)     : Atom position within mesh-cell
C integer idop(mop)     : Extended-mesh-index displacement of points
C                       : within a sphere of radius rmax
C integer indexp(nep)   : Translation from extended to normal mesh index
C integer ipa(na)       : Mesh cell in which atom is
C integer mop           : Maximum number of non-zero orbital points
c integer nmeshg(3)     : Total number of mesh points in each direction
C integer ne(3)         : Number of mesh-Extension intervals 
C                       : in each direction
C integer nem(3)        : Extended-mesh divisions in each direction
C integer nmsc(3)       : Mesh divisions of each supercell vector
C integer nsm           : Number of mesh sub-divisions in each direction
C integer nsp           : Number of sub-points of each mesh point
C real*8  rcmesh(3,3)   : Reciprocal mesh-cell vectors
C                       : (WITHOUT 2*pi factor)
C real*8  xdop(3,mop)   : Vector to mesh points within rmax
C real*8  xdsp(3,nsp)   : Vector to mesh sub-points
C ----------------------------------------------------------------------

      integer, dimension(:), allocatable, save :: idop
      integer, dimension(:), allocatable, save :: indexp
      integer, dimension(:), allocatable, save :: ipa

      real(dp),  dimension(:,:), allocatable, save :: dxa
      real(dp),  dimension(:,:), allocatable, save :: xdop
      real(dp),  dimension(:,:), pointer, save :: xdsp ! Allocated in dhscf

      integer,                save :: mop
      integer,  dimension(3), save :: ne
      integer,  dimension(3), save :: nem
      integer,  dimension(3), save :: nmsc
      integer,  dimension(3), save :: nmeshg
      integer,                save :: nsm
      integer,                save :: nsp

      real(dp),  dimension(3,3),   save :: cmesh
      real(dp),  dimension(3,3),   save :: rcmesh

      end module mesh

      module meshphi

      use precision, only : phi_grid_p
C
C Stores quantities that are connected with phi on the mesh
C
      implicit none

C ----------------------------------------------------------------------
C Phi-mesh variables
C ----------------------------------------------------------------------
C logical DirectPhi     : If true the phi is calculated on the fly
C integer endpht(0:nmpl): Last position occupied by points in lstpht
C integer lstpht(ntopl) : List of non-zero orbitals at point
C integer listp2(ntopl) : Maps orbital-mesh point to iop
C integer nphi          : Length of phi array second dimension
C real    phi(nsp,ntopl): Basis orbitals at mesh points (sparse)
C ----------------------------------------------------------------------

      logical, save :: DirectPhi

      integer, save :: nphi

      integer, dimension(:), allocatable, save :: endpht
      integer, dimension(:), allocatable, save :: lstpht
      integer, dimension(:), allocatable, save :: listp2

      real(phi_grid_p), dimension(:,:), allocatable, save :: phi

      end module meshphi

      subroutine InitMesh( na, cell, norb, iaorb, iphorb, isa,
     .                     rmax, G2max, G2mesh, nsc, nmpl, 
     .                     nm, nml, ntm, ntml, ntpl, ntopl, dvol)
C
C Initialises the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C real*8  cell(3,3)     : Auxillary lattice vectors
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer isa(na)       : Species index of all atoms in supercell
C real*8  rmax          : Maximum orbital radius
C integer nsc(3)        : Number of unit-cells in each supercell direct.
C ----------------------------------------------------------------------
C Input and output :
C ----------------------------------------------------------------------
C real*8  G2max         : Effective planewave cutoff (Ry) determines
C                       : mesh density and the precision of integrals
C                       : On input : Value required
C                       : On output: Value used, which may be larger
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C real*8  G2mesh        : Effective planewave cutoff of mesh used
C integer nmpl          : Number of mesh points in unit cell locally
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer nml(3)        : Local form of nm
C integer ntm(3)        : Total number of mesh points
C integer ntml(3)       : Total number of mesh points stored locally
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer ntopl         : Total number of nonzero orbital points locally
C real*8  dvol          : Mesh-cell volume
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C real*8  dx(3)         : Vector from atom to mesh sub-point
C real*8  dxp(3)        : Vector from atom to mesh point
C integer i             : General-purpose index
C integer ia            : Looping variable for number of atoms
C integer i1,i2,i3      : Mesh indexes in each mesh direction
C integer is            : Species index
C integer isp           : Sub-Point index
C integer ity           : Orbital-type index
C integer j             : General-porpose index
C integer j1,j2,j3      : Mesh indexes in each mesh direction
C real*8  k0(3)         : Zero-vector argument for routine chkgmx
C integer nep           : Number of extended-mesh points
C integer nmp           : Number of mesh points in unit cell
C integer nty           : Number of orbital types
C integer noty(nty)     : Number of orbitals of each type
C real*8  pldist        : Distance between mesh planes
C real*8  r             : Distance between atom and mesh point
C real*8  rcty(nty)     : Radius of each orbital type
C logical same          : Indicates whether orbital types are the same
C real*8  vecmod        : Vector modulus
C real*8  volume        : Unit cell volume
C logical within        : Is a mesh point within orbital range?
C ----------------------------------------------------------------------
C Units :
C ----------------------------------------------------------------------
C
C Energies in Rydbergs
C Distances in Bohr
C

C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : HowManyMeshPerNode, GlobalToLocalMesh
      use atmfuncs,      only : rcut
      use mesh
      use meshphi,       only : endpht
      use siesta_cml

      implicit none

C
C Passed arguments
C
      integer
     .  na, ntpl, ntopl, nm(3), nml(3), nmpl, 
     .  norb, nsc(3), ntm(3), ntml(3), iaorb(norb),
     .  iphorb(norb), isa(na)

      real(dp)
     .  cell(3,3), dvol, G2max, G2mesh, rmax

C
C Local variables
C
      integer
     .  i, ia, i1, i2, i3, indi, io, iphi, is, isp, ity, j, j1, j2, j3, 
     .  ncells, nep, ntop, ntp, nty

      integer, dimension(:), allocatable, save ::
     .  iphty, isty, noty

      real(dp)
     .  dx(3), dxp(3), pldist, r, rcell(3,3), vecmod, volume

      real(dp), save ::
     .  k0(3)

      real(dp), dimension(:), allocatable, save ::
     .  rcty

      logical
     .  same, within
C
C Functions
C
      real(dp)
     .  dismin, volcel

      external
     .  chkgmx, dismin, nfft, reclat, 
     .  timer, volcel, memory

      data
     .  k0 / 3*0.0_dp /

C ----------------------------------------------------------------------
C Orbital type initialization 
C ----------------------------------------------------------------------

C Find number of orbital types
      allocate(isty(norb))
      call memory('A','I',norb,'dhscf')
      allocate(iphty(norb))
      call memory('A','I',norb,'dhscf')
      nty = 0
      do io = 1,norb
        ia = iaorb(io)
        is = isa(ia)
        iphi = iphorb(io)
        ity = 0
        same = .false.
        do while (ity .lt. nty .and. .not.same)
          ity = ity + 1
          same = (is.eq.isty(ity) .and. iphty(ity).eq.iphi)
        enddo
        if (.not.same) then
          nty = nty + 1
          isty(nty) = is
          iphty(nty) = iphi
        endif
      enddo
      call memory('D','I',size(isty),'dhscf')
      deallocate(isty)
      call memory('D','I',size(iphty),'dhscf')
      deallocate(iphty)

C Allocate local arrays dependent on nty
      allocate(iphty(nty))
      call memory('A','I',nty,'dhscf')
      allocate(isty(nty))
      call memory('A','I',nty,'dhscf')
      allocate(noty(nty))
      call memory('A','I',nty,'dhscf')
      allocate(rcty(nty))
      call memory('A','D',nty,'dhscf')

C Find number of orbital types, number of orbitals of each type,
C and cutoff radius of each type.
      noty(1:nty) = 0
      nty = 0
      do io = 1,norb
        ia = iaorb(io)
        is = isa(ia)
        iphi = iphorb(io)
        do ity = 1,nty
          if (is.eq.isty(ity) .and. iphty(ity).eq.iphi) then
            noty(ity) = noty(ity) + 1
            goto 160
          endif
        enddo
        nty = nty + 1
        isty(nty) = is
        iphty(nty) = iphi
        noty(nty) = 1
        rcty(nty) = rcut( is, iphi )
  160   continue
      enddo

      if (allocated(iphty)) then
        call memory('D','I',size(iphty),'dhscf')
        deallocate(iphty)
      endif
      if (allocated(isty)) then
        call memory('D','I',size(isty),'dhscf')
        deallocate(isty)
      endif

C ----------------------------------------------------------------------
C Mesh initialization 
C ----------------------------------------------------------------------

C Find reciprocal cell vectors (multiplied by 2*pi)
      call reclat( cell, rcell, 1 )

C Find number of mesh intervals for each cell vector.

C Loop over cell vectors
      do i = 1,3

C The reciprocal vectors of the mesh unit cell (cell/ntm)
C are rcell*ntm, and must be larger than 2*G2max

        vecmod = sqrt(dot_product(rcell(:,i),rcell(:,i)))
        ntm(i) = 2 * sqrt(G2max) / vecmod + 1

      enddo

C Return here until mesh cut-off is large enough
   50 continue

      do i = 1,3

C NFFT selects appropriate number of points for fft
   55   call nfft( ntm(i) )

C Require that ntm(i) to be a multiple of nsm
        if ( mod( ntm(i), nsm ) .ne. 0 ) then
          ntm(i) = ntm(i) + 1
          goto 55
        endif

        nm(i) = ntm(i) / nsm
        nmsc(i) = nm(i) * nsc(i)
      enddo

C Check that effective cut-off is large enough as for non-right angled
C unit cells this is not guaranteed to be the case. If cut-off needs
C to be larger, increase ntm and try again.
      G2mesh = 1.0d6
      call chkgmx( k0, rcell, ntm, G2mesh )
      if (G2mesh .lt. G2max) then
        ntm(1:3) = ntm(1:3) + 1
        goto 50
      endif

C Store number of mesh points in module array
      nmeshg(1:3) = ntm(1:3)

C Find number of mesh points in unit cell.
      ntp = ntm(1) * ntm(2) * ntm(3)

C Find local number of Mesh points of each kind
      call HowManyMeshPerNode(nm,Node,Nodes,nmpl,nml)
      ntpl = nmpl*nsp
      ntml(1:3) = nml(1:3)*nsm

C Find volume of unit cell and of mesh cell
      volume = volcel( cell )
      dvol = volume / ntp

C Output current mesh dimensions and cut-off
      if (Node.eq.0) then
        write(6,'(/,a,3(i6,a),i12)') 'InitMesh: MESH =',
     .  ntm(1),' x',ntm(2),' x',ntm(3),' =', ntp
        write(6,'(a,2f10.3,a)')
     . 'InitMesh: Mesh cutoff (required, used) =', G2max, G2mesh, ' Ry'
      endif
      if (cml_p) then
        call cmlStartPropertyList(mainXML)
        call cmlAddProperty(xf=mainXML, value=ntm,
     .       dictref='siesta:ntm', title='Mesh',
     .       units='cmlUnits:countable')
        call cmlAddProperty(xf=mainXML, value=G2max,
     .       units='siestaUnits:Ry',
     .       dictref='siesta:g2max', title='Requested Cut-Off')     
        call cmlAddProperty(xf=mainXML, value=G2mesh,
     .       units='siestaUnits:Ry',
     .       dictref='siesta:g2mesh', title='Actual Cut-Off')
        call cmlEndPropertyList(mainXML)
      endif
      G2max = G2mesh

C Find mesh-cell vectors
      do i = 1,3
        do j = 1,3
          cmesh(j,i) = cell(j,i) / nm(i)
        enddo
      enddo

C Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
      call reclat( cmesh, rcmesh, 0 )

C Find number of extended-mesh intervals for each cell vector.
C Loop over mesh directions
      do i = 1,3

C pldist is the distance between mesh planes
        pldist = 1.0_dp / sqrt(dot_product(rcmesh(:,i),rcmesh(:,i)))

C Find number of planes spanned by rmax
        ne(i) = rmax / pldist + real(nsm-1,kind=dp) / real(nsm,kind=dp)
        ne(i) = max( ne(i), nsm-2 )

C Add ne(i) points to the left and ne(I)+1 points to the
C right, to cover the spilling rmax from an atom at any
C possible place within the unit cell.
        nem(i) = nmsc(i) + 2 * ne(i) + 1

      enddo

C Find total number of extended-mesh points.
      nep = nem(1) * nem(2) * nem(3)

C Allocate local memory relating to nep
      if (allocated(indexp)) then
        call memory('D','I',size(indexp),'mesh')
        deallocate(indexp)
      endif
      allocate(indexp(nep))
      call memory('A','I',nep,'mesh')

C Find relationship between extended and unit-cell mesh points
C Loop over extended-mesh points
      do i3 = 0, nem(3)-1
        do i2 = 0, nem(2)-1
          do i1 = 0, nem(1)-1

C Find ext-mesh indexes in range [-ne(i),nm(i)+ne(i)]
            j1 = i1 - ne(1)
            j2 = i2 - ne(2)
            j3 = i3 - ne(3)

C Find normal-mesh indexes in range [0,nm(I)]
C   1000*nm(I) is added to avoid negative numbers
C   in the argument of mod
            j1 = mod( j1 + 1000 * nmsc(1), nmsc(1) )
            j2 = mod( j2 + 1000 * nmsc(2), nmsc(2) )
            j3 = mod( j3 + 1000 * nmsc(3), nmsc(3) )

C I = combined extended-mesh index.
            I = 1 + i1 + nem(1) * i2 + nem(1) * nem(2) * i3

            if (j1.lt.nm(1) .and. j2.lt.nm(2) .and. j3.lt.nm(3)) then
C indexp(i) is the equivalent point within the unit cell
              indi = 1 + j1 + nm(1) * j2 + nm(1) * nm(2) * j3
              call GlobalToLocalMesh(indi,nm,Node,Nodes,indexp(i))
            else
              indexp(i) = 0
            endif

          enddo
        enddo
      enddo

C Find sub-points
      isp = 0
      do i3 = 0, nsm-1
        do i2 = 0, nsm-1
          do i1 = 0, nsm-1
            isp = isp + 1
            do i = 1,3
              xdsp(i,isp) = ( cmesh(i,1) * i1 +
     .                        cmesh(i,2) * i2 +
     .                        cmesh(i,3) * i3 ) / nsm
            enddo
          enddo
        enddo
      enddo

C Find number of points only within rmax (orbital points)
      mop = 0
      ntop = 0

C Loop over possible points within rmax
      do i3 = -ne(3), ne(3)+1
        do i2 = -ne(2), ne(2)+1
          do i1 = -ne(1), ne(1)+1

C Find point coordinates
            do i = 1,3
              dxp(i) = cmesh(i,1) * i1 +
     .                 cmesh(i,2) * i2 +
     .                 cmesh(i,3) * i3 
            enddo

C Loop over sub-points
            within = .false.
            do isp = 1,nsp

C Find point coordinates
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 

C Find distance from point to mesh cell
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.

            enddo
            if ( within ) then
C mop is the number of mesh points within rmax
              mop = mop + 1
              do ity = 1,nty
                if ( r .lt. rcty(ity) ) ntop = ntop + noty(ity)
              enddo
            endif

          enddo
        enddo
      enddo

C Calculate ntopl
      ncells = nsc(1) * nsc(2) * nsc(3)
      ntop = ntop / ncells
      ntopl = ((ntop/Nodes) + 1)

C Allocate local memory that depends on mop
      if (allocated(idop)) then
        call memory('D','I',size(idop),'mesh')
        deallocate(idop)
      endif
      if (allocated(xdop)) then
        call memory('D','D',size(xdop),'mesh')
        deallocate(xdop)
      endif
      allocate(idop(mop))
      call memory('A','I',mop,'mesh')
      allocate(xdop(3,mop))
      call memory('A','D',3*mop,'mesh')

C Find points within rmax (orbital points)
      mop = 0
C Loop over possible points within rmax
      do i3 = -ne(3), ne(3)+1
        do i2 = -ne(2), ne(2)+1
          do i1 = -ne(1), ne(1)+1
C Find point coordinates
            do i = 1,3
              dxp(i) = cmesh(i,1) * i1 +
     .                 cmesh(i,2) * i2 +
     .                 cmesh(i,3) * i3 
            enddo
C Loop over sub-points
            within = .false.
            do isp = 1,nsp
C Find point coordinates
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
C Find distance from point to mesh cell
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
              mop = mop + 1
C Store index-distance and vector-distance to point.
              idop(mop) = i1 + nem(1) * i2 + nem(1) * nem(2) * i3
              xdop(1:3,mop) = dxp(1:3)
            endif
          enddo
        enddo
      enddo

C Deallocate orbital related memory 
      if (allocated(noty)) then
        call memory('D','I',size(noty),'dhscf')
        deallocate(noty)
      endif
      if (allocated(rcty)) then
        call memory('D','D',size(rcty),'dhscf')
        deallocate(rcty)
      endif

C Allocate memory related to nmpl
      if (allocated(endpht)) then
        call memory('D','I',size(endpht),'meshphi')
        deallocate(endpht)
      endif
      allocate(endpht(0:nmpl))
      call memory('A','I',nmpl+1,'meshphi')

      end subroutine InitMesh

      subroutine SameMeshAndAtoms(na, xa, ucell, rmax, G2max, G2mesh,
     .                            samesh, samexa)
C
C Checks whether anything has changed that requires the mesh to be
C reinitialised or quantities relating to the atoms relative to
C the mesh.
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na           : Number of atoms in the supercell
C real*8  xa(3,na)     : Coordinates of the atoms in the supercell
C real*8  ucell(3,3)   : Current unit cell vectors
C real*8  rmax         : Maximum orbital radius
C real*8  G2max        : Requested mesh cut-off
C real*8  G2mesh       : Mesh cut-off from last call
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C logical samesh       : If .true. then the mesh must be reinitialised
C logical samexa       : If .true. then atom related quantities need
C                      : to be recalculated
C ----------------------------------------------------------------------
C Internal :
C ----------------------------------------------------------------------
C integer lastna       : Number of atoms from last call
C real*8  lastxa(3,na) : Position of atoms from last call
C real*8  lastra       : Rmax from last call
C real*8  lastc(3,3)   : Unit cell from last call
C ----------------------------------------------------------------------

      use precision,     only : dp
      implicit none

C Passed arguments
      integer,   intent(in)  :: na
      real(dp),  intent(in)  :: G2max, G2mesh
      real(dp),  intent(in)  :: rmax
      real(dp),  intent(in)  :: ucell(3,3)
      real(dp),  intent(in)  :: xa(3,na)
      logical,   intent(out) :: samesh, samexa

C Local variables
! Saved
      integer, save :: lastna = 0
      real(dp), save :: tiny       = 1.0e-12_dp,
     .                  lastc(3,3) = 0.743978657912656e50_dp,
     .                  lastra     = 0.0_dp

! non-saved
      integer
     .  i, ia, j
      real(dp), dimension(:,:), allocatable, save ::
     .  lastxa

C If number of atoms has changed then deallocate lastxa
      if (na.ne.lastna) then
        if (allocated(lastxa)) then
          call memory('D','D',size(lastxa),'SameMeshAndAtoms')
          deallocate(lastxa)
        endif
      endif

C Allocate lastxa and initialise it to an unlikely number
      if (.not.allocated(lastxa)) then
        allocate(lastxa(3,na))
        call memory('A','D',3*na,'SameMeshAndAtoms')
        lastxa(1:3,1:na) = 0.0_dp
        lastxa(1,1) = 0.987654321273567e50_dp
      endif

C Find if mesh has to be changed due to unit cell
      samesh = .true.
      do i = 1,3
        do j = 1,3
          if ( ucell(j,i) .ne. lastc(j,i) ) samesh = .false.
          lastc(j,i) = ucell(j,i)
        enddo
      enddo

C Find if mesh has to be changed due to unit cell
      if ( G2max .gt. G2mesh * (1.0_dp + tiny) ) samesh = .false.

C Find if mesh has to be changed due to rmax
      if (rmax .ne. lastra) samesh = .false.
      lastra = rmax

C Find if atoms have moved having checked the number of atoms first
      samexa = (na.eq.lastna)
      if (samexa) then
        do ia = 1,na
          do i = 1,3
            if ( xa(i,ia) .ne. lastxa(i,ia) ) samexa = .false.
          enddo
        enddo
      endif

C Copy the number of atoms and coordinates to save for next call
      lastna = na
      lastxa(1:3,1:na) = xa(1:3,1:na)

C If cell has changed then it is necessary to reinitialise coordinates
      if (.not.samesh) samexa = .false.

      end subroutine SameMeshAndAtoms

      subroutine InitAtomMesh( na, xa)
C
C Initialises the atom information relating to the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C real*8  xa(3,na)      : Atomic positions of all atoms in supercell
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C all output quantities are in the mesh module
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C integer lastna        : Number of atoms on last call
C ----------------------------------------------------------------------

      use precision,     only : dp
      use mesh

      implicit none

C
C Passed arguments
C
      integer,  intent(in) :: na
      real(dp),  intent(in) :: xa(3,na)

C
C Internal variables
C
      integer
     .  i, ia, ix(3)

      integer, save ::
     .  lastna = 0

      real(dp)
     .  dx(3)

C Check size of atom-mesh arrays nad reallocate if necessary
      if (na.ne.lastna) then
        if (allocated(ipa)) then
          call memory('D','I',size(ipa),'mesh')
          deallocate(ipa)
        endif
        if (allocated(dxa)) then
          call memory('D','D',size(dxa),'mesh')
          deallocate(dxa)
        endif
      endif
      if (.not.allocated(ipa)) then
        allocate(ipa(na))
        call memory('A','I',na,'mesh')
      endif
      if (.not.allocated(dxa)) then
        allocate(dxa(3,na))
        call memory('A','D',3*na,'mesh')
      endif
      lastna = na

C Find atomic positions relative to mesh
      do ia = 1,na

C Find index of extended-mesh cell in which atom is
        do i = 1,3
          dx(i) = dot_product(xa(:,ia),rcmesh(:,i))
          ix(i) = int( dx(i) + 100000 ) - 100000
          dx(i) = dx(i) - ix(i)
          ix(i) = mod( ix(i) + 1000*nmsc(i), nmsc(i) )
          ix(i) = ix(i) + ne(i)
        enddo
        ipa(ia) = 1 + ix(1) + nem(1) * ix(2) +
     .            nem(1) * nem(2) * ix(3)

C Find atom position within mesh cell
        do i = 1,3
          dxa(i,ia) = cmesh(i,1) * dx(1) +
     .                cmesh(i,2) * dx(2) +
     .                cmesh(i,3) * dx(3)
        enddo
      enddo

      end subroutine initatommesh

      subroutine PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, 
     .                              nsd, dvol, volume, Vscf, Vaux, 
     .                              fal, stressl, Forces, Stress )
C
C Finds the partial-core-correction energy density on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C integer nsd           : Number of diagonal spin values (1 or 2)
C real*8  dvol          : Mesh-cell volume
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of SCF density
C real*4  Vaux(ntpl)    : Auxiliary potential array
C logical Forces        : Are the forces to be calculated?
C logical Stress        : Are the stresses to be calculated?
C ----------------------------------------------------------------------
C Output : (non-gradient case)
C ----------------------------------------------------------------------
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C ----------------------------------------------------------------------
C Output : (gradient case)
C ----------------------------------------------------------------------
C real*8  fal(3,:)     : Local copy of atomic forces
C real*8  stressl(3,3)  : Local copy of stress tensor
C ----------------------------------------------------------------------

      use precision,     only : dp, grid_p
      use atmfuncs, only: rcore, chcore_sub
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp

      implicit none

C
C Passed arguments
C
      integer,  intent(in)    :: na, isa(na), ntpl
      integer,  intent(in)    :: indxua(*)
      integer,  intent(in)    :: nsd
      real(dp), intent(in)    :: dvol, volume
      real(grid_p),     intent(in)    :: Vscf(ntpl,*)
      real(grid_p),     intent(in)    :: Vaux(ntpl)
      logical,  intent(in)    :: Forces
      logical,  intent(in)    :: Stress
      real(grid_p),     intent(out)   :: rhopcc(ntpl)
      real(dp), intent(inout) :: fal(3,*)
      real(dp), intent(inout) :: stressl(3,3)

C
C Internal variables
C
      integer
     .  i, ia, iop, ip, ip0, is, isp, ispin, iua

      real(dp)
     .  dfa(3), dx(3), grrho(3), r, ra, rhop, vxc

      real(dp), save ::
     .  tiny

      logical
     .  Gradients

      data tiny / 1.0e-12_dp /

C Find out whether this is a gradient run based on the arguments passed
      Gradients = ( Forces .or. Stress )

C Initialise array to zero
      if (.not.Gradients) rhopcc(1:ntpl) = 0.0_dp

      do ia = 1,na
        iua = indxua(ia)
        is = isa(ia)
        ra = rcore( is )
        if (ra .gt. tiny) then

C Loop over mesh points inside rmax
          do iop = 1,mop
            ip0 = indexp( ipa(ia) + idop(iop) )
            if (ip0 .gt. 0) then

C Loop over sub-points
              do isp = 1,nsp
                dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
                r = sqrt(dot_product(dx,dx))
c                if (r .lt. ra .and. r .gt. tiny) then
                if (r .lt. ra) then
                  ip = isp + nsp * (ip0 - 1)
                  call chcore_sub( is, dx, rhop, grrho )
                  if (Gradients) then
C Calculate gradients of PCC
                    do ispin = 1,nsd
                      vxc = Vscf(ip,ispin) - Vaux(ip)
                      do i = 1,3
                        dfa(i) = dvol * vxc * grrho(i) / nsd
                        if (Forces) fal(i,iua) = fal(i,iua) + dfa(i)
                        if (Stress) stressl(1:3,i) = stressl(1:3,i) + 
     .                    dx(1:3)*dfa(i)/volume
                      enddo
                    enddo
                  else
C Calculate density due to PCC
                    rhopcc(ip) = rhopcc(ip) + rhop
                  endif
                endif
              enddo

            endif
          enddo

        endif
      enddo

      end subroutine PartialCoreOnMesh

      subroutine NeutralAtomOnMesh( na, isa, ntpl, vna,
     .                              indxua, dvol, volume,
     .                              drho, fal, stressl,
     .                              Forces, Stress )
C
C Finds the potential due to the neutral atoms on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C real*8  dvol          : Mesh-cell volume
C real*8  volume        : Unit cell volume
C real    drho(ntpl)    : SCF density difference
C logical Forces        : Should the forces be calculated?
C logical Stress        : Should the stress be calculated?
C ----------------------------------------------------------------------
C Input/Output : (Output for non-gradient call, must be kept otherwise)
C ----------------------------------------------------------------------
C real    vna(ntpl)     : Sum of neutral-atom potentials
C ----------------------------------------------------------------------
C Output : (gradient call)
C ----------------------------------------------------------------------
C real*8  fal(3,:)     : Local copy of atomic forces
C real*8  stressl(3,3)  : Local copy of stress tensor
C ----------------------------------------------------------------------

      use precision,     only : dp, grid_p
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp

      implicit none

C
C Passed arguments
C
      integer,  intent(in)    :: na, isa(na), ntpl
      integer,  intent(in)    :: indxua(*)
      real(dp), intent(in)    :: dvol, volume
      real(grid_p),     intent(in)    :: drho(*)
      logical,  intent(in)    :: Forces, Stress
      real(grid_p),     intent(inout)   :: vna(ntpl)
      real(dp), intent(inout) :: fal(3,*)
      real(dp), intent(inout) :: stressl(3,3)

C
C Internal variables
C
      integer
     .  i, ia, iop, ip, ip0, is, isp, iua

      real(dp)
     .  dx(3), grva(3), r, ra, va

      logical 
     .  Gradients

C Check whether forces and/or stress has been requested based on arguments
      Gradients = (Forces.or.Stress)

C Initialise array to zero if not computing gradients
C If we are computing gradients, the current value of vna
C must be kept. This is particularly important for grid-cell sampling...

      if (.not.Gradients) vna(1:ntpl) = 0.0_dp

      do ia = 1,na
        iua = indxua(ia)
        is = isa(ia)
        ra = rcut( is, 0 )

C Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then

C Loop over sub-points
            do isp = 1,nsp
              dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
              r = sqrt(dot_product(dx,dx))
              if (r .lt. ra) then
                ip = isp + nsp * (ip0 - 1)
                call phiatm( is, 0, dx, va, grva )
                if (Gradients) then
                  do i = 1,3
                    grva(i) = dvol * grva(i) * drho(ip)
                    if (Forces) fal(i,iua) = fal(i,iua) + grva(i)
                    if (Stress) stressl(1:3,i) = stressl(1:3,i) +
     .                dx(1:3) * grva(i) / volume
                  enddo
                else
                  vna(ip) = vna(ip) + va
                endif
              endif
            enddo

          endif
        enddo

      enddo

      end subroutine NeutralAtomOnMesh

      subroutine LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
C
C Finds the diffuse ionic charge, whose electrostatic potential is equal
C to Vlocal on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C ----------------------------------------------------------------------
C Output : 
C ----------------------------------------------------------------------
C real    Chlocal(ntpl)     : Sum of diffuse ionic charge densities
C ----------------------------------------------------------------------
 
      use precision,     only : dp, grid_p
      use atmfuncs, only: rcut, psch
      use atm_types, only : species
      use radial, only : rad_func
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
 
      implicit none
 
C
C Passed arguments
C
      integer, intent(in)  :: na, isa(na), ntpl
      integer, intent(in)  :: indxua(*)
      real(grid_p),    intent(out) :: Chlocal(ntpl)
 
C
C Internal variables
C
      integer
     .  ia, iop, ip, ip0, is, isp, iua
 
      real(dp)
     .  dx(3), grpsch(3), r, ra, pschloc

      type(rad_func), pointer            :: func   
 
C Initialise array to zero 
      Chlocal(1:ntpl) = 0.0_grid_p
 
      do ia = 1,na
        iua = indxua(ia)
        is = isa(ia)
        func => species(is)%chlocal
        ra = func%cutoff
 
C Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
 
C Loop over sub-points
            do isp = 1,nsp
              dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
              r = sqrt(dot_product(dx,dx))
              if (r .lt. ra) then
                ip = isp + nsp * (ip0 - 1)
                call psch( is, dx, pschloc, grpsch )
                Chlocal(ip) = Chlocal(ip) + pschloc
              endif
            enddo
 
          endif
        enddo
 
      enddo
 
      end subroutine LocalChargeOnMesh

      subroutine ModelCoreChargeOnMesh( na, isa, ntpl,
     $                                  ModelCore, indxua )
C
C     Uses the "local charge" information to construct a model total
C     core charge, useful for Bader analysis
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C ----------------------------------------------------------------------
C Output : 
C ----------------------------------------------------------------------
C real    ModelCore(ntpl)     : Sum of model core charge densities
C ----------------------------------------------------------------------
 
      use precision,     only : dp, grid_p
      use atmfuncs, only: rcut, psch, izofis, zvalfis
      use atm_types, only : species
      use radial, only : rad_func
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
 
      implicit none
 
C
C Passed arguments
C
      integer, intent(in)  :: na, isa(na), ntpl
      integer, intent(in)  :: indxua(*)
      real(grid_p),    intent(out) :: ModelCore(ntpl)
 
C
C Internal variables
C
      integer
     .  ia, iop, ip, ip0, is, isp, iua
 
      real(dp) dx(3), grpsch(3), r, ra, charge_ratio
      real(dp) core_charge, beta, pschloc

      type(rad_func), pointer            :: func   
 
C Initialise array to zero 
      ModelCore(1:ntpl) = 0.0_grid_p
 
      do ia = 1,na
        iua = indxua(ia)
        is = isa(ia)
        func => species(is)%chlocal
        ra = func%cutoff

        ! We will take advantage of the existing "local charge"
        ! by rescaling.

        ! The total local charge is Z_val
        ! We want a total extra core charge of (Z-Zval),
        ! and with the opposite sign, localized in a
        ! region of around 1.0 bohr for "normal" atoms and
        ! 0.4 bohr (= 0.21 Ang) for Hydrogen.
        ! Note that this is quite a small confinement region, so
        ! you might need a high cutoff (of the order of 300-500 Ry)
        ! to represent it on the grid. It might be better to use
        ! Denchar for this, rather than performing a Siesta calculation
        ! (even if it reads the converged DM)

        ! --- For Hydrogen, add an artificial core with 1 electron
        ! Else choose the normal core charge
        !
        ! The cutoff for chlocal is divided by beta to obtain the
        ! tighter localization radius. Note the factor of beta^3
        ! for normalization of the new function.

        if (izofis(is) == 1) then
           beta = ra / 0.4_dp
           charge_ratio = -1.0_dp * beta**3
        else
           beta = ra / 1.0_dp
           charge_ratio = - beta**3 *
     $              (izofis(is)-zvalfis(is))/zvalfis(is)
        endif
 
C Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
 
C Loop over sub-points
            do isp = 1,nsp
              dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
              r = sqrt(dot_product(dx,dx))
              if (r .lt. ra) then
                ip = isp + nsp * (ip0 - 1)
                dx = beta*dx
                call psch( is, dx, pschloc, grpsch )
                core_charge = charge_ratio*pschloc
                ModelCore(ip) = ModelCore(ip) + core_charge
              endif
            enddo
 
          endif
        enddo
 
      enddo
 
      end subroutine ModelCoreChargeOnMesh

      subroutine PhiOnMesh( nmpl, norb, iaorb, iphorb, isa)
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom)
C integer isa(na)       : Species index of all atoms in supercell
C integer nmpl          : Number of mesh points in unit cell locally
C integer norb          : Total number of basis orbitals in supercell
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshphi
C ----------------------------------------------------------------------

C
C Modules
C
      use precision,     only : dp, grid_p
      use atmfuncs, only : rcut, phiatm
      use mesh,     only : dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      use meshphi,  only : directphi, nphi, phi, lstpht, listp2, endpht
      use parallel, only : Node
      use fdf,      only : fdf_boolean
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C
C Passed arguments
C
      integer, intent(in) :: nmpl
      integer, intent(in) :: norb
      integer, intent(in) :: iaorb(norb)
      integer, intent(in) :: iphorb(norb)
      integer, intent(in) :: isa(*)

C
C Local variables
C
      integer
     .  ia, io, iop, ip, ip0, iphi, is, isp, n, nlist, nliste

#ifdef MPI
      integer
     .  MPIerror
#endif

      logical
     .  within

      logical, save :: firsttime = .true.

      integer, dimension(:), allocatable, save :: numphi

      real(dp)
     .  dxsp(3,nsp), grphi(3), phip, r2o, r2sp(nsp)

C On first call set the logical DirectPhi
      if (firsttime) then
        if (Node.eq.0) then
          DirectPhi  = fdf_boolean('DirectPhi',.false.)
        endif
#ifdef MPI
        call MPI_Bcast(DirectPhi,1,MPI_logical,0,MPI_Comm_World,
     *    MPIerror)
#endif
        firsttime = .false.
      endif

C Allocate local scratch array
      allocate( numphi(nmpl) )
      call memory('A','I',nmpl,'PhiOnMesh')
      numphi(1:nmpl) = 0

C Find number of atomic orbitals at mesh points
      do io = 1,norb
        ia = iaorb(io)
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2

C Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then

C Loop over sub-points to find if point is within range
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp) = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo

C If within range, add one to number of point orbitals
            if (within) numphi(ip0) = numphi(ip0) + 1

          endif
        enddo
      enddo

C Initialise pointer array
      endpht(0) = 0
      do ip = 1,nmpl
        endpht(ip) = endpht(ip-1) + numphi(ip)
      enddo

C Allocate phi if this is not a direct calculation
      nlist = endpht(nmpl)
      if (DirectPhi) then
        nphi = 1
      else
        nphi = nlist
      endif

C Add an extra margin of error to nlist to minimise reallocations
      nliste = 1.01*nlist

C Adjust dimensions of phi if necessary
      if (allocated(phi)) then
        if (nphi.gt.(size(phi,2))) then
          call memory('D','X',size(phi),'meshphi')
          deallocate(phi)
          if (DirectPhi) then
            allocate(phi(nsp,nphi))
            call memory('A','X',nsp*nphi,'meshphi')
          else
            allocate(phi(nsp,nliste))
            call memory('A','X',nsp*nliste,'meshphi')                       
          endif                                                             
        endif                                                               
      else                                                                  
        if (DirectPhi) then                                                 
          allocate(phi(nsp,nphi))                                           
          call memory('A','X',nsp*nphi,'meshphi')                           
        else                                                                
          allocate(phi(nsp,nliste))                                         
          call memory('A','X',nsp*nliste,'meshphi')                         
        endif                                                               
      endif

C Adjust dimensions of list arrays if necessary
      if (allocated(lstpht)) then
        if (nlist.gt.size(lstpht)) then
          call memory('D','I',size(lstpht),'meshphi')
          deallocate(lstpht)
          allocate(lstpht(nliste))
          call memory('A','I',nlist,'meshphi')
        endif
      else
        allocate(lstpht(nliste))
        call memory('A','I',nliste,'meshphi')
      endif
      if (allocated(listp2)) then
        if (nlist.gt.size(listp2)) then
          call memory('D','I',size(listp2),'meshphi')
          deallocate(listp2)
          allocate(listp2(nliste))
          call memory('A','I',nliste,'meshphi')
        endif
      else
        allocate(listp2(nliste))
        call memory('A','I',nliste,'meshphi')
      endif

C Find indexes and values of atomic orbitals at mesh points
      numphi = 0
      do io = 1,norb
        ia = iaorb(io)
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp) = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo
            if (within) then
              numphi(ip0) = numphi(ip0) + 1
              n = endpht(ip0-1) + numphi(ip0)
              lstpht(n) = io
              listp2(n) = iop
              if (.not.DirectPhi) then
                do isp = 1,nsp
                  if (r2sp(isp) .lt. r2o) then
                    call phiatm( is, iphi, dxsp(1,isp),
     .                           phip, grphi )
                    phi(isp,n) = phip
                  else
                    phi(isp,n) = 0.0_grid_p
                  endif
                enddo
              endif
            endif
          endif
        enddo
      enddo

C Free local memory
      call memory('D','I',size(numphi),'PhiOnMesh')
      deallocate( numphi )

      end subroutine PhiOnMesh
