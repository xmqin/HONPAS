! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!-------------------------------------------------------------
!     
      module meshsubs

      implicit none

      public :: initMesh      ! initialises the mesh
      public :: initAtomMesh  !   " atom info on the mesh
      public :: sameMeshAndAtoms  ! have the mesh or atoms changed?
      public :: PartialCoreOnMesh ! partial core density on the mesh
      public :: NeutralAtomOnMesh ! neutral atom potential on mesh
      public :: LocalChargeOnMesh ! local (positive ion) charge
      public :: ModelCoreChargeOnMesh ! model core charge for Bader
      public :: PhiOnMesh     ! Orbital values on the mesh
      public :: setupExtMesh  ! Sets up extended mesh
      public :: distriPhiOnMesh  ! Sets up data distribution

      private

      CONTAINS

      subroutine InitMesh( na, cell, norb, iaorb, iphorb, isa,
     &                     rmax, G2max, G2mesh, nsc, nmpl, 
     &                     nm, nml, ntm, ntml, ntpl, dvol)

!  Modules
      use precision,     only : dp, i8b
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : HowManyMeshPerNode, GlobalToLocalMesh
      use moreMeshSubs,  only : initMeshDistr, setMeshDistr
      use moreMeshSubs,  only : UNIFORM
      use alloc,         only : re_alloc, de_alloc
      use mesh,          only : ne, mop, nmeshg, nmsc, nsm, nsp
      use mesh,          only : cmesh, rcmesh, xdsp, nmuc, nusc
      use siesta_cml
      use cellsubs,      only : reclat  ! Finds reciprocal unit cell
      use cellsubs,      only : volcel  ! Finds unit cell volume
#ifdef MPI
      use mpi_siesta
#endif

! Passed arguments
      integer, intent(in):: na           ! Number of atoms in supercell
      real(dp),intent(in):: cell(3,3)    ! Unit cell vectors
      integer, intent(in):: norb         ! Total number of basis
                                         ! orbitals in supercell
      integer, intent(in):: iaorb(norb)  ! Atom to which orbitals belong
      integer, intent(in):: iphorb(norb) ! Orbital index (within atom)
                                         ! of each orbital
      integer, intent(in):: isa(na)      ! Species index of each atom
      real(dp),intent(in):: rmax         ! Maximum orbital radius
      integer, intent(in):: nsc(3)       ! Number of unit-cells in each
                                         ! supercell direction
      real(dp),intent(inout):: G2max     ! Effective planewave cutoff
                                         ! On input : Value required
                                         ! On output: Value used, which
                                         ! may be larger
      real(dp),intent(out):: G2mesh      ! Planewave cutoff of mesh used
                                         ! (same as G2max on output)
      integer, intent(out):: nmpl        ! Number of mesh points stored
                                         ! by my processor
      integer, intent(out):: nm(3)       ! Number of mesh divisions of
                                         ! each unit cell vector
      integer, intent(out):: nml(3)      ! Sizes of my processor's box
                                         ! of mesh points
      integer, intent(out):: ntm(3)      ! Mesh divisions of each unit
                                         ! cell vector, incl. subpoints
      integer, intent(out):: ntml(3)     ! Sizes of my processor's mesh
                                         ! box, including subpoints
      integer, intent(out):: ntpl        ! Mesh points stored by my
                                         ! processor, incl. subpoints
!!OLD      integer, intent(out):: ntopl       ! Total number of nonzero
                                         ! orbital points stored by my
                                         ! processor
      real(dp),intent(out):: dvol        ! Volume per mesh (super)point

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
C integer j             : General-porpose index
C integer j1,j2,j3      : Mesh indexes in each mesh direction
C integer nep           : Number of extended-mesh points
C integer nmp           : Number of mesh points in unit cell
C real*8  pldist        : Distance between mesh planes
C real*8  r             : Distance between atom and mesh point
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

C     Local variables
      integer                 :: i, ia, i1, i2, i3, io, iphi, is, isp,
     &                           ity, j, nep, j1, j2, j3, indi
      integer(i8b)            :: ntp
      real(dp)                :: dx(3), dxp(3), pldist, r,
     &                           volume
      logical                 :: within
C     Functions
      real(dp), external      :: dismin
!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE InitMesh' )
#endif

! ----------------------------------------------------------------------
! Mesh initialization
! ----------------------------------------------------------------------
      ! determines ntm
      call get_mesh_sizes(cell,G2max,nsm,ntm,G2mesh)

! Save some numbers 
!     nmeshg in module array
      nm(:) = ntm(:) / nsm
      nmsc(:) = nm(:) * nsc(:)
      nmeshg(:) = ntm(:) ! Total number of mesh points, incl. subpoints
      nmuc(:) = nm(:)    ! Mesh points in each unit cell vector
      nusc(:) = nsc(:)   ! Unit cells in each supercell direction

! Find number of mesh points in unit cell.
! Note needed integer*8 !
!!AG** Check proper kind name

      ntp = product(INT(ntm,kind=i8b))

C     Create the first mesh distribution
      call initMeshDistr( oDistr=UNIFORM, nm=nm )

C     Find and sets local number of Mesh points of each kind
!     Sets nml, nmpl, ntml and ntpl
      call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )

C     Find volume of unit cell and of mesh cell
      volume = volcel( cell )
      dvol = volume / ntp

C     Output current mesh dimensions and cut-off
      if (Node.eq.0) then
        write(6,'(/,a,3(i6,a),i12)') 'InitMesh: MESH =',
     &        ntm(1),' x',ntm(2),' x',ntm(3),' =', ntp
        write(6,'(a,3(i6,a),i12)') 'InitMesh: (bp) =',
     &        nm(1),' x',nm(2),' x',nm(3),' =', product(INT(nm,8))
        write(6,'(a,2f10.3,a)')
     &        'InitMesh: Mesh cutoff (required, used) =',
     &        G2max, G2mesh, ' Ry'
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

C     Find mesh-cell vectors
      do i = 1,3
        do j = 1,3
          cmesh(j,i) = cell(j,i) / nm(i)
        enddo
      enddo

C     Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
      call reclat( cmesh, rcmesh, 0 )

! Find number ne(:) of extended-mesh intervals for each cell vector.
      do i = 1,3
        ! pldist is the distance between mesh planes
        pldist = 1.0_dp / sqrt(dot_product(rcmesh(:,i),rcmesh(:,i)))
        ! Find number of planes spanned by rmax
        ne(i) = rmax / pldist
      enddo ! i

! For an atom at x=0, ne=rmax/pldist is the last mesh point within rmax.
! But an atom at mesh cell ix=0 can be almost at x=dx. Therefore we have
! to go up to ne+1. And subpoints in mesh cell ix=-ne-1 can be almost at
! x=-ne*dx. Therefore we have to go up to -ne-1 to the left. Thus, we 
! just increase ne and forget about these two effects from now on.
      ne(:) = ne(:) + 1

C     Find sub-points
!$OMP parallel default(shared), private(i1,i2,i3,i,isp)

!$OMP do
      do i3 = 0, nsm-1
        isp = i3 * nsm ** 2
        do i2 = 0, nsm-1
          do i1 = 0, nsm-1
            isp = isp + 1
            do i = 1,3
              xdsp(i,isp) = ( cmesh(i,1) * i1 +
     &                        cmesh(i,2) * i2 +
     &                        cmesh(i,3) * i3 ) / nsm
            enddo
          enddo
        enddo
      enddo
!$OMP end do

! Find number of points within rmax (orbital points)
!$OMP single ! Only used for OpenMP
      mop = 0
!$OMP end single 
!$OMP do private(dxp,within,dx,r), reduction(+:mop)
      do i3 = -ne(3), ne(3)    ! Loop over neighbor (super) points
        do i2 = -ne(2), ne(2)
          do i1 = -ne(1), ne(1)
            dxp(:) = cmesh(:,1) * i1 +  ! (Super) point coordinates
     .               cmesh(:,2) * i2 +
     .               cmesh(:,3) * i3 
            ! Find if any subpoint can be within rmax of an atom that
            ! is in the origin's mesh cell
            within = .false.
            do isp = 1,nsp  ! Loop over sub-points
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp)  ! Subpoint coords.
              ! Find distance from point to the mesh cell at the origin
              ! (since the atom might be at any point of the mesh cell)
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo ! isp
            if ( within ) then
              mop = mop + 1      ! Total number of points within rmax
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3
!$OMP end do nowait

!$OMP end parallel

C     Create extended mesh arrays for the first data distribution
      call setupExtMesh( UNIFORM, rmax )

#ifdef DEBUG
      call write_debug( '      POS InitMesh' )
#endif
!----------------------------------------------------------------------- END
      end subroutine InitMesh
!------------------------------------

      subroutine get_mesh_sizes(cell,G2max,nsm,ntm,RealCutoff)
      ! Computes the mesh dimensions ntm for a given cell 
      ! and mesh cutoff, making them multiples of the big-point
      ! sampling size nsm.

      use precision,     only : dp
      use cellsubs,      only : reclat  ! Finds reciprocal unit cell
      use m_chkgmx,      only : chkgmx  ! Checks planewave cutoff of a mesh
      use fft1d,         only : nfft    ! Finds allowed value for 1-D FFT
      use fdf
      use sys,           only : die, message

      real(dp), intent(in)  :: cell(3,3)
      real(dp), intent(in)  :: G2max
      real(dp), intent(out) :: RealCutoff
      integer, intent(in)   :: nsm
      integer, intent(out)  :: ntm(3)
      
      real(dp)  :: rcell(3,3), vecmod
      integer   :: i, itmp
      real(dp), parameter     :: k0(3) = (/ 0.0, 0.0, 0.0 /)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      ! Find reciprocal cell vectors (multiplied by 2*pi)
      call reclat( cell, rcell, 1 )

! Find number of mesh intervals for each cell vector.
! The reciprocal vectors of the mesh unit cell (cell/ntm)
! are rcell*ntm, and must be larger than 2*G2max
      do i = 1,3
        vecmod = sqrt(dot_product(rcell(:,i),rcell(:,i)))
        ntm(i) = 2 * sqrt(G2max) / vecmod + 1
      enddo

! Impose that mesh cut-off is large enough
      impose_cutoff: do

        ! Require that ntm is suitable for FFTs and a multiple of nsm
        do i = 1,3
          impose_fft: do
            ! nfft increases ntm to the next integer suitable for FFTs
            call nfft( ntm(i) )
            if ( mod( ntm(i), nsm )==0 ) exit impose_fft
            ntm(i) = ntm(i) + 1
          end do impose_fft
        enddo ! i

        ! Check that effective cut-off is large enough as for non-right
        ! angled unit cells this is not guaranteed to be the case.
        ! If cut-off needs to be larger, increase ntm and try again.
        ! chkgmx decreases G2mesh to the right cutoff value of the mesh
        RealCutoff = huge(1.0_dp)
        call chkgmx( k0, rcell, ntm, RealCutoff )
        if (RealCutoff >= G2max) then
          exit impose_cutoff
        else
          ntm(1:3) = ntm(1:3) + 1
        end if

      end do impose_cutoff


!     Let the user manually decide the mesh-box size
!     This option may be useful when testing Mesh.SubDivisions
!     for performance.
!     NOTE/TODO: manual specification is not checked with respect
!                to nfft (allowed integers). I don't know if this should
!                forced or not?
      if ( fdf_block('Mesh.Sizes',bfdf) ) then
        
        if ( fdf_bline(bfdf, pline) ) then
           ntm(1) = fdf_bintegers(pline,1)
           ntm(2) = fdf_bintegers(pline,2)
           ntm(3) = fdf_bintegers(pline,3)
        else
           call die('initmesh: ERROR in Mesh.Sizes block')
        endif
        call fdf_bclose(bfdf)
        
        ! Calculate the actual meshcutoff from the grid
        RealCutoff = huge(1.0_dp)
        call chkgmx( k0, rcell, ntm, RealCutoff )

      else if ( fdf_islist('Mesh.Sizes') ) then

        i = -1
        call fdf_list('Mesh.Sizes', i, ntm)
        if ( i /= 3 )
     &      call die('initmesh: ERROR in Mesh.Sizes list')
        ! Do the actual read of the mesh sizes
        call fdf_list('Mesh.Sizes', i, ntm)

        ! Calculate the actual meshcutoff from the grid
        RealCutoff = huge(1.0_dp)
        call chkgmx( k0, rcell, ntm, RealCutoff )
        
      end if

!     Issue a warning if the mesh sizes are smaller than
!     the requested mesh cutoff
      if ( RealCutoff < G2max ) then
        call message("WARNING","Read mesh sizes give a smaller"
     $      // " cutoff than requested")
      end if

!     Check whether the resulting mesh sizes are multiples
!     of nsm.
      if ( any(mod(ntm,nsm) /= 0) ) then
        call die("initmesh: Read mesh sizes are not multiples of nsm")
      end if

!     Check that the numbers are divisible with the FFT algorithms
      do i = 1, 3
        itmp = ntm(i)
        call nfft( itmp )
        if ( itmp /= ntm(i) ) then
          call die("initmesh: Read mesh sizes are not compatible "
     &        // "with the FFT algorithm: % [2, 3, 5] == 0")
        end if
      end do

      end subroutine get_mesh_sizes

!--------------------------

      subroutine SameMeshAndAtoms(na, xa, ucell, rmax, G2max, G2mesh,
     &                            samesh, samexa)
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
      use alloc,         only : re_alloc, de_alloc
      implicit none
C     Passed arguments
      integer,    intent(in)  :: na
      real(dp),   intent(in)  :: G2max, G2mesh
      real(dp),   intent(in)  :: rmax
      real(dp),   intent(in)  :: ucell(3,3)
      real(dp),   intent(in)  :: xa(3,na)
      logical,    intent(out) :: samesh, samexa
C     Local variables
C     Saved
      integer,           save :: lastna     = 0
      real(dp),          save :: tiny       = 1.0e-12_dp,
     &                           lastc(3,3) = 0.743978657912656e50_dp,
     &                           lastra     = 0.0_dp
C     non-saved
      integer                 :: i, ia, j
      real(dp), pointer, save :: lastxa(:,:)

!------------------------------------------------------------------------- BEGIN
C     If number of atoms has changed then deallocate lastxa
      if (na.ne.lastna) then
        call re_alloc(lastxa, 1, 3, 1, na, 'lastxa', 'SameMeshAndAtoms')
        lastxa(1:3,1:na) = 0.0_dp
        lastxa(1,1) = 0.987654321273567e50_dp
      endif

C     Find if mesh has to be changed due to unit cell
      samesh = .true.
      do i = 1,3
        do j = 1,3
          if ( ucell(j,i) .ne. lastc(j,i) ) samesh = .false.
          lastc(j,i) = ucell(j,i)
        enddo
      enddo

C     Find if mesh has to be changed due to unit cell
      if ( G2max .gt. G2mesh * (1.0_dp + tiny) ) samesh = .false.

C     Find if mesh has to be changed due to rmax
      if (rmax .ne. lastra) samesh = .false.
      lastra = rmax

C     Find if atoms have moved having checked the number of atoms first
      samexa = (na.eq.lastna)
      if (samexa) then
        do ia = 1,na
          do i = 1,3
            if ( xa(i,ia) .ne. lastxa(i,ia) ) samexa = .false.
          enddo
        enddo
      endif

C     Copy the number of atoms and coordinates to save for next call
      lastna = na
      lastxa(1:3,1:na) = xa(1:3,1:na)

C     If cell has changed then it is necessary to reinitialise coordinates
      if (.not.samesh) samexa = .false.

!--------------------------------------------------------------------------- END
      end subroutine SameMeshAndAtoms

! subroutine InitAtomMesh
! ----------------------------------------------------------------------
! Initialises the atom information relating to the mesh for a given data
! distribution. It creates and initialise ipa. This arrays is allocated
! inside a meshDisType structure and can be accesed by pointers of the
! mesh module.
!
! Only necesary for those data distributions who needs to use the
! extended mesh.
!
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer  distr        : Used data distribution
! integer  na           : number of atoms
! integer  xa(3,na)     : Atomic coordinates
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! All output quantities are in the mesh and moreMeshSubs modules.
! The ipa array is stored in moreMeshSubs but accessed using the mesh
! module pointer.
!
!***********************************************************************
      subroutine InitAtomMesh( distr, na, xa )
      use precision, only: dp
      use mesh, only: ipa      ! atoms' mesh points within my Extended mesh
      use mesh, only: dxa      ! atom's position relative to mesh point
      use mesh, only: cmesh    ! Mesh cell vectors
      use mesh, only: meshLim  ! My processor's box of mesh points
      use mesh, only: ne       ! Points in rmax along each lat. vector
      use mesh, only: nmsc     ! Mesh points in each supercell vector
      use mesh, only: rcmesh   ! Reciprocal mesh cell vectors
      use mesh, only: iatfold  ! Supercell vector that keeps track of the
                               ! of the folding of the atomic coordinates
                               ! in the mesh
      use alloc,only: re_alloc ! (Re)allocation routine
      use moreMeshSubs, only : allocIpaDistr
      implicit none
!     Input parameters
      integer,       intent(in) :: na        ! Number of atoms in supercell
      real(dp),      intent(in) :: xa(3,na)  ! Atomic coordinates
      integer,       intent(in) :: distr     ! Used data distribution
!     Local variables
      character(len=*),parameter:: myName = 'InitAtomMesh '
      integer,             save :: lastna = 0  ! Number of atoms on last call
      integer                   :: ia, ixa(3), jsc(3), jxa(3), nem(3),
     &                             j1, j2, j3, myBox(2,3), myExtBox(2,3)
      integer                   :: ixabeffold(3)
      real(dp)                  :: cxa(3)

#ifdef DEBUG
      call write_debug( '      PRE InitAtomMesh' )
#endif
!     (Re)allocate atom-mesh arrays
      call allocIpaDistr( distr, na )
      if (na.ne.lastna) then
        call re_alloc( dxa, 1,3, 1,na, myName//'dxa' )
        call re_alloc( iatfold, 1,3, 1,na, myName//'iatfold' )
      endif
      lastna = na

      myBox(:,:) = meshLim(:,:) - 1

!     Add 'wings extensions' to myBox of mesh points, containing all points 
!     within the orbital spheres that may intersect myBox. Thus, wings must
!     be one diameter wide.
      myExtBox(1,:) = myBox(1,:) - 2*ne(:)
      myExtBox(2,:) = myBox(2,:) + 2*ne(:)

!     Find number of extended-box points.
      nem(:) = myExtBox(2,:) - myExtBox(1,:) + 1

!     Find atomic positions relative to mesh
      do ia = 1,na

!       Find atomic coordinates in mesh basis vectors
        cxa(:) = matmul( xa(:,ia), rcmesh(:,:) )
        ixabeffold(:) = floor(cxa(:))
        ! fix: Make sure we consider equivalent positions
        !      inside the cell, even if some atoms are outside.
        cxa(:) = modulo( cxa(:), real(nmsc(:),dp) )

!       Find indexes of mesh cell in which atom is
        ixa(:) = floor(cxa(:))
        iatfold(:,ia) = ( ixa(:) - ixabeffold(:) ) / nmsc(:)

!       Find atom position within mesh cell
        cxa(:) = cxa(:) - ixa(:)                  ! Mesh coordinates
        dxa(:,ia) = matmul( cmesh(:,:), cxa(:) )  ! Cartesian coords

!       Find atom's mesh index within myExtBox, but only if the
!       (periodic) atomic sphere (of width ne) can intersect myBox.
!       Otherwise, make ipa=0. nem is the width of myExtBox.
!       nmsc is the width of the supercell.
!       Notice that, even if several periodic images of the same 
!       atom may fit within myExtBox, only one must be considered,
!       since the points within its sphere, once folded into myBox,
!       will be equivalent to those of any other periodic image.
        ipa(ia) = 0
        sc: do j3 = -1,1   ! Loop on neighbor supercells
            do j2 = -1,1
            do j1 = -1,1
              jsc(:) = (/j1,j2,j3/)
              jxa(:) = ixa(:) + jsc(:)*nmsc(:)
              if ( all(jxa(:)>=myBox(1,:)-ne(:)) .and. 
     &             all(jxa(:)<=myBox(2,:)+ne(:)) ) then
                jxa(:) = jxa(:) - myExtBox(1,:)
                ipa(ia) = 1 + jxa(1) + nem(1)*jxa(2)
     &                      + nem(1)*nem(2)*jxa(3)

                iatfold(1,ia) = iatfold(1,ia) + j1
                iatfold(2,ia) = iatfold(2,ia) + j2
                iatfold(3,ia) = iatfold(3,ia) + j3
                
                exit sc
              end if
            end do    ! j1
          end do      ! j2
        end do sc     ! j3
      enddo ! ia

#ifdef DEBUG
      call write_debug( '      POS InitAtomMesh' )
#endif
      end subroutine InitAtomMesh

      subroutine PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, 
     &                              nsd, dvol, volume, Vscf, Vaux, 
     &                              fal, stressl, Forces, Stress )
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
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcore, chcore_sub
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C     Passed arguments
      integer,      intent(in)  :: na, isa(na), ntpl
      integer,      intent(in)  :: indxua(*)
      integer,      intent(in)  :: nsd
      real(dp),     intent(in)  :: dvol, volume
      real(grid_p), intent(in)  :: Vscf(ntpl,*)
      real(grid_p), intent(in)  :: Vaux(ntpl)
      logical,      intent(in)  :: Forces
      logical,      intent(in)  :: Stress
      real(grid_p), intent(out) :: rhopcc(ntpl)
      real(dp),   intent(inout) :: fal(3,*)
      real(dp), intent(inout)   :: stressl(3,3)
C     Internal variables
      integer                   :: i, ia, iop, ip, ip0, is, isp,
     &                             ispin, iua
      real(dp)                  :: dfa(3), dx(3), grrho(3), r, ra,
     &                             rhop, vxc
      real(dp),            save :: tiny
      logical                   :: Gradients
      data tiny / 1.0e-12_dp /

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE PartialCoreOnMesh' )
#endif
C     Find out whether this is a gradient run based on the arguments passed
      Gradients = ( Forces .or. Stress )

C     Initialise array to zero
      if (.not.Gradients) rhopcc(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        ra = rcore( is )
        if (ra .gt. tiny) then
C         Loop over mesh points inside rmax
          do iop = 1,mop
            ip0 = indexp( ipa(ia) + idop(iop) )
            if (ip0 .gt. 0) then
C             Loop over sub-points
              do isp = 1,nsp
                dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
                r = sqrt(dot_product(dx,dx))
C               if (r .lt. ra .and. r .gt. tiny) then
                if (r .lt. ra) then
                  ip = isp + nsp * (ip0 - 1)
                  call chcore_sub( is, dx, rhop, grrho )
                  if (Gradients) then
C                   Calculate gradients of PCC
                    do ispin = 1,nsd
                      vxc = Vscf(ip,ispin) - Vaux(ip)
                      do i = 1,3
                        dfa(i) = dvol * vxc * grrho(i) / nsd
                        if (Forces) fal(i,iua) = fal(i,iua) + dfa(i)
                        if (Stress) stressl(1:3,i) = stressl(1:3,i) + 
     &                    dx(1:3)*dfa(i)/volume
                      enddo
                    enddo
                  else
C                   Calculate density due to PCC
                    rhopcc(ip) = rhopcc(ip) + rhop
                  endif
                endif
              enddo
            endif
          enddo
        endif
      enddo
#ifdef DEBUG
      call write_debug( '      POS PartialCoreOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine PartialCoreOnMesh

      subroutine NeutralAtomOnMesh( na, isa, ntpl, vna,
     &                              indxua, dvol, volume,
     &                              drho, fal, stressl,
     &                              Forces, Stress )
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
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcut, phiatm
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C    Passed arguments
      integer,         intent(in) :: na, isa(na), ntpl
      integer,         intent(in) :: indxua(*)
      real(dp),        intent(in) :: dvol, volume
      real(grid_p),    intent(in) :: drho(*)
      logical,         intent(in) :: Forces, Stress
      real(grid_p), intent(inout) :: vna(ntpl)
      real(dp),     intent(inout) :: fal(3,*)
      real(dp),     intent(inout) :: stressl(3,3)
C     Internal variables
      integer                     :: i, ia, iop, ip, ip0, is, isp, iua
      real(dp)                    :: dx(3), grva(3), r, ra, va
      logical                     :: Gradients

#ifdef DEBUG
      call write_debug( '      PRE NeutralAtomOnMesh' )
#endif
C     Check whether forces and/or stress has been
C     requested based on arguments
      Gradients = (Forces.or.Stress)

C Initialise array to zero if not computing gradients
C If we are computing gradients, the current value of vna
C must be kept. This is particularly important for grid-cell sampling...

      if (.not.Gradients) vna(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        ra = rcut( is, 0 )

C       Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points
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
     &                dx(1:3) * grva(i) / volume
                  enddo
                else
                  vna(ip) = vna(ip) + va
                endif
              endif
            enddo
          endif
        enddo
      enddo
#ifdef DEBUG
      call write_debug( '      POS NeutralAtomOnMesh' )
#endif
!--------------------------------------------------------------------------- END
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
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcut, psch
      use atm_types, only: species
      use radial,    only: rad_func
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C     Passed arguments
      integer, intent(in)       :: na, isa(na), ntpl
      integer, intent(in)       :: indxua(*)
      real(grid_p), intent(out) :: Chlocal(ntpl)
C     Internal variables
      integer                   :: ia, iop, ip, ip0, is, isp, iua
      real(dp)                  :: dx(3), grpsch(3), r, ra, pschloc
      type(rad_func), pointer  :: func   
 
!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE LocalChargeOnMesh' )
#endif
C     Initialise array to zero 
      Chlocal(1:ntpl) = 0.0_grid_p
 
      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        func => species(is)%chlocal
        ra = func%cutoff
C       Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points
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
#ifdef DEBUG
      call write_debug( '      POS LocalChargeOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine LocalChargeOnMesh

      subroutine ModelCoreChargeOnMesh( na, isa, ntpl,
     $                                  ModelCore, indxua )
C
C     Constructs a model total core charge, useful for Bader analysis
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
      use atmfuncs, only: izofis, zvalfis
      use atm_types, only : species
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      use fdf
      use parallel, only: ionode

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

      real(dp) dx(3), ra, charge_ratio, r
      real(dp) core_charge
      real(dp), parameter :: pi = 3.141592653589793_dp

        ! We want a total core charge of (Z-Zval), localized in a
        ! region of around 1.0 bohr.
        ! Note that we put 1 electron (extra!) for Hydrogen atoms,
        ! to provide the necessary "cuspy" behavior for the reference
        ! file for Bader analysis. For H, the radius is 0.6 bohr
       real(dp) :: radius_standard, radius_h

        ! The model function is N*f(r/ra), where 
        !  f(x) = exp( - (sinh(2.0_dp*x)/sinh(1.0_dp))**2 )
        !  f(1) is zero to within 10^-5
        ! For normalization, we need \int_{0,1} { x*x*f(x) }, which is:
        real(dp), parameter :: Integral = 0.0411548_dp

        ! so that N above is Q/(4*pi*Integral*ra**3)
        
        ! Note that this is quite a small confinement region, so
        ! you might need a high cutoff (of the order of 300-500 Ry)
        ! to represent it on the grid. It might be better to use
        ! Denchar for this, rather than performing a Siesta calculation
        ! (even if it reads the converged DM)


      ! Initialise array to zero
      ModelCore(1:ntpl) = 0.0_grid_p
        
      radius_standard = fdf_get("bader-core-radius-standard",
     $                      1.0_dp,"Bohr")
      radius_h = fdf_get("bader-core-radius-hydrogen",
     $                      0.6_dp,"Bohr")
      
      if (ionode) then
         write(6,"(a,2f6.3)") "Bader Analysis core-charge setup. " //
     $                        "Radii (standard, H): ",
     $                        radius_standard, radius_h
      endif
      
      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersect my Box
        is = isa(ia)
        if (izofis(is) == 1) then
           ra = radius_h
           charge_ratio =  1.0_dp /
     $             (4.0_dp*pi*Integral*ra**3.0_dp)
        else
           ra = radius_standard
           charge_ratio =  (izofis(is)-zvalfis(is)) /
     $             (4.0_dp*pi*Integral*ra**3.0_dp)
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
                core_charge = charge_ratio*model(r/ra)
                ModelCore(ip) = ModelCore(ip) + core_charge
              endif
            enddo

          endif
        enddo

      enddo

      CONTAINS

      function model(x)
      real(dp), intent(in) :: x
      real(dp)             :: model
      
      model = exp( - (sinh(2.0_dp*x)/sinh(1.0_dp))**2 )
      end function model
      
      end subroutine ModelCoreChargeOnMesh
!
      subroutine PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi )
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nmpl          : Number of mesh points in unit cell locally
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom)
C integer isa(na)       : Species index of all atoms in supercell
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshphi
C ----------------------------------------------------------------------

C Modules
      use precision, only : dp, grid_p
      use atmfuncs,  only : rcut, phiatm
      use mesh,      only : dxa, idop, indexp, ipa, mop, nsp, xdop,
     &                      xdsp, meshLim
      use meshphi,   only : directphi, nphi, phi, lstpht, listp2, endpht
      use meshphi,   only : resetMeshPhi
      use fdf,       only : fdf_get
      use parallel,  only : IONode
      use alloc,     only : re_alloc

      implicit none
C Passed arguments
      integer,    intent(in) :: nmpl, norb, iaorb(norb),
     &                          iphorb(norb), isa(*)

      ! This array was set in distriPhiOnMesh, and holds
      ! the pre-computed count of orbitals touching a given
      ! point. It is re-used and overwritten here.
      integer, intent(inout) :: numphi(nmpl)

C     Local variables
      integer                :: ia, io, iop, ip, ip0, iphi, is, isp, n,
     &                          nlist, nliste, nmax

      logical                :: within
      logical,          save :: firsttime = .true.
      real(dp)               :: dxsp(3,nsp), grphi(3), phip, r2o,
     &                          r2sp(nsp)

!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE PhiOnMesh' )
#endif

      call timer( 'PHION', 1 )
C     On first call set the logical DirectPhi
      if (firsttime) then
        DirectPhi = fdf_get('DirectPhi', .false. )
        firsttime = .false.
      else
        ! reset the data structures in module meshphi
        call resetMeshPhi( )
      endif

C     Allocate memory related to nmpl
      nullify( endpht, phi, lstpht, listp2 )
      call re_alloc( endpht, 0, nmpl, 'endpht', 'PhiOnMesh' )

C     Initialise pointer array
      endpht(0) = 0
      do ip = 1,nmpl
        endpht(ip) = endpht(ip-1) + numphi(ip)
      enddo

C     Allocate phi if this is not a direct calculation
      nlist = endpht(nmpl)
      if (DirectPhi) then
        nphi = 1
      else
        nphi = nlist
        if (IONode) then
           write(6,"(a,i20)")
     $          "PhiOnMesh: Number of (b)points on node 0 = ", nmpl
           write(6,"(a,i20)")
     $          "PhiOnMesh: nlist on node 0 = ", nlist
        endif
      endif

C     Add an extra margin of error to nlist to minimise reallocations
      nliste = 1.01*nlist

C     Adjust dimensions of phi if necessary
      if (associated(phi)) then
        nmax = size(phi,2)
      else
        nmax = 0
      end if
      if (nphi.gt.nmax) then
        if (DirectPhi) then
          call re_alloc( phi, 1, nsp, 1, nphi, 'phi', 'PhiOnMesh' )
       else
          call re_alloc( phi, 1, nsp, 1, nliste, 'phi', 'PhiOnMesh' )
       endif
      endif

C     Adjust dimensions of list arrays if necessary
      if (associated(lstpht) .and. associated(listp2)) then
        nmax = min( size(lstpht), size(listp2) )
      else
        nmax = 0
      end if
      if (nlist.gt.nmax) then
        call re_alloc( lstpht, 1, nliste, 'lstpht', 'PhiOnMesh' )
        call re_alloc( listp2, 1, nliste, 'listp2', 'PhiOnMesh' )
      endif

C     Find indexes and values of atomic orbitals at mesh points
      numphi = 0
      do io = 1,norb
        ia = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersect my Box
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2
        do iop = 1, mop
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
     &                           phip, grphi )
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

      call timer( 'PHION', 2 )
#ifdef DEBUG
      call write_debug( '      POS PhiOnMesh' )
#endif
!------------------------------------------------------------------------ END
      end subroutine PhiOnMesh

      subroutine setupExtMesh( distr, rmax )
! ----------------------------------------------------------------------
! Setup the extended mesh variables for a given data distribution.
! It creates and initialise indexp, xdop and idop. These arrays are
! allocated inside a meshDisType structure and can be accesed by pointers
! of the mesh module.
!
! Only necesary for those data distributions who needs to use the
! extended mesh.
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer  distr        : distribution number
! real(dp) rmax         : Maximum orbital radius
!
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! All output quantities are in the mesh and moreMeshSubs modules.
! Data is stored in moreMeshSubs but accessed using the mesh module
! pointers.
!
!***********************************************************************
      use precision, only : dp
      use parallel,  only : IONode
      use mesh,      only : meshLim  ! My processor's box of mesh points
      use mesh,      only : ne       ! Points in rmax along each lat. vector
      use mesh,      only : nem      ! Extended-mesh divisions in each direction
      use mesh,      only : nmsc     ! Mesh points in each supercell vector
      use mesh,      only : mop      ! Accumulated num. of orbital points
      use mesh,      only : nsp      ! Number of sub-points of each mesh point
      use mesh,      only : cmesh    ! Mesh cell vectors
      use mesh,      only : xdsp     ! Vector to mesh sub-points
      use mesh,      only : indexp   ! ranslation from extended to
                                     ! normal mesh index
      use mesh,      only : idop     ! Mesh-index span of points 
                                     ! within an atomic sphere
      use mesh,      only : xdop     ! Coordinates of (super)points 
                                     ! within an atomic sphere
      use moreMeshSubs, only : allocExtMeshDistr

      implicit none
C     Input parameters
      integer              :: distr
      real(dp)             :: rmax
C     Local parameters
      integer              :: i1, i2, i3, j1, j2, j3, k1, k2, k3,
     &                        i, j, k, nep, boxWidth(3), extWidth(3),
     &                        isp, myBox(2,3), myExtBox(2,3)
      real(dp)             :: r, dxp(3), dx(3)
      logical              :: within
C     Functions
      real(dp)             :: dismin

      ! Correspondence between JMS and RG box descriptors
      myBox(:,:) = meshLim(:,:) - 1

! Add 'wings extensions' to myBox of mesh points, containing all points 
! within the orbital spheres that may intersect myBox. Thus, wings must
! be one diameter wide.
      myExtBox(1,:) = myBox(1,:) - 2*ne(:)
      myExtBox(2,:) = myBox(2,:) + 2*ne(:)

! Find number of extended-box points.
      nem(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      nep = nem(1) * nem(2) * nem(3)

      if (IONode) then
         write(6,'(a,3(i6,a),i12)') 'ExtMesh (bp) on 0 =',
     &        nem(1),' x',nem(2),' x',nem(3),' =', nep
      endif

!     Allocate memory for indexp(nep), idop(mop) and xdop(3,mop) of
!     the current data distribution
      call allocExtMeshDistr( distr, nep, mop )

! Find relationship between extended and normal box points
      boxWidth(:) =    myBox(2,:)    - myBox(1,:) + 1
      extWidth(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      do i3 = myExtBox(1,3),myExtBox(2,3)
        do i2 = myExtBox(1,2),myExtBox(2,2)
          do i1 = myExtBox(1,1),myExtBox(2,1)

            ! Find periodic indexes within first supercell
            j1 = modulo( i1, nmsc(1) )
            j2 = modulo( i2, nmsc(2) )
            j3 = modulo( i3, nmsc(3) )

            ! Find indexes relative to box origins
            j1 = j1 -    myBox(1,1)
            j2 = j2 -    myBox(1,2)
            j3 = j3 -    myBox(1,3)
            k1 = i1 - myExtBox(1,1)
            k2 = i2 - myExtBox(1,2)
            k3 = i3 - myExtBox(1,3)

            ! Find combined mesh indexes.
            j = 1 + j1 + boxWidth(1)*j2 + boxWidth(1)*boxWidth(2)*j3
            k = 1 + k1 + extWidth(1)*k2 + extWidth(1)*extWidth(2)*k3

            ! Find myExtBox -> myBox index translation
            if (j1>=0 .and. j1<boxWidth(1) .and.
     .          j2>=0 .and. j2<boxWidth(2) .and.
     .          j3>=0 .and. j3<boxWidth(3)) then
              indexp(k) = j   ! Point is within myBox
            else
              indexp(k) = 0   ! Point is outside myBox
            endif ! (j>=myBox(1) and j<=myBox(2))

          enddo ! i1
        enddo ! i2
      enddo ! i3

! This loop is done for each distribution, and replicates
! the one in InitMesh, except that now we actually fill
! the data in the arrays. Note that mop does not depend
! on the distribution.

! Find points within rmax (orbital points)
      mop = 0
      do i3 = -ne(3), ne(3)
        do i2 = -ne(2), ne(2)
          do i1 = -ne(1), ne(1)
            dxp(:) = cmesh(:,1) * i1 +
     &               cmesh(:,2) * i2 +
     &               cmesh(:,3) * i3 
            within = .false.
            do isp = 1,nsp
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
              mop = mop + 1
              ! Store index-distance and vector-distance to point.
              idop(mop)     = i1 + nem(1)*i2 + nem(1)*nem(2)*i3
              xdop(1:3,mop) = dxp(1:3)
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3

      end subroutine setupExtMesh

!------------------------------------------------------------------
      subroutine distriPhiOnMesh( nm, nmpl, norb, iaorb,
     &                            iphorb, isa, numphi,
     &                            need_gradients_in_xc )

C Computes the number of orbitals that intersects every mesh point
C and calls initMeshDistr in order to compute a new data distribution
C based in the number of orbital intersections.
C ==================================================================
C
C INPUT:
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer nmpl          : Local number of mesh points
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer isa(na)       : Species index of all atoms in supercell
C logical need_gradients_in_xc : 
C OUTPUT:
C The output values are in the module moreMeshSubs (New limits of the mesh
C and the information needed to transfer data between data
C distributions)
C
C BEHAVIOR:
C The amount of computation in every point of the mesh depends on the
C number of orbitals that intersecs in that mesh point. This can lead
C to unbalanced distributions in parallel executions. We want to
C compute a new data distribution, where the load is balanced among
C the several processes.
C
C In order to compute the number of orbitals that intersects in every
C mesh point, we loop in all orbitals of every atom and we increase the
C value of numphi in every point of the mesh that is inside of the
C radius of the current orbital.
C
C Finally, we call initMeshDistr that computes the new data distributions
C that will be used inside setup_hamiltonian. We are going to have 3
C different data distributions:
C  - Uniform:   Same amount of mesh in every process. Used in Poison
C  - Quadratic: The load in every point of the mesh is proporcional
C               to the square of the number of orbitals that intersects
C               in that point. Used in rhoofd and vmat.
C  - Lineal:    We try to distribute those points of the mesh that have
C               at least one orbital intersection. Used in cellxc.
C
C ==================================================================
      use mesh,         only: mop, ipa, idop, nsp, xdop, xdsp, dxa,
     &                        indexp, meshLim
      use precision,    only: dp, grid_p
      use atmfuncs,     only: rcut
      use parallel,     only: Nodes, node
#ifdef MPI
      use moreMeshSubs, only: initMeshExtencil
#endif
      use moreMeshSubs, only: initMeshDistr
      use moreMeshSubs, only: allocASynBuffer
      use moreMeshSubs, only: UNIFORM, QUADRATIC, LINEAR
      use alloc,        only: re_alloc, de_alloc
      
      implicit none
C     Passed arguments
      integer, intent(in)  :: nm(3), nmpl, norb, iaorb(norb),
     &                        iphorb(norb), isa(*)
      integer, intent(out) :: numphi(nmpl)
      logical, intent(in)  :: need_gradients_in_xc
C     Local variables
      integer              :: ii, io, ia, iphi, is, iop, ip0, isp
      real(dp)             :: r2o, dxsp(3,nsp), r2sp(nsp)
      logical              :: within
      integer,     pointer :: wload(:)
!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE distriPhiOnMesh' )
#endif
      call timer( 'REMESH', 1 )

      do ii= 1, nmpl
        numphi(ii) = 0
      enddo

C     Find number of atomic orbitals at mesh points
      do io = 1,norb
        ia   = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iphi = iphorb(io)
        is   = isa(ia)
        r2o  = rcut(is,iphi)**2
C       Loop over mesh points inside rmax
        do iop = 1, mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points to find if point is within range
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp)     = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo
C           If within range, add one to number of point orbitals
            if (within) numphi(ip0) = numphi(ip0) + 1
          endif
        enddo
      enddo

#ifdef MPI
      if (nodes.gt.1) then
        nullify( wload )
        call re_alloc( wload, 1, nmpl, 'wload', 'distriPhiOnMesh' )
        do ip0= 1, nmpl
          wload(ip0) = numphi(ip0)*(numphi(ip0)+1) + 1
        enddo
        call initMeshDistr( UNIFORM, QUADRATIC, nm, wload )

        do ip0= 1, nmpl
!         The computational work of a mesh point in cellxc is almost 0 when
!         when the charge density is zero. We use a proportion of 400/1.
          if (numphi(ip0).ne.0) then
            wload(ip0) = 400
          else
            wload(ip0) = 1
          endif
        enddo
        call initMeshDistr( UNIFORM, LINEAR, nm, wload )
        if (need_gradients_in_xc) call initMeshExtencil( LINEAR, nm )

C       Free local memory
        call de_alloc( wload, 'wload', 'distriPhiOnMesh' )

C       Allocate memory for asynchronous communications.
C       It does nothing if we use synchronous communications.
        call allocASynBuffer( LINEAR )
      endif

#endif
      call timer( 'REMESH', 2 )
#ifdef DEBUG
      call write_debug( '      POS distriPhiOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine distriPhiOnMesh

      end module meshsubs
