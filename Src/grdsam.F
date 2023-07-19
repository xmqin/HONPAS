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
      subroutine grdsam(nspin, norb, iaorb, iphorb, 
     .                  nuo, nuotot, nua, na, isa, xa, indxua,
     .                  cell, mscell, g2max, ntm, ifa, istr, maxnd,
     .                  maxnh, numh, listhptr, listh, Dscf, Datm, Hmat,
     .                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                  Exc, Dxc, dipol, fa, stress, fal, stressl )

c ***************************************************************************
c Final call to dhscf of the scf cycle + grid-cell sampling
c
c   After a first call to DHSCF, there are more calls for rigidly shifted
c   atomic coordinates, sampling the small cell defined by the grid, to get
c   average energy, forces, stress, and dipole. It does it with fix density 
c   matrix. It can be regarded as a (discrete sampling) symmetrization
c   to restore the homogeneity of space, which was lost with the 
c   grid summations. 
c
c   The first call to DHSCF is the usual one (the point in which 
c   selfconsistency has been performed). It corresponds to the 
c   shift (0,0,0) in the grid cell.
c
c   The sampling is done on the additional points given in block
c   GridCellSampling, which are given in fractional coordinates of
c   the small grid cell. 
c
c   Even if the sum is large the symmetrization is never complete 
c   since the density matrix is the one obtained for one position of
c   the grid (no new selfconsistency for each shift).
c      
c Written by E. Artacho. January 1998. 
c *********************** INPUT TOWARDS DHSCF *******************************
c integer nspin         : number of spins considered 
c integer norb          : total number of basis orbitals
c integer iaorb(norb)   : atom to which each orbital belongs
c integer iphorb(norb)  : orbital index (within atom) of each orbital
c integer nuo           : number of orbitals in a unit cell (local)
c integer nuotot        : number of orbitals in a unit cell (global)
c integer nua           : number of atoms in unit cell
c integer na            : number of atoms in supercell
c integer isa(na)       : species indexes
c real*8  xa(3,na)      : atomic positions
c real*8  cell(3,3)     : unit cell vectors
c integer mscell(3,3)   : supercell vectors in units of ucell
c integer ifa           : scf contrib to forces calculated or not
c integer istr          : scf contrib to stress calculated or not
c integer maxnd         : first dimension of dscf
c integer maxnh         : first dimension of listh, hmat
c real*8  dscf(maxnd,nspin): scf DM elements
c real*8  datm(norb)    : Harris DM diagonal elements
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON WHETHER MESH IS CALCULATED)*******
c integer ntm(3)        : number of mesh divisions of each cell vector
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON ARGUMENT ILH) ********************
c integer numh(nuo)     : number of nonzero H elements for each row
c integer listhptr(nuo) : pointer to each row (-1) of H
c integer listh(maxnh)  : nonzero-H-element column indexes for each row
c real*8  Hmat(maxnh,nspin): Hamiltonian matrix in sparse form
c ************************* DHSCF OUTPUT ******i*****************************
c real*8  enaatm : integral of vna * rhoatm
c real*8  enascf : integral of vna * rhoscf
c real*8  uatm   : Harris hartree electron-interaction energy
c real*8  uscf   : scf hartree electron-interaction energy
c real*8  duscf  : electrostatic (hartree) energy of rhoscf-rhoatm density
c real*8  duext  : interaction energy with external electric field
c real*8  exc    : scf XC energy
c real*8  dxc    : scf double-counting correction to exc
c real*8  dipol(3): electric dipole (in a.u.)
c *********************** DHSCF INPUT AND OUTPUT ****************************
c real*8  g2max       : effective planewave cutoff in Ry 
c real*8  fa(3,na)    : atomic forces
c real*8  stress(3,3) : stress tensor
c real*8  fal(3,na)   : atomic forces (local to Node)
c real*8  stressl(3,3): stress tensor (local to Node)
c ***************************************************************************

C
C  Modules
C
      use precision
      use fdf
      use parsing
      use parallel,  only : Node
      use sys,       only : die
      use siesta_options, only: hirshpop, voropop
      use siesta_options, only: partial_charges_at_every_geometry
      use m_partial_charges, only: want_partial_charges

#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           maxnh, maxnd, na, norb, nspin, nuo, nuotot, 
     .                  nua, iaorb(norb), ifa, indxua(na), 
     .                  iphorb(norb), isa(na), istr, listh(maxnh),
     .                  listhptr(nuo), mscell(3,3), 
     .                  numh(nuo), ntm(3)

      real(dp)          cell(3,3), datm(norb), dipol(3), 
     .                  dscf(maxnd,nspin), DUscf, DUext, Dxc, 
     .                  enaatm, enascf, exc, fa(3,nua), g2max,
     .                  Hmat(maxnh,nspin), stress(3,3),
     .                  uatm, uscf, xa(3,na), fal(3,nua), stressl(3,3)

      external          memory

C Internal variables and arrays
 
      character         line*130, names*80
      logical           samesh
      integer           ni, nn, nr, nv, ipt, maxpt, iu, ia,
     .                  ix, iv, integs(4), lastc, lc(0:3)
#ifdef MPI
      integer           MPIerror
#endif
      real(dp)          reals(4), values(4), avexc,
     .                  strold(3,3), avdipo(3), avstre(3,3), avdxc,
     .                  avenaa, avenas, avuatm, avuscf, avdusc, anpt,
     .                  strold2(3,3), avstre2(3,3)

      parameter         (maxpt = 100)

      real(dp), dimension(:,:), allocatable, save ::
     .                  pt, avfa, faold, xanew, avfa2, faold2
      real(dp), dimension(:,:), allocatable, save ::
     .                  dpt


      !Saved variables:
      logical, save :: blread = .false.
      logical, save :: samgrd = .false.

      real(dp), save :: lstcll(3,3) = 0.743978657912656D50
      integer, save  :: lstntm(3) = (/10000, 1, 1/) 
      integer, save  :: npt

c ----------------------------------------------------------------------------

C Allocate local memory
      if (.not.allocated(dpt)) then
        allocate(pt(3,maxpt+1))
        call memory('A','D',3*(maxpt+1),'grdsam')
        allocate(dpt(3,maxpt+1))
        call memory('A','D',3*(maxpt+1),'grdsam')
      endif
      allocate(avfa(3,na))
      call memory('A','D',3*na,'grdsam')
      allocate(avfa2(3,na))
      call memory('A','D',3*na,'grdsam')
      allocate(faold(3,na))
      call memory('A','D',3*na,'grdsam')
      allocate(faold2(3,na))
      call memory('A','D',3*na,'grdsam')
      allocate(xanew(3,na))
      call memory('A','D',3*na,'grdsam')

C Check whether GridCellSampling block has been looked for ------------------

      if ( .not. blread ) then

C Look for block and read it if found ---------------------------------------
 
        if (Node.eq.0) then
          samgrd = fdf_block('GridCellSampling',iu)
        endif
#ifdef MPI
        call MPI_Bcast(samgrd,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

        if (samgrd) then
          if (Node.eq.0) then
            write(6,'(a)') 'grdsam: Reading %block GridCellSampling'
          endif
          do ipt = 1, maxpt + 1
            if (Node.eq.0) then
              read(iu,'(a)', end=50) line
            endif
#ifdef MPI
            call MPI_Bcast(line,130,MPI_character,0,MPI_Comm_World,
     .        MPIerror)
#endif
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (nr .ge. 3) then
              pt(1,ipt) = reals(1)
              pt(2,ipt) = reals(2)
              pt(3,ipt) = reals(3)
            else
              npt = ipt - 1
              goto 50
            endif
          enddo
          call die('grdsam: ERROR: Number of points' //
     $             ' in block larger than MAXPT')
   50     continue
        endif
        blread = .true.
      endif

c if sampling, store the forces and stresses prior to dhscf -----------------

      if (samgrd) then
        faold(1:3,1:nua) = fa(1:3,1:nua)
        faold2(1:3,1:nua) = fal(1:3,1:nua)
        strold(1:3,1:3) = stress(1:3,1:3)
        strold2(1:3,1:3) = stressl(1:3,1:3)
        if (Node.eq.0) then
          write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', 0
        endif
      endif


c first (normal) call to DHSCF, ihmat=1 --------------------------------------

!
!     At this point we have finished the scf loop for this geometry,
!     and can compute the partial charges if desired in this call to
!     dhscf.
!
      if ((hirshpop .or. voropop)
     $    .and. partial_charges_at_every_geometry) then
         want_partial_charges = .true.
      endif
      call dhscf( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     .            na, isa, xa, indxua, cell, mscell, g2max, ntm,
     .            ifa, istr, 1, ' ', ' ', ' ', ' ', ' ', ' ',
     .            maxnd, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, numh, listhptr, listh, Hmat,  
     .            Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, fa, stress, fal, stressl)
      want_partial_charges = .false.

c Sampling ------------------------------------------------------------------

      if (samgrd) then

c find if mesh is new -------------------------------------------------------

        samesh = .true.
        do iv = 1, 3
          if ( ntm(iv) .ne. lstntm(iv) ) samesh = .false.
          lstntm(iv) = ntm(iv)
          do ix = 1, 3
            if ( cell(ix,iv) .ne. lstcll(ix,iv) ) samesh = .false.
            lstcll(ix,iv) = cell(ix,iv)
          enddo
        enddo

c generate displacements (in Bohr cartesian) if first time or new mesh ------

        if (.not. samesh) then
          if (Node.eq.0) then
            write(6,'(a)') 
     .       'grdsam: Generating displacements for grid-cell sampling'
          endif
          do ipt = 1, npt
            do ix = 1, 3
              dpt(ix,ipt) = 0.0d0
              do iv = 1, 3
                dpt(ix,ipt) = dpt(ix,ipt) + 
     .                          cell(ix,iv)*pt(iv,ipt)/dble(ntm(iv))
              enddo
            enddo
          enddo
        endif

c initialize averages with output of regular DHSCF run -----------------------

        avfa(1:3,1:nua) = fa(1:3,1:nua)
        avfa2(1:3,1:nua) = fal(1:3,1:nua)
        avstre(1:3,1:3) = stress(1:3,1:3)
        avstre2(1:3,1:3) = stressl(1:3,1:3)
        avdipo(1:3) = dipol(1:3)
        avenaa = enaatm
        avenas = enascf
        avuatm = uatm
        avuscf = uscf
        avdusc = duscf
        avexc  = exc
        avdxc  = dxc

c loop sampling on displacements making averages ----------------------------

        do ipt = 1, npt

          do ia = 1, na
            do ix = 1, 3
              xanew(ix,ia) = xa(ix,ia) + dpt(ix,ipt)
            enddo
          enddo
          fa(1:3,1:nua) = faold(1:3,1:nua)
          fal(1:3,1:nua) = faold2(1:3,1:nua)
          stress(1:3,1:3) = strold(1:3,1:3)
          stressl(1:3,1:3) = strold2(1:3,1:3)
        
          if (Node.eq.0) then
          write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', ipt
          endif

          call dhscf( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     .                na, isa, xanew, indxua, cell, mscell, g2max, ntm,
     .                ifa, istr, 0, ' ', ' ', ' ', ' ', ' ', ' ',
     .                maxnd, numh, listhptr, listh, Dscf, Datm,
     .                maxnh, numh, listhptr, listh, Hmat, 
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, fal, stressl)

          avfa(1:3,1:nua) = avfa(1:3,1:nua) + fa(1:3,1:nua)
          avfa2(1:3,1:nua) = avfa2(1:3,1:nua) + fal(1:3,1:nua)
          avstre(1:3,1:3) = avstre(1:3,1:3) + stress(1:3,1:3)
          avstre2(1:3,1:3) = avstre2(1:3,1:3) + stressl(1:3,1:3)
          avdipo(1:3) = avdipo(1:3) + dipol(1:3)
          avenaa = avenaa + enaatm
          avenas = avenas + enascf
          avuatm = avuatm + uatm
          avuscf = avuscf + uscf
          avdusc = avdusc + duscf
          avexc  = avexc  + exc
          avdxc  = avdxc  + dxc
        enddo

c final averages ------------------------------------------------------------

        anpt = dble(npt+1)

        fa(1:3,1:nua) = avfa(1:3,1:nua)/anpt
        fal(1:3,1:nua) = avfa2(1:3,1:nua)/anpt
        stress(1:3,1:3) = avstre(1:3,1:3)/anpt
        stressl(1:3,1:3) = avstre2(1:3,1:3)/anpt
        dipol(1:3) = avdipo(1:3)/anpt
        enaatm = avenaa/anpt
        enascf = avenas/anpt
        uatm   = avuatm/anpt
        uscf   = avuscf/anpt
        duscf  = avdusc/anpt
        exc    = avexc/anpt
        dxc    = avdxc/anpt

      endif

C Deallocate local memory
      call memory('D','D',size(xanew),'grdsam')
      deallocate(xanew)
      call memory('D','D',size(faold2),'grdsam')
      deallocate(faold2)
      call memory('D','D',size(faold),'grdsam')
      deallocate(faold)
      call memory('D','D',size(avfa2),'grdsam')
      deallocate(avfa2)
      call memory('D','D',size(avfa),'grdsam')
      deallocate(avfa)

      return
      end

