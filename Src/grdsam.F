! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_grdsam
      public :: grdsam
      CONTAINS
      subroutine grdsam(nspin, norb, iaorb, iphorb, 
     &                  nuo, nuotot, nua, na, isa, xa, indxua,
     &                  cell, ntm, ifa, istr, maxnd,
     &                  maxnh, numh, listhptr, listh, Dscf, Datm, Hmat,
     &                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &                  Exc, Dxc, dipol, stress, fal, stressl )

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
c integer ifa           : scf contrib to forces calculated or not
c integer istr          : scf contrib to stress calculated or not
c integer maxnd         : first dimension of Dscf
c integer maxnh         : first dimension of listh, hmat
c real*8  Dscf(maxnd,h_spin_dim): scf DM elements
c real*8  datm(norb)    : Harris DM diagonal elements
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON WHETHER MESH IS CALCULATED)*******
c integer ntm(3)        : number of mesh divisions of each cell vector
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON ARGUMENT ILH) ********************
c integer numh(nuo)     : number of nonzero H elements for each row
c integer listhptr(nuo) : pointer to each row (-1) of H
c integer listh(maxnh)  : nonzero-H-element column indexes for each row
c real*8  Hmat(maxnh,h_spin_dim): Hamiltonian matrix in sparse form
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
c real*8  stress(3,3) : stress tensor
c real*8  fal(3,na)   : atomic forces (local to Node)
c real*8  stressl(3,3): stress tensor (local to Node)
c ***************************************************************************

C
C  Modules
C
      use precision
      use fdf
      use files,     only : filesOut_t    ! derived type for output file names
      use parallel,  only : IONode
      use sys,       only : die
      use units, only : eV, Ang, kbar
      use m_dhscf,   only : dhscf, dhscf_init
      use alloc,     only : re_alloc, de_alloc
      use siesta_options, only: g2cut
      use siesta_geom,    only: mscell
      use siesta_options, only: hirshpop, voropop
      use siesta_options, only: partial_charges_at_every_geometry
      use m_partial_charges, only: want_partial_charges
      
      use m_spin, only: h_spin_dim

#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      
      implicit          none

      integer           maxnh, maxnd, na, norb, nspin,
     &                  nuo, nuotot, nua, iaorb(norb), ifa, indxua(na), 
     &                  iphorb(norb), isa(na), istr, listh(maxnh),
     &                  listhptr(nuo), numh(nuo), ntm(3)

      real(dp)          cell(3,3), Datm(norb), dipol(3),
     &                  Dscf(:,:),
     &                  Hmat(:,:), stress(3,3), stressl(3,3),
     &                  xa(3,na), fal(3,nua)

      ! In/out energies
      real(dp) :: DUscf, DUext, Dxc, Enaatm, Enascf, Exc, Uatm, Uscf

C Internal variables and arrays
 
      logical           samesh
      integer           nr, ipt, ia,
     &                  iv
      real(dp) :: strold(3,3,2)
      ! Average quantities
      real(dp) :: avdipol(3)
      real(dp) :: avDUscf, avDUext, avDxc
      real(dp) :: avEnaatm, avEnascf, avExc
      real(dp) :: avUatm, avUscf
      real(dp) :: anpt

      ! Used grid-cutoff
      real(dp) :: G2max

      real(dp), pointer :: allfa(:,:,:) => null()
      real(dp), pointer :: allstr(:,:,:,:) => null()
      real(dp), pointer :: faold(:,:) => null()
      real(dp), pointer :: xanew(:,:) => null()

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline
      type(filesOut_t)           :: filesOut  ! blank output file names

      ! Saved variables:
      integer, save  :: npt = -1
      real(dp), pointer, save :: pt(:,:) => null()

      real(dp) :: d, dx(3)
      real(dp) :: famax(3), fastd(3)
      real(dp) :: strmax
      integer :: nABC(3), ix, iy, iz

c ----------------------------------------------------------------------------
#ifdef DEBUG
      call write_debug( '    PRE grdsam' )
#endif

      if ( size(Dscf,2) /= h_spin_dim ) then
         call die('Spin components is not equal to options.')
      end if
      if ( size(Hmat,2) /= h_spin_dim ) then
         call die('Spin components is not equal to options.')
      end if

C Check whether GridCellSampling block has been looked for ------------------

      if ( npt == -1 ) then

        ! Initialize number of sampled points...
        npt = 0

        ! The Grid.CellSampling may be given in a list
        ! fashion such that one can do a "MP" style
        ! sampling.

        if ( fdf_islist('Grid.CellSampling') ) then
           
           ! There are two methods of using the 
           ! Cell sampling, the diagonal one (via a list)

           call fdf_list('Grid.CellSampling', 3, nABC)
           ! all zeros are equivalent to 1
           if ( nABC(1) < 1 ) nABC(1) = 1
           if ( nABC(2) < 1 ) nABC(2) = 1
           if ( nABC(3) < 1 ) nABC(3) = 1

           ! Get total number of points
           npt = product(nABC) - 1 ! we do not have 0. 0. 0.

           if ( npt == 0 ) then
              ! This is the special case where the user
              ! says [1 1 1] / [0 0 0]
              ! in that case we sample the half diagonal voxel
              npt = 1
           end if

           if ( IONode ) then
            write(6,'(/a)') 'grdsam: Reading Grid.CellSampling list'
           end if

           ! Allocate the sampled points
           call re_alloc( pt, 1, 3, 1, npt, 'pt', 'grdsam' )
           
           ipt = 0
           do iz = 0 , nABC(3) - 1
            do iy = 0 , nABC(2) - 1
              do ix = 0 , nABC(1) - 1
                  
              ! This is 0. 0. 0., we should skip it
              if ( ix + iy + iz == 0 ) cycle

              ipt = ipt + 1
              pt(1,ipt) = real(ix,dp) / real(nABC(1),dp)
              pt(2,ipt) = real(iy,dp) / real(nABC(2),dp)
              pt(3,ipt) = real(iz,dp) / real(nABC(3),dp)

             end do
            end do
           end do

           if ( npt == 1 ) then
              ipt = ipt + 1
              pt(:,1) = 0.5_dp
           end if

           if ( ipt /= npt ) then
              call die('Error in creating the equi-grid-cell-sampling.')
           end if
           
        else if ( fdf_block('Grid.CellSampling', bfdf) ) then
           
           ! The user-defined one
           
           ! First read the number of sampled points
           npt = 0
           do while ( fdf_bline(bfdf,pline) )
              nr = fdf_bnreals(pline)
              if ( nr >= 3 ) then
                 npt = npt + 1
              end if
           end do

           ! Rewind to re-read the points
           call fdf_brewind(bfdf)

           ! Allocate arrays
           call re_alloc( pt, 1, 3, 1, npt, 'pt', 'grdsam' )
           
          if ( IONode ) then
            write(6,'(/a)') 'grdsam: Reading %block Grid.CellSampling'
          end if
          ipt = 0
          do while ( fdf_bline(bfdf,pline) )
             nr = fdf_bnreals(pline)
             if (nr .ge. 3) then
                ipt = ipt + 1
                pt(1,ipt) = fdf_breals(pline,1)
                pt(2,ipt) = fdf_breals(pline,2)
                pt(3,ipt) = fdf_breals(pline,3)
             else
                write(6,'(a,/,a)')
     &             'grdsam: syntax ERROR in %block Grid.CellSampling:',
     &             trim(fdf_getline(bfdf%mark))
             end if
          end do

          if ( ipt /= npt ) then
            call die('Error in reading the grid-cell-sampling.')
          end if

        end if
       
      end if

c if sampling, store the forces and stresses prior to dhscf -----------------

      if ( npt > 0 ) then

        ! Allocate the old forces which contain the pre-mesh local forces
        call re_alloc( faold, 1, 3, 1, nua, 'faold', 'grdsam' )
        faold(:,:) = fal(:,:)
        ! also store the pre-mesh global and local stress
        strold(:,:,1) = stress(:,:)
        strold(:,:,2) = stressl(:,:)
        if ( IONode ) then
          write(*,'(2(a,i0))') 'grdsam: Grid-cell sampling: ', 0,
     &        ' / ',npt
        end if
       
      end if


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
     &            na, isa, xa, indxua, ntm,
     &            ifa, istr, 1, filesOut,
     &            maxnd, numh, listhptr, listh, Dscf, Datm,
     &            maxnh, Hmat, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &            Exc, Dxc, dipol, stress, fal, stressl)
      want_partial_charges = .false.

c Sampling ------------------------------------------------------------------

      if ( npt > 0 ) then

c initialize averages with output of regular DHSCF run -----------------------

        ! The average part...
        anpt = dble(npt+1)

        ! Allocate the used sampling arrays
        call re_alloc( xanew, 1, 3, 1, na, 'xanew', 'grdsam' )
        call re_alloc( allfa, 1, 3, 1, nua, 0, npt, 'allfa', 'grdsam')
        call re_alloc( allstr, 1, 3, 1, 3, 1, 2, 0, npt, 
     &       'allstr', 'grdsam')

        ! It is rather instructive to save them all 
        ! and get the difference between the consecutive
        ! runs to get a magnitude of the difference.
        ! I.e. the stored forces are now _only_ the grid-cell
        ! sampled differences
        allfa(:,:,0) = fal(:,:) - faold(:,:)
        allstr(:,:,1,0) = stress(:,:) - strold(:,:,1)
        allstr(:,:,2,0) = stressl(:,:) - strold(:,:,2)

        ! We do not print out anything related to the
        ! energies and dipoles...
        avdipol(:) = dipol(:) / anpt
        avDUscf = DUscf / anpt
        avDUext = DUext / anpt
        avDxc = Dxc / anpt
        avEnaatm = Enaatm / anpt
        avEnascf = Enascf / anpt
        avExc = Exc / anpt
        avUatm = Uatm / anpt
        avUscf = Uscf / anpt

c loop sampling on displacements making averages ----------------------------

        do ipt = 1, npt

          ! Calculate the displacement for this mesh sample
          do ix = 1, 3
            dx(ix) = 0._dp
            do iy = 1, 3
              dx(ix) = dx(ix) + cell(ix,iy)*pt(iy,ipt)/dble(ntm(iy))
            end do
          end do

          do ia = 1, na
            xanew(:,ia) = xa(:,ia) + dx(:)
          end do

          ! Copy back the initial (pre-mesh) forces so we may
          ! correctly calculate the averaged forces
          fal(:,:) = faold(:,:)
          stress(:,:) = strold(:,:,1)
          stressl(:,:) = strold(:,:,2)
        
          if ( IONode ) then
           write(*,'(2(a,i0))') 'grdsam: Grid-cell sampling: ',ipt,
     &          ' / ',npt
          end if

          ! We have to call dhscf_init, among other things, to
          ! initialize the mesh structures
          G2max = g2cut
          call dhscf_init(nspin, norb, iaorb, iphorb,
     &                      nuo, nuotot, nua, na,
     &                      isa, xanew, indxua, cell,
     &                      mscell, G2max, ntm,
     &                      maxnd, numh, listhptr, listh, datm,
     &                      Fal, stressl)

          call dhscf( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     &                na, isa, xanew, indxua, ntm,
     &                ifa, istr, 0, filesOut,
     &                maxnd, numh, listhptr, listh, Dscf, Datm,
     &                maxnh, Hmat, Enaatm, Enascf, Uatm, Uscf, DUscf,
     &                DUext, Exc, Dxc, dipol, stress, fal, stressl)

          ! Collect the differences
          allfa(:,:,ipt) = fal(:,:) - faold(:,:)
          allstr(:,:,1,ipt) = stress(:,:) - strold(:,:,1)
          allstr(:,:,2,ipt) = stressl(:,:) - strold(:,:,2)

          ! Perform the average of the energies
          avdipol(:) = avdipol(:) + dipol(:) / anpt
          avDUscf = avDUscf + DUscf / anpt
          avDUext = avDUext + DUext / anpt
          avDxc = avDxc + Dxc / anpt
          avEnaatm = avEnaatm + Enaatm / anpt
          avEnascf = avEnascf + Enascf / anpt
          avExc = avExc + Exc / anpt
          avUatm = avUatm + Uatm / anpt
          avUscf = avUscf + Uscf / anpt

        end do

        ! Calculate the averaged forces

        ! Reduce the forces, local-node forces
        ! and stress so that we only compare a single part of the
        ! forces.
        fal(:,:) = allfa(:,:,0) / anpt
        stress(:,:) = allstr(:,:,1,0) / anpt
        stressl(:,:) = allstr(:,:,2,0) / anpt
        do ipt = 1 , npt
           ! node-local forces
           fal(:,:) = fal(:,:) + allfa(:,:,ipt) / anpt
           ! global stress
           stress(:,:) = stress(:,:) + allstr(:,:,1,ipt) / anpt
           ! local stress
           stressl(:,:) = stressl(:,:) + allstr(:,:,2,ipt) / anpt
        end do

        ! Add the original forces...
        fal(:,:) = fal(:,:) + faold(:,:)
        stress(:,:) = stress(:,:) + strold(:,:,1)
        stressl(:,:) = stressl(:,:) + strold(:,:,2)
        ! In fact we are now done... But we calculate
        ! the differences for informing the user.


        ! Now perform the averages to be able to compare
        ! differences.
#ifdef MPI
        do ipt = 0 , npt
           call globalize_sum(allfa(:,:,ipt), faold(:,:))
           ! After this, allfa(:,:,ipt) are the globalized
           ! mesh-dependent forces for the ipt mesh-sample
           allfa(:,:,ipt) = faold(:,:)
           ! After this, allstr(:,:,1,ipt) is the globalized
           ! mesh-dependent stress for the ipt mesh-sample
           call globalize_sum(allstr(:,:,2,ipt), strold(:,:,1))
           allstr(:,:,1,ipt) = allstr(:,:,1,ipt) + strold(:,:,1)
        end do
#else
        do ipt = 0 , npt
           ! After this, allstr(:,:,1,ipt) is the globalized
           ! mesh-dependent stress for the ipt mesh-sample
           allstr(:,:,1,ipt) = allstr(:,:,1,ipt) + allstr(:,:,2,ipt)
        end do
#endif

        ! First we need to calculate the average sampled forces
        faold(:,:) = allfa(:,:,0) / anpt
        strold(:,:,1) = allstr(:,:,1,0) / anpt
        do ipt = 1 , npt
           faold(:,:) = faold(:,:) + allfa(:,:,ipt) / anpt
           strold(:,:,1) = strold(:,:,1) + allstr(:,:,1,ipt) / anpt
        end do

        ! Calculate the maximum relative change in
        ! the sampling
        famax = 0._dp
        fastd(:) = 0._dp
        strmax = 0._dp
        do ipt = 0 , npt
         do ia = 1 , nua
          do ix = 1, 3
           d = allfa(ix,ia,ipt) - faold(ix,ia)
           if ( abs(d) > abs(famax(ix)) ) then
             ! updated maximum difference from the average
             famax(ix) = d
           end if
           ! update standard deviation calculation
           fastd(ix) = fastd(ix) + d ** 2 / (anpt * nua)
          end do
         end do
         do iv = 1 , 3
          do ix = 1, 3
           d = allstr(ix,iv,1,ipt) - strold(ix,iv,1)
           if ( abs(d) > abs(strmax) ) then
              strmax = d
           end if
          end do
         end do
        end do
        ! Calculate std. 
        fastd = sqrt(fastd)
        
        ! Print out stuff
        if ( IONode ) then
           write(*,'(/a,3(tr1,f10.6),tr1,a)')
     &      'grdsam: Max (sampled force - avg. force):',famax/eV*Ang,
     &      'eV/Ang'
           write(*,'(a,3(tr1,f10.6),tr1,a)')
     &      'grdsam: Standard deviation of sampled force:',fastd/eV*Ang,
     &      'eV/Ang'
           write(*,'(a,tr1,f10.6,tr1,a/)')
     &      'grdsam: Max (sampled stress - avg. stress):',strmax/kbar,
     &      'kbar'
        end if
      
        ! Copy back the resulting averaged energies
        ! and dipole moments
        dipol(:) = avdipol(:)
        DUscf = avDUscf
        DUext = avDUext
        Dxc = avDxc
        Enaatm = avEnaatm
        Enascf = avEnascf
        Exc = avExc
        Uatm = avUatm
        Uscf = avUscf

        call de_alloc( xanew, 'xanew', 'grdsam' )
        call de_alloc( allfa, 'allfa', 'grdsam' )
        call de_alloc( allstr, 'allstr', 'grdsam' )

        ! Since we have displaced the atoms we need to recreate
        ! the lists such that later analysis is matching the coordinates
        G2max = g2cut
        call dhscf_init(nspin, norb, iaorb, iphorb,
     &      nuo, nuotot, nua, na,
     &      isa, xa, indxua, cell,
     &      mscell, G2max, ntm,
     &      maxnd, numh, listhptr, listh, datm,
     &      Faold, strold(:,:,1))

        ! Deallocate local memory
        call de_alloc( faold, 'faold', 'grdsam' )

      end if

#ifdef DEBUG
      call write_debug( '    POS grdsam' )
#endif

      end subroutine grdsam

      end module m_grdsam
