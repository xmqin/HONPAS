! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
  subroutine overkkneig( kvector, bvector, nuo, nuotot, psiatk, psinei,   &
 &                       maxnh, delkmat, nbandsocc_loc, nbandsocc, Mkb )

   
  use precision,            only: dp             ! Real double precision type
  use parallel,             only: Nodes          ! Total number of Nodes
  use parallel,             only: Node           ! Local Node
  use parallel,             only: IONode         ! Input/output node
  use atomlist,             only: iaorb          ! Pointer to atom to which
                                                 !   orbital belongs
  use atomlist,             only: indxuo         ! Index of equivalent orbital
                                                 !   in unit cell
  use sparse_matrices,      only: numh           ! Number of nonzero element of
                                                 !   each row of the 
                                                 !   hamiltonian matrix 
  use sparse_matrices,      only: listh          ! Nonzero hamiltonian-matrix 
                                                 !   elements
  use sparse_matrices,      only: listhptr       ! Pointer to start of each row 
                                                 !   of the hamiltonian matrix
  use sparse_matrices,      only: xijo           ! Vectors between orbital
                                                 !   centers (sparse)
  use siesta_geom,          only: xa             ! Atomic positions
  use alloc,                only: re_alloc       ! Allocatation routines
  use alloc,                only: de_alloc       ! Deallocatation routines

#ifdef MPI
  use parallelsubs,         only: LocalToGlobalOrb ! Converts an orbital index
                                                 !   in the local frame 
                                                 !   to the global frame
  use parallelsubs,         only: GetNodeOrbs    ! Calculates the number of
                                                 !   orbitals stored on the 
                                                 !   local Node.
  use m_orderbands,       only: which_band_in_node  ! Given a node and a 
                                                    !   local index,
                                                    !   this array gives the
                                                    !   global index of the band
                                                    !   stored there
  use m_orderbands,       only: sequential_index_included_bands
                                                    ! Sequential number of the
                                                    !   bands included for
                                                    !   wannierization
                                                    !   (the bands are listed
                                                    !   in order of incremental
                                                    !   energy)

  use mpi_siesta
#endif

  implicit none

  real(dp),    intent(in)  :: kvector(3)         ! Wave vector for which the 
                                                 !   Overlap matrix between 
                                                 !   the periodic part of the 
                                                 !   wave functions will be 
                                                 !   computed
  real(dp),    intent(in)  :: bvector(3)         ! Vector that connects a 
                                                 !   given k-point with its 
                                                 !   neighbour (in Bohr^-1)
  integer,     intent(in)  :: nuo                ! Number of orbitals in local
                                                 !   node
                                                 ! NOTE: Running in parallel
                                                 !   this is core dependent
                                                 !   Sum_{cores} no_l = no_u
  integer,     intent(in)  :: nuotot             ! Number of orbitals in the 
                                                 !   unit cell.
  complex(dp), intent(in)  :: psiatk(nuotot,nuo) ! Coefficients of the 
                                                 !    eigenvector at the k-point
                                                 !    First  index: orbital
                                                 !   Second index: band
  complex(dp), intent(in)  :: psinei(nuotot,nuo) ! Coefficients of the 
                                                 !   eigenvector at the
                                                 !   neighbour k-point
                                                 !   First  index: orbital
                                                 !   Second index: band
! Note that, when running in parallel, the information known by each local node
! for psiatk, and psinei is: each local node knows all the coefficients
! for only a few nuo bands.
  integer,     intent(in)  :: maxnh              ! Maximum number of orbitals
                                                 !   interacting
                                                 !   NOTE: In parallel runs,
                                                 !   maxnh changes from node to 
                                                 !   node
  complex(dp), intent(in)  :: delkmat(maxnh)     ! Matrix elements of a plane
                                                 !   wave
  integer,     intent(in)  :: nbandsocc_loc      ! Number of ocuppied bands 
                                                 !   in the local node
  integer,     intent(in)  :: nbandsocc          ! Number of occupied bands
  complex(dp), intent(out) :: Mkb(nbandsocc,nbandsocc) 
                                                 ! Overlap matrix between the 
                                                 !   the periodic part of the

! Internal variables
  integer     :: imu                             ! Counter for the loop on orbit
  integer     :: inu                             ! Counter for the loop on orbit
  integer     :: nband                           ! Counter for the loop on bands
  integer     :: mband                           ! Counter for the loop on bands
  integer     :: jneig                           ! Counter for the loop on neigb
  integer     :: ia                              ! Atomic index
  integer     :: ind                             ! Index for the neighbour
                                                 !   orbital in the list
  integer     :: jo                              ! Neighbour orbital in the
                                                 !   supercell
  real(dp)    :: kxmu                            ! Dot product between the 
                                                 !   b vector and the atomic
                                                 !   position. 
                                                 !   See the term in the
                                                 !   bracket in the right hand
                                                 !   side of Eq. (5) of the
                                                 !   paper by 
                                                 !   D. Sanchez-Portal et al.
                                                 !   Fundamental Physics for 
                                                 !   Ferroelectrics 
                                                 !   (AIP Conf. Proc. Vol 535) 
                                                 !   ed R. Cohen (Melville, AIP)
                                                 !   pp 111-120 (2000).
  real(dp)    :: kxij                            ! Dot product between the  
                                                 !   wave vector k and the 
                                                 !   relative position of the 
                                                 !   two orbitals (see the  
                                                 !   exponential right after the
                                                 !   summatory in Eq. (5) of 
                                                 !   the paper by 
                                                 !   D. Sanchez-Portal et al.
                                                 !   Fundamental Physics for 
                                                 !   Ferroelectrics 
                                                 !   (AIP Conf. Proc. Vol 535) 
                                                 !   ed R. Cohen (Melville, AIP)
                                                 !   pp 111-120 (2000).

  complex(dp) :: eibr                            ! Exponential exp(i kxmu )
  complex(dp) :: eikr                            ! Exponential exp(i kxij )
  complex(dp) :: pipj                            ! Product of the coefficients
                                                 !   of the wave functions
  complex(dp), dimension(:,:), pointer :: aux => null()
  complex(dp), dimension(:,:), pointer :: aux2 => null()

  integer     :: imu_global                      ! Global index of the atom orbi
#ifdef MPI
  integer     :: MPIerror
  integer     :: inode                           ! Counter for the loop on nodes
  integer     :: norb_max_loc                    ! Maximum number of atomic 
                                                 !   orbitals that will be 
                                                 !   stored on any node
  integer     :: norb_loc                        ! Number of atomic orbitals
                                                 !   on the local node
  integer     :: moccband_global                 ! Global index of the 
                                                 !   occupied band
  integer     :: moccband_sequential             ! Global index of the 
                                                 !   occupied band in sequential
                                                 !   notation
  integer     :: noccband_global                 ! Global index of the 
                                                 !   occupied band
  integer     :: noccband_sequential             ! Global index of the 
                                                 !   occupied band in sequential
                                                 !   notation
  complex(dp), dimension(:,:), pointer :: auxtmp => null() ! Temporal arrays used to  
  complex(dp), dimension(:,:), pointer :: aux2loc => null() !   broadcast auxiliary
                                                 !   matrices
#endif


! Start time counter
  call timer( 'overkkneig', 1 )

! Initialize Mkb
  Mkb(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)

!!     For debugging
!  write(6,'(a,i5,3f12.5)')    &           
! &  'overkkneig: Node, kvector =', Node, kvector
!  write(6,'(a,i5,3f12.5)')    &           
! &  'overkkneig: Node, bvector =', Node, bvector
!  write(6,'(a,2i5)')    &           
! &  'overkkneig: Node, maxnh =', Node, maxnh
!      if( Node .eq. 1 ) then
!!      if( IOnode ) then
!        do nband = 1, nuo
!          do imu = 1, nuotot
!            write(6,'(a,2i5,2f12.5)')         &
! &            'overkkneig: nband, imu, coeff = ', imu, nband, psiatk(imu,nband)
!          enddo 
!        enddo 
!      endif
!!     End debugging


! Compute the auxiliary matrix aux, defined as:
! \sum_{\vec{R}_{l}} 
! exp^{i \vec{k} \cdot \left(\vec{R}_{\mu}-\vec{R}_{\nu}-\vec{R}_{l} \right) }
! \int d \vec{r} \phi_{\nu} \left(\vec{r}-\vec{R}_{\nu}-\vec{R}_{l} \right)
! exp^{- i \vec{b} \left(\vec{r}-\vec{R}_{\mu} \right) }
! \phi_{\mu} \left(\vec{r}-\vec{R}_{\mu} \right)
! See, for instance, Eq. (5) of the paper by D. Sanchez-Portal et al.
! Fundamental Physics for Ferroelectrics (AIP Conf. Proc. Vol 535) 
! ed R. Cohen (Melville, AIP) pp 111-120
! In aux, the second index refer to the orbital in the unit cell
! (\mu in the previous notation),
! while the first index refer to the neighbour orbital (\nu) in the 
! previous notation.
  
! Allocate local memory.
! The second dimension of aux is equal to the number of orbitals
! in the local node, nuo, since that is the maximum number of 
! orbitals mu (in the previous equation) that can be computed locally
  call re_alloc( aux, 1, nuotot, 1, nuo, name='aux', routine='overkkneig')

  aux(:,:)  = cmplx(0.0_dp,0.0_dp,kind=dp)

! Loop on all the atomic orbitals known by the local node
  do imu = 1, nuo
!   Identify the global index of the orbital
#ifdef MPI
    call LocalToGlobalOrb( imu, Node, Nodes, imu_global)
!!   For debugging
!    write(6,'(a,3i5)')' overkkneig: Node, imu, imu_global = ',  &
! &                                  Node, imu, imu_global
!!   End debugging
#else
    imu_global = imu
#endif
!   For each atomic orbital, find the atom to which that orbital belongs...
    ia = iaorb( imu_global )
!   ... the position of the atom and the dot product between the vector 
!   connecting neighbour k-points, b, and the position of the atom.
!   See the term in the bracket in the right hand side of Eq. (5) of the paper
!   D. Sanchez-Portal et al. 
!   Fundamental Physics for Ferroelectrics (AIP Conf. Proc. Vol 535) 
!   ed R. Cohen (Melville, AIP) pp 111-120, 2000
    kxmu = bvector(1) * xa(1,ia) +      &
 &         bvector(2) * xa(2,ia) +      &
 &         bvector(3) * xa(3,ia)
    eibr = cmplx( dcos(kxmu), dsin(kxmu), kind=dp )
    do jneig = 1, numh( imu ) 
      ind = listhptr(imu) + jneig
      jo  = listh(ind)
      inu = indxuo(jo)
!     Compute the dot product between the wave vector k and the relative
!     position of the two orbitals (see the exponential right after the
!     summatory in Eq. (5) of the paper by D. Sanchez-Portal et al.
!     Fundamental Physics for Ferroelectrics (AIP Conf. Proc. Vol 535) 
!     ed R. Cohen (Melville, AIP) pp 111-120 (2000).
      kxij = kvector(1) * xijo(1,ind) + &
 &           kvector(2) * xijo(2,ind) + &
 &           kvector(3) * xijo(3,ind)
!     The origin for the relative position vector is taken on
!     the atom in the unit cell.
!     In the formula for overlap matrix at neighbor k-points,
!     See, for instance, Eq. (5) of the paper by D. Sanchez-Portal et al.
!     Fundamental Physics for Ferroelectrics (AIP Conf. Proc. Vol 535) 
!     ed R. Cohen (Melville, AIP) pp 111-120 (2000).
!     just the opposite vector appears. We have to change the sign
!     of the phase factor.
      kxij = -1.0_dp * kxij
      eikr = cmplx( dcos(kxij), dsin(kxij), kind=dp )
!     imu runs on the local orbitals on this node,
!     so it has to be the second index in aux 
!     (remember that the second dimension in aux has been allocated with 
!     nuo = no_l)
      aux(inu,imu) = aux(inu,imu) + eikr * delkmat(ind) * eibr
    enddo   ! End loop on neighbour orbitals
  enddo     ! End loop on orbitals in the unit cell local to this node

!! For debugging
!!  if( IOnode ) then
!  if( Node .eq. 1 ) then
!    do imu = 1, nuotot
!      do inu = 1, nuotot
!        write(6,'(a,3i5,2f12.5)')                   &
! &        'overkkneig: Node, imu, inu, aux = ',     &
! &                     Node, imu, inu, aux(inu,imu)
!      enddo
!    enddo
!  endif
!! End debugging

! Compute \sum_{\mu} c_{\mu,m} \left( \vec{k} + \vec{b} \right) aux_{\nu,\mu}
! and store it in an auxiliary matrix

! Allocate local memory
  call re_alloc( aux2, 1, nbandsocc, 1, nuotot,      &
 &               name='aux2', routine='overkkneig' )

  aux2(:,:)  = cmplx(0.0_dp,0.0_dp,kind=dp)

#ifdef MPI

! Compute the number of atomic orbitals stored on Node 0
! We assume that this is the maximum number of orbitals that will be stored
! on any node.
  call GetNodeOrbs( nuotot, 0, Nodes, norb_max_loc)

!! For debugging
!  write(6,'(a,4i5)')' overkkneig: Node, nuotot, Nodes, norb_max_loc = ', &
! &                                Node, nuotot, Nodes, norb_max_loc
!! End debugging

! Allocate the temporal variable that will be used to broadcast
! the auxiliary matrix aux (see comments for its definition above)
! to all the other nodes
! auxtmp follows the same structure as aux:
! First index: all the atomic orbitals in the unit cell (nuotot)
! Second index: maximum number of atomic orbitals stored in a local node
  call re_alloc( auxtmp, 1, nuotot, 1, norb_max_loc, &
 &               name='auxtmp', routine='overkkneig' )
  auxtmp(:,:)  = cmplx(0.0_dp,0.0_dp,kind=dp)

! Loop on all the nodes:
! This is required because a given Node, for instance Node = 0,
! knows the coefficients of the wave function at a neighbour k-point
! for all the atomic orbitals (mu = 1, nuotot),
! while it only knows the matrix aux for some of the matrix orbitals (mu=1,nuo).
! To compute the sum on mu we have to take to the Node = 0 the full
! aux matrix from the rest of the nodes.

  do inode = 0, Nodes-1

!   Compute the number of occupied bands stored on Node inode 
    call GetNodeOrbs( nuotot, inode, Nodes, norb_loc )

!   Copy the auxiliary matrix for inode to a temporal variable
    if( Node .eq. inode ) then
      do imu = 1, norb_loc
        do inu = 1, nuotot
          auxtmp(inu,imu) = aux(inu,imu)
        enddo
      enddo
    endif

!   Broadcast the auxiliary matrix from node inode to all the other nodes
    call MPI_Bcast( auxtmp(1,1), nuotot*norb_loc, &
 &                  MPI_double_complex, inode, MPI_Comm_World, MPIerror )

!   Loop on the occupied bands stored on the local node (Node)
    do mband = 1, nbandsocc_loc

!     Identify the global index of the occupied band
       moccband_global     = which_band_in_node(Node,mband)
       moccband_sequential = sequential_index_included_bands(moccband_global)
!!      For debugging
!        write(6,'(a,4i5)')                                            &
! &        ' overkkisig: Node, mband, moccband_global, sequential = ', &
! &                      Node, mband, moccband_global, moccband_sequential
!!      End debugging

!     Loop on all the atomic orbitals of the unit cell
!     As parallelized now, a given node knows
!     the coefficients of psinei for all the atomic orbitals 
!     of nuo = no_l bands

      do inu = 1, nuotot
!       Perform the sum on mu, following the notation of Eq. (5) in the 
!       paper by D. Sanchez-Portal et al.
!       Fundamental Physics for Ferroelectrics (AIP Conf. Proc. Vol 535) (2000)
        do imu = 1, norb_loc
!         Identify the global index of the occupied band
          call LocalToGlobalOrb( imu, inode, Nodes, imu_global )
!!         For debugging
!          write(6,'(a,4i5)')' overkkneig: Node, inu, imu, imu_global = ',  &
! &                                        Node, inu, imu, imu_global
!!         End debugging
          aux2(moccband_sequential,inu) = aux2(moccband_sequential,inu) +  &
 &            psinei(imu_global,mband) * auxtmp(inu,imu)
        enddo   ! End loop on the sumatory on mu
      enddo     ! End loop on atomic orbitals in the unit cell (nuotot)
    enddo       ! End loop on local occupied bands (mband)
  enddo         ! End loop on nodes (inode)

! Allocate workspace array for global reduction
  call re_alloc( aux2loc, 1, nbandsocc, 1, nuotot,      &
 &               name='aux2loc', routine='overkkneig' )
! Global reduction of aux2loc matrix
  aux2loc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
  call MPI_AllReduce( aux2(1,1), aux2loc(1,1), nbandsocc*nuotot, &
 &                    MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
  aux2(:,:) = aux2loc(:,:)
#else
  do mband = 1, nbandsocc
    do inu = 1, nuotot
      do imu = 1, nuotot
        aux2(mband,inu) = aux2(mband,inu) +    &
 &                        psinei(imu,mband) * aux(inu,imu)
      enddo
    enddo
  enddo
#endif

!! For debugging
!  if( IOnode ) then
!!  if( Node .eq. 1 ) then
!    do nband = 1, nuo
!      do inu = 1, nuotot
!        write(6,'(a,3i5,2f12.5)')                   &
! &        'overkkneig: Node, inu, nband, psinei = ',     &
! &                     Node, inu, nband, psinei(inu,nband)
!      enddo
!    enddo
!
!    do inu = 1, nuotot
!      do nband = 1, nbandsocc
!        write(6,'(a,3i5,2f12.5)')                   &
! &        'overkkneig: Node, inu, nband, aux2 = ',     &
! &                     Node, inu, nband, aux2(nband,inu)
!      enddo
!    enddo
!  endif
!! End debugging

  do nband = 1, nbandsocc_loc
    do mband = 1, nbandsocc
      do inu = 1, nuotot
#ifdef MPI
!       Identify the global index of the occupied band
        noccband_global     = which_band_in_node(Node,nband)
        noccband_sequential = sequential_index_included_bands(noccband_global)
        Mkb(noccband_sequential,mband) = Mkb(noccband_sequential,mband) +    &
          conjg(psiatk(inu,nband)) * aux2(mband,inu)
#else
        Mkb(nband,mband) = Mkb(nband,mband) +    &
          conjg(psiatk(inu,nband)) * aux2(mband,inu)
#endif
      enddo
    enddo
  enddo

#ifdef MPI
! Allocate workspace array for global reduction
  call re_alloc( aux2loc, 1, nbandsocc, 1, nbandsocc,      &
 &               name='aux2loc', routine='overkkneig' )
! Global reduction of aux2loc matrix
  aux2loc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
  call MPI_AllReduce( Mkb(1,1), aux2loc(1,1), nbandsocc*nbandsocc, &
 &                    MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
  Mkb(:,:) = aux2loc(:,:)
#endif

!! For debugging
!  if( IOnode ) then
!!  if( Node .eq. 1 ) then
!    do inu = 1, nuotot
!      do nband = 1, nuo
!        write(6,'(a,3i5,2f12.5)')                   &
! &        'overkkneig: Node, inu, nband, psiatk = ',     &
! &                     Node, inu, nband, psiatk(inu,nband)
!      enddo
!    enddo
!
!    do nband = 1, nbandsocc
!      do mband = 1, nbandsocc
!        write(6,'(a,3i5,2f12.5)')                       &
! &        'overkkneig: Node, nband, mband, Mkb = ',     &
! &                     Node, nband, mband, Mkb(nband,mband)
!      enddo
!    enddo
!  endif
!! End debugging

  call de_alloc( aux,    name='aux', routine="overkkneig"    )
  call de_alloc( aux2,   name='aux2', routine="overkkneig"   )
#ifdef MPI
  call de_alloc( auxtmp, name='auxtmp', routine="overkkneig"   )
  call de_alloc( aux2loc, name='aux2loc', routine="overkkneig"   )
#endif

! End time counter
  call timer( 'overkkneig', 2 )

  end subroutine overkkneig

