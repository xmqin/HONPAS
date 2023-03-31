! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

subroutine Mmn( ispin ) 
!
!     In this subroutine we compute the overlaps between Bloch orbitals
!     at neighboring k points:
!
!     $M_{m n}^{(\vec{k}, \vec{b})} = 
!        \langle u_{m \vec{k}} \vbar u_{n \vec{k} + \vec{b} \rangle$
!
!     Eq. (27) of the paper by N. Marzari et al. 
!     Review of Modern Physics 84, 1419 (2012)
!
!     In the previous formula, only the periodic part of the wave functions
!     at neighbour k-points enter in the Equation.
!     We have to adapt this equation to the input provided by Siesta
!     (the coefficients of the whole wave function, not only of the periodic
!     part). 
!     This is done following Eq. (5) of the paper by 
!     D. Sanchez-Portal et al., Fundamental Physics for Ferroelectrics 
!     (AIP Conf. Proc. Vol 535) ed R. Cohen (Melville, AIP) pp 111-120 (2000).
!
!     OUTPUT: 
!     File called seedname.mmn, where the overlap matrices are written in 
!     the format required by Wannier90
!
!     Implemented by J. Junquera and R. Korytar, July 2013
!

  use precision,          only: dp           ! Real double precision type
  use m_siesta2wannier90, only: numkpoints   ! Total number of k-points
                                             !   for which the overlap of the
                                             !   periodic part of the wavefunct
                                             !   with a neighbour k-point will
                                             !   be computed
  use m_siesta2wannier90, only: kpointsfrac  ! List of k points relative to the 
                                             !   reciprocal lattice vectors.
                                             !   First  index: component
                                             !   Second index: k-point index 
                                             !      in the list
  use m_siesta2wannier90, only: nncount      ! Number of nearest k-pnt neighbors
  use m_siesta2wannier90, only: nnlist       ! nnlist(ikp,inn) is the index of
                                             !   the inn-neighbour of ikp-point
                                             !   in the Monkhorst-Pack grid 
                                             !   folded to first Brillouin zone
  use m_siesta2wannier90, only: nnfolding    ! nnfolding(i,ikp,inn) is the 
                                             !   i-component of the reciprocal 
                                             !   lattice vector 
                                             !   (in reduced units) that brings
                                             !   the inn-neighbour specified 
                                             !   in nnlist
                                             !   (which is in the first BZ) to 
                                             !   the actual \vec{k} + \vec{b} 
                                             !   that we need.
  use m_siesta2wannier90, only: numbands     ! Number of bands for wannierizatio
                                             !   before excluding bands         
  use m_siesta2wannier90, only: numincbands  ! Number of bands for wannierizatio
                                             !   after excluding bands         
  use m_siesta2wannier90, only: nincbands_loc! Number of bands for wannierizatio
                                             !   after excluding bands         
                                             !   in the local node
  use m_siesta2wannier90, only: numproj      ! Total number of projectors
  use m_siesta2wannier90, only: bvectorsfrac ! Vectors that connect each 
                                             !   mesh k-point to its 
                                             !   nearest neighbours.
  use m_siesta2wannier90, only: coeffs       ! Coefficients of the
                                             !   wavefunctions.
                                             !   First  index: orbital
                                             !   Second index: band
                                             !   Third  index: k-point
  use m_siesta2wannier90, only: Mmnkb        ! Matrix of the overlaps of 
                                             !   periodic parts of Bloch waves.
                                             !   <u_{n,k}|u_{m,k+b}>
  use atomlist,           only: no_s         ! Number of orbitals in supercell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: no_l         ! Number of orbitals in local node
                                             ! NOTE: When running in parallel,
                                             !   this is core dependent
                                             !   Sum_{cores} no_l = no_u
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: iaorb        ! Atomic index of each orbital
  use siesta_geom,        only: xa           ! Atomic positions
  use sparse_matrices,    only: maxnh        ! Maximum number of orbitals
                                             !   interacting
                                             ! NOTE: While running in parallel,
                                             !   maxnh changes from one core to 
                                             !   the other

  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines
  use parallel,           only: IOnode       ! Input/output node
  use m_digest_nnkp,      only: getdelkmatgenhandle
  use m_noccbands,        only: noccupied    ! Number of occupied bands for a 
                                             !   given spin direction
  use m_planewavematrixvar, only: delkmat    ! Matrix elements of a plane wave
                                             !   Only for one \vec{b} vector
  use m_planewavematrixvar, only: delkmatgen ! Matrix elements of a plane wave
                                             !   (array that contains the 
                                             !   matrices for all the \vec{b}
                                             !   vectors
! For debugging
  use parallel,           only: Node, Nodes
! End debugging

  implicit none

  integer,  intent(in) :: ispin              ! Spin component
! Internal variables

  integer  :: ik             ! Counter for the k-point loop
  integer  :: io             ! Counter for the orbital loop
  integer  :: jo             ! Counter for the orbital loop
  integer  :: iuo            ! Counter for the orbital loop
  integer  :: iband          ! Counter for the bands loop
  integer  :: inn            ! Counter for the k-point neighbor loop
  integer  :: indexneig      ! Index of the neighbour k-point in the list
                             !   of k-points
  integer  :: fold           ! Should we fold the coefficients of the 
                             !   wavefunction at the neighbour k-point
                             !   out of the first Brillouin zone?
  integer  :: ia             ! Index of the atom to which an atomic orbital
                             !   belongs
  integer  :: handle         ! Given a k-point vector and a neighbor, 
                             !   separated by a vector \vec{b},
                             !   it gives the position of the delkmatgen 
                             !   array where exp^(i \vec{b} \cdot \vec{r})
                             !   will be stored
  integer  :: nincbands      ! Number of occupied bands
  real(dp) :: kvector(3)     ! k-point vector for which the Overlap matrix 
                             !   between the periodic part of the 
                             !   wave functions will be computed
  real(dp) :: kvectorneig(3) ! Wave vector of the neighbor k-point
  real(dp) :: gfold(3)       ! Reciprocal lattice vector that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
  real(dp) :: bvectoraux(3)  ! Auxiliary vector
  real(dp) :: bvector(3)     ! Vector that connects a given k-point with its 
                             !   neighbour (in Bohr^-1)
  real(dp) :: gxij           ! Dot product of the reciprocal lattice vector
                             !   gfold with an atomic position
  real(dp) :: foldfrac(3)    ! Auxiliar vector to compute the folding vector
  complex(dp) :: eigx        ! Exponential exp^( i * gxij )

  complex(dp), pointer   :: coeffs2(:,:) => null()! Auxiliary array to store the 
                                         !    coefficients of the wave function


! Start time counter
  call timer( 'Mmn', 1 )

! Allocate memory related with the overlap matrix between periodic parts
! of Bloch functions at neighbour k-points.
! These matrix will depend on four indices (see Eq. (27) of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012):
! $M_{m n}^{(\vec{k}, \vec{b})} = 
!    \langle u_{m \vec{k} \vbar u_{n \vec{k} + \vec{b}} \rangle
! where $m$ and $n$ run between 1 and the number of occupied bands
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed and 
! $\vec{b}$ runs over all the neighbours of the k-point.
! These last two variables are read from the .nnkp file.
  nincbands = numincbands( ispin )
  call re_alloc( Mmnkb,          &
 &               1, nincbands,   &
 &               1, nincbands,   &
 &               1, numkpoints,  &
 &               1, nncount,     &
 &               'Mmnkb',        &
 &               'Mmn' )

! Allocate the variable to store the coefficients of the wavefunction
  call re_alloc( coeffs2,           &
 &               1, no_u,           &
 &               1, nincbands_loc,  &
 &               'coeffs2',         &
 &               'Mmn' )

!! For debugging
!  write(6,'(a,3i5)')' Mmn, no_u, nincbands_loc = ',  &
! &  Node, no_u, nincbands_loc
!! End debugging

kpoints:                         &
  do ik = 1, numkpoints

!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!   Loop on the neighbour k-points for a given k.
kneighbour:                      &
    do inn = 1, nncount

!     Get the coordinates of the neighbor k-point.
      indexneig = nnlist(ik,inn)
      call getkvector( kpointsfrac(:,indexneig), kvectorneig )
!!     For debugging
!      if( IOnode ) then
!        write(6,'(a,3f12.5)')         &
! &        'mmn: kvectorneig = ', kvectorneig
!      endif
!!     End debugging

!     Find the coefficients of the wave function for the neighbour k-point
!     Here, we obtain $\psi_{m} (\vec{k} + \vec{b})$, 
!     where m runs between 1 and the number of bands included for wannierization
!     in the local node (nincbands_loc)
      coeffs2(:,:) = coeffs(:,:,indexneig)

!     The neighbour k-point, as specified in the nnlist, 
!     is always in the first Brillouin zone.
!     To find the actual k-point, we might have to add a vector of the
!     reciprocal lattice, as specified in the nnfolding matrix
      fold =  nnfolding(1,ik,inn)**2 +  &
 &            nnfolding(2,ik,inn)**2 +  &
 &            nnfolding(3,ik,inn)**2

!     The coeffs of the wave function, as obtained in diagpol, 
!     are for a wave vector in the first Brillouin zone. 
!     If the actual neighbour is out,
!     we apply a simple transformation that holds in the periodic gauge
!     c_{i \mu} (k+b) = c_{i\mu}(k+b-G) exp(-iG\cdot r_\mu)
!     where -G brings the k+b to the first BZ
      if ( fold .gt. 0 ) then
        foldfrac(:) = nnfolding(:,ik,inn) * 1.0_dp
        call getkvector( foldfrac, gfold )  
!!       For debugging
!        if ( IOnode )
!          write(6,'(a,i5,3i5,3f12.5)')         &
! &          'mmn: inn, gfold = ', inn, nnfolding(:,ik,inn), gfold
!        endif
!!       End debugging
        do iband = 1, nincbands_loc
          do jo = 1, no_u
!           Localize the position where the atom is centered
            ia = iaorb(jo)
!           Compute exp( i G tau_{mu} )
            gxij = dot_product( gfold,xa(:,ia) )
            eigx = cmplx( dcos(gxij), dsin(gxij), kind=dp )
!           Multiply the coefficient times the gauge
            coeffs2(jo,iband) = coeffs2(jo,iband) * conjg(eigx)
          enddo
        enddo
      endif

!!     For debugging
!!      if( Node .eq. 1 ) then
!      if( IOnode ) then
!        write(6,'(a,2i5)')         &
! &        'mmn: inn, fold = ', inn, fold
!        do iband = 1, nincbands_loc
!          write(6,'(a,i5,f12.5)')         &
! &          'mmn: iband, epsilon = ', iband, epsilon(iband)
!          do io = 1, no_u
!            write(6,'(a,i5,2f12.5)')         &
! &            'mmn: io, coeff = ', io, coeffs2(io,iband)
!          enddo 
!        enddo 
!      endif
!!     End debugging

!     Now we have to determine the position in the delkmatgen array,
!     where exp^(i \vec{b} \cdot \vec{r}) is stored,
!     where \vec{b} is the vector connecting kvector and kvectorneig.
      bvectoraux(:) = kpointsfrac(:,nnlist(ik,inn)) +               &
 &                     nnFolding(:,ik,inn) - kPointsFrac(:,ik)
      handle = getdelkmatgenhandle( bvectoraux, nncount, bvectorsfrac )
      call getkvector( bvectoraux, bvector )

      delkmat(1:maxnh) = delkmatgen(handle,1:maxnh)

!! For debugging
!      if ( IOnode ) then
!        write(6,'(a,2i5)')         &
! &        'mmn: inn, handle  = ', inn, handle
!        write(6,'(a,i5,3f12.5)')   &
! &        'mmn: inn, bvector = ', inn, bvector
!!        do io = 1, maxnh
!!          write(6,'(a,i5,2f12.5)')         &
!! &          'mmn: io, delkmat  = ', io, delkmat(io)
!!        enddo
!      endif 
!! End debugging

!     Initialize Mmnkb
      Mmnkb(:,:,ik,inn) = cmplx(0.0_dp, 0.0_dp, kind=dp)

!     Compute Mmnkb, following Eq. (5) of the paper by 
!     D. Sanchez-Portal et al., Fundamental Physics for Ferroelectrics 
!     (AIP Conf. Proc. Vol 535) ed R. Cohen (Melville, AIP) pp 111-120 (2000).
      call overkkneig( kvector, bvector, no_l, no_u,                &
 &                     coeffs(:,:,ik), coeffs2,                     &
 &                     maxnh, delkmat, nincbands_loc, nincbands,    &
 &                     Mmnkb(:,:,ik,inn) )

    enddo kneighbour ! End loop on neighbour k-points (inn)

  enddo kpoints      ! End loop on the number of k-points (ik)


! Write the Mmn overlap matrices in a file, in the format required
! by Wannier90
  if( IOnode ) call writemmn( ispin )

  call de_alloc(coeffs2, 'coeffs2', 'Mmn')

! End time counter
  call timer( 'Mmn', 2 )
 
end subroutine Mmn
