! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_siesta2wannier90

! OUTPUT: 
! File called seedname.mmn, where the overlap matrices are written in 
! the format required by Wannier90
!
! File called seedname.eig, where the eigenvalues of the Hamiltonian
! at the required k-points are written according to the format required
! by Wannier90.

  use precision, only : dp                   ! Real double precision type
  use parallel,  only : IOnode               ! Input/output node
  use parallel,  only : Node                 ! This process node
  use parallel,  only : Nodes                ! Total number of processor nodes
  use sys,       only : die                  ! Termination routine
  use files,     only : label_length         ! Number of characters in slabel
  use siesta_options, only: w90_write_mmn       ! Write the Mmn matrix for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_amn       ! Write the Amn matrix for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_eig       ! Write the eigenvalues for the
                                             !   interface with Wannier
  use siesta_options, only: w90_write_unk       ! Write the unks for the
                                             !   interface with Wannier
  use TrialOrbitalClass

!
! Variables related with the atomic structure of the system
!

  real(dp) :: latvec(3,3)    ! Real lattice vectors.
                             !   Cartesian coordinates. 
                             !   Readed in Angstroms and transformed to Bohr
                             !   internally
                             !   First  index: component
                             !   Second index: vector   
                             !   This is consistent with the unit cell read
                             !   in Siesta, but it has changed with respect
                             !   the first version implemented by R. Korytar
  real(dp) :: reclatvec(3,3) ! Reciprocal lattice vectors.
                             !   Cartesian coordinates. 
                             !   Readed in Angstroms^-1 and transformed 
                             !   to Bohr^-1 internally
                             !   First  index: component 
                             !   Second index: vector
                             !   This is consistent with the reciprocal 
                             !   lattice read in Siesta
                             !   in Siesta, but it has changed with respect
                             !   the first version implemented by R. Korytar
!
! Variables related with the k-point list for which the overlap
! matrices Mmn between a k-point and its neighbor will be computed
!
  integer  :: numkpoints     ! Total number of k-points
                             !   for which the overlap of the
                             !   periodic part of the wavefunct
                             !   with a neighbour k-point will
                             !   be computed
  real(dp), pointer :: kpointsfrac(:,:) 
                             ! List of k points relative
                             !   to the reciprocal lattice vectors.
                             !   First  index: component
                             !   Second index: k-point index in the list
  real(dp), pointer :: bvectorsfrac(:,:)
                             ! The vectors b that connect each mesh-point k
                             !   to its nearest neighbours
!
! Variables related with the neighbours of the k-points
!
  integer           :: nncount  
                             ! The number of nearest neighbours belonging to 
                             !   each k-point of the Monkhorst-Pack mesh
  integer, pointer  :: nnlist(:,:)
                             ! nnlist(ikp,inn) is the index of the 
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the 
                             !   first Brillouin zone
  integer, pointer  :: nnfolding(:,:,:)
                             ! nnfolding(i,ikp,inn) is the i-component 
                             !   of the reciprocal lattice vector 
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units.
!
! Variables related with the projections with trial functions,
! initial approximations to the MLWF
!
  integer  :: numproj        ! Total number of projection centers,
                             !   equal to the number of MLWF
  type(trialorbital), target, allocatable  :: projections(:)       
!
! Variables related with the number of bands considered for wannierization
!
  integer          :: numbands(2) ! Number of bands for wannierization
                                  !    before excluding bands
  integer          :: numexcluded 
                             ! Number of bands to exclude from the calculation
                             !   of the overlap and projection matrices.
                             !   This variable is read from the .nnkp file

  integer          :: numincbands(2) 
                             ! Number of included bands in the calc. 
                             !   of the overlap and projection matrices.
                             !   after excluding the bands
  integer          :: nincbands_loc
                             ! Number of included bands in the calc. 
                             !   of the overlap and projection matrices.
                             !   in the local Node
  integer          :: blocksizeincbands ! Maximum number of bands
                            !    considered for wannierization per node

  integer, pointer :: excludedbands(:)
                             ! Bands to be excluded
                             !   This variable is read from the .nnkp file
  logical, pointer :: isexcluded(:) ! Masks excluded bands
  integer, pointer :: isincluded(:) ! Masks included bands

!
! Variables related with the coefficients of the wavefunctions and
! eigenvalues at the Wannier90 k-point mesh
!
  complex(dp), pointer :: coeffs(:,:,:) => null() ! Coefficients of the wavefunctions.
                                         !   First  index: orbital
                                         !   Second index: band
                                         !   Third  index: k-point
  real(dp),    pointer :: eo(:,:) => null()        ! Eigenvalues of the Hamiltonian 
                                         !   at the numkpoints introduced in
                                         !   kpointsfrac 
                                         !   First  index: band index
                                         !   Second index: k-point index


!
! Output matrices
!

 complex(dp), pointer :: Mmnkb(:,:,:,:) => null()  ! Matrix of the overlaps of 
                                         !   periodic parts of Bloch waves.
                                         !   <u_{ik}|u_{jk+b}>
                                         !   The first two indices refer to 
                                         !   the number of occupied bands
                                         !   (indices m and n in standard
                                         !   notation, see for instance, 
                                         !   Eq. (27) of the paper by 
                                         !   Marzari et al., RMP 84, 1419 (2012)
                                         !   The third index refer to the kpoint
                                         !   The fourth index refer to the neig
 complex(dp), pointer :: Amnmat(:,:,:) => null()  ! Projections of a trial function
                                         !   with a Bloch orbital
                                         !   <\psi_{m k}|g_n>

!
! Variables related with the input/output files
!
 character(label_length+3)  :: seedname  ! Name of the file where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, reads or dumps the 
                                         !   information.

CONTAINS

subroutine siesta2wannier90

  use m_spin,        only: nspin         ! Number of spin components
  use files,         only: slabel        ! Short system label, 
                                         !   used to generate file names
  use m_digest_nnkp, only: read_nnkp     ! Subroutine that reads the .nnkp file
  use m_digest_nnkp, only: chosing_b_vectors ! Subroutine that computes the b
                                         ! vectors that connect each mesh 
                                         ! k-point to its nearest neighbours.
  use m_digest_nnkp, only: set_excluded_bands   ! Subroutine that chooses the 
                                                !   bands that are excluded from
                                                !   the wannierization procedure
  use m_digest_nnkp, only: number_bands_wannier ! Subroutine that computes the
                                         ! number of bands for wannierization

  implicit none

  integer                    :: ispin    ! Spin counter
  integer                    :: numbandswan(2)
                                         ! Number of bands for wannierization

  external  :: timer, mmn

  call timer("siesta2wannier90",1)

  do ispin = 1, nspin
!   Append _up or _dn to the seedname if spin polarized
    seedname = trim(getFileNameRoot(ispin,nspin,slabel))
!   A priori, the .nnkp file generated by Wannier90 used as a 
!   postprocessing tool is independent of spin.
!   However, they provide examples (example08, bulk Fe) where
!   two different .nnkp files are generated, 
!   one for each component of the spin.
!   This inspired us to read the .nnkp file for each component of the spin
!   (note that the seedname will be different in every case).
!   The information is read only by the master node,
!   and broadcast to the rest of the nodes later. 
    if (IOnode) then
       write(6,'(/,a)')  &
 &       'siesta2wannier90: Reading the ' // trim(seedname) // '.nnkp file'
    endif
    call read_nnkp( seedname, latvec, reclatvec, numkpoints,          &
                    kpointsfrac, nncount, nnlist, nnfolding,          &
                    numproj, projections, numexcluded, excludedbands, &
                    w90_write_amn )

!   Compute the vectors that connect each mesh k-point to its nearest neighbours
    call chosing_b_vectors( kpointsfrac, nncount, nnlist, nnfolding,  &
                            bvectorsfrac )

!   Compute the number of bands for wannierization
    call number_bands_wannier( numbandswan )
    numbands(ispin) = numbandswan(ispin)

!   Chose which bands are excluded from the wannierization procedure
    call set_excluded_bands(  ispin, numexcluded, excludedbands, numbands, &
 &                            isexcluded, isincluded, numincbands )

!   Compute the matrix elements of the plane wave,
!   for all the wave vectors that connect a given k-point to its nearest
!   neighbours
    call compute_pw_matrix( nncount, bvectorsfrac )

!   Compute the coefficients of the wavefunctions at the 
!   k-points of the Wannier90 mesh
    call diagonalizeHk( ispin )

!   Compute the overlap between periodic parts of the wave functions
    if( w90_write_mmn ) call mmn( ispin )

!   Compute the overlaps between Bloch states onto trial localized orbitals
    if( w90_write_amn ) call amn( ispin )

!   Write the eigenvalues in a file, in the format required by Wannier90
    if( IOnode .and. w90_write_eig ) call writeeig( ispin )

!   Write the values of the Bloch states in a box
    if( w90_write_unk ) call writeunk( ispin )

  enddo

  if (IOnode) then
    write(6,'(/,a)')  &
 &    'siesta2wannier90: All the information dumped in the corresponding files'
    write(6,'(a)')  &
 &    'siesta2wannier90: End of the interface between Siesta and Wannier90'
  endif

  call timer("siesta2wannier90",2)

end subroutine siesta2wannier90

function getFileNameRoot(ispin,nspin,root)
!
! Constructs filenames for various cases of spin:
! Nonpolarized, down and up by adding an extension to
! "root". Behavior:
! (a) unpolarized: no change
! (b) spin up: _up suffix
! (c) spin down: _dn suffix
!
  character(*),intent(in)      :: root
  integer,intent(in)           :: ispin, nspin
  character(len_trim(root)+3)  :: getFileNameRoot

  select case (nspin)
!   If only one spin component, do nothing
    case (1) 
      getFileNameRoot = trim(root)
!   If two spin components, 
!   append "_up" (for the first one) or "_dn" (for the second spin component)
!   to the seed name
    case (2) 
      select case (ispin)
         case (1)
           getFileNameRoot = trim(root)//"_up"
         case (2)
           getFileNameRoot = trim(root)//"_dn"
      end select 
!   If more than two spin components, stop the program 
!   and print and error message
    case default 
      call die("getFileNameRoot: non-collinear spin not implemented yet")
  end select 
end function getFileNameRoot

endmodule m_siesta2wannier90
