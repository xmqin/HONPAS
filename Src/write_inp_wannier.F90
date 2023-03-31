! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

subroutine writemmn( ispin )
!
! In this subroutine, the overlap matrices between the periodic parts of the
! wavefunctions at neighbour k-points are written,
! $M_{m n} ^{( \vec{k} + \vec{b} )}$, see Eq. (27) of the paper by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012).
! 
! The structure of the outfile (seedname.mmn) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname ! Seed for the name of the file 
                                         !   where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, dumps the information.
  use m_siesta2wannier90, only: numincbands! Number of included bands 
                                         !   for wannierization
  use m_siesta2wannier90, only: numkpoints ! Total number of k-points
                                         !   for which the overlap of the
                                         !   periodic part of the wavefunct
                                         !   with a neighbour k-point will
                                         !   be computed
  use m_siesta2wannier90, only: nncount  ! Number of neighbour k-points
  use m_siesta2wannier90, only: nnlist   ! nnlist(ikp,inn) is the index of the 
                                         !   inn-neighbour of ikp-point
                                         !   in the Monkhorst-Pack grid 
                                         !   folded to the 
                                         !   first Brillouin zone
  use m_siesta2wannier90, only: nnfolding ! nnfolding(i,ikp,inn) is the 
                                         !   i-component 
                                         !   of the reciprocal lattice vector 
                                         !   (in reduced units) that brings
                                         !   the inn-neighbour specified in 
                                         !   nnlist (which is in the first BZ)  
                                         !   to the actual \vec{k} + \vec{b} 
                                         !   that we need.
                                         !   In reciprocal lattice units.
  use m_siesta2wannier90, only: Mmnkb    ! Matrix of the overlaps of 
                                         !   periodic parts of Bloch waves.
                                         !   <u_{m k}|u_{n k+b}>
                                         !   The first two indices refer to 
                                         !   the number of occupied bands
                                         !   (indices m and n in standard
                                         !   notation, see for instance, 
                                         !   Eq. (27) of the paper by 
                                         !   Marzari et al., RMP 84, 1419 (2012)
                                         !   The third index refer to the kpoint
                                         !   The fourth index refer to the neig

  implicit none

  integer, intent(in)     :: ispin      ! Spin component

!
! Local variables
!

  character(len=len_trim(seedname)+4) :: mmnfilename ! Name of the file where
                                                     !   the overlap Mmn
                                                     !   matrices will be
                                                     !   written

  integer      :: ik        ! Counter for the loop on k-points
  integer      :: inn       ! Counter for the loop on beighbours
  integer      :: nband     ! Counter for the loop on bands
  integer      :: mband     ! Counter for the loop on bands
  integer      :: mmnunit   ! Logical unit assigned to the file where the 
                            !   Mmn matrices will be written
  integer      :: g(3)      ! Auxiliary variable to write the folding vector
  integer      :: eof

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit


! Asign a logical unit to the file where the Mmn matrices will be written
  call io_assign( mmnunit )

! Asign a name to the file where the Mmn matrices will be written
  mmnfilename = trim( seedname )// ".mmn"

! Open the output file where the Mmn matrices will be written
  open( unit=mmnunit, err=1984, file=mmnfilename, status="replace", &
 &      iostat=eof)

! The first line of the mmn file is a user comment
  write( unit=mmnunit, fmt="(a38)", err=1984 )                      &
 &  "siesta2wannier90: Siesta for Wannier90"

! The second line includes three integers:
!   numincbands: the number of bands for wannierization
!   numkpoints:  the number of k-points for which the overlap of the
!      periodic part of the wavefunction with a neighbour k-point will be 
!      computed
!   nncount:     the number of neighbour k-points 
  write( unit=mmnunit, fmt="(i5,1x,i5,1x,i2)", err=1984 )             &
 &  numincbands(ispin),numkpoints,nncount

! Then, there are numkpoints x nncount blocks of data.
  do ik = 1,numkpoints
    do inn = 1,nncount
   ! override type for the sake of formatting
      g(:) = nnfolding(:,ik,inn)

!     The first line of each block contains five integers:
!     The first specifies the k-points (i.e. gives the ordinal correspondence
!        to its position in the list of k-points in seedname.win).
!     The second to the fifth integers specify \vec{k} + \vec{b}
!     The second integer, in particular, points to the k-point on the list
!        that is a periodic image of \vec{k} + \vec{b}, and in particular
!        is the image that is actually mentioned in the list.
!     The last three integers specify the \vec{G} vector, 
!        in reciprocal lattice units, that brings the k-point specified by
!        the second integer, and that this lives inside the first BZ zone,
!        to the actual \vec{k} + \vec{b} that we need.

      write( unit=mmnunit, fmt="(i5,i5,3x,3i4)", err=1984 )         &
 &      ik, nnlist(ik,inn), g

!     Subsequent numincbands x numincbands lines of each block: 
!     two real numbers per line. 
!     These are the real and imaginary parts, respectively, of the actual
!     scalar product $M_{m n} ^{( \vec{k} + \vec{b} )} for
!     $m, n \in [1, numincbands].
!     The order of these elements is such that the first index $m$ is fastest

      do nband = 1,numincbands(ispin)
        do mband = 1,numincbands(ispin)
        write( unit=mmnunit, fmt="(f12.5,2x,f12.5)", err=1984 )     &
 &         real(  Mmnkb(mband,nband,ik,inn) ),                      &
 &         aimag( Mmnkb(mband,nband,ik,inn) )
        enddo  ! End loop on bands      (mband)
      enddo    ! End loop on bands      (nband)
    enddo      ! End loop on neighbours (inn)
  enddo        ! End loop on k-points   (ik) 

! Close the output file where the Mmn matrices will be written
  call io_close(mmnunit)

  write(6,'(/,a)')  &
 &  'mmn: Overlap matrices between periodic part of wavefunctions'
  write(6,'(a,/)')  &
 &  'mmn: written in ' // trim(seedname) // '.mmn file'

  return

1984 call die("writemmn: Error writing to .mmn file")

end subroutine writemmn

!
! ------------------------------------------------------------------------------
!

subroutine writeamn( ispin )
!
! In this subroutine, the overlap matrices between the trial projection
! functions and the Hamiltonian eigenstates are written,
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle 
! (see paragraph after Eq. (16) of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! 
! The structure of the outfile (seedname.amn) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname ! Seed for the name of the file 
                                         !   where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, dumps the information.
  use m_siesta2wannier90, only: numincbands! Number of included bands 
                                         !   for wannierization
  use m_siesta2wannier90, only: numkpoints ! Total number of k-points
                                         !   for which the overlap of the
                                         !   periodic part of the wavefunct
                                         !   with a neighbour k-point will
                                         !   be computed
  use m_siesta2wannier90, only: numproj  ! Total number of projection centers,
                                         !   equal to the number of MLWF
  use m_siesta2wannier90, only: Amnmat   ! Projections of a trial function
                                         !   with a Bloch orbital
                                         !   <\psi_{m k}|g_n>

  integer, intent(in) :: ispin

!
! Local variables
!

  character(len=len_trim(seedname)+4) :: amnfilename ! Name of the file where
                                                     !   the overlap Amn
                                                     !   matrices will be
                                                     !   written
  integer      :: ik        ! Counter for the loop on k-points
  integer      :: iproj     ! Counter for the loop on projectors
  integer      :: mband     ! Counter for the loop on bands
  integer      :: amnunit   ! Logical unit assigned to the file where the 
                            !   Amn matrices will be written
  integer      :: eof

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit


! Asign a logical unit to the file where the Amn matrices will be written
  call io_assign( amnunit )

! Asign a name to the file where the Amn matrices will be written
  amnfilename = trim( seedname )// ".amn"

! Open the output file where the Amn matrices will be written
  open( unit=amnunit, err=1992, file=amnfilename, status="replace", &
 &      iostat=eof)

! The first line of the amn file is a user comment
  write( unit=amnunit, fmt="(a38)", err=1992 )                      &
 &  "siesta2wannier90: Siesta for Wannier90"

! The second line includes three integers:
!   numincbands: the number of bands for wannierization
!   numkpoints:  the number of k-points for which the overlap of the
!      periodic part of the wavefunction with a neighbour k-point will be 
!      computed
!   numproj:     the number of projections
  write( unit=amnunit, fmt="(i5,1x,i5,1x,i2)", err=1992 )             &
 &  numincbands(ispin), numkpoints, numproj

! Subsequent numincbands x numproj x numkpoint lines:  
! three integers and two real numbers on each line
! The first two integers are the band indices m and n.
! The third integer specifies the \vec{k} by giving the ordinal 
! corresponding to its position in the list of k-points in seedname.win
! The real numbers are the real and the imaginary parts,
! respectively, of the actual $A_{m n} (\vec{k})$.

  do ik = 1, numkpoints
    do iproj = 1, numproj
      do mband = 1, numincbands(ispin)
        write(unit=amnunit,fmt="(3i5,1x,f12.5,2x,f12.5)",err=1992)      &
 &         mband, iproj, ik,                                           &
 &         real(Amnmat(mband,iproj,ik)),aimag(Amnmat(mband,iproj,ik))
      enddo
    enddo
  enddo

! Close the output file where the Amn matrices will be written
  call io_close(amnunit)

  write(6,'(/,a)')  &
 &  'amn: Overlap matrices between trial projection functions and wavefunctions'
  write(6,'(a)')  &
 &  'amn: written in ' // trim(seedname) // '.amn file'

  return

1992 call die("writeamn: Error writing to .amn file")

end subroutine writeamn

!
! ------------------------------------------------------------------------------
!

subroutine writeeig( ispin )
!
! In this subroutine, the Kohn-Sham eigenvalues $\epsilon_{n \vec{k}}$ 
! (in eV) at each point in the Monkhorst-Pack mesh are written
!
! This is required if any of disentanglement, plot_bands, plot_fermi_surface
! or hr_plot options are set to true in Wannier90.
! 
! NOTE: The actual extension of the file is .eigW to avoid name clashes
!       with Siesta's standard eigenvalue file in case-insensitive filesystems
!
! The structure of the outfile (seedname.eig) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname  ! Seed for the name of the file 
                                          !   where the Wannier90
                                          !   code, when used as postprocessing
                                          !   tool, dumps the information.
  use m_siesta2wannier90, only: numkpoints! Total number of k-points
                                          !   for which the overlap of the
                                          !   periodic part of the wavefunct
                                          !   with a neighbour k-point will
                                          !   be computed
  use m_siesta2wannier90, only: numincbands! Number of included bands 
                                          !   for wannierization
  use m_siesta2wannier90, only: eo        ! Eigenvalues of the Hamiltonian 
                                          !   at the numkpoints introduced in
                                          !   kpointsfrac 
                                          !   First  index: band index
                                          !   Second index: k-point index

  implicit none

  integer, intent(in)     :: ispin       ! Spin component

!
! Local variables
!

  character(len=len_trim(seedname)+5) :: eigfilename ! Name of the file where
                                                     !   the eigenvalues
                                                     !   written

  integer      :: ik        ! Counter for the loop on k-points
  integer      :: iband     ! Counter for the loop on bands
  integer      :: eigunit   ! Logical unit assigned to the file where the 
                            !   eigenvalues will be written
  integer      :: eof       ! Status of the output file

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit

! Asign a logical unit to the file where the eignvalues will be written
  call io_assign(eigunit)

! Asign a name to the file where the eigenvalues will be written
  eigfilename = trim(seedname)// ".eigW"

! Open the output file where the eigenvalues will be written
  open( unit=eigunit, err=1983, file=eigfilename, status="replace", &
 &      iostat=eof )

! Each line consists of two integers and a real number.
! The first integer is the band index, 
! The second integer gives the ordinal corresponding to the k-point in
!   the list of k-points in seedname.win.
! The real number is the eigenvalue (in eV).
  do ik = 1, numkpoints
    do iband = 1,numincbands(ispin)
      write( unit=eigunit, fmt="(i5,i5,3x,f12.5)", err=1983 )       &
 &      iband, ik, eo(iband,ik)
    enddo
  enddo

  write(6,'(/,a)')  &
 &  'eig: Eigenvalues of the Hamiltonian '
   write(6,'(a)')  &
 &  'eig: written in ' // trim(seedname) // '.eigW file'


! Close the output file where the eigenvalues will be written
  call io_close(eigunit)

  return

  1983 call die("writeeig: Error writing to .eigW file")
end subroutine writeeig

!
! ------------------------------------------------------------------------------
!

subroutine writeunk( ispin )
!
! Produces UNKXXXXX.X files which contain the periodic
! part of a Bloch function in the unit cell on a grid given by
! global unk_nx, unk_ny, unk_nz variables. 
! 
! The periodic part of the Bloch functions is defined by 
! \begin{equation}
!   u_{n \vec{k}} (\vec{r}) =
!   \sum_{\vec{R} \mu} c_{n \mu}(\vec{k})
!        e^{i \vec{k} \cdot ( \vec{r}_{\mu} + \vec{R} - \vec{r} )}
!        \phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} ) ,
!\end{equation}
!
!\noindent where $\phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} )$
! is an atomic orbital of the basis set centered on atom $\mu$ in
! the unit cell $\vec{R}$, and $c_{n \mu}(\vec{k})$ are the coefficients
! of the wave function, that must be identical that the ones used
! for wannierization in Mmn.
!
! Hamiltonian must be diagonalized
! for every k before the routine is invoked. 
! and the eigenvectors used must be the same in writeunk and in Mmn.
! Therefore, the diagonalizatinon routine cannot be called twice 
! for every k-point (an arbitrary phase factor can be added in 
! the different calls). For this reason, the subroutine writeunk
! is called with the coefficients of the wavefunctins computed in diagonHk
!
! Below is the paragraph of the Wannier90 User Guide relative to these files
! INPUT. Read if wannier_plot=.TRUE. and used to plot the MLWF. 
! Read if !transport_mode=lcr and tran_read_ht=.FALSE. for use in 
! automated lcr transportcalculations.
!
! The periodic part of the Bloch states represented on a
! regular real space grid, 
! indexed by k-point p (from 1 to num_kpts) and
!  spin s (‘1’ for ‘up’, ‘2’ for‘down’).
!
! The name of the wavefunction file is assumed to have the form:
!    write(wfnname,200) p,spin
! 200 format (’UNK’,i5.5,’.’,i1)
!
! The first line of each file should contain 5 integers: 
! the number of grid points in each direction (ngx, ngy and ngz), 
! the k-point number ik and 
! the total number of bands num_band in the file. 
!
! The full file will be read by wannier90 as:
! read(file_unit) ngx,ngy,ngz,ik,nbnd
! do loop_b=1,num_bands
!   read(file_unit) (r_wvfn(nx,loop_b),nx=1,ngx*ngy*ngz)
! end do
! The file can be in formatted or unformatted style, this is controlled by the
! logical keyword wvfn_formatted.
!

  use precision,          only: dp            ! Real double precision type
  use fdf                                     ! FDF library                
  use alloc,              only: re_alloc      ! Reallocation routines
  use alloc,              only: de_alloc      ! Deallocation routines
  use parallel,           only: IOnode        ! Input/output node
  use parallel,           only: Node          ! Working Node
  use parallel,           only: Nodes         ! Total number of nodes
  use neighbour,          only: maxnna        ! Maximum number of neighbours
  use neighbour,          only: jan           ! Atom-index of neighbours
  use neighbour,          only: r2ij          ! Squared distances to neighbors
  use neighbour,          only: xij           ! Vector from a given atom
                                              !   to its neighbours
  use neighbour,          only: x0            ! Position of the point around
                                              !   which we are going to compute
                                              !   the neighbours
  use siesta_geom,        only: xa            ! Atomic positions
  use siesta_geom,        only: na_u          ! Number of atoms in the unit cell
  use atomlist,           only: rmaxo         ! Maximum cutoff for atomic orb.
  use atomlist,           only: no_u          ! Number of orbitals in unit cell
  use m_siesta2wannier90, only: latvec        ! Lattice vectors in real 
                                              !   space

  use neighbour,          only: mneighb       ! Subroutine to compute the
                                              !   number of neighbours
  use m_siesta2wannier90, only: numincbands   ! Number of bands for 
                                              !   wannierization
                                              !   after excluding bands  
  use m_siesta2wannier90, only: nincbands_loc ! Number of bands for 
                                              !   wannierization
                                              !   after excluding bands  
                                              !   in the local node
  use m_siesta2wannier90, only: numkpoints    ! Total number of k-points
                                              !   for which the overlap of
                                              !   the periodic part of the
                                              !   wavefunct with a 
                                              !   neighbour k-point will
                                              !   be computed
  use m_siesta2wannier90, only: kpointsfrac   ! List of k points relative
                                              !   to the reciprocal 
                                              !   lattice vectors.
                                              !   First  index: component
                                              !   Second index: k-point  
                                              !      index in the list
  use m_siesta2wannier90, only: coeffs        ! Coefficients of the
                                              !   wavefunctions.
                                              !   First  index: orbital
                                              !   Second index: band
                                              !   Third  index: k-point
  use m_ntm,              only: ntm           ! Number of integration mesh
                                              !   divisions of each cell vector
#ifdef MPI
  use parallelsubs,       only: GetNodeOrbs      ! Calculates the number of
                                                 !   orbitals stored on the 
                                                 !   local Node.
  use parallelsubs,       only: LocalToGlobalOrb ! Converts an orbital index
                                                 !   in the local frame 
                                                 !   to the global frame
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
  use mpi_siesta
#endif


  implicit none

!
! Arguments of the subroutine
!
  integer,  intent(in) :: ispin      ! Spin component


!
! Variables related with the discretization of space
!
  integer      :: unk_nx, unk_nx_default        ! Number of points along x
                                !   to plot the periodic part 
                                !   of the Bloch function
  integer      :: unk_ny, unk_ny_default        ! Number of points along y
                                !   to plot the periodic part 
                                !   of the Bloch function
  integer      :: unk_nz, unk_nz_default        ! Number of points along z
                                !   to plot the periodic part 
                                !   of the Bloch function
  logical      :: unk_format, unk_format_default ! Are the UNK files written in 
                                !   binary (.true.) or ASCII (.false.)
!
  integer      :: nneig        ! Dummy variable for initialization call
!
! Variables related with the wave function
!
  real(dp)     :: kvector(3)   ! k-point vector for which the 
                               !   periodic part of the wave function
                               !   will be written in a file
                               !   between the projection function and the
                               !   eigenvector of the Hamiltonian will be
                               !   computed
  integer      :: nincbands    ! Number of included bands for wannierization

#ifdef MPI
  integer     :: MPIerror
  complex(dp), dimension(:,:,:,:), pointer :: auxloc => null()! Temporal array for the
                                             !   the global reduction of buffer
#endif

  complex(dp), pointer :: buffer(:,:,:,:) => null()  ! Variable where the 
                                             !   periodic part of the wave
                                             !   functions at the points of the
                                             !   grid will be stored
!
! Variables related with the input/output
!
  character(len=11) :: unkfilename ! Name of the file where the periodic
                                   !   part of the wave functions at the
                                   !   points of the grid will be saved
  integer      :: unkfileunit  ! Logical unit of the file
  integer      :: eof          ! Code error

!
! Counters
!
  integer      :: ik           ! Counter for the loop on kpoints
  integer      :: ix           ! Counter for the loop on points along x
  integer      :: iy           ! Counter for the loop on points along y
  integer      :: iz           ! Counter for the loop on points along z
  integer      :: iband        ! Counter for the loop on bands
  integer      :: iband_global ! Global index of the band
  integer      :: io           ! Counter for the loop on atomic orbitals

  integer, allocatable  :: iband_sequential(:)
  complex(dp), allocatable :: tmp_arr(:)

! Start time counter
  call timer('writeunk',1)

! Set number of bands for wannierization
  nincbands = numincbands( ispin )

! Set default real space discretisation same to the Siesta mesh
  unk_nx_default     = ntm(1)
  unk_ny_default     = ntm(2)
  unk_nz_default     = ntm(3)
  unk_format_default = .true. 

!
! Read FDF tags
!
  unk_nx     = fdf_get( 'Siesta2Wannier90.UnkGrid1',      unk_nx_default     )
  unk_ny     = fdf_get( 'Siesta2Wannier90.UnkGrid2',      unk_ny_default     )
  unk_nz     = fdf_get( 'Siesta2Wannier90.UnkGrid3',      unk_nz_default     )
  unk_format = fdf_get( 'Siesta2Wannier90.UnkGridBinary', unk_format_default )

    call re_alloc( buffer,              &
 &                 1, nincbands,        &
 &                 1, unk_nx,           &
 &                 1, unk_ny,           &
 &                 1, unk_nz,           &
 &                 name = 'buffer',     &
 &                 routine='writeunk' )

    allocate(tmp_arr(1:nincbands), iband_sequential(1:nincbands))

kpoints:                 &
  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )
!
!   Initialize neighbour subroutine.
!   The reallocation of the different arrays is done within neighb,
!   in the call to the subroutine sizeup_neighbour_arrays
!
    nneig = 0
    x0(:) = 0.0_dp
    call mneighb( latvec, rmaxo, na_u, xa, 0, 0, nneig )
  
!   Open the output file
    if( IOnode ) then
      write(unkfilename,"('UNK',i5.5,'.',i1)") ik, ispin
      call io_assign(unkfileunit)
      if( .not. unk_format ) then
        open( unit=unkfileunit, err=1992, file=unkfilename,          &
 &            status='replace', form='formatted', iostat=eof )
      else
        open( unit=unkfileunit, err=1992, file=unkfilename,          &
 &            status='replace', form='unformatted', iostat=eof )
      endif
    endif

!   Dump in the output file the number of points in the mesh, 
!   the index of the k-point and the number of occupied bands
    if( IOnode ) then
      if( .not. unk_format ) then
        write(unkfileunit,'(5(i5,2x))') unk_nx, unk_ny, unk_nz, ik, nincbands
      else
        write(unkfileunit) unk_nx, unk_ny, unk_nz, ik, nincbands
      endif
    endif


#ifndef MPI
!
! The serial version is easy enough
!
   do iz = 1, unk_nz
      do iy = 1, unk_ny
         do ix = 1, unk_nx
            call periodicpart(ix,iy,iz,nincbands,coeffs(:,1:nincbands,ik), &
                                       buffer(1:nincbands,ix,iy,iz))
         enddo
      enddo
   enddo

BAND_LOOP:  do iband = 1, nincbands

   if( .not. unk_format ) then
           
      do iz = 1, unk_nz
         do iy = 1, unk_ny
            do ix = 1, unk_nx
               write(unkfileunit,'(2f12.5)') &
                   real(buffer(iband,ix,iy,iz)), aimag(buffer(iband,ix,iy,iz))
            enddo 
         enddo 
      enddo  
   else
      write(unkfileunit) & 
           &          (((buffer(iband,ix,iy,iz),ix=1,unk_nx),iy=1,unk_ny),iz=1,unk_nz)
   endif

enddo  BAND_LOOP

#else
!
!   In MPI, the goal of writing directly the UNK files complicates
!   the logic and forces inefficiencies. This should be rewritten.
!   Each node should compute its bands and send the info to the master
!   node with the appropriate tags. The master node can then write
!   the information to a properly indexed netCDF file.
!
!   For each k-point, initialize the buffer in all the nodes
    buffer = cmplx(0.0_dp,0.0_dp,kind=dp)

    do iband = 1, nincbands_loc
       !     Identify the global index of the local band          
       iband_global = which_band_in_node(Node,iband)
       iband_sequential(iband) = sequential_index_included_bands(iband_global)
    enddo

   do iz = 1, unk_nz
      do iy = 1, unk_ny
         do ix = 1, unk_nx
            ! Compute all (local) bands at once
            call periodicpart(ix,iy,iz,nincbands_loc,coeffs(:,1:nincbands_loc,ik), &
                                       tmp_arr(1:nincbands_loc))
            
            ! Place in the correct global slot
            do iband = 1, nincbands_loc
               buffer(iband_sequential(iband),ix,iy,iz) = tmp_arr(iband)
            enddo
         enddo
      enddo
   enddo


!   A reduction is needed since we need all the matrix buffer in IOnode
!   (to dump it into a file), but the results for some of the bands might
!   be computed in other nodes.

!   Allocate workspace array for reduction
!   We assume that the IOnode is the root node...

    if (IOnode) then
       call re_alloc( auxloc, 1, nincbands, 1, unk_nx, 1, unk_ny, 1, unk_nz,  &
            &                 name='auxloc', routine='writeunk' )
       auxloc(:,:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    else
       auxloc => buffer(:,:,:,:)  ! Un-touched in the other nodes
                                  ! but has to be defined for checks
    endif

    call MPI_Reduce( buffer(1,1,1,1), auxloc(1,1,1,1),                  &
                     nincbands*unk_nx*unk_ny*unk_nz,                    &
                     MPI_double_complex,MPI_sum,0,MPI_Comm_World,MPIerror )

    if( IOnode ) then
      buffer(:,:,:,:) = auxloc(:,:,:,:)

      if( .not. unk_format ) then
        do iband = 1, nincbands
          do iz = 1, unk_nz
            do iy = 1, unk_ny
              do ix = 1, unk_nx
                write(unkfileunit,'(2f12.5)')   &
 &                real(buffer(iband,ix,iy,iz)), aimag(buffer(iband,ix,iy,iz))
              enddo ! Enddo in ix
            enddo ! Enddo in iy
          enddo ! Enddo in iz
        enddo ! Enddo in bands
      else
        do iband = 1, nincbands
          write(unkfileunit)                                   & 
 &          ((( buffer(iband,ix,iy,iz), ix = 1, unk_nx ),      &
 &                                      iy = 1, unk_ny ),      & 
 &                                      iz = 1, unk_nz )
        enddo 
      endif
    endif  ! IOnode
#endif

! Close the output file
    if( IOnode ) then
      call io_close(unkfileunit)
    endif

  enddo kpoints

! Deallocate some of the variables
#ifdef MPI
  if (IONode) then
     call de_alloc( auxloc, name='auxloc', routine='writeunk' )
  endif
#endif
  call de_alloc( buffer, name='buffer', routine='writeunk' )

  deallocate(tmp_arr,iband_sequential)

! End time counter
  call timer('writeunk',2)
  return

  1992 call die('swan: Error creating UNK files')

CONTAINS

  ! This subroutine will see the relevant variables by
  ! host association

  subroutine periodicpart(ix,iy,iz,nbands,psi,values) 

  use siesta_geom,        only: isa           ! Species index of each atom
  use atomlist,           only: lasto         ! Position of last orbital 
                                              !   of each atom
  use atomlist,           only: iphorb        ! Orbital index of each  orbital
                                              !   in its atom
  use atomlist,           only: indxuo        ! Index of equivalent orbital  
                                              !   in "u" cell
  use atmfuncs,           only: rcut          ! Function that determines the
                                              !   cutoff radius of a given
                                              !   orbital of a given specie
  use atmfuncs,           only: phiatm        ! Subroutine to compute the
                                              !   atomic orbital at a point

!-----------------------------------------------------------------
    integer, intent(in) :: ix, iy, iz
    integer, intent(in) :: nbands          ! Number of bands to process
    complex(dp), intent(in) :: psi(:,:)    ! Coefficients of wfns
    complex(dp), intent(out) ::  values(:) ! Output values


! Variables related with the mesh point
!
  real(dp)     :: rvector(3)    ! Coordinates of the mesh point in real space
  real(dp)     :: rvectorarg(3) ! Relative position of the orbital center with
                                !   respect to the real space point

!
! Variables related with the list of non-vanishing atoms at a given point
!
  integer      :: ineig        ! Counter for the loop on neighbours
  integer      :: iatom        ! Atomic index of the neighbour
  integer      :: ispecie      ! Atomic species of the neighbour
  integer      :: iso          ! Orbital index of each orbital in its atom
  integer      :: iorbital     ! Orbital index
  integer      :: iorbital0    ! Orbital index in the unit cell

  real(dp)     :: phi          ! Value of an atomic orbital at a point
  real(dp)     :: grphi(3)     ! Value of the gradient of an atomic orbital
                               !   at a given point
  real(dp)     :: phase        ! Phase of the exponential
!                  e^{i \vec{k} \cdot ( \vec{r}_{\mu} + \vec{R} - \vec{r} )}
  real(dp)     :: distance
  complex(dp)  :: exponential  ! Value of the previous exponential


            rvector(:) = ( latvec(:,1) * (ix-1) ) / unk_nx             &
 &                     + ( latvec(:,2) * (iy-1) ) / unk_ny             &
 &                     + ( latvec(:,3) * (iz-1) ) / unk_nz

            values(:) = 0.0_dp

!           Find atoms within rmaxo sphere and store them in jan
            x0(:) = rvector(:)
            call mneighb( latvec, rmaxo, na_u, xa, 0, 0, nneig )

            if ( nneig .gt. maxnna )                                     &
 &           call die('swan: insufficient array shapes; see neighb(..)')

!           Loop over atoms in the list of non-vanishing atoms at a given point
            do ineig = 1, nneig
!             Identify the atomic index.
!             In other subroutines, like overfsm, kinefsm, etc.
!             mneighb is called for the supercell,
!             with atom indices ranging 
!             from 1 to na_s (number atoms in supercell)
!             Here, mneighb is called for the unit cell, and the
!             indices of the atoms range 
!             from 1 to na_u (number atoms in unit cell
!             If different periodic images of the same atom do not vanish
!             at a given point, they appear in the list with the same
!             index iatom. This is irrelevant for the purpouses 
!             of this subroutine
              iatom   = jan(ineig)
              ispecie = isa(iatom)
!             The vector xij computed in mneighb has the origin at the
!             mesh point. However, in the argument of the orbital,
!             the origin is located at its center
!             \phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} )
!             We change the sign here.
              rvectorarg(:) = -xij(:,ineig)
!                 Compute the phase
!                 e^{i \vec{k} \cdot ( \vec{r}_{\mu} + \vec{R} - \vec{r} )}
!                 We have to change the sign of rvectorarg again.
                  phase = -1.0_dp * dot_product(kvector,rvectorarg)
                  exponential = exp( phase * cmplx(0.0_dp,1.0_dp,kind=dp) )

              distance = sqrt(r2ij(ineig))

!             Loop on all the orbitals of a given atom
              do iorbital = lasto(iatom-1)+1, lasto(iatom)
                iso = iphorb(iorbital)
!               If the point is within the range of the orbital
                if ( rcut(ispecie,iso) .gt. distance ) then
                  iorbital0 = indxuo(iorbital)
!                 Compute the value of the orbital at the point
!                 \phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} )
                  call phiatm( ispecie, iso, rvectorarg, phi, grphi )
!                 Compute the sum that gives the periodic part of the orbital
!                 for all bands
                  values(1:nbands) = values(1:nbands) + &
                                   exponential * psi(iorbital0,1:nbands) * phi
                endif
              enddo ! Enddo on orbitals that do not vanish
            enddo ! Enddo on atoms that have orbitals that do not vanish 
!
!           Transform Bohr^(-3/2) to Ang^(-3/2)
!
            values(1:nbands) =  2.59775721_dp * values(1:nbands)

          end subroutine periodicpart

end subroutine writeunk
