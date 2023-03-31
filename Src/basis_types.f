! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module basis_types
!
!=======================================================================
!
!     This module defines data structures to handle the specification
!     of the basis set and KB projectors in SIESTA. This specification
!     is read from the fdf file by routine 'basis_specs.basis_read', and
!     then converted to the form expected by routine 'atom.atom_main'
!
!     At present, for historical reasons, these data structures are 
!     different from those in module 'atm_types', but in principle could
!     be made to contain the same information (except for the mapping
!     of 'nl' to 'nlm' orbitals and projectors and the polarization
!     orbitals) by the use of the 'rad_func' pointers indicated within 
!     types 'shell_t' and 'kbshell_t'. (See module 'radial')
!     
!
      use atmparams, only: lmaxd, nzetmx, nsemx, nkbmx
      use pseudopotential, only: pseudopotential_t, pseudo_init_constant
      use precision, only: dp
      use sys, only : die

      Implicit None

      type, public ::  ground_state_t
          integer                   ::  lmax_valence
          integer                   ::  n(0:3)
          real(dp)                  ::  occupation(0:3)
          logical                   ::  occupied(0:4)   ! note 0..4
          real(dp)                  ::  z_valence
      end type ground_state_t

      type, public :: shell_t
          integer                   ::  n          ! n quantum number
          integer                   ::  l          ! angular momentum
          integer                   ::  nzeta      ! Number of PAOs
          logical                   ::  polarized  
          integer                   ::  nzeta_pol
          real(dp)                  ::  split_norm ! Split norm value
          logical                   ::  split_norm_specified ! in a
                                                             ! S value
                                                             ! construct
          real(dp)                  ::  rinn       ! Soft confinement
          real(dp)                  ::  vcte       ! Soft confinement
          real(dp)                  ::  filtercut  ! Filter cutoff
          real(dp)                  ::  qcoe       ! Charge confinement
          real(dp)                  ::  qyuk       ! Charge confinement
          real(dp)                  ::  qwid       ! Charge confinement
          real(dp), pointer         ::  rc(:) => null()! rc's for PAOs
          real(dp), pointer         ::  lambda(:) => null() ! Contraction factors
          !!! type(rad_func), pointer   ::  orb(:) ! Actual orbitals 
      end type shell_t

      type, public :: lshell_t
          integer                   ::  l          ! angular momentum
          integer                   ::  nn         ! number of n's for this l
          type(shell_t), pointer    ::  shell(:) => null() ! One shell for each n
      end type lshell_t

      type, public :: kbshell_t
          integer                   ::  l          ! angular momentum
          integer                   ::  nkbl       ! No. of projs for this l
          real(dp), pointer         ::  erefkb(:) => null() ! Reference energies
          !!! type(rad_func), pointer  ::  proj(:) => null() ! Actual projectors
      end type kbshell_t

      type, public :: dftushell_t
          integer                   ::  n          ! principal quantum number
                                                   !   of the atomic orbital
                                                   !   where the U correction
                                                   !   will be applied
          integer                   ::  l          ! angular quantum number 
                                                   !   of the atomic orbital
                                                   !   where the U correction
                                                   !   will be applied
          real(dp)                  ::  rinn       ! Soft confinement 
                                                   !   inner radius
          real(dp)                  ::  vcte       ! Soft confinement potential
                                                   !  prefactor of the potential
          real(dp)                  ::  rc         ! rc's for DFT+U projectors
          integer                   ::  nrc        ! Point in the log grid where
                                                   !  the DFT+U proj. vanishes
          real(dp)                  ::  lambda     ! Contraction factors
          real(dp)                  ::  dnrm_rc    ! Parameter used to determine
                                                   !   the cutoff radius of the
                                                   !   Fermi distrib. used 
                                                   !   to cut the DFT+U proj.
          real(dp)                  ::  width      ! Width of the Fermi distrib.
                                                   !   to cut the DFT+U proj.
          real(dp)                  ::  u          ! Value of the U parameter
          real(dp)                  ::  j          ! Value of the J parameter
          !!! type(rad_func), pointer  ::  dftu_proj(:) => null() ! Actual projectors
                                                   !   all these radial function
                                                   !   are now defined in the 
                                                   !   derived type "species"
                                                   !   in module atm_types
      end type dftushell_t
!
!     Main data structure
!
      type, public :: basis_def_t
          character(len=20)         ::  label      ! Long label
          integer                   ::  z          ! Atomic number
          type(ground_state_t)      ::  ground_state
          type(pseudopotential_t)   ::  pseudopotential
          integer                   ::  lmxo       ! Max l for basis
          integer                   ::  lmxkb      ! Max l for KB projs
          integer                   ::  lmxdftupj  ! Max l for DFT+U projs
          type(lshell_t), pointer   ::  lshell(:) => null() ! One shell per l 
          type(kbshell_t), pointer  ::  kbshell(:) => null() ! One KB shell per l
          real(dp)                  ::  ionic_charge
          real(dp)                  ::  mass   
          !
          ! The rest of the components are auxiliary
          ! 
          logical                   ::  floating   
          logical                   ::  bessel
          logical                   ::  synthetic
          character(len=20)         ::  basis_type
          character(len=20)         ::  basis_size
          logical                   ::  semic      ! 
          integer                   ::  nshells_tmp
          integer                   ::  nkbshells
          integer                   ::  ndftushells      ! For a given atomic
                                                         !  species, on how many
                                                         !  orbitals we are 
                                                         !  going to apply the
                                                         !  U correction
          integer                   ::  ndftuprojs_lm    ! How many projectors
                                                         !  including angular 
                                                         !  dependencies.
          integer                   ::  lmxkb_requested
          integer                   ::  lmxdftupj_requested
          type(shell_t), pointer    ::  tmp_shell(:) => null()
          type(dftushell_t), pointer::  dftushell(:) => null()
      end type basis_def_t

      integer, save, public              :: nsp  ! Number of species
      type(basis_def_t), public,
     $     allocatable, save, target     :: basis_parameters(:)

!=====================================================================
!     OLD ARRAYS
!=====================================================================
!
      logical      ,save, public, pointer :: semic(:) => null()
      integer      ,save, public, pointer :: lmxkb(:) => null()
      integer      ,save, public, pointer :: lmxo(:) => null()
      integer      ,save, public, pointer :: nsemic(:,:) => null()
      integer      ,save, public, pointer :: nkbl(:,:) => null()
      integer      ,save, public, pointer :: cnfigmx(:,:) => null()
      integer      ,save, public, pointer :: polorb(:,:,:) => null()
      integer      ,save, public, pointer :: nzeta(:,:,:) => null()
      real(dp)     ,save, public, pointer :: split_norm(:,:,:) => null()
      real(dp)     ,save, public, pointer :: vcte(:,:,:) => null()
      real(dp)     ,save, public, pointer :: rinn(:,:,:) => null()
      real(dp)     ,save, public, pointer :: qcoe(:,:,:) => null()
      real(dp)     ,save, public, pointer :: qyuk(:,:,:) => null()
      real(dp)     ,save, public, pointer :: qwid(:,:,:) => null()
      real(dp)     ,save, public, pointer :: erefkb(:,:,:) => null()
      real(dp)     ,save, public, pointer :: charge(:) => null()
      real(dp)     ,save, public, pointer :: lambda(:,:,:,:) => null()
      real(dp)     ,save, public, pointer :: filtercut(:,:,:) => null()
      real(dp)     ,save, public, pointer :: rco(:,:,:,:) => null()
      integer      ,save, public, pointer :: iz(:) => null()
      real(dp)     ,save, public, pointer :: smass(:) => null()
      character(len=10), save, public, pointer :: basistype(:) => null()
      character(len=20), save, public, pointer :: atm_label(:) => null()

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      interface destroy
        module procedure destroy_shell, 
     $                   destroy_lshell, destroy_basis_def
      end interface
      interface initialize
        module procedure init_shell, init_kbshell, init_dftushell,
     $                   init_lshell, init_basis_def
      end interface

!---------------------------------------------------------
      public  :: destroy, copy_shell, initialize
      public  :: write_basis_specs, basis_specs_transfer
      public  :: deallocate_spec_arrays
      public  :: print_dftushell
!---------------------------------------------------------

      PRIVATE

      CONTAINS

!-----------------------------------------------------------------------
      subroutine copy_shell(source,target)
!
!     This is a "deep-copy" of a structure, including the 
!     *allocated memory* associated to the internal pointer components.
!
      type(shell_t), intent(in)          :: source
      type(shell_t), intent(out)         :: target

      target%l = source%l
      target%n = source%n
      target%nzeta = source%nzeta
      target%polarized = source%polarized
      target%nzeta_pol = source%nzeta_pol
      target%rinn = source%rinn
      target%vcte = source%vcte
      target%qcoe = source%qcoe
      target%qyuk = source%qyuk
      target%qwid = source%qwid
      target%split_norm = source%split_norm
      target%filtercut = source%filtercut

      allocate(target%rc(1:size(source%rc)))
      allocate(target%lambda(1:size(source%lambda)))
      target%rc(:) = source%rc(:)
      target%lambda(:) = source%lambda(:)
      end subroutine copy_shell
!-----------------------------------------------------------------------

      subroutine init_shell(p)
      type(shell_t)          :: p

      p%l = -1
      p%n = -1
      p%nzeta = 0
      p%polarized = .false.
      p%nzeta_pol = 0
      p%rinn = 0._dp
      p%vcte = 0._dp
      p%qcoe = 0._dp
      p%qyuk = 0._dp
      p%qwid = 0.01_dp
      p%split_norm = 0.0_dp
      p%split_norm_specified = .false.
      p%filtercut = 0.0_dp
      nullify(p%rc,p%lambda)
      end subroutine init_shell


!-----------------------------------------------------------------------
      subroutine init_kbshell(p)
      type(kbshell_t)          :: p
      p%l = -1
      p%nkbl = -1
      nullify(p%erefkb)
      end subroutine init_kbshell

!-----------------------------------------------------------------------
      subroutine init_dftushell(p)
      type(dftushell_t)          :: p
      p%n       = -1
      p%l       = -1
      p%rinn    = 0.0_dp
      p%vcte    = 0.0_dp
      p%u       = 0.0_dp
      p%j       = 0.0_dp
      p%rc      = 0.0_dp
      p%nrc     = 0
      p%lambda  = 1.0_dp
      p%dnrm_rc = 0.9_dp
      p%width   = 0.5_dp
      end subroutine init_dftushell

!-----------------------------------------------------------------------
      subroutine init_lshell(p)
      type(lshell_t)          :: p

      p%l = -1
      p%nn = 0
      nullify(p%shell)
      end subroutine init_lshell

!-----------------------------------------------------------------------
      subroutine init_basis_def(p)
      type(basis_def_t)          :: p

      p%lmxo = -1
      p%lmxkb = -1
      p%lmxdftupj = -1
      p%lmxkb_requested = -1
      p%lmxdftupj_requested = -1
      p%nkbshells = -1
      p%ndftushells = -1
      p%ndftuprojs_lm = -1
      p%nshells_tmp = -1
      p%label = 'Unknown'
      p%semic = .false.
      p%ionic_charge = huge(1.0_dp) ! To signal it was not set
      call pseudo_init_constant(p%pseudopotential)
      nullify(p%lshell)
      nullify(p%kbshell)
      nullify(p%tmp_shell)
      nullify(p%dftushell)
      end subroutine init_basis_def

!-----------------------------------------------------------------------
      subroutine destroy_shell(p)
      type(shell_t), pointer   :: p(:)

      integer i
      type(shell_t), pointer   :: q

      if (.not. associated(p)) return
      do i = 1, size(p)
         q=>p(i)
         deallocate(q%rc)
         deallocate(q%lambda)
      enddo
      deallocate(p)
      end subroutine destroy_shell

!-----------------------------------------------------------------------
      subroutine destroy_lshell(p)
      type(lshell_t), pointer   :: p(:)

      integer i
      type(lshell_t), pointer   :: q

      if (.not. associated(p)) return
      do i = 1, size(p)
         q=>p(i)
         call destroy_shell(q%shell)
      enddo
      deallocate(p)
      end subroutine destroy_lshell

      
      subroutine destroy_dftushell(p)
      type(dftushell_t), pointer   :: p(:)

      if (.not. associated(p)) return
      deallocate(p)
      end subroutine destroy_dftushell

!-----------------------------------------------------------------------
      subroutine destroy_basis_def(p)
      type(basis_def_t)          :: p

      call destroy_dftushell(p%dftushell)
      call destroy_lshell(p%lshell)
      call destroy_shell(p%tmp_shell)

      end subroutine destroy_basis_def

!-----------------------------------------------------------------------
      subroutine print_shell(p)
      type(shell_t)            :: p

      integer i

      write(6,*) 'SHELL-------------------------'
      write(6,'(5x,a20,i20)') 'Angular momentum',     p%l
      write(6,'(5x,a20,i20)') 'n quantum number',     p%n
      write(6,'(5x,a20,i20)') 'Nzeta'           ,     p%nzeta
      write(6,'(5x,a20,l20)') 'Polarized?       ',    p%polarized
      write(6,'(5x,a20,i20)') 'Nzeta pol'           , p%nzeta_pol
      write(6,'(5x,a20,g20.10)') 'split_norm'     , p%split_norm
      write(6,'(5x,a20,g20.10)') 'filter cutoff'  , p%filtercut
      write(6,'(5x,a20,g20.10)') 'rinn'           , p%rinn
      write(6,'(5x,a20,g20.10)') 'vcte'           , p%vcte
      write(6,'(5x,a20,g20.10)') 'qcoe'           , p%qcoe
      write(6,'(5x,a20,g20.10)') 'qyuk'           , p%qyuk
      write(6,'(5x,a20,g20.10)') 'qwid'           , p%qwid
      write(6,'(5x,a)') 'rc and lambda for each nzeta:'
      do i = 1, p%nzeta
         write(6,'(5x,i2,2x,2g20.10)') i, p%rc(i), p%lambda(i)
      enddo
      write(6,*) '--------------------SHELL'

      end subroutine print_shell

!-----------------------------------------------------------------------
      subroutine print_kbshell(p)
      type(kbshell_t)            :: p

      integer i

      write(6,*) 'KBSHELL-------'
      write(6,'(5x,a20,i20)') 'Angular momentum',     p%l
      write(6,'(5x,a20,i20)') 'number of projs',     p%nkbl
      write(6,'(5x,a)') 'ref energy for each proj:'
      do i = 1, p%nkbl
         write(6,'(5x,i2,2x,g20.5)') i, p%erefkb(i)
      enddo
      write(6,*) '---------------------KBSHELL'

      end subroutine print_kbshell

!-----------------------------------------------------------------------
      subroutine print_dftushell(p)
      type(dftushell_t)            :: p

      write(6,*) 'DFTUSHELL-------'
      write(6,'(5x,a25,i20)')   'Principal quantum number',  p%n
      write(6,'(5x,a25,i20)')   'Angular momentum',          p%l
      write(6,'(5x,a25,g20.5)') 'U parameter:', p%u
      write(6,'(5x,a25,g20.5)') 'J parameter:', p%j
      write(6,'(5x,a25,g20.5)') 'rinn:',        p%rinn
      write(6,'(5x,a25,g20.5)') 'vcte:',        p%vcte
      write(6,'(5x,a25,g20.5)') 'lambda:',      p%lambda
      write(6,'(5x,a25,g20.5)') 'width:',       p%width
      write(6,'(5x,a25,g20.5)') 'dnrm_rc:',     p%dnrm_rc
      write(6,'(5x,a25,g20.5)') 'rc:',          p%rc
      write(6,'(5x,a25,i10)')   'nrc:',         p%nrc
      write(6,*) '---------------------DFTUSHELL'

      end subroutine print_dftushell

!-----------------------------------------------------------------------
      subroutine print_lshell(p)
      type(lshell_t)            :: p

      integer i

      write(6,*) 'LSHELL --------------'
      write(6,'(5x,a20,i20)') 'Angular momentum',     p%l
      write(6,'(5x,a20,i20)') 'Number of n shells',     p%nn
      if (.not. associated(p%shell)) return
      do i=1, p%nn
         call print_shell(p%shell(i))
      enddo
      write(6,*) '--------------LSHELL'
      end subroutine print_lshell

!-----------------------------------------------------------------------
      subroutine print_basis_def(p)
      type(basis_def_t)            :: p

      integer :: i

      write(6,*)
      write(6,*) 'SPECIES---'

      write(6,'(5x,a20,a20)') 'label',           p%label
      write(6,'(5x,a20,i20)') 'atomic number',   p%z
      write(6,'(5x,a20,a20)') 'basis type',      p%basis_type
      write(6,'(5x,a20,a20)') 'basis size',      p%basis_size
      write(6,'(5x,a20,g20.10)') 'ionic charge', p%ionic_charge
      write(6,'(5x,a20,i20)') 'lmax basis',      p%lmxo

      if ( associated(p%lshell) ) then
         do i = 0 , p%lmxo
            call print_lshell(p%lshell(i))
         end do
      else
         write(6,*) 'No L SHELLS, lmxo=', p%lmxo
      end if
      if ( associated(p%kbshell) ) then
         do i = 0 , p%lmxkb
            call print_kbshell(p%kbshell(i))
         end do
      else
         write(6,*) 'No KB SHELLS, lmxkb=', p%lmxkb
      end if
      if ( associated(p%dftushell) ) then
         do i=1, p%ndftushells
            call print_dftushell(p%dftushell(i))
         end do
      else
         write(6,*) 'No DFT+U PROJECTORS, lmxdftupj=', p%ndftushells
      end if

      write(6,*) '------------SPECIES'
      write(6,*) 

      end subroutine print_basis_def

!-----------------------------------------------------------------------
      subroutine basis_specs_transfer()
      use alloc, only: re_alloc

      integer lmax, lmaxkb, nzeta_max, nsemi_max, nkb_max
      integer l, isp, n, inz

      type(basis_def_t), pointer::   basp
      type(shell_t), pointer::  s
      type(lshell_t), pointer::  ls
      type(kbshell_t), pointer:: k

!
!     Check dimensions
!
      lmax = -1
      lmaxkb = -1
      nzeta_max = 0
      nsemi_max = 0
      nkb_max = 0

      do isp=1,nsp

         basp=>basis_parameters(isp)

         lmax = max(lmax,basp%lmxo)
         lmaxkb = max(lmaxkb,basp%lmxkb)
         do l=0,basp%lmxo
            ls=>basp%lshell(l)
            nsemi_max = max(nsemi_max,ls%nn)
            do n=1,ls%nn
               s=>ls%shell(n)
               nzeta_max = max(nzeta_max,s%nzeta)
            enddo
         enddo
         do l=0,basp%lmxkb
            k=>basp%kbshell(l)
            nkb_max = max(nkb_max,k%nkbl)
         enddo

      enddo

      lmax = max(lmax,lmaxkb)
      if (lmax .gt. lmaxd) then
         write(6,*) "Increment lmaxd to ", lmax
         call die()
      endif
      if (nzeta_max .gt. nzetmx) then
         write(6,*) "Increment nzetmx to ", nzeta_max
         call die()
      endif
      if (nsemi_max .gt. nsemx) then
         write(6,*) "Increment nsemx to ", nsemi_max
         call die()
      endif
      if (nkb_max .gt. nkbmx) then
         write(6,*) "Increment nkbmx to ", nkb_max
         call die()
      endif
!
!     ALLOCATE old arrrays
!
      nullify( semic )
      call re_alloc( semic, 1, nsp, 'semic', 'basis_types' )
      nullify( lmxkb )
      call re_alloc( lmxkb, 1, nsp, 'lmxkb', 'basis_types' )
      nullify( lmxo )
      call re_alloc( lmxo, 1, nsp, 'lmxo', 'basis_types' )
      nullify( nsemic )
      call re_alloc( nsemic, 0, lmaxd, 1, nsp, 'nsemic', 'basis_types' )
      nullify( cnfigmx )
      call re_alloc( cnfigmx, 0, lmaxd, 1, nsp,
     &               'cnfigmx',  'basis_types' )
      nullify( nkbl )
      call re_alloc( nkbl, 0, lmaxd, 1, nsp, 'nkbl', 'basis_types' )
      nullify( polorb )
      call re_alloc( polorb, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'polorb', 'basis_types' )
      nullify( nzeta )
      call re_alloc( nzeta, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'nzeta', 'basis_types' )
      nullify( split_norm )
      call re_alloc( split_norm, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'split_norm', 'basis_types' )
      nullify( filtercut )
      call re_alloc( filtercut, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'filtercut', 'basis_types' )
      nullify( vcte )
      call re_alloc( vcte, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'vcte', 'basis_types' )
      nullify( rinn )
      call re_alloc( rinn, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'rinn', 'basis_types' )
      nullify( qcoe )
      call re_alloc( qcoe, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'qcoe', 'basis_types' )
      nullify( qyuk )
      call re_alloc( qyuk, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'qyuk', 'basis_types' )
      nullify( qwid )
      call re_alloc( qwid, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'qwid', 'basis_types' )
      nullify( erefkb )
      call re_alloc( erefkb, 1, nkbmx, 0, lmaxd, 1, nsp,
     &               'erefkb', 'basis_types' )
      nullify( charge )
      call re_alloc( charge, 1, nsp, 'charge', 'basis_types' )
      nullify( lambda )
      call re_alloc( lambda, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'lambda', 'basis_types' )
      nullify( rco )
      call re_alloc( rco, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsp,
     &               'rco', 'basis_types' )
      nullify( iz )
      call re_alloc( iz, 1, nsp, 'iz', 'basis_types' )
      nullify( smass )
      call re_alloc( smass, 1, nsp, 'smass', 'basis_types' )
      nullify( basistype )
      allocate(basistype(nsp))
!      call re_alloc( basistype, 1, nsp, 'basistype', 'basis_types' )
      nullify( atm_label )
      allocate(atm_label(nsp))
!      call re_alloc( atm_label, 1, nsp, 'atm_label', 'basis_types' )

!
!     Transfer
!
      nkbl(:,:) = 0
      nzeta(:,:,:) = 0
      split_norm(:,:,:) = 0._dp
      filtercut(:,:,:) = 0._dp
      vcte(:,:,:) = 0._dp
      rinn(:,:,:) = 0._dp
      qcoe(:,:,:) = 0._dp
      qyuk(:,:,:) = 0._dp
      qwid(:,:,:) = 0.01_dp
      polorb(:,:,:) = 0
      rco(:,:,:,:) = 0._dp
      lambda(:,:,:,:) = 0._dp
      erefkb(:,:,:) = 0._dp
      semic(:) = .false.
      nsemic(:,:) = 0
      cnfigmx(:,:) = 0
      
      do isp=1,nsp

         basp=>basis_parameters(isp)

         iz(isp) = basp%z
         lmxkb(isp) = basp%lmxkb
         lmxo(isp) = basp%lmxo
         atm_label(isp) = basp%label
         charge(isp) = basp%ionic_charge
         smass(isp) = basp%mass
         basistype(isp) = basp%basis_type
         semic(isp) = basp%semic

         do l=0,basp%lmxo
            ls=>basp%lshell(l)
!
!           If there are no "non-polarizing" orbitals for
!           this l, nn=0. Set nsemic to 0 in that case.
!           (Kludge for now until future reorganization)
!
            nsemic(l,isp) = max(ls%nn -1 ,0)
            cnfigmx(l,isp) = 0
            do n=1,ls%nn
               s=>ls%shell(n)
               cnfigmx(l,isp) = max(cnfigmx(l,isp),s%n)
               nzeta(l,n,isp) = s%nzeta
               polorb(l,n,isp) = s%nzeta_pol
               split_norm(l,n,isp) = s%split_norm
               filtercut(l,n,isp) = s%filtercut
               vcte(l,n,isp) = s%vcte
               rinn(l,n,isp) = s%rinn
               qcoe(l,n,isp) = s%qcoe
               qyuk(l,n,isp) = s%qyuk
               qwid(l,n,isp) = s%qwid
!
!              This would make the code act in the same way as
!              siesta 0.X, but it does not seem to be necessary...
!
!               if (s%nzeta_pol.ne.0) then
!                  vcte(l+1,n,isp) = s%vcte
!                  rinn(l+1,n,isp) = s%rinn
!                  qcoe(l+1,n,isp) = s%qcoe
!                  qyuk(l+1,n,isp) = s%qyuk
!                  qwid(l+1,n,isp) = s%qwid
!               endif

!
!              Avoid referencing s%rc and s%lambda for
!              polarization orbitals.
!
               do inz = 1, s%nzeta
                  rco(inz,l,n,isp) = s%rc(inz)
                  lambda(inz,l,n,isp) = s%lambda(inz)
               enddo
            enddo
!
!           Fix for l's without PAOs
!
            if (cnfigmx(l,isp).eq.0)
     $           cnfigmx(l,isp) = basp%ground_state%n(l)

         enddo

         ! NOTE: cnfigmx and nsemic are only initialized for l up to lmxo
         !       in the above loop
         !       Extend them so that we can deal properly with outer polarization states

         do l=basp%lmxo+1, lmaxd
            !     gs is only setup up to l=3 (f)
            if (l <= 3) then
               cnfigmx(l,isp) = basp%ground_state%n(l)
            else
               ! g orbitals. Use the "l+1" heuristic
               ! For example, 5g pol orb associated to a 4f orb.
               cnfigmx(l,isp) = l + 1
            endif
            nsemic(l,isp) = 0
         enddo
         

         do l=0,basp%lmxkb
            k=>basp%kbshell(l)
            nkbl(l,isp) = k%nkbl
            erefkb(1:k%nkbl,l,isp) = k%erefkb(1:k%nkbl)
         enddo

      enddo

      end subroutine basis_specs_transfer

!-----------------------------------------------------------------------
      subroutine write_basis_specs(lun,is)
      integer, intent(in) :: lun
      integer, intent(in) :: is

      ! Pointer to basis specification
      type(basis_def_t), pointer :: basp
      type(dftushell_t), pointer :: dftu

      integer :: l, n, i
      integer :: nprin
      character(len=4) :: orb_id
      character(len=1), parameter   ::
     $                           sym(0:4) = (/ 's','p','d','f','g' /)

      basp => basis_parameters(is)

      write(lun,'(/a/79("="))') '<basis_specs>'
      write(lun,'(a20,1x,a2,i4,4x,a5,g12.5,4x,a7,g12.5)')
     $     atm_label(is), 'Z=',iz(is),
     $     'Mass=', smass(is), 'Charge=', charge(is)
      ! Allow a 2-char width for lmxkb (=-1 for floating and bessel orbs)
      write(lun,'(a5,i1,1x,a6,i2,4x,a10,a10,1x,a6,l1)')
     $     'Lmxo=', lmxo(is), 'Lmxkb=', lmxkb(is),
     $     'BasisType=', basistype(is), 'Semic=', semic(is)
      do l=0,lmxo(is)
         write(lun,'(a2,i1,2x,a7,i1,2x,a8,i1)')
     $        'L=', l, 'Nsemic=', nsemic(l,is),
     $        'Cnfigmx=', cnfigmx(l,is)
         do n=1,nsemic(l,is)+1
            if (nzeta(l,n,is) == 0) exit
            nprin = cnfigmx(l,is) - nsemic(l,is) + n - 1
            write(orb_id,"(a1,i1,a1,a1)") "(",nprin, sym(l), ")"
            write(lun,'(10x,a2,i1,2x,a6,i1,2x,a7,i1,2x,a4)')
     $            'i=', n, 'nzeta=',nzeta(l,n,is),
     $            'polorb=', polorb(l,n,is), orb_id
            if (basistype(is).eq.'filteret') then
               write(lun,'(10x,a10,2x,g12.5)') 
     $              'fcutoff:', filtercut(l,n,is)
            else
               write(lun,'(10x,a10,2x,g12.5)') 
     $              'splnorm:', split_norm(l,n,is)
            end if
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'vcte:', vcte(l,n,is)
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'rinn:', rinn(l,n,is)
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'qcoe:', qcoe(l,n,is)
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'qyuk:', qyuk(l,n,is)
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'qwid:', qwid(l,n,is)
            write(lun,'(10x,a10,2x,4g12.5)') 'rcs:',
     $           (rco(i,l,n,is),i=1,min(4,nzeta(l,n,is)))
            write(lun,'(10x,a10,2x,4g12.5)') 'lambdas:',
     $           (lambda(i,l,n,is),i=1,min(4,nzeta(l,n,is)))
         end do
      end do
      if ( lmxkb(is) > 0 ) then
         write(lun,'(79("-"))')
         do l=0,lmxkb(is)
            write(lun,'(a2,i1,2x,a5,i1,2x,a6,4g12.5)')
     $           'L=', l, 'Nkbl=', nkbl(l,is),
     $           'erefs:  ', (erefkb(i,l,is),i=1,nkbl(l,is))
         end do
      end if
      if ( associated(basp%dftushell) ) then
         write(lun,'(79("-"))')
         do l = 1 , basp%ndftushells
            dftu => basp%dftushell(l)
            write(lun,'(a2,i1,2x,a2,i1)')
     $           'L=', dftu%l, 'n=', dftu%n
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'U:', dftu%U
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'J:', dftu%J
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'rinn:', dftu%rinn
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'vcte:', dftu%vcte
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'lambda:', dftu%lambda
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'width:', dftu%width
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'rc:', dftu%rc
            write(lun,'(10x,a10,2x,g12.5)') 
     $           'dnrm_rc:', dftu%dnrm_rc
         end do
      end if
      write(lun,'(79("="))')
      write(lun,'(a/)') '</basis_specs>'

      end subroutine write_basis_specs

      subroutine deallocate_spec_arrays()
!
      use alloc, only: de_alloc

      call de_alloc( semic,      'semic',      'basis_types' )
      call de_alloc( lmxkb,      'lmxkb',      'basis_types' )
      call de_alloc( lmxo,       'lmxo',       'basis_types' )
      call de_alloc( nsemic,     'nsemic',     'basis_types' )
      call de_alloc( cnfigmx,    'cnfigmx',    'basis_types' )
      call de_alloc( nkbl,       'nkbl',       'basis_types' )
      call de_alloc( polorb,     'polorb',     'basis_types' )
      call de_alloc( nzeta,      'nzeta',      'basis_types' )
      call de_alloc( split_norm, 'split_norm', 'basis_types' )
      call de_alloc( filtercut,  'filtercut',  'basis_types' )
      call de_alloc( vcte,       'vcte',       'basis_types' )
      call de_alloc( rinn,       'rinn',       'basis_types' )
      call de_alloc( qcoe,       'qcoe',       'basis_types' )
      call de_alloc( qyuk,       'qyuk',       'basis_types' )
      call de_alloc( qwid,       'qwid',       'basis_types' )
      call de_alloc( erefkb,     'erefkb',     'basis_types' )
      call de_alloc( charge,     'charge',     'basis_types' )
      call de_alloc( lambda,     'lambda',     'basis_types' )
      call de_alloc( rco,        'rco',        'basis_types' )
      call de_alloc( iz,         'iz',         'basis_types' )
      call de_alloc( smass,      'smass',      'basis_types' )
      deallocate( basistype )
      deallocate( atm_label )
!      call de_alloc( basistype, 'basistype', 'basis_types' )
!      call de_alloc( atm_label, 'atm_label', 'basis_types' )
!
      end subroutine deallocate_spec_arrays

!-----------------------------------------------------------------------
      end module basis_types
