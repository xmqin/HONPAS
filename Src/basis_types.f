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
      use pseudopotential, only: pseudopotential_t
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
          real(dp), pointer         ::  rc(:)      ! rc's for PAOs
          real(dp), pointer         ::  lambda(:)  ! Contraction factors
          !!! type(rad_func), pointer   ::  orb(:) ! Actual orbitals 
      end type shell_t

      type, public :: lshell_t
          integer                   ::  l          ! angular momentum
          integer                   ::  nn         ! number of n's for this l
          type(shell_t), pointer    ::  shell(:)   ! One shell for each n
      end type lshell_t

      type, public :: kbshell_t
          integer                   ::  l          ! angular momentum
          integer                   ::  nkbl       ! No. of projs for this l
          real(dp), pointer         ::  erefkb(:)  ! Reference energies
          !!! type(rad_func), pointer  ::  proj(:) ! Actual projectors
      end type kbshell_t
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
          type(lshell_t), pointer   ::  lshell(:)  ! One shell per l 
          type(kbshell_t), pointer  ::  kbshell(:) ! One KB shell per l
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
          integer                   ::  lmxkb_requested
          type(shell_t), pointer    ::  tmp_shell(:)
      end type basis_def_t

      integer, save, public              :: nsp  ! Number of species
      type(basis_def_t), public,
     $     allocatable, save, target     :: basis_parameters(:)

!=====================================================================
!     OLD ARRAYS
!=====================================================================
!
      logical      ,save, public, pointer :: semic(:)
      integer      ,save, public, pointer :: lmxkb(:), lmxo(:)
      integer      ,save, public, pointer :: nsemic(:,:), nkbl(:,:)
      integer      ,save, public, pointer :: cnfigmx(:,:)
      integer      ,save, public, pointer :: polorb(:,:,:)
      integer      ,save, public, pointer :: nzeta(:,:,:)
      real(dp)     ,save, public, pointer :: split_norm(:,:,:)
      real(dp)     ,save, public, pointer :: vcte(:,:,:)
      real(dp)     ,save, public, pointer :: rinn(:,:,:)
      real(dp)     ,save, public, pointer :: erefkb(:,:,:)
      real(dp)     ,save, public, pointer :: charge(:)
      real(dp)     ,save, public, pointer :: lambda(:,:,:,:)
      real(dp)     ,save, public, pointer :: rco(:,:,:,:)
      integer      ,save, public, pointer :: iz(:)
      real(dp)     ,save, public, pointer :: smass(:)
      character(len=10), save, public, pointer :: basistype(:)
      character(len=20), save, public, pointer :: atm_label(:)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      interface destroy
        module procedure destroy_shell,
     $                   destroy_lshell, destroy_basis_def
      end interface
      interface initialize
        module procedure init_shell, init_kbshell,
     $                   init_lshell, init_basis_def
      end interface

!---------------------------------------------------------
      public  :: destroy, copy_shell, initialize
      public  :: write_basis_specs, basis_specs_transfer
      public  :: deallocate_spec_arrays
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
      target%split_norm = source%split_norm

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
      p%split_norm = 0.0_dp
      p%split_norm_specified = .false.
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
      p%lmxkb_requested = -1
      p%nkbshells = -1
      p%nshells_tmp = -1
      p%label = 'Unknown'
      p%semic = .false.
      p%ionic_charge = huge(1.0_dp)  ! To signal it was not set
      nullify(p%tmp_shell)
      nullify(p%lshell)
      nullify(p%kbshell)
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

!-----------------------------------------------------------------------
      subroutine destroy_basis_def(p)
      type(basis_def_t)          :: p

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
      write(6,'(5x,a20,g20.10)') 'rinn'           , p%rinn
      write(6,'(5x,a20,g20.10)') 'vcte'           , p%vcte
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

      integer i

      write(6,*) ' '
      write(6,*) 'SPECIES---'

      write(6,'(5x,a20,a20)') 'label',     p%label
      write(6,'(5x,a20,i20)') 'atomic number',     p%z
      write(6,'(5x,a20,a20)') 'basis type',     p%basis_type
      write(6,'(5x,a20,a20)') 'basis size',     p%basis_size
      write(6,'(5x,a20,g20.10)') 'ionic charge'  , p%ionic_charge
      write(6,'(5x,a20,i20)') 'lmax basis',     p%lmxo

      if (.not. associated(p%lshell)) then
         write(6,*) 'No L SHELLS, lmxo=', p%lmxo
         return
      endif
      do i=0, p%lmxo
         call print_lshell(p%lshell(i))
      enddo
      if (.not. associated(p%kbshell)) then
         write(6,*) 'No KB SHELLS, lmxkb=', p%lmxkb
         return
      endif
      do i=0, p%lmxkb
         call print_kbshell(p%kbshell(i))
      enddo

      write(6,*) '------------SPECIES'
      write(6,*) 

      end subroutine print_basis_def

!-----------------------------------------------------------------------
      subroutine basis_specs_transfer
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
      call re_alloc( semic, 1, nsp, name='semic',
     &               routine='basis_specs_transfer' )
      nullify( lmxkb )
      call re_alloc( lmxkb, 1, nsp, name='lmxkb',
     &               routine='basis_specs_transfer' )
      nullify( lmxo )
      call re_alloc( lmxo, 1, nsp, name='lmxo',
     &               routine='basis_specs_transfer' )
      nullify( nsemic )
      call re_alloc( nsemic, 0, lmaxd, 1, nsp, name='nsemic',
     &               routine='basis_specs_transfer' )
      nullify( cnfigmx )
      call re_alloc( cnfigmx, 0, lmaxd, 1, nsp, name='cnfigmx',
     &               routine='basis_specs_transfer' )
      nullify( nkbl )
      call re_alloc( nkbl, 0, lmaxd, 1, nsp, name='nkbl',
     &               routine='basis_specs_transfer' )
      nullify( polorb )
      call re_alloc( polorb, 0, lmaxd, 1, nsemx, 1, nsp, name='polorb',
     &               routine='basis_specs_transfer' )
      nullify( nzeta )
      call re_alloc( nzeta, 0, lmaxd, 1, nsemx, 1, nsp, name='nzeta',
     &               routine='basis_specs_transfer' )
      nullify( split_norm )
      call re_alloc( split_norm, 0, lmaxd, 1, nsemx, 1, nsp,
     &               name='split_norm',
     &               routine='basis_specs_transfer' )
      nullify( vcte )
      call re_alloc( vcte, 0, lmaxd, 1, nsemx, 1, nsp, name='vcte',
     &               routine='basis_specs_transfer' )
      nullify( rinn )
      call re_alloc( rinn, 0, lmaxd, 1, nsemx, 1, nsp, name='rinn',
     &               routine='basis_specs_transfer' )
      nullify( erefkb )
      call re_alloc( erefkb, 1, nkbmx, 0, lmaxd, 1, nsp, name='erefkb',
     &               routine='basis_specs_transfer' )
      nullify( charge )
      call re_alloc( charge, 1, nsp, name='charge',
     &               routine='basis_specs_transfer' )
      nullify( lambda )
      call re_alloc( lambda, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsp,
     &               name='lambda',
     &               routine='basis_specs_transfer' )
      nullify( rco )
      call re_alloc( rco, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsp,
     &               name='rco',
     &               routine='basis_specs_transfer' )
      nullify( iz )
      call re_alloc( iz, 1, nsp, name='iz',
     &               routine='basis_specs_transfer' )
      nullify( smass )
      call re_alloc( smass, 1, nsp, name='smass',
     &               routine='basis_specs_transfer' )
      nullify( basistype )
      call re_alloc( basistype, 1, nsp, name='basistype',
     &               routine='basis_specs_transfer' )
      nullify( atm_label )
      call re_alloc( atm_label, 1, nsp, name='atm_label',
     &               routine='basis_specs_transfer' )
!
!     Transfer
!

      nkbl(:,:) = 0
      nzeta(:,:,:) = 0
      split_norm(:,:,:) = 0.d0
      vcte(:,:,:) = 0.d0
      rinn(:,:,:) = 0.d0
      polorb(:,:,:) = 0
      rco(:,:,:,:) = 0.d0
      lambda(:,:,:,:) = 0.d0
      erefkb(:,:,:) = 0.d0
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
               vcte(l,n,isp) = s%vcte
               rinn(l,n,isp) = s%rinn
!
!              This would make the code act in the same way as
!              siesta 0.X, but it does not seem to be necessary...
!
!               if (s%nzeta_pol.ne.0) then
!                  vcte(l+1,n,isp) = s%vcte
!                  rinn(l+1,n,isp) = s%rinn
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
         do l=0,basp%lmxkb
            k=>basp%kbshell(l)
            nkbl(l,isp) = k%nkbl
            erefkb(1:k%nkbl,l,isp) = k%erefkb(1:k%nkbl)
         enddo

      enddo

      end subroutine basis_specs_transfer

!-----------------------------------------------------------------------
      subroutine write_basis_specs(lun,is)
      integer, intent(in)  :: lun
      integer, intent(in)  :: is

      integer l,  n, i

      write(lun,'(/a/79("="))') '<basis_specs>'
         write(lun,'(a20,1x,a2,i4,4x,a5,g12.5,4x,a7,g12.5)')
     $        atm_label(is), 'Z=',iz(is),
     $        'Mass=', smass(is), 'Charge=', charge(is)
         write(lun,'(a5,i1,1x,a6,i2,4x,a10,a10,1x,a6,l1)')
     $        'Lmxo=', lmxo(is), 'Lmxkb=', lmxkb(is),
     $        'BasisType=', basistype(is), 'Semic=', semic(is)
         do l=0,lmxo(is)
            write(lun,'(a2,i1,2x,a7,i1,2x,a8,i1)')
     $           'L=', l, 'Nsemic=', nsemic(l,is),
     $           'Cnfigmx=', cnfigmx(l,is)
            do n=1,nsemic(l,is)+1
               if (nzeta(l,n,is) == 0) exit
               write(lun,'(10x,a2,i1,2x,a6,i1,2x,a7,i1)')
     $                         'n=', n, 'nzeta=',nzeta(l,n,is),
     $                         'polorb=', polorb(l,n,is)
               write(lun,'(10x,a10,2x,g12.5)') 
     $                         'splnorm:', split_norm(l,n,is)
               write(lun,'(10x,a10,2x,g12.5)') 
     $                         'vcte:', vcte(l,n,is)
               write(lun,'(10x,a10,2x,g12.5)') 
     $                         'rinn:', rinn(l,n,is)
               write(lun,'(10x,a10,2x,4g12.5)') 'rcs:',
     $                         (rco(i,l,n,is),i=1,nzeta(l,n,is))
               write(lun,'(10x,a10,2x,4g12.5)') 'lambdas:',
     $                         (lambda(i,l,n,is),i=1,nzeta(l,n,is))
            enddo
         enddo
         write(lun,'(79("-"))')
         do l=0,lmxkb(is)
            write(lun,'(a2,i1,2x,a5,i1,2x,a6,4g12.5)')
     $           'L=', l, 'Nkbl=', nkbl(l,is),
     $           'erefs:  ', (erefkb(i,l,is),i=1,nkbl(l,is))
         enddo
         write(lun,'(79("="))')
      write(lun,'(a/)') '</basis_specs>'

      end subroutine write_basis_specs

      subroutine deallocate_spec_arrays()
!
      use alloc, only: de_alloc

      call de_alloc( semic, name='semic')
      call de_alloc( lmxkb, name='lmxkb')
      call de_alloc( lmxo, name='lmxo')
      call de_alloc( nsemic, name='nsemic')
      call de_alloc( cnfigmx, name='cnfigmx')
      call de_alloc( nkbl, name='nkbl')
      call de_alloc( polorb, name='polorb')
      call de_alloc( nzeta, name='nzeta')
      call de_alloc( split_norm, name='split_norm')
      call de_alloc( vcte, name='vcte')
      call de_alloc( rinn, name='rinn')
      call de_alloc( erefkb, name='erefkb')
      call de_alloc( charge, name='charge')
      call de_alloc( lambda, name='lambda')
      call de_alloc( rco, name='rco')
      call de_alloc( iz, name='iz')
      call de_alloc( smass, name='smass')
      call de_alloc( basistype, name='basistype')
      call de_alloc( atm_label, name='atm_label')
!
      end subroutine deallocate_spec_arrays

!-----------------------------------------------------------------------
      end module basis_types
