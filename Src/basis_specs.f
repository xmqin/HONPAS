! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module basis_specs
! 
! Alberto Garcia, August 2000, based on 'classic' redbasis.
! 
! Processes the basis information in an fdf file and populates
! the "basis specifications" data structures.
! 
! Here is a guide to the behavior of the main routine "read_basis_specs":
! 
! * Find out the number of species
! * Allocate storage for the data structures
! * Determine any "global" basis specification parameters:
!   - basis_size   (sz, dz, dzp, etc)
!   - basis_type   (split, nodes, etc) ("USER" is no longer valid)
! * Process the mandatory "Chemical-species-label" block to read the
!   species labels (not the same as chemical symbols) and atomic numbers.
!   (note that negative atomic numbers have a special meaning)
! * Assign default values to 'basis_size', 'basis_type'.
! * Assign default values to the 'mass' based on the atomic number.
! * Set up the information about the ground state of the atom.
! * Read information about the valence charge density in the
!   pseudopotential file and determine whether there are semicore
!   states. (Note that in this case  the user is explicitly asked 
!   (see below) to input the 'n' quantum numbers in the PAO.Basis block.)
! * Read the optional fdf blocks:
!   AtomicMass  (routine remass)
!   PAO.BasisSize - Overrides the default 'basis_size' on a per-species basis.
!   PAO.Basis - This is the most complex block, very flexible but in
!               need  of spelling-out the specific restrictions. Line by
!               line, these are:
! 
!   1st:   Species_label number_of_l_shells [basis_type] [ionic_charge] 
! 
!   For each l_shell:
!     [n= n] l nzeta [P [ nzeta_pol]] [E vcte rinn] [Q qcoe [qyuk [qwid]]] [F cutoff]
!   where 'n= n' is *mandatory* if the species has any semicore states,
!   and the 'P' (polarization), E (soft confinement potential), 
!   'Q' (charge confinement), and 'F' (filteret) sections are optional and can appear
!   in any order after nzeta. 
! 
!          rc_1  rc_2 .... rc_nzeta  
! 
!   are the cutoff radii in bohrs. This line is mandatory
!   If the number of rc's for a given shell is less than the number of
!   'zetas', the program will assign the last rc value to the remaining
!   zetas, rather than stopping with an error. This is particularly
!   useful for Bessel suites of orbitals.
!   A line containing contraction (scale) factors is optional, but if it
!   appears, the values *must be* real numbers, and there must be at
!   least nzeta of them.
!   --------------------------------------------------------------------
! 
!   After processing PAO.Basis (if it exists), whatever PAO information
!   which is not already available is determined in routine 'autobasis', 
!   using the following defaults:
! 
!   (There is a first check to make sure that species with semicore
!   states, or with Bessel floating orbitals,  have been included in 
!   the PAO.Basis block)
!   
!   The maximum angular momentum for the basis (lmxo) (excluding any
!   polarization orbitals) is that of the ground state, as returned by
!   routine 'lmxofz' from Z.
! 
!   Each l-shell contains just one shell of PAOS, with n set
!   to the appropriate ground state value.
! 
!   Nzeta is determined by the first two characters of the 'basis_size'
!   string: 1 for 'sz' and  2 for 'dz'.
! 
!   There are no 'per l-shell' polarization orbitals, except if the
!   third character of 'basis_size' is 'p' (as in 'dzp'), in which
!   case polarization orbitals are defined so they have the minimum      
!   angular momentum l such that there are no occupied orbitals 
!   with the same l in the valence shell of the ground-state 
!   atomic configuration. They polarize the corresponding l-1 shell.
! 
!   The Soft-Confinement parameters 'rinn' and 'vcte' are set to 0.0
!   The Charge-Confinement parameters 'qcoe', 'qyuk' and 'qwid' 
!   are set to 0.0, 0.0 and 0.01
! 
!   rc(1:nzeta) is set to 0.0
!   lambda(1:nzeta) is set to 1.0  (this is a change from old practice)
!
!  ----------------------------------
!  
!   Next come the blocks associated to the KB projectors:
! 
!   Block PS.lmax (routine relmxkb) assigns a per-species maximum l for
!   the KB projectors. Internally, it sets the 'lmxb_requested'
!   component of the data structure, which by default is -1.
!  
!   Block PS.KBprojectors (optional) sets the number of KB projectors
!   per angular momentum. The same routine which reads this block
!   (readkb) also sets any appropriate defaults:
! 
! 
!   If the species does not appear in the  PS.KBprojectors block, then
!       Setting of lmxkb:
!       If lmxkb is set in PS.lmax, 
!          set it to that
!       else
!          Use PAO information: (routine set_default_lmxkb)
!          (Set it to a meaningless value for floating orbitals)
!           Set it to lmxo+1, or to lpol+1, where lpol is the angular
!           momentum of the highest-l polarization orbital.
!       endif
!
!       There is a hard limit for lmxkb: if the pseudopotential file
!       contains semilocal pseudopotentials up to lmax_pseudo, then
!       lmxkb <= lmax_pseudo.
!     
!       The  number of KB projectors per l is set to the number of
!       n-shells in the corresponding PAO shell with the same l. For l
!       greater than lmxo, it is set to 1. The reference energies are
!       in all cases set to huge(1.d0) to mark the default.
!
!       Archaeological note: The first implementation of the basis-set generation
!       module had only non-polarization orbitals. Polarization orbitals were added
!       later as "second-class" companions. This shows in details like "lmxo" (the
!       maximum l of the basis set) not taking into account polarization orbitals.
!       Polarization orbitals were tagged at the end, without maintaining l-shell
!       ordering.
!       Even later, support for "semicore" orbitals was added. A new ('nsm') index was
!       used to distinguish the different orbitals in a l-shell. Polarization orbitals
!       were not brought into this classification. 
!       The code is strained when semicore and polarization orbitals coexist for the same l,
!       as in Ti, whose electronic structure is []3s2 3p6 3d2 4s2 (4p0*)  with 4p as the
!       polarization orbital.  
!       Here 'nsemic' (the number of semicore states) is kept at zero for l=1 
!       (the polarization orbital is not counted), with the side-effect that only one KB
!       projector is generated for l=1.
!       Special checks have been implemented to cover these cases.
! 
!       Future work should probably remove the separate treatment of polarization orbitals.
!       (Note that if *all* orbitals are specified in a PAO.Basis block, in effect turning
!       perturbative polarization orbitals into 'normal' orbitals, this problem is not present.)
!     
! =======================================================================
!
      use precision
      use basis_types, only: basis_def_t, shell_t, lshell_t, kbshell_t
      use basis_types, only: nsp, basis_parameters, ground_state_t
      use basis_types, only: destroy, copy_shell, initialize
      use pseudopotential, only: pseudo_read, pseudo_reparametrize
      use pseudopotential, only: pseudo_init_constant
      use periodic_table, only: qvlofz, lmxofz, cnfig, atmass
      use chemical
      use sys
      use fdf

      Implicit None

      type(basis_def_t), pointer :: basp => null()
      type(shell_t), pointer :: s => null()
      type(lshell_t), pointer :: ls => null()
      type(kbshell_t), pointer :: k => null()

      character(len=1), parameter   ::
     $                           sym(0:4) = (/ 's','p','d','f','g' /)

!     Default Soft-confinement parameters set by the user
!
      logical, save   :: lsoft
      real(dp), save  :: softRc, softPt

! Default norm-percentage for the automatic definition of
! multiple-zeta orbitals with the 'SPLIT' option

      real(dp), parameter       :: splnorm_default=0.15_dp
      real(dp), parameter       :: splnormH_default=-1.0_dp
      real(dp), save            :: global_splnorm, global_splnorm_H
      integer           isp  ! just an index dummy variable for the whole module

!
      logical, save, public     :: restricted_grid
      logical, parameter        :: restricted_grid_default = .true.

      real(dp), save, public    :: rmax_radial_grid
      
      public :: read_basis_specs
      public :: label2species

      private

      CONTAINS

!---
      subroutine read_basis_specs()

      character(len=15), parameter  :: basis_size_default='standard'
      character(len=10), parameter  :: basistype_default='split'

      character(len=15) :: basis_size
      character(len=10) :: basistype_generic

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      type(ground_state_t), pointer :: gs

      integer nns, noccs, i, ns_read, l
      logical synthetic_atoms, found, reparametrize_pseudos
      real(dp) :: new_a, new_b

!------------------------------------------------------------------------
      reparametrize_pseudos =
     $   fdf_boolean('ReparametrizePseudos',.false.)
!
!        r(i) = b * ( exp(a*(i-1)) - 1 )
!        spacing near zero: delta = a*b
!        If the desired spacing at r=R is delta*(1+beta), 
!        then:
!               a = beta*delta/R
!               b = R/beta
!               N at R = 1 + R*ln(1+beta)/(beta*delta)
!
      if (reparametrize_pseudos) then
         ! These default values will provide a grid spacing
         ! of 1.0e-5 near r=0, and 0.01 near r=10 (bohr units)
         ! with N at Rmax=100 on the order of 10000 points.
         new_a = fdf_double("NewAParameter",0.001_dp)
         new_b = fdf_double("NewBParameter",0.01_dp)
      endif

!
!     Whether or not to restrict integration grids to an odd number
!     of points (and hence to shift rcs accordingly)     
!
      restricted_grid =
     $     fdf_boolean("Restricted.Radial.Grid",
     $                  restricted_grid_default)
!
!     If non-zero, the value will be the maximum value of the
!     radial coordinate
!
      if (reparametrize_pseudos) then
         rmax_radial_grid = fdf_double('Rmax.Radial.Grid',50.0_dp)
      else
         rmax_radial_grid = fdf_double('Rmax.Radial.Grid',0.0_dp)
      endif

!
      basis_size = fdf_string('PAO.BasisSize',basis_size_default)
      call size_name(basis_size)
      basistype_generic = fdf_string('PAO.BasisType',basistype_default)
      call type_name(basistype_generic)

C Read information about defaults for soft confinement

      lsoft  = fdf_boolean('PAO.SoftDefault',.false.)
      softRc = fdf_double('PAO.SoftInnerRadius',0.9d0)
      softPt = fdf_double('PAO.SoftPotential',40.0d0)

C Sanity checks on values

      softRc = max(softRc,0.00d0)
      softRc = min(softRc,0.99d0)
      softPt = abs(softPt)
!
!     Read defaults for split_norm parameter

      global_splnorm = fdf_double('PAO.SplitNorm',splnorm_default)
      global_splnorm_H = fdf_double('PAO.SplitNormH',splnormH_default)
      if (global_splnorm_H < 0.0_dp) global_splnorm_H = global_splnorm

!------------------------------------------------------------------
!
!     Use standard routine in chemical module to process the
!     chemical species
!
      nsp = size(basis_parameters)

      synthetic_atoms = .false.

      do isp=1,nsp
        basp=>basis_parameters(isp)

        basp%label = species_label(isp)
        basp%z = atomic_number(isp)
        basp%floating = is_floating(isp)
        basp%bessel = is_bessel(isp)
        basp%synthetic = is_synthetic(isp)

        basp%basis_size = basis_size
        basp%basis_type = basistype_generic
        if (basp%floating) then
          basp%mass = 1.d40   ! big but not too big, as it is used
                              ! later in computations
        else if (basp%synthetic) then
          basp%mass = -1.0_dp      ! Signal -- Set later
        else
          basp%mass = atmass(abs(int(basp%z)))
        endif
        if (basp%bessel) then
          ! Initialize a constant pseudo
          call pseudo_init_constant(basp%pseudopotential)
        else if (basp%synthetic) then
          synthetic_atoms = .true.
          ! Will set gs later
          call pseudo_read(basp%label,basp%pseudopotential)
        else
          call ground_state(abs(int(basp%z)),basp%ground_state)
          call pseudo_read(basp%label,basp%pseudopotential)
        endif
        if (reparametrize_pseudos.and. .not. basp%bessel)
     .    call pseudo_reparametrize(p=basp%pseudopotential,
     .                             a=new_a, b=new_b,label=basp%label)
      enddo

      if (synthetic_atoms) then

        found = fdf_block('SyntheticAtoms',bfdf)
        if (.not. found )
     .    call die("Block SyntheticAtoms does not exist.")
        ns_read = 0
        do while(fdf_bline(bfdf, pline))

          ns_read = ns_read + 1
          if (.not. fdf_bmatch(pline,'i'))
     .      call die("Wrong format in SyntheticAtoms")
          isp = fdf_bintegers(pline,1)
          if (isp .gt. nsp .or. isp .lt. 1)
     .      call die("Wrong specnum in SyntheticAtoms")
          basp => basis_parameters(isp)
          gs   => basp%ground_state
          if (.not. fdf_bline(bfdf, pline)) call die("No n info")
          nns = fdf_bnintegers(pline)
          if (nns .lt. 4)
     .      call die("Please give all valence n's " //
     .               "in SyntheticAtoms block")
          gs%n = 0
          do i = 1, nns
            gs%n(i-1) = fdf_bintegers(pline,i)
          enddo
          if (.not. fdf_bline(bfdf, pline))
     .      call die("No occupation info")
          noccs = fdf_bnvalues(pline)
          if (noccs .lt. nns) call die("Need more occupations")
          gs%occupation(:) = 0.0_dp
          do i = 1, noccs
            gs%occupation(i-1) = fdf_bvalues(pline,i)
          enddo
          ! Determine the last occupied state in the atom
          do i = nns, 1, -1
            if (gs%occupation(i-1) /=0) then
              gs%lmax_valence = i-1
              exit
            endif
          enddo
          gs%occupied(0:3) = (gs%occupation .gt. 0.0_dp)
          gs%occupied(4) = .false.
          gs%z_valence = sum(gs%occupation(0:noccs-1))
          write(6,'(a,i2)',advance='no')
     .         'Ground state valence configuration (synthetic): '
          do l=0,3
            if (gs%occupied(l))
     .        write(6,'(2x,i1,a1,f8.5)',advance='no')
     .          gs%n(l),sym(l), gs%occupation(l)
          enddo
          write(6,'(a)') ''

        enddo
        write(6,"(a,i2)") "Number of synthetic species: ", ns_read

      endif
!
!  Defer this here in case there are synthetic atoms
      do isp=1,nsp
        basp=>basis_parameters(isp)
        if (basp%synthetic) then
          gs => basp%ground_state
          if (gs%z_valence .lt. 0.001)
     .      call die("Synthetic species not detailed")
         endif
         call semicore_check(isp)
      enddo

      call remass()
      call resizes()
      call repaobasis()
      call autobasis()
      call relmxkb()
      call readkb()

!      do isp=1,nsp
!        call print_basis_def(basis_parameters(isp))
!      enddo

      end subroutine read_basis_specs
!-----------------------------------------------------------------------

      function label2species(str) result(index)
      integer index
      character(len=*), intent(in) ::  str

      integer i

      index = 0
      do i=1,nsp
        if(leqi(basis_parameters(i)%label,str)) index = i
      enddo
      end function label2species
!-----------------------------------------------------------------------

      subroutine ground_state(z,gs)
      integer, intent(in)               ::  z
      type(ground_state_t), intent(out) :: gs
!
!     Determines ground state valence configuration from Z
!
      integer l, latm

      gs%z_valence = 0.d0
      do l=0,3
        gs%occupation(l)=0.0d0
      enddo

      call lmxofz(z,gs%lmax_valence,latm)
      call qvlofz(z,gs%occupation(:))
      do l=0,gs%lmax_valence
        gs%z_valence = gs%z_valence + gs%occupation(l)
      enddo
      call cnfig(z,gs%n(0:3))

      write(6,'(a,i2)',advance='no')
     .     'Ground state valence configuration: '
      gs%occupied(4) = .false.         !! always
      do l=0,3
        gs%occupied(l) =  (gs%occupation(l).gt.0.1d0)
        if (gs%occupied(l))
     .    write(6,'(2x,i1,a1,i2.2)',advance='no')
     .       gs%n(l),sym(l),nint(gs%occupation(l))
      enddo
      write(6,'(a)') ''

      end subroutine ground_state

!---------------------------------------------------------------
      subroutine readkb()
      integer lpol, isp, ish, i, l
      character(len=20) unitstr

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer :: lmax_pseudo
      
      lpol = 0

      if (fdf_block('PS.KBprojectors',bfdf) ) then
 
! First pass to find out about lmxkb and set any defaults.

        do while(fdf_bline(bfdf, pline))    !! over species
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.KBprojectors")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PS.KBprojectors:",
     .        trim(fdf_bnames(pline,1))
            call die()
          endif
          basp => basis_parameters(isp)
          basp%nkbshells = fdf_bintegers(pline,1)
          do ish= 1, basp%nkbshells
            if (.not. fdf_bline(bfdf, pline)) call die("No l nkbl")
            if (.not. fdf_bmatch(pline,'ii'))
     .        call die("Wrong format l nkbl")
            l = fdf_bintegers(pline,1)
            if (l .gt. basp%lmxkb) basp%lmxkb = l
            if (.not. fdf_bline(bfdf, pline)) then
               if (ish .ne. basp%nkbshells)
     .          call die("Not enough kb shells for this species...")
              ! There is no line with ref energies
            else  if (fdf_bmatch(pline,'ni')) then
                ! We are seeing the next species' section
               if (ish .ne. basp%nkbshells)
     .            call die("Not enough shells for this species...")
               if (.not. fdf_bbackspace(bfdf))
     .               call die('readkb: ERROR in PS.KBprojectors block')
            else if (fdf_bmatch(pline,'ii')) then
                ! We are seeing the next shell's section
               if (ish .gt. basp%nkbshells)
     .            call die("Too many kb shells for this species...")
               if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
            endif
          enddo       ! end of loop over shells for species isp

        enddo
      endif

      do isp=1,nsp
!
!      Fix defaults and allocate kbshells
!
        basp=>basis_parameters(isp)
        if (basp%lmxkb .eq. -1) then ! not set in KBprojectors 
          if (basp%lmxkb_requested.eq.-1) then ! not set in PS.lmax
            basp%lmxkb = set_default_lmxkb(isp) ! Use PAO info
          else
            basp%lmxkb = basp%lmxkb_requested
          endif
          allocate(basp%kbshell(0:basp%lmxkb))
          do l=0,basp%lmxkb
            call initialize(basp%kbshell(l))
            k=>basp%kbshell(l)
            k%l = l
            if (l.gt.basp%lmxo) then
              k%nkbl = 1
            else
              ! Set equal to the number of PAO shells with this l
              k%nkbl = basp%lshell(l)%nn     
              ! Should include polarization orbs (as in Ti case: 3p..4p*)
              ! (See 'archaeological note' in the header of this file)
              if (l>0) then
                 do i = 1, basp%lshell(l-1)%nn
                    if (basp%lshell(l-1)%shell(i)%polarized) then
                       k%nkbl = k%nkbl + 1
                       write(6,"(a,i1,a)") trim(basp%label) //
     $                  ': nkbl increased for l=',l,
     $                  ' due to the presence of a polarization orbital'
                    endif
                 enddo
              endif
              if (k%nkbl.eq.0) then
                write(6,*) 'Warning: Empty PAO shell. l =', l
                write(6,*) 'Will have a KB projector anyway...'
                k%nkbl = 1
              endif
            endif
            allocate(k%erefkb(1:k%nkbl))
            k%erefkb(1:k%nkbl) = huge(1.d0)
          enddo
        else          ! Set in KBprojectors
          if (basp%lmxkb_requested.ne.-1) then ! set in PS.lmax
            if (basp%lmxkb.ne.basp%lmxkb_requested) then
              call die("LmaxKB conflict between " //
     $                 "PS.Lmax and PS.KBprojectors blocks")
            endif
          endif
          !! OK, we have a genuine lmxkb
          allocate(basp%kbshell(0:basp%lmxkb))
          do l=0,basp%lmxkb
            call initialize(basp%kbshell(l))
          enddo
        endif
        if (basp%z .le. 0) then
          if (basp%lmxkb .ne. -1)
     .      call die("Floating orbs cannot have KB projectors...")
        endif
      enddo
!
!     Now re-scan the block (if it exists) and fill in as instructed
!            
      if (fdf_block('PS.KBprojectors',bfdf) ) then

        do while(fdf_bline(bfdf, pline))     !! over species
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.KBprojectors")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PS.KBprojectors:",
     .        trim(fdf_bnames(pline,1))
            call die()
          endif
          basp => basis_parameters(isp)
          basp%nkbshells = fdf_bintegers(pline,1)
          do ish=1, basp%nkbshells
            if (.not. fdf_bline(bfdf,pline)) call die("No l nkbl")
            if (.not. fdf_bmatch(pline,'ii'))
     .        call die("Wrong format l nkbl")
            l = fdf_bintegers(pline,1)
            k => basp%kbshell(l)
            k%l = l
            k%nkbl = fdf_bintegers(pline,2)
            if (k%nkbl < 0) then
               call die("nkbl < 0 in PS.KBprojectors")
            endif
            if (k%nkbl == 0) then
               call message("WARNING","nkbl=0 in PS.KBprojectors")
            endif
            allocate(k%erefkb(k%nkbl))
            if (.not. fdf_bline(bfdf,pline)) then
              if (ish .ne. basp%nkbshells)
     .          call die("Not enough shells for this species...")
              ! There is no line with ref energies
              ! Use default values
              k%erefKB(1:k%nkbl) = huge(1.d0)
            else if (fdf_bmatch(pline,'ni')) then
                ! We are seeing the next species' section
               if (ish .ne. basp%nkbshells)
     .            call die("Not enough kb shells for this species...")
                ! Use default values for ref energies
               k%erefKB(1:k%nkbl) = huge(1.d0)
               if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
            else if (fdf_bmatch(pline,'ii')) then
                ! We are seeing the next shell's section
                if (ish .gt. basp%nkbshells)
     .            call die("Too many kb shells for this species...")
                ! Use default values for ref energies
                k%erefKB(1:k%nkbl) = huge(1.d0)
                if (.not. fdf_bbackspace(bfdf))
     .            call die('readkb: ERROR in PS.KBprojectors block')
             else
                if (fdf_bnreals(pline) .ne. k%nkbl)
     .            call die("Wrong number of energies")
                unitstr = 'Ry'
                if (fdf_bnnames(pline) .eq. 1)
     .            unitstr = fdf_bnames(pline,1)
                ! Insert ref energies in erefkb
                do i= 1, k%nkbl
                  k%erefKB(i) =
     .                 fdf_breals(pline,i)*fdf_convfac(unitstr,'Ry')
                enddo
             endif
          enddo            ! end of loop over shells for species isp

          ! For those l's not specified in block, use default values
          do l=0, basp%lmxkb
            k => basp%kbshell(l)
            if (k%l.eq.-1) then
              k%l = l
              k%nkbl = 1
              allocate(k%erefkb(1))
              k%erefkb(1) = huge(1.d0)
            endif
          enddo
        enddo   !! Over species
      endif

      do isp=1,nsp
!
!      Check that we have enough semilocal components...
!
         basp=>basis_parameters(isp)
         lmax_pseudo = basp%pseudopotential%npotd - 1 
         if (basp%lmxkb > lmax_pseudo) then
            write(6,'(a,i1,a)')
     .           trim(basp%label) //
     .           " pseudopotential only contains V_ls up to l=",
     .           lmax_pseudo, " -- lmxkb reset."
            basp%lmxkb = lmax_pseudo
         endif
      enddo

      end subroutine readkb
!---------------------------------------------------------------

      subroutine repaobasis()

      integer isp, ish, nn, i, ind, l, indexp, index_splnorm
      integer nrcs_zetas

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      if (.not. fdf_block('PAO.Basis',bfdf)) RETURN

      do while(fdf_bline(bfdf, pline))     !! over species
        if (.not. fdf_bmatch(pline,'ni'))
     .    call die("Wrong format in PAO.Basis")
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(6,'(a,1x,a)')
     .      "WRONG species symbol in PAO.Basis:",
     .      trim(fdf_bnames(pline,1))
          call die()
        endif

        basp => basis_parameters(isp)
        basp%label = fdf_bnames(pline,1)
        basp%nshells_tmp = fdf_bintegers(pline,1)
        basp%lmxo = 0
        !! Check whether there are optional type and ionic charge
        if (fdf_bnnames(pline) .eq. 2)
     .    basp%basis_type = fdf_bnames(pline,2)
        if (fdf_bnvalues(pline) .eq. 2)
     .    basp%ionic_charge = fdf_bvalues(pline,2)
        allocate(basp%tmp_shell(basp%nshells_tmp))

        shells: do ish= 1, basp%nshells_tmp
          s => basp%tmp_shell(ish)
          call initialize(s)
          if (.not. fdf_bline(bfdf,pline)) call die("No l nzeta, etc")

          if (fdf_bmatch(pline,'niii')) then
            s%n = fdf_bintegers(pline,1)
            s%l = fdf_bintegers(pline,2)
            basp%lmxo = max(basp%lmxo,s%l)
            s%nzeta = fdf_bintegers(pline,3)
          elseif (fdf_bmatch(pline,'ii')) then
            !    l, nzeta

            if (basp%semic)
     .        call die("Please specify n if there are semicore states")

            s%l = fdf_bintegers(pline,1)
            s%n = basp%ground_state%n(s%l)
            s%nzeta = fdf_bintegers(pline,2)
            basp%lmxo = max(basp%lmxo,s%l)
          else
            call die("Bad format of (n), l, nzeta line in PAO.Basis")
          endif
!
! If this is a filteret basis then the number of zetas input must be one
!
          if (basp%basis_type.eq.'filteret') then
            s%nzeta = 1
          endif
!
! Optional stuff: Polarization, Soft-confinement Potential and Filteret Cutoff
!
! Split norm
!
          if (fdf_bsearch(pline,'S',index_splnorm)) then
            if (fdf_bmatch(pline,'v',after=index_splnorm)) then
              s%split_norm = fdf_bvalues(pline, ind=1,
     .                                   after=index_splnorm)
              if (s%split_norm .eq. 0.0_dp)
     .          write(6,"(a)")
     .            "WARNING: zero split_norm after S in PAO.Basis"
              s%split_norm_specified = .TRUE.
            else
              call die("Specify split_norm after S in PAO.Basis")
            endif
          else
            if (abs(basp%z) .eq. 1) then
              s%split_norm = global_splnorm_H
            else
              s%split_norm = global_splnorm
            endif
          endif
!
! Polarization functions
!
          if (fdf_bsearch(pline,'P',indexp)) then
            s%polarized = .TRUE.
            if (fdf_bmatch(pline,'i',after=indexp)) then
              s%nzeta_pol = fdf_bintegers(pline,ind=1,after=indexp)
            else
              s%nzeta_pol = 1
            endif
          endif
!
! Soft-confinement
!
          if (fdf_bsearch(pline,'E',indexp)) then
            if (fdf_bmatch(pline,'vv',after=indexp)) then
              s%vcte = fdf_bvalues(pline,ind=1,after=indexp)
              s%rinn = fdf_bvalues(pline,ind=2,after=indexp)
            else
              call die("Need vcte and rinn after E in PAO.Basis")
            endif
          elseif (lsoft) then
            s%vcte = softPt 
            s%rinn = -softRc
          else
            s%vcte = 0.0_dp
            s%rinn = 0.0_dp
          endif
!
! Charge confinement
! 
          if (fdf_bsearch(pline,'Q',indexp)) then
            if (fdf_bmatch(pline,'vvv',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
               s%qwid = fdf_bvalues(pline,ind=3,after=indexp)
            elseif (fdf_bmatch(pline,'vv',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = fdf_bvalues(pline,ind=2,after=indexp)
               s%qwid = 0.01_dp
            elseif (fdf_bmatch(pline,'v',after=indexp)) then
               s%qcoe = fdf_bvalues(pline,ind=1,after=indexp)
               s%qyuk = 0.0_dp
               s%qwid = 0.01_dp
            else
               call die("Need one, two or three real numbers after Q in 
     .                   PAO.Basis")
            endif
          else
            s%qcoe = 0.0_dp
            s%qyuk = 0.0_dp
            s%qwid = 0.01_dp
          endif
!
! Filteret cutoff
!
          if (fdf_bsearch(pline,"F",indexp)) then
            if (fdf_bmatch(pline,"v",after=indexp)) then
              s%filtercut = fdf_bvalues(pline,ind=1,after=indexp)
            else
              call die("Need cut-off after F in PAO.Basis")
            endif
          else
            s%filtercut = 0.0_dp
          endif

          allocate(s%rc(s%nzeta),s%lambda(s%nzeta))
          s%rc(:) = 0.d0
          s%lambda(:) = 1.d0
          if (.not. fdf_bline(bfdf,pline)) call die("No rc's")
          
          ! Use the last rc entered for the successive zetas
          ! if there are not enough values (useful for Bessel)
          nrcs_zetas = fdf_bnvalues(pline)
          if (nrcs_zetas < 1) then
           call die("Need at least one rc per shell in PAO.Basis block")
          endif
          do i= 1, s%nzeta
             if (i <= nrcs_zetas) then
                s%rc(i) = fdf_bvalues(pline,i)
             else
                s%rc(i) = s%rc(nrcs_zetas)
             endif
          enddo

          if (s%split_norm_specified) then
            do i = 2,s%nzeta
              if (s%rc(i) /= 0.0_dp) then
                write(6,"(/,a,i1,a,f8.4,/)")
     .            "*Warning: Per-shell split_norm parameter " //
     .            "will not apply to zeta-", i, ". rc=", s%rc(i)
              endif
            enddo
          endif

          ! Optional scale factors. They MUST be reals, or else...
          if (.not. fdf_bline(bfdf,pline)) then
            if (ish .ne. basp%nshells_tmp)
     .        call die("Not enough shells")
              ! Default values for scale factors
          else
            if (.not. fdf_bmatch(pline,'r')) then
              ! New shell or species
              ! Default values for the scale factors
              if (.not. fdf_bbackspace(bfdf)) 
     .          call die('repaobasis: ERROR in PAO.Basis block')
              cycle shells
            else
              ! Read scale factors
              ! Use the last scale factor entered for the successive zetas
              ! if there are not enough values 
               nrcs_zetas = fdf_bnreals(pline)
               if (nrcs_zetas < 1) then
                 call die("Need at least one scale factor in PAO.Basis")
               endif
               do i= 1, s%nzeta
                  if (i <= nrcs_zetas) then
                     s%lambda(i) = fdf_breals(pline,i)
                  else
                     s%lambda(i) = s%lambda(nrcs_zetas)
                  endif
               enddo
            endif
          endif

        enddo shells
        ! Clean up for this species
      enddo
!
!        OK, now classify the states by l-shells
!
      do isp = 1, nsp
        basp => basis_parameters(isp)
        if (basp%lmxo .eq. -1) cycle !! Species not in block
                                     !!
        allocate (basp%lshell(0:basp%lmxo))
        loop_l: do l= 0, basp%lmxo
          ls => basp%lshell(l)
          call initialize(ls)
          ls%l = l
          ! Search for tmp_shells with given l
          nn = 0
          do ish= 1, basp%nshells_tmp
            s => basp%tmp_shell(ish)
            if (s%l .eq. l) nn=nn+1
          enddo
          ls%nn = nn
          if (nn.eq.0) then
            !! What else do we do here?
            cycle loop_l  
          endif
                                     !! 
          allocate(ls%shell(1:nn))
          ! Collect previously allocated shells
          ind = 0
          do ish=1, basp%nshells_tmp
            s => basp%tmp_shell(ish)
            if (s%l .eq. l) then
              ind = ind+1
              call copy_shell(source=s,target=ls%shell(ind))
            endif
          enddo
          if (nn.eq.1) then
            ! If n was not specified, set it to ground state n
            if (ls%shell(1)%n.eq.-1)
     .        ls%shell(1)%n=basp%ground_state%n(l)
          endif
          !! Do we have to sort by n value????
          !!
        enddo loop_l
        !! Destroy temporary shells in basp
        !! Warning: This does seem to destroy information!!
        call destroy(basp%tmp_shell)
      enddo

      end subroutine repaobasis
!_______________________________________________________________________


         function set_default_lmxkb(is) result(lmxkb)
         integer lmxkb
         integer, intent(in) :: is

         integer lpol, l, i

         lmxkb = -1
         basp=>basis_parameters(is)
         if (basp%z .le. 0) return     ! Leave it at -1 for floating orbs.

         lmxkb = basp%lmxo + 1

         ! But watch out for polarization orbitals...
         ! We ASSUME that these do not count towards lshell%nn...

         lpol = 0
         do l = 0, basp%lmxo
            ls=>basp%lshell(l)
            do i = 1, ls%nn
               s=>ls%shell(i)
               if (s%polarized) lpol = s%l + 1
            enddo
         enddo
         if (lpol .gt. basp%lmxo) lmxkb = lpol + 1

         write(6,'(3a,i1,/,2a,/,a)') 'For ', trim(basp%label),
     .              ', standard SIESTA heuristics set lmxkb to ',
     .              lmxkb,
     .              ' (one more than the basis l,',
     .              ' including polarization orbitals).',
     .  'Use PS.lmax or PS.KBprojectors blocks to override.'
!
!        But there is an upper limit for sanity: f is the highest
!
         if (lmxkb.gt.3) then
            write(6,'(3a,i1)') 'Warning: For ', trim(basp%label),
     .           ' lmxkb would have been set to ', lmxkb
            write(6,'(a)')
     .           'Setting it to maximum value of 3 (f projector)'
            lmxkb = 3
         endif

         end function set_default_lmxkb

!-----------------------------------------------------------------------

      subroutine resizes()

c Reading atomic basis sizes for different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have the basis sizes
c given by the general input PAO.BasisSize, or its default value.

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

      if (fdf_block('PAO.BasisSizes',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'nn'))
     .      call die("Wrong format in PAO.BasisSizes")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     .        "WRONG species symbol in PAO.BasisSizes:",
     .        trim(fdf_bnames(pline,1))
            call die()
          else
            basp => basis_parameters(isp)
            basp%basis_size = fdf_bnames(pline,2)
            call size_name(basp%basis_size)   !!! DEPRECATED
            write(6,'(4a)')
     .           'resizes: Read basis size for species ',
     .            trim(fdf_bnames(pline,1)),' = ',basp%basis_size
          endif
        enddo
      endif

      end subroutine resizes

!-----------------------------------------------------------------------

      subroutine relmxkb()

c Reads the maximum angular momentum of the Kleinman-Bylander
c projectors for the different species.
c
c Reads fdf block. Not necessarily all species have to be given.

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

      if (fdf_block('PS.lmax',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'ni'))
     .      call die("Wrong format in PS.lmax")
          isp = label2species(fdf_bnames(pline,1))
          if (isp .eq. 0) then
            write(6,'(a,1x,a)') "WRONG species symbol in PS.lmax:",
     .                           trim(fdf_bnames(pline,1))
            call die()
          else
            basp => basis_parameters(isp)
            basp%lmxkb_requested = fdf_bintegers(pline,1)
            write(6,"(a, i4, 2a)")
     .            'relmxkb: Read Max KB Ang. Momentum= ',
     .             basp%lmxkb_requested,
     .            ' for species ', trim(fdf_bnames(pline,1))
          endif
        enddo
      endif

      end subroutine relmxkb
!-----------------------------------------------------------------------

      subroutine remass()


      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      integer isp

c Read atomic masses of different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have their natural mass
c (according to atmass subroutine).

      if (fdf_block('AtomicMass',bfdf)) then
        do while(fdf_bline(bfdf,pline))
          if (.not. fdf_bmatch(pline,'iv'))
     .      call die("Wrong format in AtomicMass")
          isp = fdf_bintegers(pline,1)                                 
          if (isp .gt. nsp .or. isp .lt. 1)                   
     .      call die("Wrong specnum in AtomicMass")        
          basp => basis_parameters(isp)                         
          basp%mass = fdf_bvalues(pline,2)                             
          write(6,"(a, i4, a, f12.5)")                        
     .         'remass: Read atomic mass for species ', isp,  
     .         ' as ', basp%mass                              
        enddo
      endif

      end subroutine remass
!-----------------------------------------------------------------------

      subroutine size_name(str)
      character(len=*), intent(inout)  ::  str


      if(leqi(str,'MINIMAL')) str='sz'
      if(leqi(str,'SZ'))  str='sz'
      if(leqi(str,'SZP')) str='szp'
      if(leqi(str,'SZP1')) str='szp'
      if(leqi(str,'SZSP')) str='szp'
      if(leqi(str,'SZ1P')) str='szp'
!
      if(leqi(str,'DZ')) str='dz'
      if(leqi(str,'STANDARD')) str='dzp'
      if(leqi(str,'DZP'))  str='dzp'
      if(leqi(str,'DZP1'))  str='dzp'
      if(leqi(str,'DZ1P'))  str='dzp'
      if(leqi(str,'DZSP'))  str='dzp'
      if(leqi(str,'DZP2'))  str='dzp2'
      if(leqi(str,'DZDP'))  str='dzp2'
      if(leqi(str,'DZ2P'))  str='dzp2'
!
      if(leqi(str,'TZ')) str='tz'
      if(leqi(str,'TZP')) str='tzp'
      if(leqi(str,'TZ1P')) str='tzp'
      if(leqi(str,'TZP1')) str='tzp'
      if(leqi(str,'TZSP')) str='tzp'
      if(leqi(str,'TZP2')) str='tzp2'
      if(leqi(str,'TZ2P')) str='tzp2'
      if(leqi(str,'TZDP')) str='tzp2'
      if(leqi(str,'TZP3')) str='tzp3'
      if(leqi(str,'TZ3P')) str='tzp3'
      if(leqi(str,'TZTP')) str='tzp3'
 
      if ( (str.ne.'szp').and.(str.ne.'sz').and.
     .     (str.ne.'dz') .and.(str.ne.'dzp') .and.
     .     (str.ne.'tz') .and.(str.ne.'tzp') .and.
     .     (str.ne.'dzp2') .and.
     .     (str.ne.'tzp2') .and. (str.ne.'tzp3') ) then

        write(6,'(/,2a,/,9(a,/))')
     .    'size_name: Incorrect basis-size option specified,',
     .    ' active options are:',
     .    '  SZ or MINIMAL', 
     .    '  SZP, SZSP, SZ1P, SZP1',
     .    '  DZ ',
     .    '  DZP, DZSP, DZP1, DZ1P or STANDARD',
     .    '  DZDP, DZP2, DZ2P ',
     .    '  TZ ',
     .    '  TZP, TZSP, TZP1, TZ1P',
     .    '  TZDP, TZP2, TZ2P',
     .    '  TZTP, TZP3, TZ3P'

        call die()
      endif

      end subroutine size_name
!-----------------------------------------------------------------------

      subroutine type_name(basistype)

      character basistype*(*)

      if(leqi(basistype,'NODES')) then
        basistype='nodes'
      elseif(leqi(basistype,'NONODES')) then
        basistype='nonodes'
      elseif(leqi(basistype,'SPLIT')) then
        basistype='split'
      elseif(leqi(basistype,'SPLITGAUSS')) then
        basistype='splitgauss'
      elseif (leqi(basistype,'FILTERET')) then
        basistype='filteret'
      else
        write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .    'type_name: Incorrect basis-type option specified,',
     .    ' active options are:',
     .    'NODES','SPLIT','SPLITGAUSS','NONODES','FILTERET'
         call die
      endif

      end subroutine type_name
!-----------------------------------------------------------------------
!
!     Find out whether semicore states are implied by the valence
!     charge density in the pseudopotential file.
!
      subroutine semicore_check(is)
      integer, intent(in)  :: is

      real(dp), parameter :: tiny = 1.d-5
      integer ndiff
      real(dp) zval, zval_vps, charge_loc

      basp => basis_parameters(is)

      basp%semic = .false.
      if (basp%bessel) return

      zval_vps = basp%pseudopotential%zval
      zval = basp%ground_state%z_valence
      
      if (abs(Zval-zval_vps).lt.tiny) return

      ndiff = nint(abs(Zval-zval_vps))
      if (abs(ndiff-abs(Zval-zval_vps)).gt.tiny) then
        write(6,'(2a)')
     .    'ERROR expected valence charge for species ',
     .    basp%label
        write(6,'(a)')
     .    'ERROR and the value read from the vps file'
        write(6,'(a,f6.3,a,f6.3)')
     .    'ERROR differ:  Zval(expected)= ', Zval,
     .    ' Zval(vps)= ',zval_vps
        call die()
      endif

      basp%semic = .true.
      charge_loc = Zval_vps-Zval
      write(6,'(a,i2,a)')
     .  'Semicore shell(s) with ', nint(charge_loc),
     .  ' electrons included in the valence for', trim(basp%label)

      end subroutine semicore_check
!----------------------------------------------------------------------
      subroutine autobasis()

!
!     It sets the defaults if a species has not been included
!     in the PAO.Basis block
!
      integer l, nzeta, nzeta_pol

      loop: do isp=1, nsp
         basp=>basis_parameters(isp)
         if (basp%lmxo .ne. -1) cycle loop   ! Species already set
                                             ! in PAO.Basis block
         if (basp%semic) then
            write(6,'(2a)') basp%label,
     .           ' must be in PAO.Basis (it has semicore states)'
            call die()
         endif
         if (basp%bessel) then
            write(6,'(2a)') basp%label,
     .      ' must be in PAO.Basis (it is a floating Bessel function)'
            call die()
         endif
         !
         ! Set the default max l 
         !
         basp%lmxo = basp%ground_state%lmax_valence

         allocate (basp%lshell(0:basp%lmxo))

         if (basp%basis_size(1:2) .eq. 'sz') nzeta = 1
         if (basp%basis_size(1:2) .eq. 'dz') nzeta = 2
         if (basp%basis_size(1:2) .eq. 'tz') nzeta = 3

         loop_l: do l=0, basp%lmxo
            ls=>basp%lshell(l)
            call initialize(ls)
            ls%l = l
            ls%nn = 1
            allocate(ls%shell(1:1))
            s => ls%shell(1)
            call initialize(s)
            s%l = l
            s%n = basp%ground_state%n(l)
            if (basp%ground_state%occupied(l)) then
               s%nzeta = nzeta
            else
               s%nzeta = 0
            endif
            s%polarized = .false.
            if (abs(basp%z).eq.1) then
               s%split_norm = global_splnorm_H
            else
               s%split_norm = global_splnorm
            endif
            s%nzeta_pol = 0

            if (lsoft) then
               s%vcte = softPt 
               s%rinn = -softRc
            else
               s%rinn = 0.d0
               s%vcte = 0.d0
            endif

            ! Default filteret cutoff for shell
            s%filtercut = 0.0d0

            if (s%nzeta .ne.0) then
               allocate(s%rc(1:s%nzeta))
               allocate(s%lambda(1:s%nzeta))
               s%rc(1:s%nzeta) = 0.0d0
               s%lambda(1:s%nzeta) = 1.0d0
            endif
         enddo loop_l

         if (basp%basis_size(3:3) .eq. 'p') then

         ! Polarization orbitals are defined so they have the minimum      
         ! angular momentum l such that there are not occupied orbitals 
         ! with the same l in the valence shell of the ground-state 
         ! atomic configuration. They polarize the corresponding l-1    
         ! shell.

         ! note that we go up to l=4 to make the loop simpler.
         ! (that is the reason why 'occupied' is dimensioned to 0:4)

            select case (basp%basis_size(4:4))
              case (' ', '1')
                 nzeta_pol = 1
              case ('2')
                 nzeta_pol = 2
              case ('3')
                 nzeta_pol = 3
            end select

            loop_angmom: do l=1,4
               if (.not. basp%ground_state%occupied(l)) then
                  ls=>basp%lshell(l-1)
                  s => ls%shell(1)
        
              ! Check whether shell to be polarized is occupied in the gs.
              ! (i.e., whether PAOs are going to be generated for it)
              ! If it is not, mark it for PAO generation anyway.
              ! This will happen for confs of the type s0 p0 dn

                  ! Default filter cutoff for shell
                  s%filtercut = 0.0d0

                  if (s%nzeta == 0) then
                     write(6,"(a,i2,a)") "Marking shell with l=",
     .                l-1, " for PAO generation and polarization."
                     s%nzeta = nzeta
                     allocate(s%rc(1:s%nzeta))
                     allocate(s%lambda(1:s%nzeta))
                     s%rc(1:s%nzeta) = 0.0d0
                     s%lambda(1:s%nzeta) = 1.0d0
                  endif
                  s%polarized = .true.
                  s%nzeta_pol = nzeta_pol

                  exit loop_angmom  ! Polarize only one shell!
               endif
            enddo loop_angmom

         endif

      enddo loop
            
      end subroutine autobasis

      End module basis_specs
