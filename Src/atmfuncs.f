! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module atmfuncs

C     This file contains a set of routines which provide all the information
C     about the basis set, pseudopotential, atomic mass, etc... of all the
C     chemical species present in the calculation.

C     The routines contained in this file can only be called after they
C     are initialized by calling the subroutine 'atom' for all the 
C     different chemical species in the calculation:

      use precision
      use sys, only: die
      use atm_types
      use radial, only: rad_get, rad_func
      use spher_harm, only: rlylm

      implicit none 
!
      character(len=79) message
      integer, parameter               :: max_l = 5
      integer, parameter               :: max_ilm = (max_l+1)*(max_l+1)

      real(dp), parameter              :: tiny20=1.e-20_dp
      real(dp), parameter              :: tiny12=1.e-12_dp

      private :: chk, max_l, max_ilm, message

      public  :: nofis, nkbfis, izofis, massfis
      public  :: rcore, rchlocal, rcut, chcore_sub, epskb, uion
      public  :: atmpopfio, psch, zvalfis, floating, psover
      public  :: lofio, symfio, cnfigfio, zetafio, mofio
      public  :: labelfis, lomaxfis, nztfl, rphiatm, lmxkbfis
      public  :: phiatm, all_phi
      public  :: pol   ! Added JMS Dec.2009

      public  :: orb_gindex, kbproj_gindex, vna_gindex, dftu_gindex
      private
      
      contains

      subroutine chk(name,is)
      character(len=*), intent(in) :: name

      integer, intent(in) :: is

      if ((is.lt.1).or.(is.gt.nspecies)) then 
         write(message,'(2a,i3,a,i3)')
     $           name, ": Wrong species", is, ". Have", nspecies
         call die(message)
      endif
      end subroutine chk
!
!
      function floating(is)
!       Returns .true. if the species is really a "fake" one, intended
!       to provide some floating orbitals.
      logical floating
      integer, intent(in) :: is

      floating = izofis(is) .lt. 0
      end function floating

      FUNCTION IZOFIS( IS )
      integer :: izofis ! Atomic number
      integer, intent(in) :: is ! Species index

      call chk('izofis',is)

      izofis = species(is)%z

      end function izofis

      FUNCTION ZVALFIS( IS )
      real(dp) :: zvalfis          ! Valence charge
      integer, intent(in) :: is            ! Species index

      call chk('zvalfis',is)
 
      zvalfis= species(is)%zval
      end function zvalfis
!
      FUNCTION LABELFIS (IS)
      character(len=20) ::  labelfis  ! Atomic label
      integer, intent(in) :: is            ! Species index

      call chk('labelfis',is)
      labelfis= species(is)%label
      end function labelfis

      FUNCTION LMXKBFIS (IS)
      integer :: lmxkbfis    ! Maximum ang mom of the KB projectors
      integer, intent(in) :: is            ! Species index

      call chk('lmxkbfis',is)
      lmxkbfis= species(is)%lmax_projs
      end function lmxkbfis
!
      FUNCTION LOMAXFIS (IS)
      integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
      integer, intent(in) :: is            ! Species index

      call chk('lomaxfis',is)

      lomaxfis = species(is)%lmax_basis
      end function lomaxfis
!
      FUNCTION MASSFIS(IS)
      real(dp) :: massfis            ! Mass
      integer, intent(in) :: is            ! Species index

      call chk('massfis',is)
      massfis=species(is)%mass
      end function massfis
!
      FUNCTION NKBFIS(IS)
      integer :: nkbfis    ! Total number of KB projectors
      integer, intent(in) :: is            ! Species index

      call chk('nkbfis',is)
      nkbfis = species(is)%nprojs
      end function nkbfis
!

      FUNCTION NOFIS(IS)
      integer :: nofis    ! Total number of Basis functions
      integer, intent(in) :: is            ! Species index

      call chk('nofis',is)
      nofis = species(is)%norbs
      end function nofis

      FUNCTION UION ( IS )
      real(dp) uion
      integer, intent(in) :: is    ! Species index
      call chk('uion',is)
      uion = species(is)%self_energy
      end function uion

      FUNCTION RCORE(is)
      real(dp) rcore
      integer, intent(in) :: is    ! Species index

C  Returns cutoff radius of the pseudo-core charge density for the non-linear
C   core corrections for xc potential.
C  Distances in Bohr

      call chk('rcore',is)
      rcore = species(is)%core%cutoff

      end function rcore

      FUNCTION RCHLOCAL(is)
      real(dp) rchlocal
      integer, intent(in) :: is    ! Species index

C  Returns cutoff radius of the Vlocal charge density
C  Distances in Bohr

      call chk('rchlocal',is)
      rchlocal = species(is)%Chlocal%cutoff

      end function rchlocal

!         AMENOFIS
!
!---- Global index helpers------------------------------

      FUNCTION orb_gindex (IS,IO)
      integer orb_gindex
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the global index of a basis orbital

      call chk('orb_gindex',is)
      if ( (io .gt. species(is)%norbs) .or.
     $     (io .lt. 1))   call die("orb_gindex: Wrong io")

      orb_gindex = species(is)%orb_gindex(io)
      end function orb_gindex

      FUNCTION kbproj_gindex (IS,IO)
      integer kbproj_gindex
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! KBproj index 
                                   ! (within atom, <0, for compatibility)

C Returns the global index of a KB projector

      integer :: ko

      call chk('kbproj_gindex',is)
      ko = -io

      if ( (ko .gt. species(is)%nprojs) .or.
     $     (ko .lt. 1)) then
         call die("kbproj_gindex: Wrong io")
      endif

      kbproj_gindex = species(is)%pj_gindex(ko)
      end function kbproj_gindex

      FUNCTION dftu_gindex (IS,IO)
      integer dftu_gindex
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the global index of a DFT+U projector

      call chk('dftu_gindex',is)
      if ( (io .gt. species(is)%nprojsdftu) .or.
     $     (io .lt. 1))   call die("dftu_gindex: Wrong io")

      dftu_gindex = species(is)%pjdftu_gindex(io)
      end function dftu_gindex

      FUNCTION vna_gindex (IS)
      integer vna_gindex
      integer, intent(in) :: is    ! Species index

C Returns the global index for a Vna function

      call chk('vna_gindex',is)
      vna_gindex = species(is)%vna_gindex
      end function vna_gindex
!------------------------------------------------

      FUNCTION ATMPOPFIO (IS,IO)
      real(dp) atmpopfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the population of the atomic basis orbitals in the atomic 
C ground state configuration.

      call chk('atmpopfio',is)
      if ( (io .gt. species(is)%norbs) .or.
     $     (io .lt. 1))   call die("atmpopfio: Wrong io")

      atmpopfio = species(is)%orb_pop(io)
      end function atmpopfio

      FUNCTION CNFIGFIO(IS,IO)
      integer cnfigfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the valence-shell configuration in the atomic ground state
C (i.e. the principal quatum number for orbitals of angular momentum l)

C   INTEGER CNFIGFIO: Principal quantum number of the shell to what 
C                     the orbital belongs ( for polarization orbitals
C                     the quantum number corresponds to the shell which
C                     is polarized by the orbital io) 


      call chk('cnfigfio',is)
      if ( (io .gt. species(is)%norbs) .or.
     $     (io .lt. 1))   call die("cnfigfio: Wrong io")

      cnfigfio = species(is)%orb_n(io)

      end function cnfigfio

      FUNCTION LOFIO (IS,IO)
      integer lofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns total angular momentum quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C                    IO = 0 => Local pseudopotential
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector
      type(species_info), pointer :: spp

      call chk('lofio',is)

      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("lofio: No such orbital")
         lofio = spp%orb_l(io)
      else if (io.lt.0) then
         if (-io.gt.spp%nprojs)  call die("lofio: No such projector")
         lofio = spp%pj_l(-io)
      else
         lofio = 0
      endif

      end function lofio

      FUNCTION MOFIO (IS,IO)
      integer mofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C   Returns m quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C                    IO = 0 => Local pseudopotential
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number m of orbital or KB projector
      type(species_info), pointer :: spp

      call chk('mofio',is)
      
      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("mofio: No such orbital")
         mofio = spp%orb_m(io)
      else if (io.lt.0) then
         if (-io.gt.spp%nprojs)  call die("mofio: No such projector")
         mofio = spp%pj_m(-io)
      else
         mofio = 0
      endif

      end function mofio

      FUNCTION ZETAFIO (IS,IO)
      integer zetafio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C   Returns zeta number of a
C   basis orbital 

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C************************OUTPUT*****************************************
C   INTEGER ZETAFIO  : Zeta number of orbital
      type(species_info), pointer :: spp

      call chk('mofio',is)

      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("zetafio: No such orbital")
         zetafio = spp%orbnl_z(spp%orb_index(io))
      else 
         call die('zetafio only deals with orbitals')
      endif

      end function zetafio

      function rcut(is,io)
      real(dp) rcut
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io> => basis orbitals
                                   ! io<0  => KB projectors
                                   ! io=0 : Local screened pseudopotential

C  Returns cutoff radius of Kleynman-Bylander projectors and
C  atomic basis orbitals.
C  Distances in Bohr
      type(species_info), pointer :: spp

      call chk('rcut',is)
      
      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("rcut: No such orbital")

         rcut = spp%orbnl(spp%orb_index(io))%cutoff
      else if (io.lt.0) then
         if (-io.gt.spp%nprojs)  call die("rcut: No such projector")
         rcut = spp%pjnl(spp%pj_index(-io))%cutoff
      else
         rcut = spp%vna%cutoff
      endif

      end function rcut
!
      
!
      FUNCTION SYMFIO (IS,IO)
      character(len=20) symfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns a label describing the symmetry of the
C   basis orbital or Kleynman-Bylander projector.
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors

C   INTEGER SYMFIO  : Symmetry of the orbital or KB projector
C  2) Returns 's' for IO = 0

      integer ilm, i, lorb, morb
      integer, parameter  :: lmax_sym=4

      character(len=11)  sym_label((lmax_sym+1)*(lmax_sym+1)) 

      data  sym_label(1)
     .  / 's' /
      data (sym_label(i),i=2,4)
     .  / 'py', 'pz', 'px' /
      data (sym_label(i),i=5,9)
     .  / 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2' / 
      data (sym_label(i),i=10,16)
     .  / 'fy(3x2-y2)', 'fxyz', 'fz2y', 'fz3',
     .    'fz2x', 'fz(x2-y2)', 'fx(x2-3y2)' /
      data (sym_label(i),i=17,25)
     .  / 'gxy(x2-y2)', 'gzy(3x2-y2)', 'gz2xy', 'gz3y', 'gz4',
     .    'gz3x', 'gz2(x2-y2)', 'gzx(x2-3y2)', 'gx4+y4' /
      type(species_info), pointer :: spp

      call chk('rcut',is)
      
      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("symfio: No such orbital")
      else if (io.lt.0) then
         if (-io.gt.spp%nprojs)  call die("symfio: No such projector")
      else
         symfio = 's'
      endif

      lorb=lofio(is,io)
      morb=mofio(is,io)

      if(lorb.gt.lmax_sym ) then 
         symfio=' '
      else
         ilm=lorb*lorb+lorb+morb+1  
         if(pol(is,io)) then 
            symfio='P'//sym_label(ilm)
         else
            symfio=sym_label(ilm) 
         endif 
      endif         

      end function symfio
!
!  End of FIOs ----------------------------------------------------
!
      FUNCTION POL (IS,IO)
      logical pol
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io>0 => basis orbitals

C If true, the orbital IO is a perturbative polarization orbital
      type(species_info), pointer :: spp

      spp => species(is)
      if ( (io .gt. species(is)%norbs) .or.
     $     (io .le. 0))   call die("pol: Wrong io")

      pol = spp%orbnl_ispol(spp%orb_index(io))
      
      end function pol

      FUNCTION EPSKB (IS,IO)
      real(dp) epskb
      integer, intent(in)   ::  is   ! Species index
      integer, intent(in)   ::  io   ! KB proyector index (within atom)
                                     ! May be positive or negative 
                                     ! (only ABS(IO) is used).

C  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
C       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C  where Phi_lm is returned by subroutine PHIATM.
C  Energy in Rydbergs.
      type(species_info), pointer :: spp

      integer ik

      spp => species(is)
      ik = abs(io)
      if ((ik.gt.spp%nprojs) .or.
     $    (ik .lt. 1) )  call die("epskb: No such projector")
      epskb = spp%pjnl_ekb(spp%pj_index(ik))

      end function epskb

!--------------------------------------------------------------------
      subroutine vna_sub(is,r,v,grv)
      integer, intent(in) :: is      ! Species index
      real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
      real(dp), intent(out) :: v       ! Value of local pseudopotential
      real(dp), intent(out) :: grv(3)  ! Gradient of local pseudopotential

C Returns local part of neutral-atom Kleynman-Bylander pseudopotential.
C Distances in Bohr,  Energies in Rydbergs
C  2) Returns exactly zero when |R| > RCUT(IS,0)
      type(rad_func), pointer :: func

      real(dp) rmod, dvdr

      call chk('vna_sub',is)

      v = 0.0_dp
      grv(1:3) = 0.0_dp

      if (floating(is)) return

      func => species(is)%vna

      rmod = sqrt(sum(r*r))
      if (rmod .gt. func%cutoff) return

      call rad_get(func,rmod,v,dvdr)
      rmod = rmod + tiny20
      grv(1:3) = dvdr * r(1:3)/rmod
 
      end subroutine vna_sub

      subroutine psch(is,r,ch,grch)
      integer, intent(in) :: is      ! Species index
      real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
      real(dp), intent(out) :: ch      ! Local pseudopot. charge dens.
      real(dp), intent(out) :: grch(3) ! Gradient of local ps. ch. dens.

C Returns 'local-pseudotential charge density'.
C Distances in Bohr, Energies in Rydbergs
C Density in electrons/Bohr**3
C  2) Returns exactly zero when |R| > Rchloc
      type(rad_func), pointer :: func

      real(dp) :: rmod, dchdr

      call chk('psch',is)

      ch = 0.0_dp 
      grch(1:3) = 0.0_dp 

      if (floating(is)) return

      func => species(is)%chlocal
      rmod = sqrt(sum(r*r))
      if (rmod .gt. func%cutoff) return

      call rad_get(func,rmod,ch,dchdr)
      rmod = rmod + tiny20
      grch(1:3) = dchdr * r(1:3)/rmod

      end subroutine psch

      subroutine chcore_sub(is,r,ch,grch)
      integer, intent(in) :: is      ! Species index
      real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
      real(dp), intent(out) :: ch      ! Value of pseudo-core charge dens.
      real(dp), intent(out) :: grch(3) ! Gradient of pseudo-core ch. dens.

C Returns returns pseudo-core charge density for non-linear core correction
C in the xc potential.
C Distances in Bohr, Energies in Rydbergs, Density in electrons/Bohr**3
C  2) Returns exactly zero when |R| > Rcore
      type(rad_func), pointer :: func

      real(dp) rmod, dchdr

      call chk('chcore_sub',is)

      ch = 0.0_dp
      grch(1:3) = 0.0_dp

      if (floating(is)) return

      func => species(is)%core
      rmod = sqrt(sum(r*r))
      rmod = rmod + tiny20                   ! Moved here. JMS, Dec.2012
      if (rmod .gt. func%cutoff) return

      call rad_get(func,rmod,ch,dchdr)
!      rmod=rmod+tiny20                   ! Removed. JMS, Dec.2012
      grch(1:3) = dchdr * r(1:3)/rmod
 
      end subroutine chcore_sub

      subroutine phiatm(is,io,r,phi,grphi)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real(dp), intent(in)  :: r(3)    ! Point vector, relative to atom
      real(dp), intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential
      real(dp), intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C Distances in Bohr
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS and/or IO
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) PHIATM with IO = 0 is strictly equivalent to VNA_SUB
      type(species_info), pointer :: spp
      type(rad_func), pointer :: func

      real(dp) rmod, phir, dphidr
      real(dp) rly(max_ilm), grly(3,max_ilm)
      integer l, m, ik, ilm

      phi = 0.0_dp
      grphi(1:3) = 0.0_dp

      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("phiatm: No such orbital")
         func => spp%orbnl(spp%orb_index(io))
         l = spp%orb_l(io)
         m = spp%orb_m(io)
      else if (io.lt.0) then
         if (floating(is)) return
         ik = -io
         if (ik.gt.spp%nprojs)  call die("phiatm: No such projector")
         func => spp%pjnl(spp%pj_index(ik))
         l = spp%pj_l(ik)
         m = spp%pj_m(ik)
      else     ! io=0
         if (floating(is)) return
         func => spp%vna
         l = 0
         m = 0
      endif

      rmod = sqrt(sum(r*r)) + tiny20
      if(rmod.gt.func%cutoff-tiny12) return

      call rad_get(func,rmod,phir,dphidr)

      if (io.eq.0) then
         phi=phir
         grphi(1:3)=dphidr*r(1:3)/rmod
      else

         ilm = l*l + l + m + 1
         call rlylm( l, r, rly, grly )
         phi = phir * rly(ilm)
         grphi(1)=dphidr*rly(ilm)*r(1)/rmod+phir*grly(1,ilm)
         grphi(2)=dphidr*rly(ilm)*r(2)/rmod+phir*grly(2,ilm)
         grphi(3)=dphidr*rly(ilm)*r(3)/rmod+phir*grly(3,ilm)

      endif

      end subroutine phiatm


      subroutine rphiatm(is,io,r,phi,dphidr)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real(dp), intent(in)  :: r       ! Radial distance, relative to atom
      real(dp), intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential
      real(dp), intent(out) :: dphidr  ! Radial derivative of BO, 
                                     !  KB proj, or Loc pseudopot.

C  Returns the radial component of 
C  Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their radial drivatives)
C Distances in Bohr
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) RPHIATM with ITYPE = 0 is strictly equivalent to VNA_SUB
      type(species_info), pointer :: spp
      type(rad_func), pointer :: func

      real(dp) rmod, phir
      integer l, m, ik

      phi = 0.0_dp
      dphidr = 0._dp

      spp => species(is)
      if (io.gt.0) then
         if (io.gt.spp%norbs)  call die("rphiatm: No such orbital")
         func => spp%orbnl(spp%orb_index(io))
         l = spp%orb_l(io)
         m = spp%orb_m(io)
      else if (io.lt.0) then
         if (floating(is)) return
         ik = -io
         if (ik.gt.spp%nprojs)  call die("rphiatm: No such projector")
         func => spp%pjnl(spp%pj_index(ik))
         l = spp%pj_l(ik)
         m = spp%pj_m(ik)
      else
         if (floating(is)) return
         func => spp%vna
         l = 0
         m = 0
      endif

      rmod = r + tiny20
      if(rmod.gt.func%cutoff-tiny12) return

      call rad_get(func,rmod,phir,dphidr)

      if (l.eq.0) then
         phi=phir
      elseif (l.eq.1) then
         phi=phir*r
         dphidr=dphidr*r
         dphidr=dphidr+phir 
      else
         phi=phir*r**l 
         dphidr=dphidr * r**l
         dphidr=dphidr + l * phir * r**(l-1)
      endif

      end subroutine rphiatm


      subroutine all_phi( is, it, r, maxnphi, nphi, phi, grphi )
      integer,   intent(in) :: is     ! Species index
      integer,   intent(in) :: it     ! Orbital-type switch:
                                      ! IT > 0 => Basis orbitals
                                      ! IT < 0 => KB projectors
      real(dp),  intent(in) :: r(3)   ! Point vector, relative to atom
      integer,  intent(in) :: maxnphi ! Maximum number of phi's
      integer,  intent(out) :: nphi   ! Number of phi's
      real(dp), intent(out) :: phi(maxnphi) ! Basis orbital, KB projector, or
                                      !  local pseudopotential
      real(dp), optional, intent(out) :: grphi(3,maxnphi) ! Gradient of phi

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C  Same as phiatm but returns all orbitals or KB projectors of the atom
C  Written by D.Sanchez-Portal and J.M.Soler. Jan. 2000 

C Distances in Bohr
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 8) If arrays phi or grphi are too small, returns with the required
C    value of nphi
      type(species_info), pointer :: spp

      integer i, jlm, l, lmax, m, maxlm
      double precision  rmod, phir, dphidr
      real(dp) rly(max_ilm), grly(3,max_ilm)

      integer :: ilm(maxnphi)
      double precision :: rmax(maxnphi)
      logical :: within(maxnphi)

      call chk('all_phi',is)
      spp => species(is)

!     Find number of orbitals
      if (it.gt.0) then
        nphi=spp%norbs
      elseif (it.lt.0) then
        nphi=spp%nprojs
      else
         call die("all_phi: Please use phiatm to get Vna...")
      endif
      
      if (nphi.gt.maxnphi) call die('all_phi: maxphi too small')

      if (it.gt.0) then
         do i = 1, nphi
            l = spp%orb_l(i)
            m = spp%orb_m(i)
            ilm(i) = l*(l+1)+m+1
            rmax(i) = spp%orbnl(spp%orb_index(i))%cutoff
         enddo
      else
         do i = 1, nphi
            rmax(i) = spp%pjnl(spp%pj_index(i))%cutoff
            l = spp%pj_l(i)
            m = spp%pj_m(i)
            ilm(i) = l*(l+1)+m+1
         enddo
      endif

!     Initialize orbital values
      phi(1:nphi) = 0._dp
      if (present(grphi)) grphi(:,1:nphi) = 0._dp

      if ((it.lt.0) .and. floating(is)) return

!     Find for which orbitals rmod < rmax and test for quick return
      rmod = sqrt(sum(r*r)) + tiny20
      within(1:nphi) = ( rmax(1:nphi) > rmod )
      if (.not.any(within(1:nphi))) return

!     Find spherical harmonics
      maxlm = maxval( ilm(1:nphi), mask=within(1:nphi) )
      lmax=nint(sqrt(real(maxlm,dp)))-1
      call rlylm(lmax,r,rly,grly)

!     Find values

      i_loop: do i=1,nphi

!       Check if rmod > rmax
        if (.not.within(i)) cycle i_loop
          
!       Find radial part
        if (it.gt.0) then
          call rad_get(spp%orbnl(spp%orb_index(i)),rmod,phir,dphidr)
        else
          call rad_get(spp%pjnl(spp%pj_index(i)),rmod,phir,dphidr)
        end if

!       Multiply radial and angular parts
        jlm = ilm(i)
        phi(i) = phir * rly(jlm)
        if (present(grphi))
     .    grphi(:,i) = dphidr * rly(jlm) * r(:) / rmod + 
     .                 phir * grly(:,jlm)

      enddo i_loop

      end subroutine all_phi
!
!
!     This routine takes two species arguments
!
      subroutine psover(is1,is2,r,energ,dedr)
      integer, intent(in) :: is1, is2     ! Species indexes
      real(dp), intent(in)  :: r       ! Distance between atoms
      real(dp), intent(out) :: energ   ! Value of the correction
                                     !  interaction energy
      real(dp), intent(out) :: dedr    ! Radial derivative of the correction

C Returns electrostatic correction to the ions interaction energy
C due to the overlap of the two 'local pseudopotential charge densities'
C Distances in Bohr, Energies in Rydbergs
C  2) Returns exactly zero when |R| > Rchloc
      type(rad_func), pointer :: func

      integer ismx, ismn, indx
      real(dp) r_local

      call chk('psover',is1)
      call chk('psover',is2)
      
      energ=0.0_dp 
      dedr=0.0_dp 
      
      if (floating(is1) .or. floating(is2)) return

      ismx=max(is1,is2)
      ismn=min(is1,is2)
      indx=((ismx-1)*ismx)/2 + ismn
      func => elec_corr(indx)

      if ( r .gt. func%cutoff - tiny12 ) return

      call rad_get(func,r,energ,dedr)
      r_local = r+tiny20
      energ=2.0_dp*energ/r_local
      dedr=(-energ + 2.0_dp*dedr)/r_local

      end subroutine psover

!
!     Deprecated
!
      FUNCTION NZTFL (IS,L)
      integer nztfl
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the basis funcs
C Returns the number of different basis functions
C with the same angular momentum and for a given species
      type(species_info), pointer :: spp

      integer i

      call chk('nztfl',is)
      spp => species(is)

      nztfl = 0
      do i = 1, spp%norbs
         if (spp%orb_l(i).eq.l) nztfl = nztfl+1
      enddo

      end function nztfl

      FUNCTION NKBL_FUNC (IS,L)
      integer nkbl_func
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the basis funcs

C Returns the number of different KB projectors
C with the same angular momentum and for a given species
      type(species_info), pointer :: spp

      integer i

      call chk('nkbl_func',is)

      spp => species(is)

      nkbl_func = 0
      do i = 1, spp%nprojs
         if (spp%pj_l(i).eq.l) nkbl_func = nkbl_func+1
      enddo

      end function nkbl_func

      end module atmfuncs





