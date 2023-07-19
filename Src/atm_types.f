!
! This file is part of the HOPAS package.
!
! Copyright USTC Yang group
!
! Use of this software constitutes agreement with the full conditions
! given in the HONPAS license, as signed by all legitimate users.
!
! Origin from SIESTA program,
C Modified by shanghui, 2010
C Modified by xmqin, November 2013


      module atm_types

      use precision, only: dp
      use radial, only: rad_func
!
!     Derived types for orbitals and KB projectors
!
      implicit none
!
!     Storage of orbital and projector real-space tables and other
!     characteristics
!
!     These parameters are over-dimensioned, but there is no storage
!     penalty, as the real information is packed and indexed.
!
      integer, parameter, public  :: maxnorbs = 100
!       Maximum number of nlm orbitals
!
      integer, parameter, public  :: maxn_pjnl = 10
!       Maximum number of projectors (not counting different "m" copies)
      integer, parameter, public  :: maxn_orbnl = 20
!       Maximum number of nl orbitals (not counting different "m" copies)
      integer, parameter, public  :: maxnprojs = 50
!       Maximum number of nlm projectors
!

C----------------------------add for nao2gto-----------------------
C NAO2GTO: NAO = sum RGTOs = sum CGTOs
C RGTOs :  Real Spherical Harmonic GTOs or slater-type GTOs.
C CGTOs:   Cartesian GTOs.

      integer, parameter, public :: maxnorbs_cphi = 120 ! larger than maxnorbs (PAOs of an atom)
C number of CGTOs
      integer, parameter, public :: maxn_contract = 10  ! Max number of RGTOs
      integer, parameter, public :: ncon_max = 60
c  support Max 6 GTOs for f orbitals, 10 GTOs for d orbitals. 
C  considered magnetic quantum number m for RGTOs
      integer, parameter, public :: l_max = 3
      integer, public :: nco(0:l_max),nso(0:l_max),
     .         ncosum(-1:l_max),co(0:l_max,0:l_max,0:l_max),
     .         coset(0:l_max,0:l_max,0:l_max), indco(3,20) 

C CO: CGTOS so: RGTOs
C 1s+3p+6d =10
C 1s+3p+6d+10f=20
C ncosum(lmax)=sum_l (l+1)(l+2)/2 ; l=0,  l_max

c Calculation of the spherical harmonics and the corresponding orbital
!>      transformation matrices.
!> \par Literature
!>      H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)
!> \par History

      TYPE, public :: orbtramat_type
      REAL(dp), DIMENSION(:,:), POINTER :: c2s
      END TYPE orbtramat_type

      TYPE(orbtramat_type), DIMENSION(:), POINTER, public :: orbtramat

!----------------------------end add-------------------------------


!     Species_info: Consolidate all the pieces of information in one place
!
      type, public :: species_info
         character(len=2)                ::  symbol
         character(len=20)               ::  label
         integer                         ::  z          ! Atomic number
         real(dp)                        ::  mass
         real(dp)                        ::  zval       ! Valence charge
         real(dp)                        ::  self_energy !Electrostatic
                                                         !self-energy
!
!        Orbitals
!             We keep track of just one orbital for each
!             "nl" family
!
         integer                         ::  n_orbnl    ! num of nl orbs
         integer                         ::  lmax_basis ! basis l cutoff
         integer, dimension(maxn_orbnl)  ::  orbnl_l    ! l of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_n    ! n of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_z    ! z of each nl orb
         logical, dimension(maxn_orbnl)  ::  orbnl_ispol! is it a pol. orb?

         real(dp),
     $            dimension(maxn_orbnl)  ::  orbnl_pop  ! pop. of nl orb
                                                        ! (total of 2l+1
                                                        ! components)
!

C sphi: Slater-type orbital ! spherical
c cphi : Cartesian orbital
!-----------------------shanghui add for ------------------------------
         integer, dimension(maxn_orbnl)  ::  orbnl_contract
         integer, dimension(maxn_orbnl)  ::  orbnl_index_cphi
         integer, dimension(maxn_orbnl)  ::  orbnl_index_sphi
         real(dp),dimension(maxn_contract,maxn_orbnl) :: orbnl_zeta
         real(dp),
     $   dimension(maxn_contract,maxn_orbnl) :: orbnl_coefficient

         real(dp),
     $   dimension(maxn_contract,maxn_orbnl) :: pgf_radius
         real(dp),dimension(maxn_orbnl) ::   shell_radius

         real(dp) :: kind_radius

! added by xmqin, zeta and coeff of the most diffuse GTO. Normalized adjoint GTO
! NAO = sum Cm*GTOm

         real(dp), dimension(maxn_orbnl) ::  orbnl_adjoined_zeta ! min zeta, the most diffuse gto         
         !real(dp), dimension(maxn_orbnl) ::  orbnl_diff_coeff ! the corresponding coefficient

! added by xmqin, max contraction coeffecient:  C: cartesian, R: spherical
! shell : same n and l,  
! CGTO (nco*k)  ------------->    RGTO (nso*k) ------------------------------->  PAO (nso)
!                c2s(nco, nso)                 Fit: Sum_k D_k*RGTO_k,  k = 0, M 
!     RGTO (nso) = T[c2s(nco, nso)] * CGTO (nco),  matrix * vector, here T is "Transpose" 
!     RGTO (io)   = sum c2s(1:nco, io) * CGTO(1:nco)
!
!     ERI (nsoa*nsob*nsoc*nsod),  nso PAO once !
! Contribution of each CGTO shell (nco) to PAO shell (nso) ERIs : 
!
!     Sum_k D_k* c2s_k(nco,nso)*GTO_k(nco) ), k=1, M
!
! Actually, both D_k and normorlized coeff have been involved in c2s_k (nco, nso), 
! We calc the trans matrix sphi(nco*M, nso) !
!
! So max contribution coeff of echo CGTO shell is:
!  max( sum [sphi(1: nco, i)],  i = 1, nso ) !
         real(dp), dimension(maxn_orbnl) ::  orbnl_contraction_coeff
!
!------------------------------------------------------------------------!       
         !----------------cphi's data------------------------!
         integer                         ::  norbs_cphi
         integer,dimension(maxnorbs_cphi)  :: orb_cphi_lx
         integer,dimension(maxnorbs_cphi)  :: orb_cphi_ly
         integer,dimension(maxnorbs_cphi)  :: orb_cphi_lz
         real(dp),dimension(maxnorbs_cphi) :: norm_cphi
         integer, dimension(maxnorbs_cphi) :: orb_index_cphi
         integer, dimension(maxnorbs_cphi) ::  orb_n_cphi
         integer, dimension(maxnorbs_cphi) ::  orb_l_cphi
         real(dp),dimension(ncon_max,maxnorbs_cphi) :: cphi
         !---------------end cphi's data---------------------!

         !---------------sphi's data-------------------------!
         !-----In sesta, in fact all orb is sphi,------------!
         !-----here we add sphi(:,:)-------------------------!
         real(dp),dimension(ncon_max,maxnorbs) :: sphi
         !---------------end cphi's data---------------------!
!--------------------- end add---------------------------------------




!        Projectors
!             For each l, there can be several projectors. Formally, we 
!             can can use the "nl" terminology for them. n will run from
!             1 to the total number of projectors at that l.
!             
!
         integer                         ::  n_pjnl     ! num of "nl" projs
         integer                         ::  lmax_projs ! l cutoff for projs
         integer, dimension(maxn_pjnl)   ::  pjnl_l     ! l of each nl proj
         integer, dimension(maxn_pjnl)   ::  pjnl_n     ! n of each nl proj
         real(dp), dimension(maxn_pjnl)
     $                                   ::  pjnl_ekb   ! energy of
                                                         ! each nl proj
!
!                        ------------------------------
!
!        Aggregate numbers of orbitals and projectors (including 2l+1
!        copies for each "nl"), and index arrays keeping track of
!        which "nl" family they belong to, and their n, l, and m (to avoid
!        a further dereference)
!
         integer                         ::  norbs
         integer, dimension(maxnorbs)    ::  orb_index
         integer, dimension(maxnorbs)    ::  orb_n
         integer, dimension(maxnorbs)    ::  orb_l
         integer, dimension(maxnorbs)    ::  orb_m
         real(dp),
     $            dimension(maxnorbs)    ::  orb_pop   ! pop. of nl orb

         integer                         ::  nprojs
         integer, dimension(maxnprojs)   ::  pj_index
         integer, dimension(maxnprojs)   ::  pj_n
         integer, dimension(maxnprojs)   ::  pj_l
         integer, dimension(maxnprojs)   ::  pj_m
!
         type(rad_func), dimension(:), pointer       ::  orbnl
         type(rad_func), dimension(:), pointer       ::  pjnl
         type(rad_func)                              ::  vna
         type(rad_func)                              ::  chlocal
         type(rad_func)                              ::  reduced_vlocal
         logical                                     ::  there_is_core
         type(rad_func)                              ::  core

         logical                        :: read_from_file
      end type species_info

!
      integer, save, public             :: nspecies
      integer, save, public             :: npairs

      type(species_info), target, allocatable,
     $                            save, public   ::  species(:)
      type(rad_func), allocatable, target,
     $                            save, public   ::  elec_corr(:)
!

      private

      end module atm_types
