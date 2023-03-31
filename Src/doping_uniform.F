! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_doping_uniform

! Adds a background charge density to simulate doping.
! This routine calculates the net charge of the system, and
! adds a compensating background charge that makes the system 
! neutral. This background charge is constant at points of
! the mesh near the atoms, and zero at points far from the atoms.
! This simulates situations like doped slabs, where the extra
! electrons (holes) are compensated by oposite charges at the
! material (the ionized dopant impurities), but not at the vacuum.

      use precision, only : dp

      implicit none

      private
      integer, save   ::  ntp_nonzero
      logical, save   ::  doping_active
      real(dp), save  ::  qtot
      real(dp), save  ::  threshold_charge_per_spin

      public :: doping_active
      public :: initialize_doping_uniform,
     $          compute_doping_structs_uniform
      public :: doping_uniform

      CONTAINS

      subroutine initialize_doping_uniform()
      use fdf,      only: fdf_get
      use parallel, only: ionode
      use sys,      only: die

      doping_active=fdf_get('SimulateDoping',.false.)
      if (doping_active) then
         if (ionode) then
            write(6,'(/,(a))')
     .        'doping: SimulateDoping = .true. in input file',
     .        'doping: Neutralizing background will be added to points',
     .        'doping: other than vacuum, to simulate doping'
         endif
         qtot = fdf_get('NetCharge',0.0_dp)
      endif
      end subroutine initialize_doping_uniform
!---------------------------------------------------------------

      subroutine compute_doping_structs_uniform(ntpl,rhoatm,nsd)

!     Determines the number of fine mesh points in which a compensating
!     charge will be placed.  The points in which the background charge
!     is removed are those for which the absolute value of the atomic
!     charge rhoatm is larger than a threshold value, defined in this
!     subroutine as thres=0.005 (empirically found to work well). The
!     user can experiment with this at will...

      use alloc,        only : re_alloc
      use precision,    only : dp, grid_p
#ifdef MPI
      use mpi_siesta
#endif

      integer      ::  ntpl         ! Number of local fine points
      integer      ::  nsd         ! Spinor size
      real(grid_p) ::  rhoatm(ntpl)

C     Threshold charge for points where background is added
      real(dp), parameter  ::   thres = 0.005_dp

      integer      :: ntpl_nonzero

#ifdef MPI
      integer           MPIerror
#endif
      ! rhoatm is the charge per spin, hence rhoatm*nsd is the total
      ! atomic charge density at a point
      ! (note that there is no provision for spin
      ! polarization of the atomic charge)

      threshold_charge_per_spin = thres/nsd
      ntpl_nonzero = count(rhoatm > threshold_charge_per_spin)
   
#ifdef MPI
      ! Globalize
      call MPI_AllReduce(ntpl_nonzero,ntp_nonzero,1,MPI_integer,MPI_sum,
     .                   MPI_Comm_World,MPIerror)
      ! We will get the total number of fine points from the mesh modules
#else
      ntp_nonzero = ntpl_nonzero
#endif
      
      end subroutine compute_doping_structs_uniform
!-------------------------------------------------------------------

      subroutine doping_uniform(cell,ntpl,task,rho,rhoatm,add_q)

! Adds a background charge density to simulate doping.  This routine
! adds a compensating background charge that makes the system
! neutral. This background charge is constant at points of the mesh near
! the atoms, and zero at points far from the atoms.  This simulates
! situations like doped slabs, where the extra electrons (holes) are
! compensated by oposite charges at the material (the ionized dopant
! impurities), but not at the vacuum.
! The routine must be first called with 'task=0', and later with
! 'task=1' 
!
! Written by P. Ordejon, July 2009
! Edited by Nick P. Andersen, January 2015
!
! Charges in electrons/bohr**3
! Energies in Rydbergs

      use precision,    only : dp, grid_p
      use mesh,         only : nsp
      use sys,          only : die
      use m_ntm,        only : ntm

!
!     The second is the "main-point" index. There is no spin index as
!     the charges involved are the total ones.
!
      ! Number of fine mesh points in this MPI process
      integer, intent(in)          ::  ntpl

      ! Unit cell vectors
      real(dp), intent(in)         ::  cell(3,3)

!     ! Task :   0 = add background
!                1 = remove background added in a previous call
      integer, intent(in)          ::  task

!     Charge density, to which the background charge is added.
      real(grid_p), intent(inout)  ::  rho(ntpl)
!     Atomic charge density, which determines where to put the charge
      real(grid_p), intent(in)     ::  rhoatm(ntpl)
!     Any additional charge added to the simulation
      real(dp), intent(in), optional :: add_q

      real(dp)          fact, q
      real(dp), external  ::          volcel

      integer             :: ntp   !  total number of fine points in the unit cell

      ntp = product(ntm(1:3))
      q = qtot
      if ( present(add_q) ) q = q + add_q
      if (task .eq. 0) then
        fact = q*(ntp)/(ntp_nonzero*volcel(cell))
      else if (task .eq. 1) then
        fact = -q*(ntp)/(ntp_nonzero*volcel(cell))
      else 
        call die('ERROR: wrong task in doping_uniform.F')
      endif

      where (rhoatm > threshold_charge_per_spin) rho = rho + fact

      end subroutine doping_uniform

      end module m_doping_uniform
