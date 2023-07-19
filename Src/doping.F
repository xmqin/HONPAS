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
      module m_doping
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

      integer, save, dimension(:), pointer  ::  indip => null() 
      integer, save, dimension(:), pointer  ::  indisp => null() 

      integer, save   ::  ntloc, nt, nptot
      logical, save   ::  doping_active
      real(dp), save  ::  qtot

      private
      public :: doping_active
      public :: initialize_doping, compute_doping_structs
      public :: doping

      integer, parameter, public :: ADD_BACKGROUND = 0
      integer, parameter, public :: REMOVE_BACKGROUND = 1


      CONTAINS

      subroutine initialize_doping()
      use parallel, only: ionode
      use m_fdf_global, only: fdf_global_get

      call fdf_global_get(doping_active,'SimulateDoping',.false.)

      if (doping_active) then
         if (ionode) write(6,'(/,(a))')
     .        'doping: SimulateDoping = .true. in input file',
     .        'doping: Neutralizing background will be added to points',
     .        'doping: other than vacuum, to simulate doping'
         call fdf_global_get(qtot,'NetCharge',0.0_dp)
      endif
      end subroutine initialize_doping
!---------------------------------------------------------------

      subroutine compute_doping_structs(np,rhoatm,nsd)

!     Determines in which mesh points the compensating charge will be
!     placed.  The points in which the background charge is removed are
!     those for which the absolute value of the atomic charge rhoatm is
!     larger than a threshold value, defined in this subroutine as
!     thres=0.005 (empirically found to work well). The user can
!     experiment with this at will...

      use alloc,        only : re_alloc
      use precision,    only : dp, grid_p
      use mesh,         only : nsp
      use meshphi,      only : endpht
#ifdef MPI
      use mpi_siesta
#endif

      integer      ::  np, nsd
      real(grid_p) ::  rhoatm(nsp,np)

!     Threshold charge for points where background is added
      real(dp), parameter  ::   thres = 0.005_dp

      integer :: ip, nc, isp
      logical :: it

      real(dp) ::  charge
#ifdef MPI
      integer           MPIerror, npl
#endif

      call re_alloc(indip,1,np*nsp, name="indip",routine="doping")
      call re_alloc(indisp,1,np*nsp, name="indisp",routine="doping")

      ntloc = 0
!  Loop over grid points
      do ip = 1,np
!  Find number of nonzero orbitals at this point
         nc = endpht(ip) - endpht(ip-1)
!  Check if there is a nonzero orbital at that point
         if (nc .ne. 0) then
            do isp = 1,nsp
               it = .false.
!  Check if charge at that point is larger than threshold
               charge = rhoatm(isp,ip)*nsd
               if(abs(charge).ge.thres) it = .true.
               if (it) then
                  ntloc = ntloc+1
                  indip(ntloc) = ip
                  indisp(ntloc) = isp
               endif
            enddo
         endif
      enddo

#ifdef MPI
      call MPI_AllReduce(ntloc,nt,1,MPI_integer,MPI_sum,
     .                   MPI_Comm_World,MPIerror)
      call MPI_AllReduce(np,nptot,1,MPI_integer,MPI_sum,
     .                   MPI_Comm_World,MPIerror)
#else
      nt = ntloc
      nptot = np
#endif
      
      end subroutine compute_doping_structs
!-------------------------------------------------------------------

      subroutine doping(cell,np,task,rho)

! Adds a background charge density to simulate doping.  This routine
! adds a compensating background charge that makes the system
! neutral. This background charge is constant at points of the mesh near
! the atoms, and zero at points far from the atoms.  This simulates
! situations like doped slabs, where the extra electrons (holes) are
! compensated by opposite charges at the material (the ionized dopant
! impurities), but not at the vacuum.
! The routine must be first called with 'task=ADD_BACKGROUND', and later with
! 'task=REMOVE_BACKGROUND' 
!
! Written by P. Ordejon, July 2009
!
! Charges in electrons/borh**3
! Energies in Rydbergs

      use precision,    only : dp, grid_p
      use mesh,         only : nsp
      use sys,          only : die
!
!     Note that the rho array is shaped in such a way that the
!     first index is the "fine point" corresponding to a given main point
!     The second is the "main-point" index. There is no spin index as
!     the charges involved are the total ones.
!
!     Unit cell vectors
      real(dp), intent(in)         ::  cell(3,3)
      ! Number of "big" mesh points in this MPI process
      integer, intent(in)          ::  np

!     ! Task :   0 = add background
!                1 = remove background added in a previous call
!     (Symbolic names are now used and exported for the use of calling routines)
      integer, intent(in)          ::  task
!     Charge density, to which the background charge is added.
      real(grid_p), intent(inout)  ::  rho(nsp, np)

      integer           i,ip, isp
      real(dp)          fact
      real(dp), external  ::          volcel

      if (task .eq. ADD_BACKGROUND) then
        fact = qtot*(nsp*nptot)/(nt*volcel(cell))
      else if (task .eq. REMOVE_BACKGROUND) then
        fact = -qtot*(nsp*nptot)/(nt*volcel(cell))
      else 
        call die('ERROR: wrong task in doping.F')
      endif

!  Loop over grid points where charge has to be added/removed
      do i=1,ntloc
        ip = indip(i)
        isp = indisp(i)
        rho(isp,ip) = rho(isp,ip) + fact
      enddo

      end subroutine doping

      end module m_doping
