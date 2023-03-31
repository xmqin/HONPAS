! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.

subroutine vacuum_level( np, nspin, rho, V, npvac, Vmax, Vmean )

  ! Finds the maximum and mean potential in the vacuum region, defined as
  ! the set of grid points (in all processors) with zero electron density. 
  ! Returns zero if there are no empty points. J.M.Soler. Dec.2009

  use precision,        only: dp            ! Double precision real kind
  use precision,        only: gp => grid_p  ! Real kind of grid arrays
  use moreParallelSubs, only: miscAllReduce ! Parallel reduction of misc. arrays
  use moreParallelSubs, only: miscAllReduceInt ! Integer version

  implicit none
  integer, intent(in) :: np            ! Number of grid points in my processor
  integer, intent(in) :: nspin         ! Number of spin components
  real(gp),intent(in) :: rho(np,nspin) ! Electron (spin) density at grid points
  real(gp),intent(in) :: V(np,nspin)   ! Effective (spin) potential
  integer, intent(out):: npvac         ! Number of points in vacuum region
                                       ! (in all processors)
  real(dp),intent(out):: Vmax          ! Max. eff. potential in vaccuum region
  real(dp),intent(out):: Vmean         ! Mean eff. potential in vaccuum region

  real(dp),parameter:: rho_min = 1.e-30_dp  ! Min. density of non-empty points

  integer :: ip, ns
  real(dp):: rho_tot, Vsum

  ns = min( nspin, 2 )  ! Number of diagonal components of local dens. matrix
  npvac = 0             ! Number of points in vacuum (empty) region
  Vmax = -huge(Vmax)    ! Max. potential in vacuum region
  Vsum = 0              ! Sum of potential at vacuum points
  do ip = 1,np                  ! Loop over all my processor points
    rho_tot = sum(rho(ip,1:ns)) ! Total electron density at this point
    if (rho_tot < rho_min) then ! Vacuum (empty) region
      npvac = npvac + 1
      Vmax = max( Vmax, maxval(V(ip,1:ns)) )
      Vsum = Vsum + sum(V(ip,1:ns))
    end if
  end do ! ip

  call miscAllReduceInt( 'sum', npvac )  ! Add vacuum points of all processors
  call miscAllReduce( 'sum', Vsum )   ! Add vacuum potential of all processors
  call miscAllReduce( 'max', Vmax )   ! Get max. vac. pot. of all processors

  if (npvac>0) then             ! Some vacuum points
    Vmean = Vsum / npvac / ns   ! Mean potential over vacuum points and spins
  else                          ! No vacuum region
    Vmax = 0
    Vmean = 0
  end if

end subroutine vacuum_level

