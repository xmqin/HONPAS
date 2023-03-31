! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_born_charge

      private
      public :: born_charge

      CONTAINS
      
      subroutine born_charge()
      use precision
      use siesta_options, only:dx
      use sparse_matrices
      use parallel,       only: Node, Nodes, IOnode
      USE m_ksvinit,      only: nkpol, kpol, wgthpol
      use siesta_geom,    only: na_u, na_s, xa, scell, ucell, isa
      use atomlist,       only: no_u, indxuo, iphorb, iaorb, lasto
      use atomlist,       only: rmaxo, no_s, no_l
      use m_spin,         only: nspin, SPpol, NonCol, SpOrb ! CC RC Added
      use m_ksv,          only: KSV_pol
      use siesta_geom,    only: shape
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none

      real(dp) :: polxyz(3,nspin)  ! Automatic, small
      real(dp) :: polR(3,nspin)    ! Automatic, small

      real(dp) :: qspin(4)         ! Local variable

      integer  :: j, ind, io, ispin
#ifdef MPI
      real(dp):: qtmp(4)   ! Temporary for globalized spin charge
#endif

!     Computes Born-effective charges by finite differences,
!     using the "force constant calculation" mode of operation.
!     It is called only for those steps (even) in which the dx
!     is positive.
!
      if (NonCol .or. SpOrb) then
       if ( IONode ) then
          write(*,'(/,a)')
     .         'born_charge: implemented only for spin-polarized '
          write(*,'(a)')
     .         'born_charge: or paramagnetic calculations'
       end if
       
       return
       
      endif

      if (nkpol.lt.1) then
         if (IOnode) write(6,'(/,a,f12.6)')
     .        'siesta: specify polarization grid for BC calculation'
         if (IOnode) write(6,'(a,f12.6)')
     .        'siesta: The Born charge matrix will not be calculated'
         RETURN
      endif
      if (IOnode) write(6,'(/,a,f12.6)')
     .     'siesta: Calculating polarization. '

!     Find total population of spin up and down
      if (nspin .ge. 2) then
         do ispin = 1,nspin
            qspin(ispin) = 0.0_dp
            do io = 1,no_l
               do j = 1,numh(io)
                  ind = listhptr(io) + j
                  qspin(ispin) = qspin(ispin)
     .                 + Dscf(ind,ispin)*S(ind)
               enddo
            enddo
         enddo
#ifdef MPI
!     Global reduction of spin components
         call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
         qspin(1:nspin) = qtmp(1:nspin)
#endif
      endif
      if (nkpol.gt.0) then
         ! Note use of xa instead of xa_last, since this routine
         ! is called at every geometry step, before moving the atoms.
         !
         call KSV_pol(na_u, na_s, xa, rmaxo, scell, ucell,
     .        no_u, no_l, no_s, nspin, qspin,
     .        maxnh, nkpol, numh, listhptr, listh,
     .        H, S, xijo, indxuo, isa, iphorb,
     .        iaorb, lasto, shape,
     .        nkpol,kpol,wgthpol, polR, polxyz)
      endif
      if (nkpol.gt.0.and.IOnode) then
         call obc( polR, ucell, dx, nspin)
      endif

      end subroutine born_charge

      end module m_born_charge
