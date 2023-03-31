! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine savepsi(psiprev, psi, nuo, nuotot, nocc)
C *********************************************************************
C Copies psi into psiprev allowing for the fact that the parallel form
C requires a re-distribution of the data.
C Written by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)     : Wavefunctions in current k point
C integer nuo                   : Number of (local) orbitals in the unit cell
C integer nuotot                : Number of orbitals in the unit cell
C integer nocc                  : number of occupied states
C *************************** OUTPUT **********************************
C real*8  psiprev(2,nuo,nuotot) : Wavefunctions in previous k point
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
 
      integer nuo, nuotot, nocc

      real(dp)
     .  psiprev(2,nuo,nuotot), 
     .  psi(2,nuotot,nuo)

C**** Internal variables ***********************************************

      integer 
     .  iuo, juo
#ifdef MPI
      integer 
     .  MPIerror, iuog, juog, n, noccloc
      real(dp), dimension(:,:,:), pointer ::  psitmp
#endif

#ifdef MPI
! AG
! Allocate as in detover, using the number of orbitals on the first node,
! as some of the nodes might have zero orbitals!

      call GetNodeOrbs(nocc,0,Nodes,noccloc)
      nullify( psitmp )
      call re_alloc( psitmp, 1, 2, 1, nuotot, 1, noccloc,
     &                 name='psitmp', routine='savepsi' )

      do n = 0,Nodes-1

C Broadcast copy of psi on node n to all other nodes
        call GetNodeOrbs(nocc,n,Nodes,noccloc)

        if (Node .eq. n) then
           psitmp(1:2,1:nuotot,1:noccloc) = psi(1:2,1:nuotot,1:noccloc)
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccloc,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

C Save local part of psiprev
        do iuo = 1,noccloc
          call LocalToGlobalOrb(iuo,n,Nodes,iuog)
          do juo = 1,nuo
            call LocalToGlobalOrb(juo,Node,Nodes,juog)
            psiprev(1,juo,iuog) = psitmp(1,juog,iuo)
            psiprev(2,juo,iuog) = psitmp(2,juog,iuo)
          enddo
        enddo


      enddo
      call de_alloc( psitmp,  name='psitmp' , routine='savepsi')

#else
C Straight serial copy

      psiprev(1:2,1:nuotot,1:nocc) = psi(1:2,1:nuotot,1:nocc)

#endif

      end
