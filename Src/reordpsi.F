! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine reordpsi(psisave, psi, nuo, nuotot, nocc, ninc)
C *********************************************************************
C Reorders the coefficient of the bands, so the bands whose bands index
C range between 1 and numincbands are the bands included for
C wannierization.
C Written by J. Junquera, October 2013, based on subroutine savepsi,
C by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)      : Wavefunctions in current k point
C integer nuo                    : Number of (local) orbitals in the unit cell
C integer nuotot                 : Number of orbitals in the unit cell
C integer nocc                   : Number of bands considered for
C                                  Wannierization (before excluding bands)
C integer ninc                   : Number of maximum bands whose
C                                  coefficients will be stored in a node
C *************************** OUTPUT **********************************
C complex(dp)  psisave(nuotot,ninc) : Wavefunctions saved in the matrix of
C                                  coefficients
C                                  In psi (and last two in psisave):
C                                  First  index: real or complex value
C                                  Second index: atomic orbital
C                                  Third  index: band index
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision, only: dp
      use m_siesta2wannier90, only : isexcluded
#ifdef MPI
      use alloc,              only : re_alloc, de_alloc
      use parallel,           only : Node, Nodes, BlockSize
      use parallelsubs,       only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,       only : GlobalToLocalOrb
      use mpi_siesta
      use m_orderbands,       only : index_included_band
      use m_orderbands,       only : node_included_band
      use m_orderbands,       only : band_index_in_node
      use m_orderbands,       only : which_band_in_node
#endif

      implicit none
 
      integer, intent(in) ::  nuo, nuotot, nocc, ninc
      real(dp), intent(in) ::   psi(2,nuotot,nuo)
      complex(dp), intent(out) :: psisave(nuotot,ninc)



C**** Internal variables ***********************************************

      integer :: jband, jbandloc, iuo, juo, iband

#ifdef MPI
      integer 
     .  MPIerror, iuog, juog, n, noccloc
      real(dp), dimension(:,:,:), pointer ::  psitmp
#endif

#ifdef MPI
! AG
! Allocate as in detover, using the number of orbitals on the first node,
! as some of the nodes might have zero orbitals!
!
! Actually, this works because the first node handles the largest number
! of orbitals when using the block-cyclic distribution.

      call GetNodeOrbs(nocc,0,Nodes,noccloc)

      nullify( psitmp )
      call re_alloc( psitmp, 1, 2, 1, nuotot, 1, noccloc,
     &                 name='psitmp', routine='reordpsi' )

      jband = 0
      do n = 0, Nodes-1

!       Broadcast copy of psi on node n to all other nodes
!       Note how not all the bands are broadcasted. 
!       Only the "occupied" bands in the corresponding node. 
        call GetNodeOrbs(nocc,n,Nodes,noccloc)

        if (Node .eq. n) then
           psitmp(1:2,1:nuotot,1:noccloc) = psi(1:2,1:nuotot,1:noccloc)
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccloc,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

!       Save local part of psisave
        do iuo = 1, noccloc
          call LocalToGlobalOrb(iuo,n,Nodes,iuog)

!         Select here if the local occupied band is included or not.
!         If it is included, then a new index (jband) that ranges
!         between 1 and numincbands(ispin) is assigned 
!
!         This is the same logic as in m_orderbands, so the arrays
!         can be re-used

          if( .not.  isexcluded(iuog) ) then
            jband    = jband + 1
!           The coefficients of the included band will be stored
!           in Node (node_included_band(jband)),
!           and will be the band_index_in_node(jband)-th band
!           stored on that Node
            jbandloc = band_index_in_node(jband)
            if( node_included_band(jband) .eq. Node ) then
              do juo = 1, nuotot
                psisave(juo,jbandloc) =
     $             cmplx(psitmp(1,juo,iuo),psitmp(2,juo,iuo),kind=dp)
              enddo ! atomic orbitals
            endif   ! If the local Node is the same where the included
                    !    band should be stored
          endif     ! Is the band included for wannierization?
        enddo       ! local bands on a given node
      enddo         ! nodes
      call de_alloc( psitmp,  name='psitmp' , routine='reordpsi')
#else
!     Straight serial copy
      jband = 0
      do iband = 1, nocc
        if (.not. isexcluded(iband)) then
          jband = jband + 1
          psisave(:,jband) =
     $        cmplx(psi(1,:,iband),psi(2,:,iband),kind=dp)
        endif
      enddo
#endif


      end subroutine reordpsi

