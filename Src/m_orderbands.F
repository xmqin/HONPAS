! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      module m_orderbands

      integer, dimension(:),   save, public, pointer :: 
     .   index_included_band
!     This array has the length of the total number of included bands 
!     for wannierization.
!     It points to the global index of the band that is included.
!     Depending on the number of included bands and the number of nodes
!     on which the calculation is run, 
!     it might not be ordered sequentially.
!     This means that the 
!     index_included_bands(1) is usually the lowest energy band for a
!     given k-point,
!     but 
!     index_included_bands(2) can be other band in the list
!     (not necesarilly the next lower energy band).
      integer, dimension(:),   save, public, pointer ::
     .   node_included_band
!     This array has the length of the total number of included bands 
!     for wannierization.
!     It points to the node where the coefficients of the included band
!     will be stored
      integer, dimension(:),   save, public, pointer ::  
     .   node_diagon_all
!     This array has the length of all the bands in the unit cell.
!     It tells you on which node a given band is stored
!     when all the bands are considered!
!     (i. e. before the eventual exclusion of some bands)
      integer, dimension(:),   save, public, pointer :: 
     .   band_index_in_node
!     This array has the length of the total number of included bands 
!     for wannierization.
!     For each included band, it tells you what is the local index of
!     the band in the node.     
      integer, dimension(:,:), save, public, pointer::
     .   which_band_in_node
!     First dimension of the array: number of nodes
!     Second dimension of the array: number of bands stored on this node
!     Given a node and a local index of the band within the node,
!     this array tells you the global index of the band stored there
      integer, dimension(:),   save, public, pointer ::
     .   sequential_index_included_bands
!     This array has the length of the total number of bands, i.e.,
!     the number of atomic orbitals in the unit cell.
!     If the band is not included for wannierization, a zero is
!     assigned.
!     If the band is included for wannierization, then a sequential
!     index is assigned, from 1 to the number of included bands.
!     1 is for the lowest band that will be included for wannierization.
!     2 for the second lowest band
!     and so on

#ifdef MPI
      CONTAINS

      subroutine order_index( nuo, nuotot, numincbands )

C *********************************************************************
C Populates some arrays that identify in which node the different bands
C are stored.
C See explanation of the different arrays above
C Written by J. Junquera, May 2015
C **************************** INPUT ********************************** 
C integer nuo                    : Number of (local) orbitals in the unit cell
C integer nuotot                 : Number of orbitals in the unit cell
C integer numincbands            : Number of bands considered for
C                                  Wannierization (after excluding bands)
C *************************** OUTPUT **********************************
C The different arrays declared above in the module
C *********************************************************************

      use m_siesta2wannier90, only : isexcluded
      use m_siesta2wannier90, only : blocksizeincbands
      use m_siesta2wannier90, only : nincbands_loc

      use parallel,           only : Node, Nodes
      use parallelsubs,       only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,       only : GlobalToLocalOrb
      use alloc,              only : re_alloc, de_alloc

!     For debugging
      use sys,                only : die
!     End debugging

      implicit none 

      integer, intent(in)  :: nuo           ! Number of bands stored in
                                            !    the local node
      integer, intent(in)  :: nuotot        ! Total number of bands
      integer, intent(in)  :: numincbands   ! Number of bands considered
                                            !    for wannierization

      integer 
     .  jband, nodejband

      integer :: n           ! Counter for loops on nodes
      integer :: noccloc     ! Number of bands in the local node
                             !    (before excluding bands for diagonalization)
      integer :: ibandl      ! Counter for loops on the local bands 
                             !    on a given node
      integer :: iband_included ! Counter for loops on included bands
      integer :: Nodeibandg  ! On which node will be stored the
                             !    coefficients of a particular band after 
                             !    the exclusion of all not considered bands
      integer :: icounter    ! Counter to count the number of included
                             !    bands for wannierization
      integer :: ibandg      ! Counter for loops on all the bands (global)
      integer :: maxilocalb  ! Maximum number of ilocal bands
     

      integer, allocatable, dimension(:) ::  bandspernode ! Counts how many 
                                                          !   bands are stored
                                                          !   on each node
      integer 
     .  MPIerror
      integer, external :: indxg2p


!     Define the maximum number of local bands to be stored in a Node
!     after wannierization.
      maxilocalb = numincbands/Nodes + 1

!     Allocate the arrays where the indices of the bands included for 
!     wannierization, and the node that will contain its coefficients 
!     will be references
      nullify( index_included_band )
      call re_alloc( index_included_band, 1, numincbands, 
     &               name='index_included_band', routine='order_index' )

      nullify( node_included_band )
      call re_alloc( node_included_band, 1, numincbands, 
     &               name='node_included_band', routine='order_index' )

      nullify( node_diagon_all )
      call re_alloc( node_diagon_all, 1, nuotot, 
     &               name='node_diagon_all', routine='order_index' )

      nullify( band_index_in_node )
      call re_alloc( band_index_in_node, 1, numincbands,
     &               name='band_index_in_node', routine='order_index' )

      nullify( which_band_in_node )
      call re_alloc( which_band_in_node, 0, Nodes-1, 1, maxilocalb,
     &               name='which_band_in_node', routine='order_index')

      nullify( sequential_index_included_bands )
      call re_alloc( sequential_index_included_bands, 1, nuotot,
     &               name='sequential_index_included_bands', 
     &               routine='order_index' )

      if (.not. allocated(bandspernode)) then
        allocate(bandspernode(0:Nodes-1))
      endif

      do n = 0, Nodes - 1
        call GetNodeOrbs( nuotot, n, Nodes, noccloc )
!! For debugging
!        write(6,'(a,4i5)')' order_index: Node, n, nuotot, noccloc = ', 
!     .                       Node, n, nuotot, noccloc
!! End debugging
        do ibandg = 1, nuotot
          call GlobalToLocalOrb(ibandg,n,Nodes,ibandl)
          if( ibandl .ne. 0 ) then
            node_diagon_all(ibandg) = n
          endif
        enddo
      enddo 

!! For debugging
!        do ibandg = 1, nuotot
!          write(6,'(a,4i5)')
!     .     ' order_index: Node, ibandg, node_diagon_all=', 
!     .                    Node, ibandg, node_diagon_all(ibandg)
!        enddo 
!! End debugging

!     Set up the indices of the included bands
!
!     Note that, in parallel operation, the Scalapack distribution
!     results in each node handling all the coefficients of no_l
!     eigenvectors, where no_l is the number of locally handled
!     orbitals. (When Scalapack is used, the orbital distribution
!     has to be block-cyclic).

!     In the loop below, no_l is called noccloc.
!     The global index of each of those noccloc eigenvectors can
!     be computed with the LocalToGlobalOrb routine. These helper
!     routines use implicitly the Blocksize determined at the
!     beginning of the program.

      jband = 0 
      do n = 0, Nodes-1 
        call GetNodeOrbs( nuotot, n, Nodes, noccloc )
        do ibandl = 1, noccloc
          call LocalToGlobalOrb( ibandl, n, Nodes, ibandg )
          if( .not.  isexcluded(ibandg) ) then
            jband = jband + 1
!           Identify the global band index of the included band
            index_included_band(jband) = ibandg

!           Identify in which node the coefficients of this band will be stored
            ! Use scalapack helper routine with the appropriate blocksize

            Nodeibandg = indxg2p(jband,blocksizeincbands,0,0,Nodes)
            node_included_band(jband) = Nodeibandg
          endif
        enddo
      enddo

!     Identify the local index of the included band inside the node
      bandspernode(:) = 0
      do jband = 1, numincbands
        nodejband = node_included_band(jband)
        bandspernode(nodejband) = bandspernode(nodejband) + 1
        band_index_in_node(jband) = bandspernode(nodejband)
!       Given a band and a local index, identiy the global index 
!       of the band that is stored there
        which_band_in_node(nodejband,band_index_in_node(jband)) = 
     .       index_included_band(jband)
      enddo

!     Set up the sequential index for the included band for wannierization
      icounter = 0
      sequential_index_included_bands(:) = 0
      do ibandg = 1, nuotot
        if( .not.  isexcluded(ibandg) ) then
          icounter = icounter + 1
          sequential_index_included_bands(ibandg) = icounter
        endif
      enddo 

!!     For debugging
!      do ibandg = 1, nuotot
!        write(6,'(a,3i5)')
!     .    ' Node, iband, sequential_index_included_bands = ',
!     .      Node, ibandg, sequential_index_included_bands(ibandg)
!      enddo
!
!      do iband_included = 1, numincbands
!        jband = index_included_band(iband_included)
!        write(6,'(a,6i5)')
!     .   ' Node, iband_included, index, node, local index, sequential=',
!     .     Node, iband_included, index_included_band(iband_included), 
!     .     node_included_band(iband_included),
!     .     band_index_in_node(iband_included),
!     .     sequential_index_included_bands(jband)
!      enddo
!
!      do n = 0, Nodes-1 
!        do ibandl = 1, nincbands_loc
!          write(6,'(a,4i5)')
!     .      ' Node, n, ibandl, orbital =',
!     .        Node, n, ibandl, which_band_in_node(n,ibandl)
!        enddo 
!      enddo 
!
!      call die()
!!     End debugging

      deallocate(bandspernode)

      end subroutine order_index
#endif
      end module m_orderbands
