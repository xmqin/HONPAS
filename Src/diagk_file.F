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
#ifndef CDF
      subroutine dummy_diagk_file()
      end subroutine dummy_diagk_file
#else
      subroutine diagk_file( nspin, nuo, no, maxspn,
     $                  maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, listdptr, 
     .                  listd, H, S, getD, getPSI, fixspin, qtot, qs, 
     .                  temp, e1, e2, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     .                  Haux, Saux, psi, Dk, Ek, aux, nuotot, occtol,
     .                  iscf, n_eigenvectors )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version.
C Writen by J.Soler, August 1998.
C Modified by A. Garcia, June 2008, to trade file space for speed (!)
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxo                : First dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd and DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical getPSI              : Find and print wave functions?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer n_eigenvectors      : Number of eigenvectors to be saved
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C *************************** OUTPUT **********************************
C real*8 eo(maxo,maxspn,nk)   : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuotot,nuo) : Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,nuotot,nuo) : Auxiliary space for the overlap matrix
C real*8 psi(2,nuotot,nuo)  : Auxiliary space for the eigenvectors
C real*8 aux(2,nuotot*5)    : Extra auxiliary space
C real*8 Dk(2,nuotot,nuo)   : Aux. space that may be the same as Haux 
C real*8 Ek(2,nuotot,nuo)   : Aux. space that may be the same as Saux 
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetric and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave,     only : writew
      use m_fermid,      only : fermid, fermispin, stepf
#ifdef MPI
      use mpi_siesta
#endif

      use iowfs_netcdf,  only : setup_wfs_netcdf_file, write_wfs_netcdf
      use iowfs_netcdf,  only : get_wfs_netcdf, get_wfs_block_netcdf
      use iowfs_netcdf,  only : open_wfs_netcdf_file,
     $                          close_wfs_netcdf_file

      implicit          none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer           maxnd, maxnh, maxspn, maxuo, maxo, nk, no,
     .                  nspin, nuo, nuotot, iscf, n_eigenvectors
      integer           indxuo(no), listh(maxnh), numh(nuo),
     .                  listd(maxnd), numd(nuo), listhptr(nuo),
     .                  listdptr(nuo)
      real(dp)          Dnew(maxnd,nspin), 
     .                  e1, e2, ef, efs(nspin), Enew(maxnd,nspin),
     .                  Entropy, eo(maxo,maxspn,nk), H(maxnh,nspin),
     .                  kpoint(3,nk), qo(maxo,maxspn,nk), qtot, 
     .                  S(maxnh), temp, wk(nk), occtol,
     .                  xij(3,maxnh), qs(nspin)
      real(dp)          Dk(2,nuotot,nuo), Ek(2,nuotot,nuo),
     .                  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo),
     .                  psi(2,nuotot,nuo), aux(2,nuotot*5)
      logical           getD, getPSI, fixspin
      external          cdiag

C  Internal variables .............................................
      integer
     .  BNode, BTest, ie, ierror, iie, ik, ind, io, iio,
     .  ispin, iuo, j, jo, juo, nd, neigneeded, iu

      real(dp)
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
      real(dp) :: qpipj1, qpipj2, epipj1, epipj2

      integer :: ifirst, ilast, nvecs_in_block, nvecs_read, i
      integer :: max_block_size
      real(dp), dimension(:,:,:), allocatable, target   :: psi_block
      real(dp), dimension(:,:), pointer                 :: psi_aux
C  ....................

C Find eigenvalues and eigenvectors

      ! Setup netCDF file
      !
      if (Node == 0) then
         call setup_wfs_netcdf_file(nuotot, nk, nspin, n_eigenvectors)
      endif
      do ik = 1,nk
        do ispin = 1,nspin
          call timer( 'c-eigval', 1 )

          call build_dense_HS(ik,ispin)
          ! Compute all eigenvectors for now
          call cdiag(Haux,Saux,nuotot,nuo,nuotot,eo(1,ispin,ik),psi,
     .      nuotot,iscf,ierror)

          if (ierror.gt.0) then
            call die('Terminating due to failed diagonalisation')
          elseif (ierror.lt.0) then
             ! Repeat diagonalisation with increased memory to handle clustering

             call build_dense_HS(ik,ispin)

             call cdiag(Haux,Saux,nuotot,nuo,nuotot,eo(1,ispin,ik),psi,
     .            nuotot,iscf,ierror)
          endif
          call timer( 'c-eigval', 2 )
!
!         SAVE eigenvectors to file  (only up to those specified)
!
          call timer( 'c-save', 1 )
          call write_wfs_netcdf(nuotot,nuo,ik,ispin, psi,
     $                 n_eigenvectors,eo(1:n_eigenvectors,ispin,ik))
          call timer( 'c-save', 2 )

        enddo
      enddo

      call close_wfs_netcdf_file()

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, maxo, nuotot, eo,
     .               temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef, Entropy )
      endif

      if (Node .eq. 0) then
         call io_assign(iu)
         open(iu,file="OCCS",form="formatted",status="unknown",
     $           action="write",position="rewind")
      endif

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.e-6_dp )
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuotot
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf((eo(io,ispin,ik)-e2)/t) -
     .               stepf((eo(io,ispin,ik)-e1)/t)) * 2.0_dp/nspin
            enddo
          enddo
        enddo
      endif

C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.0_dp
      Enew(1:nd,1:nspin) = 0.0_dp

      call open_wfs_netcdf_file()

      max_block_size = n_eigenvectors / Nodes
      allocate(psi_block(2,nuotot,max_block_size))

C Loop over k points
      do ik = 1,nk
        do ispin = 1,nspin

            neigneeded = 0
            ie = nuotot
            do while (ie.gt.0.and.neigneeded.eq.0)
              qe = qo(ie,ispin,ik)
              if (abs(qe).gt.occtol) neigneeded = ie
              ie = ie - 1
            enddo
            if (Node .eq. 0) then
             write(iu,"(i6,i3,f12.8,i8)")
     $              ik, ispin, wk(ik), neigneeded
            endif
            if (neigneeded > n_eigenvectors) then
               if (Node == 0) then
                  write(6,"(a)") "ERROR: More eigenvectors are needed!"
                  write(6,"(a,2i6)") "Saved, needed: ", n_eigenvectors,
     $                 neigneeded
               endif
               call die()
            endif

            call timer( 'c-buildD', 1 )

C Global operation to form new density matrix

!           Get the eigenvectors in blocks
            nvecs_read = 0
            do while (nvecs_read < neigneeded)

               ifirst = nvecs_read + 1
               ilast = ifirst + max_block_size - 1
               ilast = min(ilast,neigneeded)
               nvecs_in_block = ilast - ifirst + 1

               ! Get eigenvectors (will broadcast to all nodes)
               call get_wfs_block_netcdf(ifirst,nvecs_in_block,
     $                                   ik,ispin,nuotot,psi_block)
               nvecs_read = nvecs_read + nvecs_in_block

            ! Loop over slots in density-matrix in this node

            do iuo = 1,nuo
               call LocalToGlobalOrb(iuo,Node,Nodes,iio)
               do j = 1,numd(iuo)
                ind = listdptr(iuo) + j
                jo = listd(ind)
                juo = indxuo(jo)

                kxij = kpoint(1,ik) * xij(1,ind) +
     .                 kpoint(2,ik) * xij(2,ind) +
     .                 kpoint(3,ik) * xij(3,ind)
                ckxij = cos(kxij)
                skxij = sin(kxij)

                qpipj1 = 0.0_dp
                qpipj2 = 0.0_dp
                epipj1 = 0.0_dp
                epipj2 = 0.0_dp
                do i = 1, nvecs_in_block
                   ie = ifirst + (i-1)
                   psi_aux => psi_block(:,:,i)
                   qe = qo(ie,ispin,ik)
                   ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
                   pipj1 = psi_aux(1,iio) * psi_aux(1,juo) +
     .                  psi_aux(2,iio) * psi_aux(2,juo)
                   pipj2 = psi_aux(1,iio) * psi_aux(2,juo) -
     .                  psi_aux(2,iio) * psi_aux(1,juo)
                   qpipj1 = qpipj1 + qe * pipj1
                   qpipj2 = qpipj2 + qe * pipj2
                   epipj1 = epipj1 + ee * pipj1
                   epipj2 = epipj2 + ee * pipj2
                enddo

                Dnew(ind,ispin)=Dnew(ind,ispin)+ qpipj1*ckxij -
     .                                           qpipj2*skxij
                Enew(ind,ispin)=Enew(ind,ispin)+ epipj1*ckxij -
     .                                           epipj2*skxij
              enddo
            enddo

         enddo                  ! while


         call timer( 'c-buildD', 2 )

        enddo   ! ispin
      enddo     ! ik

      deallocate(psi_block)
      if (Node .eq. 0) then
         call io_close(iu)
      endif
      call close_wfs_netcdf_file()

      CONTAINS

      subroutine build_dense_HS(ik,ispin)
      integer, intent(in) :: ik, ispin

      ! All other variables taken from parent scope

          call timer( 'c-buildHS', 1 )
          Saux = 0.0_dp
          Haux = 0.0_dp
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
            enddo
          enddo
          call timer( 'c-buildHS', 2 )
          end subroutine build_dense_HS

      end subroutine diagk_file
#endif
