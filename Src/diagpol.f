! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diagpol( ispin, nspin, nuo, no, nuotot,
     .                    maxnh, numh, listhptr, listh, H, S,
     .                    xij, indxuo, kpoint, eo, psi, ng, 
     .                    Haux, Saux )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, 
C for given Hamiltonian and Overlap matrices (including
C spin polarization), and for a given k_point
C Written by DSP, March 1999. From the routine diagk of J.Soler.
C Modified for parallel execution by J.D.Gale, March 2000.
C **************************** INPUT **********************************
C integer ispin               : Spin component which will be calculated
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuotot              : First dimension of eo, qo, last of xij
C integer maxnh               : Maximum number of orbitals interacting  
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xij(3,*)            : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  kpoint(3)           : k point vectors
C real*8  psi                 : Workspace array
C integer ng                  : first dimension of Haux, and Saux
C real*8  Haux(ng,nuotot,nuo) : Workspace for dense H
C real*8  Saux(ng,nuotot,nuo) : Workspace for dense S
C *************************** OUTPUT **********************************
C real*8 eo(nuotot)           : Eigenvalues
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C eo H.
C *********************************************************************

      use precision
      use parallel,      only : BlockSize
      use sys

      implicit          none

      integer
     .  maxnh, nuotot, no, nspin, nuo, indxuo(no), listh(maxnh), 
     .  listhptr(nuo), numh(nuo), ng

      real(dp)
     .  eo(nuotot), H(maxnh,nspin), kpoint(3), S(maxnh), 
     .  xij(3,*), psi(ng,nuotot,nuo), Haux(ng,nuotot,nuo),
     .  Saux(ng,nuotot,nuo)
      external          rdiag, cdiag

C  Internal variables .............................................
      integer
     .  ierror, ind, ispin, iuo, j, jo, juo
      real(dp)
     .  ckxij, kxij, skxij

C Solve eigenvalue problem .........................................
      Saux = 0.0d0
      Haux = 0.0d0
      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          if(ng.eq.2) then 
           kxij = kpoint(1) * xij(1,ind) +
     .            kpoint(2) * xij(2,ind) +
     .            kpoint(3) * xij(3,ind)
           ckxij = cos(kxij)
           skxij = sin(kxij)
           Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
           Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
           Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
           Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
          else 
           Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)
           Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)
          endif
        enddo
      enddo
      if(ng.eq.2) then 
       call cdiag( Haux, Saux, nuotot, nuo, nuotot, eo, psi,
     .            nuotot, 1, ierror, BlockSize)
      else
       call rdiag( Haux, Saux, nuotot, nuo, nuotot, eo, psi,
     .            nuotot, 1, ierror, BlockSize)
      endif
C Check error flag and take appropriate action
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
        Saux = 0.0d0
        Haux = 0.0d0
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
          if(ng.eq.2) then
           kxij = kpoint(1) * xij(1,ind) +
     .            kpoint(2) * xij(2,ind) +
     .            kpoint(3) * xij(3,ind)
           ckxij = cos(kxij)
           skxij = sin(kxij)
           Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
           Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
           Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
           Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
          else
           Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)
           Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)
          endif
        enddo
      enddo
      if(ng.eq.2) then
       call cdiag( Haux, Saux, nuotot, nuo, nuotot, eo, psi,
     .            nuotot, 1, ierror, BlockSize)
      else
       call rdiag( Haux, Saux, nuotot, nuo, nuotot, eo, psi,
     .            nuotot, 1, ierror, BlockSize)
      endif

      endif

      end
