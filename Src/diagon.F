! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_diagon

      private
      public :: diagon

      contains 

      subroutine diagon(no, nspin, maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, 
     .                  listdptr, listd, H, S, qtot, fixspin, 
     .                  qs, temp, e1, e2, xij, indxuo, gamma, nk, 
     .                  kpoint, wk, eo, qo, Dnew, Enew, ef, efs, 
     .                  Entropy, nuotot, occtol, iscf, neigwanted,
     .                  dealloc_psi)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization).
C Written by J.Soler and P.Ordejon, August 1997.
C Now includes the option to parallelise run over K points.
C Modified by J.Gale, November 1999.
C Added new diagk_file routine, A. Garcia, June 2008.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Spin polarization, spinor (1 or 2)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnh)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,spin%H)     : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  qtot                : Number of electrons in unit cell
C logical fixspin             : Fix the spin of the system?
C real*8  qs(nspin)          : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
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
C logical Gamma               : Whether only the Gamma point is sampled.
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer nuotot              : total number of orbitals in unit cell 
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer iscf                : SCF cycle number
C integer neigwanted          : Number of eigenvalues wanted
C logical, dealloc_psi        : [OPTIONAL] whether or not to deallocate
C                               the densematrices used [[densematrix]]
C *************************** OUTPUT **********************************
C real*8 eo(maxo,nspin,nk)  : Eigenvalues 
C real*8 qo(maxo,nspin,nk)  : Occupations of eigenstates
C real*8 Dnew(maxnd,spin%DM)  : Output Density Matrix
C real*8 Enew(maxnd,spin%EDM)  : Output Energy-Density Matrix
C real*8 ef                    : Fermi energy
C real*8 efs(nspin)         : Fermi energy for each spin
C                                   (for fixed spin calculations)
C real*8 Entropy                 : Electronic entropy
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** Parallel ********************************
C When running in parallel some of the dimensions are now the 
C maximum per node and the corresponding number passed in as
C an argument is the number of locally stored values. The
C variables for which this is the case are:
C
C maxuo/no
C
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use densematrix,  only : allocDenseMatrix, resetDenseMatrix
      use densematrix,  only : Haux, Saux, psi
      use sys,          only : die

      use fdf
      use alloc
      use m_memory
      use m_spin,       only: Spiral, qSpiral
      use m_spin,       only: spin
#ifdef MPI
      use m_diag_option, only: ParallelOverK
      use mpi_siesta
#endif
      use siesta_options, only: diag_wfs_cache

      implicit none

      real(dp), intent(in) :: H(:,:)
      logical, intent(in) :: gamma

      integer
     .  iscf, maxnd, maxnh, maxuo, maxo, nk, no, nuotot,
     .  neigwanted, nspin

      integer 
     .  indxuo(no), listh(maxnh), numh(*), listd(maxnd), numd(*),
     .  listhptr(*), listdptr(*)

      real(dp)
     .  Dnew(maxnd,size(H,2)), e1, e2, ef, Enew(:,:), 
     .  Entropy, eo(maxo,nspin,nk),  
     .  kpoint(3,nk), 
     .  qo(maxo,nspin,nk), qtot, S(maxnh), temp, wk(nk), occtol,
     .  xij(3,maxnh), qs(:), efs(nspin)
     
      logical
     .  fixspin, getD, getPSI

      logical, intent(in), optional :: dealloc_psi
      
      external
     .  diag2g, diag2k, diag2kspiral
#ifdef CDF
      external diagk_file
#endif

#ifdef MPI
      external           :: diagkp
#endif
      integer            :: io, iuo, naux, nhs, npsi, nuo
      real(dp),  pointer :: aux(:)
C  ....................

C Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

C Start time counter ................................................
      call timer( 'diagon', 1 )

C     Check internal dimensions ..........................................
      
      if ( spin%none .or. spin%Col ) then
       if (gamma) then
        nhs  = nuotot * nuo
        npsi = nuotot * maxuo * nspin
        naux = nuotot * 5
       else
        nhs  = 2 * nuotot * nuo
        npsi = 2 * nuotot * nuo
        naux = 2 * nuotot * 5            !!!! AG: Put factor 2 !!
#ifdef MPI
        if ( ParallelOverK ) then
          nhs  = 2 * nuotot * nuotot
          npsi = 2 * nuotot * nuotot
        endif
#endif
       endif
      elseif ( spin%NCol .or. spin%SO ) then
        nhs  = 2 * (2*nuotot) * (2*nuo)
        npsi = 2 * (2*nuotot) * (2*nuo)
        naux = 2 * 2 * nuotot ! COMPLEX, also for SPIRAL
#ifdef MPI
        if ( ParallelOverK ) then
          nhs  = 2 * 2 * nuotot * 2 * nuotot
          npsi = 2 * 2 * nuotot * 2 * nuotot
          naux = 2 * nuotot     ! REAL
          if ( Spiral ) call die('diagon: ParallelOverK and Spiral &
     &are not compatible.')
        end if
#endif
      else
         call die('diagon: ERROR: incorrect value of nspin')
      endif

      call allocDenseMatrix(nhs, nhs, npsi)

C Allocate local arrays
      nullify(aux)
      call re_alloc( aux,  1, naux, 'aux',  'diagon' )


C Call apropriate routine
      getD = .true.
      getPSI = .false.

      if ( spin%none .or. spin%Col ) then
       if (gamma) then
        call diagg( spin%H, nuo, maxuo, maxnh, maxnd, maxo,
     .              numh, listhptr, listh, numd, listdptr, listd, 
     .              H, S, getD, getPSI, fixspin, qtot, qs, temp, e1, e2,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, nuotot, occtol, iscf,
     .              neigwanted )
       else
#ifdef MPI
        if (ParallelOverK) then
          call diagkp( spin%H, nuo, no, nspin, maxnh,
     .              maxnd, maxo, numh,
     .              listhptr, listh, H, S, getD, fixspin, qtot, qs,
     .              temp, e1, e2, xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux, 
     .              nuotot, occtol, iscf, neigwanted )
        else
#endif
#ifndef CDF
           if ( diag_wfs_cache == 1 ) then
              if (Node == 0) then
                 write(6,"(a)")
     $           "** Cannot use new diagk without netCDF support"
                 write(6,"(a)")
     $           "** Falling back to standard diagk routine"
              endif
              diag_wfs_cache = 0
           endif
#endif
           if ( diag_wfs_cache == 1 ) then
#ifdef CDF
              ! Use new routine with file storage
              call diagk_file(spin%H, nuo, no, nspin,
     $              maxnh, maxnd, 
     .              maxo, numh, listhptr, listh, numd, listdptr,
     .              listd, H, S, getD, getPSI, fixspin, qtot, qs, temp, 
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux,
     .              nuotot, occtol, iscf, neigwanted )
#endif
           else
              call diagk( spin%H, nuo, no, nspin,
     .              maxnh, maxnd, 
     .              maxo, numh, listhptr, listh, numd, listdptr,
     .              listd, H, S, getD, getPSI, fixspin, qtot, qs, temp, 
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux, 
     .              nuotot, occtol, iscf, neigwanted )

           endif
#ifdef MPI
        endif
#endif
       endif
      
      elseif ( spin%NCol ) then
       if (gamma) then
        call diag2g( nuo, no, maxnh, maxnd, maxo, numh,
     .               listhptr, listh, numd, listdptr, listd,
     .               H, S, getD, getPSI, qtot, temp, e1, e2, eo, qo, 
     .               Dnew, Enew, ef, Entropy,
     .               Haux, Saux, psi, aux,
     .               nuotot, occtol, iscf, neigwanted )
      elseif ( Spiral ) then
        call diag2kspiral( nuo, no, maxnh, maxnd, maxo,
     .               numh, listhptr, listh, numd, listdptr,
     .               listd, H, S, getD, qtot, temp, e1, e2,
     .               xij, indxuo, nk, kpoint, wk, eo, qo,
     .               Dnew, Enew, ef, Entropy, qspiral, Haux,
     .               Saux, psi, Haux, Saux, aux, nuotot, occtol,
     .               iscf, neigwanted)
#ifdef MPI
       else if ( ParallelOverK ) then
         call diag2kp( spin, maxuo, nuotot, no,
     .       maxnh, numh, listhptr, listh, H, S, getD,
     .       qtot, temp, e1, e2, xij,
     .       nk, kpoint, wk,
     .       eo, qo, Dnew, Enew, Ef, Entropy,
     .       Haux, Saux, psi, aux, occtol, iscf, neigwanted)
#endif
       else
         call diag2k( nuo, no, maxnh, maxnd, maxo, numh,
     .       listhptr, listh, numd, listdptr, listd,
     .       H, S, getD, getPsi, qtot, temp, e1, e2, xij,
     .       indxuo, nk, kpoint, wk, eo, qo, Dnew,
     .       Enew, ef, Entropy,
     .       Haux, Saux, psi, Haux, Saux, aux,
     .       nuotot, occtol, iscf, neigwanted )
       endif
      elseif ( spin%SO ) then
        if (gamma) then
          call diag3g( nuo, no, maxnh, maxnd, maxo, numh,
     .                 listhptr, listh, numd, listdptr, listd,
     .                 H, S, getD, getPSI, qtot, temp, e1, e2, eo, qo,
     .                 Dnew, Enew, ef, Entropy, 
     .                 Haux, Saux, psi, aux,
     .                 nuotot, occtol, iscf, neigwanted )
#ifdef MPI
        else if ( ParallelOverK ) then
          call diag3kp( spin, maxuo, nuotot, no,
     .        maxnh, numh, listhptr, listh, H, S, getD,
     .        qtot, temp, e1, e2, xij,
     .        nk, kpoint, wk,
     .        eo, qo, Dnew, Enew, Ef, Entropy,
     .        Haux, Saux, psi, aux, occtol, iscf, neigwanted)
#endif
        else
          call diag3k( nuo, no, maxnh, maxnd, maxo, numh,
     .        listhptr, listh, numd, listdptr, listd,
     .        H, S, getD, getPSI, qtot, temp, e1, e2, xij,
     .        indxuo, nk, kpoint, wk, eo, qo, Dnew,
     .        Enew, ef, Entropy,
     .        Haux, Saux, psi, Haux, Saux, aux,
     .        nuotot, occtol, iscf, neigwanted )
        endif
      endif      


C Free local arrays
      call de_alloc( aux, name='aux', routine='diagon' )

      call resetDenseMatrix(dealloc_psi=dealloc_psi)

C Stop time counter
      call timer( 'diagon', 2 )

      end subroutine diagon
   
      end module m_diagon
