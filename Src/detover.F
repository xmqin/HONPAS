! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine detover(psiprev, psi, S, Sr, 
     .               numh, listhptr, listh, indxuo,
     .               no, nuo, xij, maxnh, nuotot, nocc,
     .               kpoint, dk, detr, deti )
C *********************************************************************
C Finds the determinant of the overlap matrix
C between the periodic Bloch functions corresponding to neighboring  
C k points
C Written by DSP. March 1999
C Modified for parallel execution by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)   : Wavefunctions in current k point
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  Sr(maxnh)           : Position operator matrix elements (sparse)
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row 
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuo                 : Number of basis orbitals in the unit cell
C integer maxnh               : Maximum number of orbitals interacting
C integer nuotot              : Third dimension of xij
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)  
C integer nocc                : number of occupied states
C real*8  kpoint(3)           : Current kpoint
C real*8  dk(3)               : Vector joining the previous and current 
C                               kpoint
C *************************** INPUT/OUTPUT ****************************
C real*8  psiprev(2,nuo,nuotot) : Wavefunctions in previous k point
C real*8  detr                  : Real part of the determinant
C real*8  deti                  : Imaginary part of the determinant
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
      use alloc,        only : re_alloc, de_alloc
      use m_linpack,    only : zgedi, zgefa
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
 
      integer 
     .  nuo, maxnh, nuotot, no, nocc,
     .  listh(maxnh), listhptr(nuo), numh(nuo), indxuo(no), 
     .  info, job
      parameter(job=10)

      real(dp)
     .  psiprev(2,nuo,nuotot), dk(3),detr, deti,
     .  xij(3,maxnh), S(maxnh), Sr(maxnh),
     .  psi(2,nuotot,nuo), kpoint(3)

      complex(dp) ::   det(2), pipj, determ

C**** Internal variables ***********************************************

      integer :: iuo, juo, j, ie, je, jo, ind
      integer, dimension(:), pointer ::  Auxint
#ifdef MPI
      integer :: MPIerror, jeg, n, noccloc, noccmax
      real(dp), dimension(:,:,:), pointer :: psitmp
#endif
      real(dp) ::  kxij, skxij, ckxij,pipj1, pipj2

      complex(dp), dimension(:,:), pointer :: Aux, Aux2

C Allocate local memory
      nullify( Aux )
      call re_alloc( Aux, 1, nocc, 1, nocc, name='Aux',
     &               routine='detover' )
      nullify( Aux2 )
      call re_alloc( Aux2, 1, nuotot, 1, nuo, name='Aux2',
     &               routine='detover' )
      nullify( Auxint )
      call re_alloc( Auxint, 1, nocc, name='Auxint', routine='detover' )

      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo)+j
          jo = listh(ind)
          juo = indxuo(jo)
          kxij = (kpoint(1)-0.5d0*dk(1)) * xij(1,ind) +
     .           (kpoint(2)-0.5d0*dk(2)) * xij(2,ind) +
     .           (kpoint(3)-0.5d0*dk(3)) * xij(3,ind) 
          ckxij = dcos(kxij)
          skxij = dsin(kxij)
          Aux2(juo,iuo)=Aux2(juo,iuo)+
     .      cmplx(  S(ind)*ckxij + Sr(ind)*skxij, 
     .      S(ind)*skxij - Sr(ind)*ckxij , dp )  
 
        enddo 
      enddo 
 
#ifdef MPI
      call GetNodeOrbs(nuotot,0,Nodes,noccmax)
      nullify( psitmp )
      call re_alloc( psitmp, 1, 2, 1, nuotot, 1, noccmax, name='psitmp',
     &               routine='detover' )

C Ultimately this needs modifying so that Aux is distributed
      do n = 0,Nodes-1

C Broadcast copy of psi on node n to all other nodes
        call GetNodeOrbs(nocc,n,Nodes,noccloc)
        call GetNodeOrbs(nuotot,n,Nodes,noccmax)
        if (Node .eq. n) then
          do iuo = 1,noccmax
            do juo = 1,nuotot
              psitmp(1,juo,iuo) = psi(1,juo,iuo)
              psitmp(2,juo,iuo) = psi(2,juo,iuo)
            enddo
          enddo
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccmax,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

        do ie = 1,nocc
          do je = 1, noccloc
            call LocalToGlobalOrb(je,n,Nodes,jeg)

            do iuo=1,nuo
              do juo=1,nuotot

                pipj1 = psiprev(1,iuo,ie) * psitmp(1,juo,je) +
     .                  psiprev(2,iuo,ie) * psitmp(2,juo,je)
                pipj2 = psiprev(1,iuo,ie) * psitmp(2,juo,je) -
     .                  psiprev(2,iuo,ie) * psitmp(1,juo,je)
                pipj=cmplx(pipj1,pipj2, dp)
  
                Aux(jeg,ie) = Aux(jeg,ie) + pipj*Aux2(juo,iuo)

              enddo
            enddo

          enddo
        enddo

      enddo

      call de_alloc( psitmp, name='psitmp' , routine='detover')
#else
      do ie = 1,nocc
        do je = 1, nocc
              
          do iuo=1,nuo
            do juo=1,nuotot

              pipj1 = psiprev(1,iuo,ie) * psi(1,juo,je) +
     .                psiprev(2,iuo,ie) * psi(2,juo,je)
              pipj2 = psiprev(1,iuo,ie) * psi(2,juo,je) -
     .                psiprev(2,iuo,ie) * psi(1,juo,je) 
              pipj=cmplx(pipj1,pipj2, dp)

              Aux(je,ie) = Aux(je,ie) + pipj*Aux2(juo,iuo)

            enddo 
          enddo  

        enddo 
      enddo 
#endif

C Resize Aux2 for re-use
      call re_alloc( Aux2, 1, nocc, 1, nocc, name='Aux2',
     &               routine='detover', copy=.false. )

#ifdef MPI
      call MPI_AllReduce(Aux(1,1),Aux2(1,1),nocc*nocc,
     .  MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror)
      Aux = Aux2
#endif

      call zgefa(Aux,nocc,nocc,Auxint,info)
      call zgedi(Aux,nocc,nocc,Auxint,det,Aux2,job)

      determ=det(1)
      deti=aimag(determ)
      detr=real(determ)   

C Deallocate local memory
      call de_alloc( Auxint, name='Auxint' )
      call de_alloc( Aux2, name='Aux2' )
      call de_alloc( Aux, name='Aux' )

      end
