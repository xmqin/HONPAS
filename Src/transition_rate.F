! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine transition_rate( ng, psi, ek, efermi, 
     .                           ekloc, efermiloc, temp, 
     .                           smooth, wmin, wmax,
     .                           Sr, Aux, Aux2, numh, listhptr, listh, 
     .                           indxuo, no, nuo, nuotot, xij, 
     .                           maxnh, nbands, kpoint, matrix,
     .                           intraband )
C *********************************************************************
C Finds the matrix element for the dipolar transition between to states
C Written by DSP. August 1999
C Restyled for f90 version by JDG. June 2004
C Modified by DSP. April 2010
C **************************** INPUT ********************************** 
C integer ng                  : first dimension of psi, Aux and Aux2
C real*8  psi(ng,nuotot,nuo)  : Wavefunctions in current k point
C real*8  ek(nuotot)          : Eigenvalues
C real*8  temp                : Electronic temperature
C real*8  efermi              : Fermi level 
C real*8  ekloc(nuotot)       : Eigenvalues (possibly modified by escissor)
C real*8  efermiloc           : Fermi level (possibly modified by escissor)
C real*8  wmin                : minimum transition energy required
C real*8  wmax                : maximum transition energy required
C real*8  smooth              : artificial width given to transitions
C real*8  Sr(maxnh)           : Position operator matrix elements (sparse)
C real*8  Aux(ng,nuotot,nuo)  : Auxiliary space
C real*8  Aux2(ng,nuotot,nuo)  : Auxiliary space
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to start of row in listh
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row 
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuo                 : Number of orbitals in the cell (locally)
C integer nuotot              : Number of orbitals in the cell (globally)
C integer maxnh               : Maximum dimension of listh
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)  
C integer nocc                : number of occupied states
C real*8  kpoint(3)           : Current kpoint
C real*8  dk(3)               : Vector joining the previous and current 
C                               kpoint
C character matrix*1          : 'R' or 'P' for position or momentum operator
C *************************** OUTPUT ****************************
C real*8  Aux(2,nuotot,nuo)   : matrix elements of the dipolar transition
C real*8  intraband(2,nuo)    : weight of intraband transitions (if a band
C                               crosses the Fermi energy
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************

      use precision
      use parallel,     only : BlockSize, Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb
      use m_fermid,     only : stepf
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none
 
C Passed variables
      integer nuo, nuotot, maxnh, no,
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  indxuo(no), nbands, ng

      real(dp)
     .  xij(3,maxnh), Sr(maxnh), 
     .  psi(ng,nuotot,nuo), kpoint(3), Aux(ng,nuotot,nuo), 
     .  Aux2(ng,nuotot,nuo), ek(nuotot), temp, efermi,
     .  intraband(ng,nuo), smooth, wmin, wmax, ekloc(nuotot),
     .  efermiloc, ediff

      character
     .  matrix*1

C Internal variables
      integer 
     .  ind, iuo, juo, j, ie, iie, iio, jje, je, jo, 
     .  BNodei, Bnodej, BTest, ig 
      real(dp)
     .  kxij, skxij, ckxij, pipj1, pipj2,
     .  f1,  f2, intrb(2), intrb2(2)

#ifdef MPI
      integer                     :: MPIerror
      real(dp), pointer           :: AuxLocal2(:,:)
#endif
      real(dp), pointer           :: AuxLocal(:,:)
      real(dp), pointer           :: psibandi(:,:)
      real(dp), pointer           :: psibandj(:,:)
      real(dp),              save :: ediffmin = 1.0d-3
      real(dp),              save :: tiny = 1.0d-9


C Start timer
      call timer('transrate',1)

C Check input matrix
      if(matrix.ne.'P'.and.matrix.ne.'R')
     $  call die('transrate: matrix only can take values R or P')

C Initialise matrix elements to zero
      do iuo = 1,nuo
        do juo = 1,nuotot
         do ig=1,ng
          Aux(ig,juo,iuo) = 0.0_dp
          Aux2(ig,juo,iuo) = 0.0_dp
         enddo 
        enddo 
        do ig=1,ng
           intraband(ig,iuo)=0.0_dp
        enddo 
      enddo 
          
C Compute matrix elements
      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          if(ng.gt.1) then 
           kxij = kpoint(1) * xij(1,ind) +
     .            kpoint(2) * xij(2,ind) +
     .            kpoint(3) * xij(3,ind) 
           ckxij = dcos(kxij)
           skxij = dsin(kxij)
           Aux2(1,juo,iuo) = Aux2(1,juo,iuo) + Sr(ind)*ckxij
           Aux2(2,juo,iuo) = Aux2(2,juo,iuo) + Sr(ind)*skxij
          else
           Aux2(1,juo,iuo) = Aux2(1,juo,iuo) + Sr(ind)
          endif
        enddo 
      enddo 

C Allocate workspace array
#ifdef MPI
      nullify( AuxLocal2 )
      call re_alloc( AuxLocal2, 1, ng, 1, nuotot,
     &               name    = 'AuxLocal2',
     &               routine = 'transition_rate' )
#endif
      nullify( AuxLocal )
      call re_alloc( AuxLocal, 1, ng, 1, nuotot,
     &               name    = 'AuxLocal',
     &               routine = 'transition_rate' )
      nullify( psibandi )
      call re_alloc( psibandi, 1, ng, 1, nuotot,
     &               name    = 'psibandi',
     &               routine = 'transition_rate' )
      nullify( psibandj )
      call re_alloc( psibandj, 1, ng, 1, nuotot,
     &               name    = 'psibandj',
     &               routine = 'transition_rate' )
 
      BNodei = 0
      iie = 0
      do ie = 1,nbands
        f1 = 2.0d0*stepf((ekloc(ie)-efermiloc)/temp)
        if (Node.eq.BNodei) then
          iie = iie + 1 
          do j = 1,nuotot 
            do ig=1,ng
              psibandi(ig,j) = psi(ig,j,iie)
            enddo 
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(psibandi(1,1),ng*nuotot,MPI_double_precision,
     .    BNodei,MPI_Comm_World,MPIerror)
#endif
            
C Compute the intra-band contribution to include a Drude-like term
C for metals
       if(matrix.eq.'P') then 
             intrb(:)=0.0_dp
             ediff=ek(ie)-efermi
             if(dabs(ediff).lt.2.0d0*smooth) then 
              do iuo = 1,nuo
                call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                do juo = 1,nuotot
                 if(ng.eq.2) then
                   pipj1 = psibandi(1,iio)*psibandi(1,juo) +
     .                     psibandi(2,iio)*psibandi(2,juo)
                   pipj2 = psibandi(1,iio)*psibandi(2,juo) -
     .                     psibandi(2,iio)*psibandi(1,juo)

C The factor 1/2 is necessary because the in subroutine phirphi_opt
C we have added a factor of 2, to be canceled by the energy denominator
C in Ry (rather than in Ha). For the calculation of the 
C Drude term we do not have the energy denominator 
                   intrb(1) = intrb(1)
     .                + 0.5_dp * (pipj1*Aux2(1,juo,iuo)
     .                - pipj2*Aux2(2,juo,iuo))
                                                                                
                   intrb(2) = intrb(2)
     .                + 0.5_dp * (pipj1*Aux2(2,juo,iuo)
     .                + pipj2*Aux2(1,juo,iuo))
                  else
                   pipj1 = psibandi(1,iio)*psibandi(1,juo) 
                                                                                
                   intrb(1) = intrb(1)
     .                + 0.5_dp* (pipj1*Aux2(1,juo,iuo))
                   intrb(2) = 0.0_dp
                  endif
                 enddo
               enddo
              endif
#ifdef MPI
        call MPI_Reduce(intrb,intrb2,2,
     .    MPI_double_precision,MPI_sum,BNodei,MPI_Comm_World,MPIerror)
        if (Node.eq.BNodei) then
          intraband(1:ng,iie) = intrb2(1:ng)
        endif
#else
        if (Node.eq.BNodei) then
          intraband(1:ng,iie) = intrb(1:ng)
        endif
#endif
        endif

        AuxLocal(1:ng,1:nuotot) = 0.0d0

        BNodej = 0
        jje = 0
        do je = 1,nbands
          if (Node.eq.BNodej) then
            jje = jje + 1 
          endif
          if (dabs(ek(ie)-ek(je)).gt.ediffmin) then   
            f2 = 2.0d0*stepf((ekloc(je)-efermiloc)/temp)
                  
            if (f1*(2.0d0-f2).gt.tiny) then 

             ediff = ekloc(je) - ekloc(ie)
             if ((ediff.ge.max(wmin-2.0d0*smooth,0.0d0)).and.
     .               (ediff.le.wmax+2.0d0*smooth)) then

              if (Node.eq.BNodej) then
                do j = 1,nuotot 
                  do ig=1,ng
                    psibandj(ig,j) = psi(ig,j,jje)
                  enddo 
                enddo
              endif
#ifdef MPI
              call MPI_Bcast(psibandj(1,1),ng*nuotot,
     $        MPI_double_precision,BNodej,MPI_Comm_World,MPIerror)
#endif

              do iuo = 1,nuo
                call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                do juo = 1,nuotot
                 if(ng.eq.2) then  
                   pipj1 = psibandi(1,iio)*psibandj(1,juo) +
     .                     psibandi(2,iio)*psibandj(2,juo)
                   pipj2 = psibandi(1,iio)*psibandj(2,juo) -
     .                     psibandi(2,iio)*psibandj(1,juo) 

                   AuxLocal(1,je) = AuxLocal(1,je)
     .                + pipj1*Aux2(1,juo,iuo)
     .                - pipj2*Aux2(2,juo,iuo)

                   AuxLocal(2,je) = AuxLocal(2,je)
     .                + pipj1*Aux2(2,juo,iuo)
     .                + pipj2*Aux2(1,juo,iuo)
                  else 
                   pipj1 = psibandi(1,iio)*psibandj(1,juo) 
                      
                   AuxLocal(1,je) = AuxLocal(1,je)
     .                + pipj1*Aux2(1,juo,iuo)
                  endif

                enddo 
              enddo  
              if (matrix.eq.'P') then 
                do ig=1,ng
                 AuxLocal(ig,je) = AuxLocal(ig,je)/(ek(je)-ek(ie))
                enddo 
              endif 
             endif
            endif 
          endif 
          BTest = je/BlockSize
          if (BTest*BlockSize.eq.je) then
            BNodej = BNodej + 1
            if (BNodej .gt. Nodes-1) BNodej = 0
          endif
        enddo 
#ifdef MPI
        call MPI_Reduce(AuxLocal(1,1),AuxLocal2(1,1),ng*nuotot,
     .    MPI_double_precision,MPI_sum,BNodei,MPI_Comm_World,MPIerror)
        if (Node.eq.BNodei) then
          Aux(1:ng,1:nuotot,iie) = AuxLocal2(1:ng,1:nuotot)
        endif
#else
        if (Node.eq.BNodei) then
          Aux(1:ng,1:nuotot,iie) = AuxLocal(1:ng,1:nuotot)
        endif
#endif
        BTest = ie/BlockSize
        if (BTest*BlockSize.eq.ie) then
          BNodei = BNodei + 1
          if (BNodei .gt. Nodes-1) BNodei = 0
        endif
      enddo 

C Free workspace array
      call de_alloc( psibandj, 'psibandj', 'transition_rate')
      call de_alloc( psibandi, 'psibandi', 'transition_rate')
      call de_alloc( AuxLocal, 'AuxLocal', 'transition_rate')
#ifdef MPI
      call de_alloc( AuxLocal2, 'AuxLocal2', 'transition_rate')
#endif
         
C Stop timer
      call timer('transrate',2)

      end
