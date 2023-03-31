! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
c
c     Put dummy subroutine here so that compilers do not complain
c     about empty file in serial compilation.
c
      subroutine dummy_sub()
      end
c
#ifdef MPI
      subroutine diagkp( nspin, nuo, no, maxspn, maxnh, maxnd, 
     &                   maxo, numh, listhptr, listh, H, S, getD,
     &                   fixspin, qtot, qs, temp, e1, e2, xij, indxuo,
     &                   nk, kpoint, wk, eo, qo, Dnew, Enew, ef, efs,
     &                   Entropy, Haux, Saux, psi, Dk, Ek, aux, nuotot,
     &                   occtol, iscf, neigwanted)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version. 
C Created from diagk, written by J.Soler.
C Uses parallelisation over K points instead of parallelisation 
C within them.
C Parallel modifications by J.Gale, November 1999.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
C integer maxo                : First dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix locally
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix locally
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
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
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer neigwanted          : maximum number of eigenvalues wanted
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
C real*8 aux(2,nuotot)      : Extra auxiliary space
C real*8 Dk(2,nuotot,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,nuotot,nuo)   : Aux. space that may be the same as Saux
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetric and so the order
C of referencing has been changed in several places to reflect this.
C Note : It is assumed in a couple of places that the sparsity of
C H/S and Dscf are the same (which is the case in the code at present)
C and therefore the only one set of pointer arrays are globalised.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
      use mpi_siesta
      use m_fermid,     only : fermid, fermispin, stepf
      use alloc,        only : re_alloc, de_alloc

      implicit          none

      integer 
     &  MPIerror

      integer           maxnd, maxnh, maxspn, maxo, nk, no,
     &                  nspin, nuo, nuotot, iscf, neigwanted
      integer           indxuo(no), listh(maxnh), numh(nuo),
     &                  listhptr(nuo)
      real(dp)          Dnew(maxnd,nspin), 
     &                  e1, e2, ef, efs(nspin), Enew(maxnd,nspin),
     &                  Entropy, eo(maxo,maxspn,nk), H(maxnh,nspin),
     &                  kpoint(3,nk), qo(maxo,maxspn,nk), qtot, 
     &                  S(maxnh), temp, wk(nk), occtol,
     &                  xij(3,maxnh), qs(nspin)
      real(dp)          Dk(2,nuotot,nuotot), Ek(2,nuotot,nuotot),
     &                  Haux(2,nuotot,nuotot), Saux(2,nuotot,nuotot),
     &                  psi(2,nuotot,nuotot), aux(2,nuotot)
      logical           getD, fixspin
      external          cdiag, memory

C  Internal variables .............................................
      integer           :: BNode, ie, ierror, ik, io, iio, ind, is,
     &                     ispin, iuo, j, jo, juo, maxnhg, nuog,
     &                     maxnumh, neigneeded
      integer,  pointer :: numhg(:) => null(),
     $                     listhptrg(:) => null(),
     $                     listhg(:) => null()
      real(dp)          :: ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
      real(dp), pointer :: Snew(:) => null(),
     $                     Dloc(:) => null(),
     $                     Hnew(:,:) => null(),
     $                     Dnewloc(:,:) => null(),
     $                     Enewloc(:,:) => null(),
     $                     xijloc(:,:) => null()

#ifdef DEBUG
      call write_debug( '    PRE diagkp' )
#endif

C Globalise list arrays - assumes listh and listd are the same

      nullify(numhg, listhptrg)
C Allocate local memory for global list arrays
      call re_alloc( numhg, 1, nuotot, name='numhg',
     &               routine= 'diagkp' )
      call re_alloc( listhptrg, 1, nuotot, name='listhptrg',
     &               routine= 'diagkp' )

C Find maximum value in numh and create local storage
      maxnumh = 0
      do io = 1,nuo
        maxnumh = max(maxnumh,numh(io))
      enddo
      nullify(Dloc)
      call re_alloc( Dloc, 1, maxnumh, name='Dloc',
     &               routine= 'diagkp' )

C Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     &    MPI_Comm_World,MPIerror)
      enddo

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Globalse listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      nullify(listhg)
      call re_alloc( listhg, 1, maxnhg, name='listhg',
     &               routine= 'diagkp' )
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = 
     &        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     &    BNode,MPI_Comm_World,MPIerror)
      enddo

C     Create new distribution of H and S
      nuog = nuotot
      nullify(Snew,Hnew,xijloc)
      call re_alloc( Snew, 1, maxnhg, name='listhg',
     &               routine= 'diagkp' )
      call re_alloc( Hnew, 1, maxnhg, 1, nspin, name='Hnew',
     &               routine= 'diagkp' )
      call re_alloc( xijloc, 1, 3, 1, maxnhg, name='xijloc',
     &               routine= 'diagkp' )

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(iio)
              Hnew(listhptrg(io)+jo,is) = H(listhptr(iio)+jo,is)
            enddo
          enddo
          do jo = 1,numh(iio)
            Snew(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          enddo
          do jo = 1,numh(iio)
            xijloc(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
          enddo
        endif
        do is = 1,nspin
          call MPI_Bcast(Hnew(listhptrg(io)+1,is),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        enddo
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),
     &    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(xijloc(1,listhptrg(io)+1),3*numhg(io),
     &    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
      enddo

C Find eigenvalues .........................................
      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin
!$OMP parallel default(shared),
!$OMP&private(iuo,juo,io,j,ind,jo,kxij,ckxij,skxij)
!$OMP do
          do iuo = 1,nuog
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
              Haux(1,juo,iuo) = 0.0d0
              Haux(2,juo,iuo) = 0.0d0
            enddo
          enddo
!$OMP end do
!$OMP do
          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     &               kpoint(2,ik) * xijloc(2,ind) +
     &               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind)*ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind)*skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,ispin)*ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,ispin)*skxij
            enddo
          enddo
!$OMP end do nowait
!$OMP end parallel
          call cdiag(Haux,Saux,nuotot,nuotot,nuog,eo(1,ispin,ik),psi,
     &               -neigwanted,iscf,ierror, -1) ! dummy blocksize
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif
        enddo
      enddo

C Globalise eigenvalues
      BNode = -1
      do ik = 1,nk
        BNode = BNode + 1
        if (BNode.eq.Nodes) BNode = 0
        call MPI_Bcast(eo(1,1,ik),nuotot*nspin,MPI_double_precision,
     &    BNode,MPI_Comm_World,MPIerror)
      enddo

C Check if we are done ................................................
      if (.not.getD) goto 999

      nullify(Dnewloc,Enewloc)
C Allocate local copies of Dk and Ek for building matrices
      call re_alloc( Dnewloc, 1, maxnhg, 1, 2, name='Dnewloc',
     &               routine= 'diagkp' )
      call re_alloc( Enewloc, 1, maxnhg, 1, 2, name='Enewloc',
     &               routine= 'diagkp' )

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
         call fermispin( nspin, nspin, nk, wk, maxo,
     &        neigwanted, eo, temp, qs, qo, efs, Entropy )
      else
         call fermid(nspin, maxspn, nk, wk, maxo,
     &        neigwanted, eo, temp, qtot, qo, ef, Entropy )
      endif


C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
!$OMP parallel do default(shared), private(ik,ispin,io), firstprivate(t)
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuotot
              qo(io,ispin,ik) = wk(ik) * 
     &             ( stepf((eo(io,ispin,ik)-e2)/t) -
     &               stepf((eo(io,ispin,ik)-e1)/t)) * 2.0d0 / nspin
            enddo
          enddo
        enddo
!$OMP end parallel do
      endif
      
C New density and energy-density matrices of unit-cell orbitals .......

      Dnewloc(1:maxnhg,1:nspin) = 0.0d0
      Enewloc(1:maxnhg,1:nspin) = 0.0d0

      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin

C Find maximum eigenvector that is required for this k point and spin
          neigneeded = 0
          ie = nuog
          do while (ie.gt.0.and.neigneeded.eq.0)
            qe = qo(ie,ispin,ik)
            if (abs(qe).gt.occtol) neigneeded = ie
            ie = ie - 1
          enddo

C Find eigenvectors
          Saux = 0.0d0
          Haux = 0.0d0

!$OMP parallel do default(shared),
!$OMP&private(io,j,ind,jo,iuo,juo,kxij,ckxij,skxij)
          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     &               kpoint(2,ik) * xijloc(2,ind) +
     &               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind)*ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind)*skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,ispin)*ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,ispin)*skxij
            enddo
          enddo
!$OMP end parallel do
          call cdiag(Haux,Saux,nuotot,nuotot,nuog,aux,psi,
     &               neigneeded,iscf,ierror, -1) !dummy blocksize

C Check error flag and take appropriate action
          if (ierror.gt.0) then
            call die('Terminating due to failed diagonalisation')
          elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
            Saux = 0.0d0
            Haux = 0.0d0
!$OMP parallel do default(shared),
!$OMP&private(io,j,ind,jo,iuo,juo,kxij,ckxij,skxij)
            do io = 1,nuog
              do j = 1,numhg(io)
                ind = listhptrg(io) + j
                jo = listhg(ind)
                iuo = indxuo(io)
                juo = indxuo(jo)
                kxij = kpoint(1,ik) * xijloc(1,ind) +
     &                 kpoint(2,ik) * xijloc(2,ind) +
     &                 kpoint(3,ik) * xijloc(3,ind)
                ckxij = cos(kxij)
                skxij = sin(kxij)
                Saux(1,juo,iuo) = Saux(1,juo,iuo) +Snew(ind)*ckxij
                Saux(2,juo,iuo) = Saux(2,juo,iuo) -Snew(ind)*skxij
                Haux(1,juo,iuo) = Haux(1,juo,iuo) +Hnew(ind,ispin)*ckxij
                Haux(2,juo,iuo) = Haux(2,juo,iuo) -Hnew(ind,ispin)*skxij
              enddo
            enddo
!$OMP end parallel do
            call cdiag(Haux,Saux,nuotot,nuotot,nuog,aux,psi,
     &                 neigneeded,iscf,ierror, -1) ! dummy blocksize
          endif
!
!         Here we could write the wavefunctions, taking care to
!         use a suitable version of writew. There are also problems
!         of ordering of the different functions...
!
!          if (getPSI) then
!             call writew_one_node(nuotot,nuog,ik,kpoint(1,ik),ispin,
!     .                 aux,psi,.false.)
!          endif

          Dk = 0.0d0
          Ek = 0.0d0

!$OMP parallel default(shared),
!$OMP&private(ie,qe,j,ee,iuo,juo,io,ind,jo,kxij,ckxij,skxij),
!$OMP&private(pipj1,pipj2)

C Add contribution to density matrices of unit-cell orbitals

C Global operation to form new density matrix
          do ie = 1,nuotot
            qe = qo(ie,ispin,ik)
            if (abs(qe).gt.occtol) then
!$OMP do
              do j = 1,nuotot
                aux(1,j) = psi(1,j,ie)
                aux(2,j) = psi(2,j,ie)
              enddo
!$OMP end do
              ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
!$OMP do
              do iuo = 1,nuog
                do juo = 1,nuotot
                  pipj1 = aux(1,iuo) * aux(1,juo) +
     &                    aux(2,iuo) * aux(2,juo)
                  pipj2 = aux(1,iuo) * aux(2,juo) -
     &                    aux(2,iuo) * aux(1,juo)
                  Dk(1,juo,iuo) = Dk(1,juo,iuo) + qe * pipj1
                  Dk(2,juo,iuo) = Dk(2,juo,iuo) + qe * pipj2
                  Ek(1,juo,iuo) = Ek(1,juo,iuo) + ee * pipj1
                  Ek(2,juo,iuo) = Ek(2,juo,iuo) + ee * pipj2
                enddo
              enddo
!$OMP end do
            endif
          enddo
!$OMP do
          do io = 1,nuog
            iuo = indxuo(io)
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     &               kpoint(2,ik) * xijloc(2,ind) +
     &               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Dnewloc(ind,ispin) = Dnewloc(ind,ispin) + 
     &          Dk(1,juo,iuo)*ckxij - Dk(2,juo,iuo)*skxij
              Enewloc(ind,ispin) = Enewloc(ind,ispin) + 
     &          Ek(1,juo,iuo)*ckxij - Ek(2,juo,iuo)*skxij
            enddo
          enddo
!$OMP end do nowait

!$OMP end parallel

        enddo
      enddo


      do io = 1,nuog
        call WhichNodeOrb(io,Nodes,BNode)
        call GlobalToLocalOrb(io,BNode,Nodes,iio)
        do ispin = 1,nspin
          call MPI_Reduce(Dnewloc(listhptrg(io)+1,ispin),
     &      Dloc(1),numhg(io),MPI_double_precision,
     &      MPI_sum,BNode,MPI_Comm_World,MPIerror)
          if (Node.eq.BNode) then
            do j = 1,numh(iio)
              Dnew(listhptr(iio)+j,ispin) = Dloc(j)
            enddo
          endif
          call MPI_Reduce(Enewloc(listhptrg(io)+1,ispin),
     &      Dloc(1),numhg(io),MPI_double_precision,
     &      MPI_sum,BNode,MPI_Comm_World,MPIerror)
          if (Node.eq.BNode) then
            do j = 1,numh(iio)
              Enew(listhptr(iio)+j,ispin) = Dloc(j)
            enddo
          endif
        enddo
      enddo

C Deallocate local memory needed only if getD is true
      call de_alloc( Enewloc, name='Enewloc', routine= 'diagkp' )
      call de_alloc( Dnewloc, name='Dnewloc', routine= 'diagkp' )

C Exit point 
  999 continue

C Free local memory
      call de_alloc( xijloc, name='xijloc', routine= 'diagkp' )
      call de_alloc( Hnew, name='Hnew', routine= 'diagkp' )
      call de_alloc( Snew, name='Snew', routine= 'diagkp' )
      call de_alloc( Dloc, name='Dloc', routine= 'diagkp' )
      call de_alloc( listhg, name='listhg', routine= 'diagkp' )
      call de_alloc( listhptrg, name='listhptrg', routine= 'diagkp' )
      call de_alloc( numhg, name='numhg', routine= 'diagkp' )

#ifdef DEBUG
      call write_debug( '    POS diagkp' )
#endif

      end subroutine diagkp
#endif
