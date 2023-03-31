! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diagk( h_spin_dim, nuo, no, maxspn, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, listdptr, 
     .                  listd, H, S, getD, getPSI, fixspin, qtot, qs, 
     .                  temp, e1, e2, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     .                  Haux, Saux, psi, Dk, Ek, aux, nuotot, occtol,
     .                  iscf, neigwanted )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer maxspn=spinor_dim   : # of spin components of eo and qo, qs, efs
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
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
C real*8  H(maxnh,h_spin_dim) : Hamiltonian in sparse form
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
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer neigwanted          : maximum number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo,spinor_dim,nk)   : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,spinor_dim,nk)    : Occupations of eigenstates
C real*8 Dnew(maxnd,h_spin_dim): Output Density Matrix
C real*8 Enew(maxnd,e_spin_dim): Output Energy-Density Matrix
C real*8 ef                    : Fermi energy
C real*8 Entropy               : Electronic entropy
C *************************** AUXILIARY *******************************
C integer h_spin_dim          : Number of spin components of H and D
C integer e_spin_dim          : Number of spin components of E_dm  
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
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C Optimization of Dscf building by Alberto Garcia (August 2009)
C
C NOTE that the diagk_file routine is still much more efficient in terms of
C CPU time.
C
C  Modules
C
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave,     only : writew
      use m_fermid,      only : fermid, fermispin, stepf
      use m_spin,        only : spinor_dim, e_spin_dim
      use iso_c_binding, only : c_loc, c_f_pointer
#ifdef MPI
      use mpi_siesta
#endif
      use m_norm

      implicit          none


#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer           maxnd, maxnh, maxspn, maxo, nk, no,
     .                  nuo, nuotot, iscf, neigwanted
      integer           h_spin_dim
      integer           indxuo(no), listh(maxnh), numh(nuo),
     .                  listd(maxnd), numd(nuo), listhptr(nuo),
     .                  listdptr(nuo)
      real(dp)          Dnew(maxnd,h_spin_dim), 
     .                  e1, e2, ef, efs(spinor_dim),
     .                  Enew(maxnd,e_spin_dim), Entropy, 
     .                  eo(maxo,spinor_dim,nk), H(maxnh,h_spin_dim), 
     .                  kpoint(3,nk), qo(maxo,spinor_dim,nk), qtot, 
     .                  S(maxnh), temp, wk(nk), occtol,
     .                  xij(3,maxnh), qs(spinor_dim)
      real(dp)          Dk(2,nuotot,nuo), Ek(2,nuotot,nuo),
     .                  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo)
      real(dp), target :: psi(2,nuotot,nuo), aux(2,nuotot)
      logical           getD, getPSI, fixspin
      external          cdiag

C  Internal variables .............................................
      integer
     .  BNode, BTest, ie, ierror, iie, ik, ind, io, iio,
     .  ispin, iuo, j, jo, juo, nd, neigneeded

      real(dp)
     .  ckxij, kxij, qe, ee, skxij, t, rt
      real(dp) :: qp1, qp2, eqp1, eqp2

      real(dp), pointer :: lpsi(:,:)
      real(dp), pointer :: psi_real_1d(:)

      logical, allocatable :: done_juo(:)
C     ....................

#ifdef DEBUG
      call write_debug( '    POS diagk' )
#endif

      if (h_spin_dim /= spinor_dim) then
         call die("Spin size mismatch in diagk")
      endif


      allocate(done_juo(nuotot))

C Skip eigenvalues if DM is not required and WF's are required ....
      if (getPSI .and. .not. getD) goto 10

C Find eigenvalues ................................................
      do ik = 1,nk
        do ispin = 1,spinor_dim
          call timer( 'c-eigval', 1 )
          call timer( 'c-buildHS', 1 )
!$OMP parallel do default(shared), 
!$OMP&private(iuo,j,ind,jo,juo,kxij,ckxij,skxij)
          do iuo = 1,nuo
            Saux(:,:,iuo) = 0._dp
            Haux(:,:,iuo) = 0._dp
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
!$OMP end parallel do
          call timer( 'c-buildHS', 2 )
          call cdiag( Haux, Saux, nuotot, nuo, nuotot, eo(1,ispin,ik),
     &                psi, -neigwanted, iscf, ierror, BlockSize )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif
          call timer( 'c-eigval', 2 )
        enddo
      enddo

C Check if we are done ................................................
      if (.not.getD .and. .not. getPSI) then
#ifdef DEBUG
         call write_debug( '    PRE diagk' )
#endif
         return
      end if

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
         call fermispin( spinor_dim, spinor_dim, nk, wk, maxo,
     &        neigwanted, eo, temp, qs, qo, efs, Entropy )
      else
         call fermid(spinor_dim, spinor_dim, nk, wk, maxo,
     &        neigwanted, eo, temp, qtot, qo, ef, Entropy )
      endif

      nd = listdptr(nuo) + numd(nuo)

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.e-6_dp )
        rt = 1._dp / t
!$OMP parallel do default(shared), private(ik,ispin,io)
        do ik = 1,nk
          do ispin = 1,spinor_dim
            do io = 1,nuotot
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf((eo(io,ispin,ik)-e2)*rt) -
     .               stepf((eo(io,ispin,ik)-e1)*rt)) * 2.0_dp/spinor_dim
            enddo
          enddo
        enddo
!$OMP end parallel do
      endif

C New density and energy-density matrices of unit-cell orbitals .......
      if (nuo.gt.0) then
        Dnew(1:nd,1:h_spin_dim) = 0.0_dp
        Enew(1:nd,1:e_spin_dim) = 0.0_dp
      endif


  10  continue

C Loop over k points
      do ik = 1,nk

        do ispin = 1, spinor_dim
C Find maximum eigenvector that is required for this k point and spin
C Compute all eigenvectors if getPSI = .true.

          if (getPSI) then
             ! We do not have info about occupation
             ! In the current implementation, neigneeded
             ! is ALWAYS nuotot.
             ! This is because writewave expects the full spectrum
             neigneeded = neigwanted
          else
            neigneeded = 0
            ie = nuotot
            do while (ie.gt.0.and.neigneeded.eq.0)
              qe = qo(ie,ispin,ik)
              if (abs(qe).gt.occtol) neigneeded = ie
              ie = ie - 1
            enddo
         endif

          call timer( 'c-eigvec', 1 )
C Build dense H and S - note that only eigenvectors required are now
C generated by the diagonaliser and that the scratch array aux is 
C used to hold the eigenvalues to prevent corruption of the full set
C determined on the first call.
!$OMP parallel do default(shared), 
!$OMP&private(iuo,j,ind,jo,juo,kxij,ckxij,skxij)
          do iuo = 1,nuo
            Saux(:,:,iuo) = 0._dp
            Haux(:,:,iuo) = 0._dp
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
            enddo
          enddo
!$OMP end parallel do

C Find eigenvectors
          call cdiag(Haux,Saux,nuotot,nuo,nuotot,aux,psi,
     .      neigneeded,iscf,ierror,BlockSize)

C Check error flag and take appropriate action
          if (ierror.gt.0) then
            call die('Terminating due to failed diagonalisation')
          elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
!$OMP parallel do default(shared), 
!$OMP&private(iuo,j,ind,jo,juo,kxij,ckxij,skxij)
            do iuo = 1,nuo
              Saux(:,:,iuo) = 0._dp
              Haux(:,:,iuo) = 0._dp
              do j = 1,numh(iuo)
                ind = listhptr(iuo) + j
                jo = listh(ind)
                juo = indxuo(jo)
                kxij = kpoint(1,ik) * xij(1,ind) +
     .                 kpoint(2,ik) * xij(2,ind) +
     .                 kpoint(3,ik) * xij(3,ind)
                ckxij = cos(kxij)
                skxij = sin(kxij)
                Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
                Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
                Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
                Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
              enddo
            enddo
!$OMP end parallel do
            call cdiag(Haux,Saux,nuotot,nuo,nuotot,aux,psi,
     .        neigneeded,iscf,ierror,BlockSize)
          endif

          call timer( 'c-eigvec', 2 )

          if (getPSI) then
             ! We do not have info about occupation
             ! In the current implementation, neigneeded
             ! is ALWAYS nuotot.
             ! This is because writewave expects the full spectrum
             if (.not. getD) then
                ! Fill in the missing eigenvalue information
                ! This is useful if we are requesting WFS in a "bands" setting
               call dcopy(nuotot, aux(1,1), 1, eo(1,ispin,ik), 1)
             endif
             call c_f_pointer(c_loc(psi),psi_real_1d,[size(psi)])
             call writew(nuotot,nuo,ik,kpoint(1,ik),ispin,
     .                 aux,psi_real_1d,gamma=.false.,non_coll=.false.,
     $                 blocksize=BlockSize)
          endif

          if (getD) then
             call timer( 'c-buildD', 1 )

C Add contribution to density matrices of unit-cell orbitals
C WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
!$OMP parallel do default(shared), private(iuo,juo)
            do iuo = 1,nuo
              do juo = 1,nuotot
                Dk(1,juo,iuo) = 0.0_dp
                DK(2,juo,iuo) = 0.0_dp
                Ek(1,juo,iuo) = 0.0_dp
                Ek(2,juo,iuo) = 0.0_dp
              enddo
            enddo
!$OMP end parallel do

C Global operation to form new density matrix
            BNode = 0
            iie = 0
            do ie = 1,nuotot
              if ( Node == BNode ) iie = iie + 1
              qe = qo(ie,ispin,ik)
              if (abs(qe).gt.occtol) then
                 
                if ( Node == BNode ) then
                   lpsi => psi(:,:,iie)
                else
                   lpsi => aux(:,:)
                endif
#ifdef MPI
                call MPI_Bcast(lpsi(1,1),2*nuotot,MPI_double_precision,
     &               BNode,MPI_Comm_World,MPIerror)
#endif

                ee = eo(ie,ispin,ik)
                
!$OMP parallel do default(shared), 
!$OMP&private(iuo,iio,done_juo,ind,juo,qp1,qp2,eqp1,eqp2)
                do iuo = 1,nuo
                   call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                   
                   qp1  = qe * lpsi(1,iio)
                   qp2  = qe * lpsi(2,iio)

                   ! Process only the elements really needed
                   ! Dk and Ek are re-used storage, so
                   ! no memory is really wasted  [AG, Aug 2009]
                   done_juo(1:nuotot) = .false.
                   do ind = listdptr(iuo) + 1 ,
     &                  listdptr(iuo) + numd(iuo)
                      juo = indxuo(listd(ind))
                      if (done_juo(juo)) cycle
                      eqp1 = qp1 * lpsi(1,juo) + qp2 * lpsi(2,juo)
                      eqp2 = qp1 * lpsi(2,juo) - qp2 * lpsi(1,juo)

                      Dk(1,juo,iuo) = Dk(1,juo,iuo) + eqp1
                      Dk(2,juo,iuo) = Dk(2,juo,iuo) + eqp2
                      Ek(1,juo,iuo) = Ek(1,juo,iuo) + eqp1 * ee
                      Ek(2,juo,iuo) = Ek(2,juo,iuo) + eqp2 * ee
                      done_juo(juo) = .true.
                   enddo
                enddo
!$OMP end parallel do
              endif
              BTest = ie/BlockSize
              if (BTest*BlockSize.eq.ie) then
                BNode = BNode + 1
                if (BNode .gt. Nodes-1) BNode = 0
              endif
            enddo

!$OMP parallel do default(shared),
!$OMP&private(iuo,ind,juo,kxij,ckxij,skxij)
            do iuo = 1,nuo
               do ind = listdptr(iuo) + 1 ,
     &              listdptr(iuo) + numd(iuo)
                juo = indxuo(listd(ind))
                kxij = kpoint(1,ik) * xij(1,ind) +
     .                 kpoint(2,ik) * xij(2,ind) +
     .                 kpoint(3,ik) * xij(3,ind)
                ckxij = cos(kxij)
                skxij = sin(kxij)
                Dnew(ind,ispin)=Dnew(ind,ispin)+ Dk(1,juo,iuo)*ckxij -
     .                                           Dk(2,juo,iuo)*skxij
                Enew(ind,ispin)=Enew(ind,ispin)+ Ek(1,juo,iuo)*ckxij -
     .                                           Ek(2,juo,iuo)*skxij
              enddo
            enddo
!$OMP end parallel do

            call timer( 'c-buildD', 2 )
          endif

        enddo
      enddo
      deallocate(done_juo)

      end subroutine diagk
