! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_chempot
      use precision, only: dp
      use alloc,     only: re_alloc, de_alloc

      implicit none

      private
      public :: chempot

        ! Formerly in numbvect
        integer   p, nb
        real(dp)    qtot
        real(dp), dimension(:), pointer ::   c(:), rr(:,:)


        ! Formerly in common
        real(dp)   :: betap

      CONTAINS

      subroutine chempot(h,s,numh,listhptr,listh,rcoor,pmax,beta,
     .                   lasto,cell,xa,enum,nbasis,nbasisloc,
     .                   natoms,maxnh,chpot,emax,emin)

C     .,gap,homo,lumo)
C *****************************************************************
C Calculates the maximum and minimum eigenvalues, the chemical 
C potential and the HOMO-LUMO gap.
C The calculation of the max. and min. eigenvalues is done with
C the Lanczos method. The chemical potential is calculated 
C with the projection method of Goedecker, and the HOMO-LUMO
C gap is obtained by the Folded-Spectrum method (combined with
C Lanczos)
C
C NOTE: In this version, the calculation of the HOMO, LUMO
C       and gap is disabled because it doesnot work properly.
C
C NOTE: Parallelism has only been introduced to the extent of
C       handling the distribution of H and S. Needs more
C       extensive work to distribute local data to save
C       memory for large problems.
C
C Written by Maider Machado and P.Ordejon, June'98 
C ***************************** INPUT ***************************** 
C real*8 h(maxnh)                 : Hamiltonian in sparse form
C real*8 s(maxnh)                 : Overlap in sparse form
C integer numh(nbasisloc)         : Control vector of sparse hamilt.
C integer listhptr(nbasisloc)     : Control vector of sparse hamilt.
C integer listh(maxnh)            : Control vector of sparse hamilt.
C real*8 rcoor                    : Cutoff range for the projection 
C                                   vectors
C integer pmax                    : Maximum order of Chebishev expansion
C real*8 beta                     : Inverse Temperature for Chebi expansion
C integer lasto(0:natoms)         : Index vector of last orbital of 
C                                   each atom
C real*8 cell(3,3)                : Lattice vectors
C real*8 xa(3,natoms)             : Atomic positions
C real*8 enum                     : Total number of electrons
C integer nbasis                  : Number of basis orbitals (global)
C integer nbasisloc               : Number of basis orbitals (local)
C integer natoms                  : Number of atoms
C integer maxnh                   : Maximum number of non-zero elements
C                                   of sparse Hamiltonian
C **************************** OUTPUT *******************************
C real*8 chpot                    : Chemical potential
C real*8 emax                     : Maximum eigenvalue
C real*8 emin                     : Minimum eigenvalue
C *************************** DISABLED ******************************
C real*8 gap                      : Energy gap
C real*8 homo                     : Highest occ. molecular orbital energy
C real*8 lumo                     : Lowest unocc. molecular orbital energy
C *********************************************************************

C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use sys,           only : die
      use neighbour,     only : jan, r2ij, xij, mneighb, maxnna,
     &                          reset_neighbour_arrays
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer 
     .  natoms, nbasis, nbasisloc, maxnh, pmax

      integer
     .  lasto(0:natoms), listh(maxnh), numh(nbasisloc), 
     .  listhptr(nbasisloc)

      real(dp)
     .  beta, cell(3,3), chpot, emin, emax, enum, 
     .  h(maxnh), rcoor, s(maxnh), xa(3,natoms)
C    .  homo, lumo, gap

C Internal variables...

      integer ::  ind, nhmax, nvmaxnew

C nvmax = maximum number of non-zero elements within a sparse vector v
      integer, save ::   nvmax = 1000

      real(dp),  parameter :: mu1=-1.2d0 , mu2=1.2d0

      integer
     .  i, ia, ii, iil, imu, iorb, j, ja, ji, jj, jorb, 
     .  k, m, mu, nna, norb, nu, num, numloc, numv

      integer, dimension(:), pointer :: listvt, numhp, indexloc,
     &                                  ibuffer, listv

      real(dp)  chpotsh, deltae, emean,
     .  rmax, rrmod, ri(3), tol

      integer,  dimension(:,:), pointer :: listhp, listhpp
      real(dp), dimension(:),   pointer :: paux, vec
      real(dp), dimension(:,:), pointer :: hbar, Hdense, Sdense
      real(dp), dimension(:,:), pointer :: v

#ifdef MPI
      integer ::  MPIerror
      real(dp), dimension(:),   pointer :: dpbuffer1
#endif
C ...

C Start timer
      call timer( 'chempot', 1 )
      print *, "CHEMPOT CALLED"

C Assign information for Chebyshev expansion (module variables)
      qtot = enum
      nb = nbasis
      p = pmax

C Find size of maxnh
      nhmax = 0
      do i = 1,nbasisloc
        nhmax = max(nhmax,numh(i))
      enddo

C Allocate arrays that depend on maxnh

      nullify( listhp )
      call re_alloc( listhp, 1, nhmax, 1, nbasis, name='listhp',
     &               routine='chempot' )
      nullify( listhpp )
      call re_alloc( listhpp, 1, nhmax, 1, nbasis, name='listhpp',
     &               routine='chempot' )
      nullify( paux )
      call re_alloc( paux, 1, nhmax, name='paux', routine='chempot' )
      nullify( vec )
      call re_alloc( vec, 1, nhmax, name='vec', routine='chempot' )
      nullify( hbar )
      call re_alloc( hbar, 1, nhmax, 1, nbasis, name='hbar',
     &               routine='chempot' )
      nullify( Hdense )
      call re_alloc( Hdense, 1, nhmax, 1, nhmax, name='Hdense',
     &               routine='chempot' )
      nullify( Sdense )
      call re_alloc( Sdense, 1, nhmax, 1, nhmax, name='Sdense',
     &               routine='chempot' )

C Allocate arrays that depend on nbasis

      nullify( listvt )
      call re_alloc( listvt, 1, nbasis, name='listvt',
     &               routine='chempot' )
      nullify( numhp )
      call re_alloc( numhp, 1, nbasis, name='numhp',
     &               routine='chempot' )
      nullify( rr )
      call re_alloc( rr, 1, nbasis, 0, pmax, name='rr',
     &               routine='chempot' )

C Allocate arrays that depend on pmax

      nullify( c )
      call re_alloc( c, 1, pmax, name='c', routine='chempot' )

C Allocate arrays that depend on maxnna

      nullify( indexloc )
      call re_alloc( indexloc, 1, max(maxnna,natoms),
     $           name='indexloc',  routine='chempot' )

C Allocate arrays that depend on nvmax

      nullify( listv )
      call re_alloc( listv, 1, nvmax, name='listv',
     &               routine='chempot' )
      nullify( v )
      call re_alloc( v, 1, nvmax, 1, 3, name='v', routine='chempot' )

C *** Calculate H' = S(-1)*H  (See Gibson et al, PRB 47, 9229 (92)) ***
C This is done by Choleski decomposition of S, and solving the
C linear system S H' = H, in the subspace of orbitals which overlap
C with those of a given atom. The advantage is that the Choleski
C decomposition can be done only once for each atom, and use
C the result for all the orbitals.
C Since different orbitals in the same atom have differen
C cutoff radii, this must be done carefully, since the list of
C neighbors is not the same for all orbitals.


C  Loop over atoms.............
      do ia = 1,natoms

C   form S and S-1 matrices in reduced space...

C   first determine which is the longer range orbital of atom ia
C   (which will determine the reduce space) -
        num = 0
        mu = 0
        do i = lasto(ia-1)+1,lasto(ia)
          call GlobalToLocalOrb(i,Node,Nodes,ii)
          if (numh(ii) .gt. num) then
            mu = ii
            num = numh(ii)
          endif
        enddo
        if (mu .eq. 0) then
          if (Node.eq.0) then
            write(6,*) 'chempot: ERROR: zero neighbors for orbital ',i
          endif
          CALL DIE()
        endif
        nu = numh(mu) 

C Initialize overlap and interaction in reduced space -
        do i = 1,nu
          do j = 1,nu
            Sdense(j,i) = 0.0d0
            Hdense(j,i) = 0.0d0
          enddo
        enddo

C Construct S and H in reduced space  -
C  loop over orbitals ii in reduced space
        do i = 1,nu
          ii = listh(listhptr(mu)+i)
          call GlobalToLocalOrb(ii,Node,Nodes,iil)
          if (iil.gt.0) then
C  loop over orbitals jj in reduced space
            do j = i,nu
              jj = listh(listhptr(mu)+j)
C  see if orbitals ii and jj interact
              do k = 1,numh(iil)
                ind = listhptr(iil) + k
                if (listh(ind) .eq. jj) then
                  Sdense(i,j) = s(ind)
                  Sdense(j,i) = Sdense(i,j)
                  Hdense(i,j) = h(ind)
                  Hdense(j,i) = Hdense(i,j)
                endif
              enddo
            enddo
          endif
        enddo

#ifdef MPI
C  Globalise Hdense and Sdense for now
        nullify( dpbuffer1 )
        call re_alloc( dpbuffer1, 1, nu, name='dpbuffer1',
     &                 routine='chempot' )
        do i = 1,nu
          call MPI_AllReduce(Sdense(1,i),dpbuffer1,nu,
     .      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
          do j = 1,nu
            Sdense(j,i) = dpbuffer1(j)
          enddo
          call MPI_AllReduce(Hdense(1,i),dpbuffer1,nu,
     .      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
          do j = 1,nu
            Hdense(j,i) = dpbuffer1(j)
          enddo
        enddo
        call de_alloc( dpbuffer1, name='dpbuffer1',routine="chempot")
#endif

C  Cholesky factorization:
        call cholDcmp(Sdense,nu,nhmax,paux)

C Loop over orbitals of atom ia ....
        do i = 1,nu
          ii = listh(listhptr(mu)+i)
C  Check if ii is in atom ia
          if (ii .ge. (lasto(ia-1)+1) .and. ii .le. lasto(ia)) then
            do j = 1,nu
              vec(j) = Hdense(j,i)
            enddo
            call cholLinSys(Sdense,nu,nhmax,paux,vec,vec)

C  vec contains the elements of H' in reduced space.
C  now H' must be formed in sparse format -

c  global indes of orbital i of reduced space
            ii = listh(listhptr(mu)+i)
            call GlobalToLocalOrb(ii,Node,Nodes,iil)
            do j = 1,nu
c  global indes of orbital j of reduced space
              jj = listh(listhptr(mu)+j)
              do k = 1, numh(iil)
                if (listh(listhptr(iil)+k) .eq. jj) hbar(k,ii) = vec(j)
              enddo
            enddo

          endif
        enddo

      enddo
C ............

C *** Compute smallest and largest eigenvalues (Lanczos Method) ***

      call lanc1(2,hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .  emin,Node)
      call lanc1(1,hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .  emax,Node)

      emean = 0.5d0*(emax+emin)
      deltae = 0.55d0*(emax-emin)

C *** Calculate Chemical Potential using the Projection Method of
C              Goedecker (PRB 51,9455 (95)). ***
C 
C rr(in,ip) stores the in-th element of the vector resulting from 
C the application of the ip-th Chebyshev polynomial to the in-th atomic
C orbital. This is all what is needed to calculate the number of 
C electrons.

C  First scale and shift the hamiltonian ...

      do j = 1,nbasis
        do i = 1,numh(j)
          if (listh(listhptr(j)+i).eq.j) then
            hbar(i,j) = (hbar(i,j)-emean)/deltae
          else
            hbar(i,j) = hbar(i,j)/deltae
          endif
        enddo
      enddo     

C Calculate maximum length in unit cell ...
  ! i.e, maximum distance among (-1:1) lattice points
!AG: This seems to be wrong, as the cell vectors are given
!    by columns...
      rmax = 0.0d0
      do i = -1,1
        do j = -1,1
          do k = -1,1
            !AG: ri(:) = i*cell(:,1) + j*cell(:,2) + k*cell(:,3)
            ri(1) = i*cell(1,1) + j*cell(2,1) + k*cell(3,1)
            ri(2) = i*cell(1,2) + j*cell(2,2) + k*cell(3,2)
            ri(3) = i*cell(1,3) + j*cell(2,3) + k*cell(3,3)
            rrmod = sqrt( ri(1)**2 + ri(2)**2 + ri(3)**2 )
            if (rrmod .gt. rmax) rmax = rrmod
          enddo
        enddo
      enddo

C Initialize routine for neighbour search
      if (2.*rcoor .lt. rmax) then
        call mneighb(cell,rcoor,natoms,xa,0,0,nna)
      endif

C initialize control vectors to zero 
      listvt(1:nbasis) = 0

C Loop over atoms ...............
      do ia = 1,natoms

        if (2.0*rcoor .lt. rmax) then
C  look for neighbors of atom ia
          call mneighb(cell,rcoor,natoms,xa,ia,0,nna)
            ! No need to copy old contents in indexloc
          call re_alloc( indexloc, 1, maxnna, name='indexloc',
     &                   routine='chempot')
        else    
          ! Make sure we have enough space
          ! (Cleaner alternative to initial allocation with  
          ! max(maxnna,natoms)
          !call sizeup_neighbour_arrays(natoms)
          !call re_alloc( indexloc, 1, natoms, name='indexloc')
          nna = natoms
          do jj = 1,natoms
            jan(jj) = jj
          enddo
        endif

C Build structure of sparse vector v ...

C Clear list ot atoms considered within loc. range ...
        indexloc(1:nna) = 0
        numloc = 0

        numv = 0
        do 30 j = 1,nna
          ja = jan(j)

C Check if ja has already been included in current vector ...
          do jj = 1,numloc
            if (ja .eq. indexloc(jj)) goto 30
          enddo
          numloc = numloc + 1
          indexloc(numloc) = ja

          do jorb = 1,lasto(ja) - lasto(ja-1)
            nu = jorb + lasto(ja-1)
            numv = numv + 1
            if (numv .gt. nvmax) then
              nvmaxnew = numv + nint(0.1*numv)
C
              call re_alloc( listv, 1, nvmaxnew, name='listv',
     &                       routine='chempot', copy=.true. )
              call re_alloc( v, 1, nvmaxnew, 1, 3, name='v',
     &                       routine='chempot', copy=.true. )
C
              nvmax = nvmaxnew
            endif
            listv(numv) = nu
            listvt(nu) = numv
          enddo
30      continue

c number of orbitals of atom ia
        norb = lasto(ia) - lasto(ia-1)

c loop over orbitals of atom ia ...
        do iorb = 1,norb
          mu = iorb + lasto(ia-1)

          v(1:numv,1:2) = 0.0d0

          imu = listvt(mu)

          v(imu,1) = 1.0d0
          rr(mu,0) = v(imu,1)

          do j = 1,numv
            ji = listv(j)
            numhp(ji) = 0
            do i = 1,numh(ji)
              m = listh(listhptr(ji)+i)
              jj = listvt(m)
              if (jj .ne. 0) then
                numhp(ji) = numhp(ji) + 1
                listhp(numhp(ji),ji) = jj
                listhpp(numhp(ji),ji) = i
                v(j,2) = v(j,2) + hbar(i,ji)*v(jj,1)
              endif
            enddo
          enddo

          rr(mu,1) = v(imu,2)

          do 40 k = 2,p-1

            do j = 1,numv
              v(j,3) = - v(j,1)
              ji = listv(j)
              do i = 1,numhp(ji)
                jj = listhp(i,ji)
                ii = listhpp(i,ji)
                v(j,3) = v(j,3) + 2.0*hbar(ii,ji)*v(jj,2)
              enddo
            enddo

            rr(mu,k) = v(imu,3)

            j = 0
            do i = 1,numv
              v(i,1) = v(i,2)
              v(i,2) = v(i,3)
            enddo

 40       continue  

        enddo

C Reset control vectors to cero
        do i = 1,numv
          imu = listv(i)
          listvt(imu) = 0
        enddo

      enddo
C ............

      tol = 0.0001d0
C Calculate chemical potential as the root of Nel - Tr(rho) = 0 ...

C Inverse temperature (in units of energy scaled so that the spectrum
C  lays between (-1,+1)

      betap = beta * deltae

C      do p=1,pmax
         chpotsh=root(numb,mu1,mu2,tol)
         ! Shift and scale the result to absolute energy scale ...
         chpot=chpotsh*deltae+emean

C      enddo

C Deallocate local arrays 
      if (2.*rcoor .lt. rmax) then
        call reset_neighbour_arrays( )
      endif
      call de_alloc( listhp, name='listhp' )
      call de_alloc( listhpp, name='listhpp' )
      call de_alloc( paux, name='paux' )
      call de_alloc( vec, name='vec' )
      call de_alloc( hbar, name='hbar' )
      call de_alloc( Hdense, name='Hdense' )
      call de_alloc( Sdense, name='Sdense' )
      call de_alloc( listvt, name='listvt' )
      call de_alloc( numhp, name='numhp' )
      call de_alloc( rr, name='rr' )
      call de_alloc( c, name='c' )
      call de_alloc( indexloc, name='indexloc' )
      call de_alloc( listv, name='listv' )
      call de_alloc( v, name='v' )

C  CALCULATION OF THE GAP DISABLED; IT DOES NOT WORK PROPERLY

C Stop timer
      call timer( 'chempot', 2 )

      return

CC Scale and shift the hamiltonian to absolute energy scale ...
C
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(listhptr(j)+i).eq.j) then
C            hbar(i,j)=deltae*hbar(i,j)+emean
C          else
C            hbar(i,j)=deltae*hbar(i,j)
C          endif
C        enddo
C      enddo
C
CC ...
C
C
CC  *** The gap is obtained by applying the Lanczos Method to
CC        the Folded Spectrum Method. See:
CC        Capaz-Koiler, J. Appl. Phys, 74, 5531 (93)
CC        Grosso et al, Nuovo Cimento D 15, 269 (93)
CC        Wang-Zunger, J. Chem. Phys. 100, 2394 (94) ***
C
C      eref1=chpot
C
CC Shift Hamiltonian...
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(listhptr(j)+i).eq.j) then
C            hbar(i,j)=hbar(i,j)-eref1
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref1)**2 by Lanczos...
C      call lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,eig1,Node)
CC ...
C
C      eig1 = eig1 + eref1
C
C      delta=eig1-eref1
C      eref2=eref1-delta
C
C50    continue
C
CC Shift Hamiltonian ...
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(listhptr(j)+i).eq.j) then
C            hbar(i,j)=hbar(i,j)+eref1-eref2
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref2)**2 by Lanczos ...
C      call lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,eig2,Node)
CC ...
C
C      eig2=eig2+eref2
C      gap=abs(eig1-eig2)
C
CC Check that levels are above and below the Fermi Level ...
C      if ((eig1 .gt. chpot .and. eig2 .gt. chpot) .or.
C     .    (eig1 .lt. chpot .and. eig2 .lt. chpot)) then
C        eref1=eref2
C        eref2=eref2-delta
C        goto 50
C      endif
CC ...
C
CC Convert to absolute energy scale ...
Cc      eig1 = eig1*deltae + emean
Cc      eig2 = eig2*deltae + emean
Cc      gap = gap*deltae
CC ...
C
CC Assign HOMO and LUMO ...
C      if (eig1. gt. eig2) then
C        homo=eig2
C        lumo=eig1
C      else
C        homo=eig1
C        lumo=eig2
C      endif
CC ...
C
C      return
      end subroutine chempot


      function numb(mu)
C **********************************************************************
C This function calculates the difference between the true number of
C electrons of the system, and the output number of electrons for a
C given Chemical Potential mu.
C
C Written by Maider Machado and P.Ordejon, June'98
C **********************************************************************
      use precision, only : dp

      implicit none

      real(dp) , intent(in) :: mu
      real(dp) :: Ne,numb

! Host association:  p, c, rr, nb, qtot

      integer k,n

      call chebfd(p,mu,c) 
c      write(16,*) qtot,nb
c      do ix=1,1001
c      x = -1. + 2.*(ix-1)/1000.
      Ne=0.0d0

      do n=1,nb
        Ne=Ne+0.5*c(1)*rr(n,0)+c(2)*rr(n,1)
        do k=1,p-2
           Ne=Ne+c(k+2)*rr(n,k+1) 
        enddo
      enddo

      numb=qtot-2.0d0*Ne

c      Ne=Ne+0.5d0*c(1)+c(2)*x
c      txm1 = 1
c      tx   = x
c      do k=1,p-2
c         txp1 = 2*x*tx - txm1
c         Ne=Ne+c(k+2)*txp1
c         txm1 = tx
c         tx   = txp1
c      enddo
c      write(7,*) x,Ne
c      enddo

      end function numb
     

      subroutine chebfd(n,mu,cof)
C ***********************************************************************
C Calculates the coefficients of the Chebyshev polynomials
C expansion of the Fermi-Dirac function.
C Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
C Wrtten by Maider Machado and P. Ordejon, June'98
C Restyled by J.M.Soler, May 2015
C ****************************** INPUT **********************************
C integer n                    : order of the expansion
C real*8 mu                    : chemical potential for F-D function
C ****************************** OUTPUT *********************************
C real*8 cof(n)                : expansion coefficients
C ***********************************************************************
  
      use precision, only: dp

!     Host association: betap

      integer, intent(in)   ::  n
      real(dp), intent(in)  ::  mu
      real(dp), intent(out) ::  cof(n)

      real(dp),parameter:: ymax = 50._dp
      real(dp):: f(n), pi, y(n)
      integer :: j, k

      ! Find Fermi-Dirac function
      pi = acos(-1.0_dp)
      do k = 1,n
        y(k) = cos(pi*(k-0.5_dp)/n)
        y(k) = betap*(y(k)-mu)
        if (abs(y(k))<ymax) then
          f(k)=1.0_dp/(exp(y(k))+1)
        elseif (y(k)<0._dp) then
          f(k) = 1
        else
          f(k) = 0
        endif
      enddo

      ! Find expansion coefficients
      do j=1,n
        forall(k=1:n) y(k) = cos(pi*(j-1)*(k-0.5_dp)/n)
        cof(j) = 2*sum(f*y)/n
      enddo
 
      END subroutine chebfd


      function root(func,xmin,xmax,tol) result(x0)
c  ***************************************************
c Returns the root of function func(x) bracketed by x1 and x2,
c with accuracy tol. Uses the bisection method.
c J.M.Soler, May.2015
c  ***************************************************     
      use precision, only: dp
      use sys, only      : die

      ! Arguments
      implicit none
      real(dp),intent(in):: xmin, xmax  ! interval bracketing root
      real(dp),intent(in):: tol     ! precision required for root
      real(dp)           :: x0      ! root value
      interface
         function func(x) result(y)
         use precision, only: dp
         real(dp), intent(in) :: x
         real(dp)             :: y
         end function func
      end interface

      real(dp):: x1, x2, y1, y2, y0

      ! Check input arguments
      x1 = xmin
      x2 = xmax
      y1 = func(x1)
      y2 = func(x2)
      if (y1*y2>0._dp) then
        call die('chempot/root: xmin,xmax not a root bracket')
      elseif (tol<abs(x2-x1)*1.e-12_dp) then
       call die('chempot/root: tol too small')
      endif

      ! Use bisection to reduce interval
      do while (abs(x2-x1)>tol)
        x0 = (x1+x2)/2
        y0 = func(x0)
        if (y0*y1>0.0_dp) then
          x1 = x0
          y1 = y0
        else
          x2 = x0
          y2 = y0
        endif
      enddo

      ! Linear approximation within last interval
      y0 = func(x0)
      x0 = (x1*(y2-y0)+x2*(y0-y1))/(y2-y1)

      end function root


      subroutine cholDcmp(a,n,np,p)
C Cholesky decomposition of a symmetric matrix a
C Ref: "Numerical Recipes", W.Press et al, Cambridge U.P.
C Written by J.M.Soler, May.2015

      use precision, only: dp
      use sys,       only: die

      implicit none
      integer, intent(in)   :: n        ! true dimension of matrix a
      integer, intent(in)   :: np       ! physical size of array a
      real(dp),intent(inout):: a(np,np) ! matrix to be decomposed
      real(dp),intent(out)  :: p(n)     ! auxiliary vector

      integer :: i,j
      real(dp):: a2

      do i=1,n
        a2 = sum(a(i,1:i-1)**2)
        if (a2>a(i,i)) call die('chempot/cholDcmp failed')
        p(i) = sqrt(a(i,i)-a2)        
        do j=i+1,n
          a2 = sum( a(i,1:i-1) * a(j,1:i-1) )
          a(j,i) = (a(i,j)-a2)/p(i)
        enddo
      enddo

      end subroutine cholDcmp


      subroutine cholLinSys(a,n,np,p,b,x)
C Solves linear system a*x=b for a symmetric matrix a
C (in Cholesky form, as output from output of cholDcmp)
C Ref: "Numerical Recipes", W.Press et al, Cambridge U.P.
C Written by J.M.Soler, May.2015
      implicit none
      integer, intent(in) :: n        ! true dimension of matrix a
      integer, intent(in) :: np       ! physical size of array a
      real(dp),intent(in) :: a(np,np) ! matrix decomposed by cholDcmp
      real(dp),intent(in) :: p(n)     ! auxiliary vector from cholDcmp
      real(dp),intent(in) :: b(n)     ! vector of linear system
      real(dp),intent(out):: x(n)     ! solution vector

      integer i,k
      real(dp) sum

      do i=1,n
        x(i) = ( b(i) - sum(a(i,1:i-1)*x(1:i-1)) ) / p(i)
      enddo
      do i=n,1,-1
        x(i) = ( x(i) - sum(x(i+1:n)*a(i+1:n,i)) ) / p(i)
      enddo

      end subroutine cholLinSys


      subroutine lanc1(opt,hbar,nhmax,numh,listhptr,listh,
     .                   maxnh,nbasis,ener,Node)
C *********************************************************************
C Routine to calculate the mimimum or maximum eigenvalues of
C a given sparse Hamiltonian, by (2nd order) the Lanczos Method.
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C integer opt                     : 1 = compute minimun eigenval.
C                                   2 = compute maximum eigenval.
C real*8 hbar(maxnh,nbasis)       : Sparse Hamiltonian
C integer nhmax                   : Lower dimension of hbar
C integer numh(nhmax)             : control vector of hbar
C integer listhptr(nhmax)         : control vector of hbar
C integer listh(maxnh)            : control vector of hbar
C integer maxnh                 : Maximum number of nonzero elements of
C                                   the hamiltonian
C integer nbasis                  : number of basis orbitals
C integer Node                    : local node number
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue
C **********************************************************************

      use precision, only : dp

      implicit none

      integer
     .  nhmax, nbasis, opt, Node, maxnh

      integer 
     .  listh(maxnh), numh(nbasis), listhptr(nbasis)

      real(dp)
     .  ener, hbar(nhmax,nbasis)

C  Internal variables ...

      integer
     .  itmax
      parameter (itmax=200)

      real(dp)
     .  tol
      parameter (tol=0.0001d0)

      integer
     .  i, ii, j, k

      real(dp)
     .  a0, b1, b2, c1, diff, ecomp, eivec(2), 
     .  mod, norm, randomg, wr(2)

      real(dp), dimension(:), pointer ::  hv0, u0, v0

      external  randomg

C ...

C  Allocate local memory
      nullify( hv0 )
      call re_alloc( hv0, 1, nbasis, name='hv0', routine='lanc1' )
      nullify( v0 )
      call re_alloc( v0, 1, nbasis, name='v0', routine='lanc1' )
      nullify( u0 )
      call re_alloc( u0, 1, nbasis, name='u0', routine='lanc1' )

C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.0d0
      do i=1,nbasis
        u0(i)=2.0d0*randomg(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      v0(1:nbasis)=0.0d0

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
           enddo
        enddo
        a0=0.0d0
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.0d0
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        hv0(1:nbasis)=0.0d0

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              hv0(ii)=hv0(ii)+hbar(i,j)*v0(j)
           enddo
        enddo
        b2=0.0d0
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.0d0
        do i=1,nbasis
           c1=c1+v0(i)*hv0(i)/norm
        enddo


C eigenvalues and eigenvectors ...

        wr(1) = 0.5d0*((a0+c1) 
     .          + sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
        wr(2) = 0.5d0*((a0+c1) 
     .          - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))

        eivec(1)=1/sqrt(1+((a0+b1-wr(opt))/(b2+c1-wr(opt)))**2)
        eivec(2)=-eivec(1)*(a0+b1-wr(opt))/(b2+c1-wr(opt))

        ener=wr(opt)
C ...

        norm=0.0d0
        diff=abs((wr(opt)-ecomp)/wr(opt))
        if (diff.gt.tol) then
          do i=1,nbasis 
             u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
             norm=norm+u0(i)**2
             v0(i)=0.0d0
          enddo
          do j=1,nbasis
            u0(j)=u0(j)/sqrt(norm)
          enddo
          ecomp=wr(opt)
        else 
          goto 20
        endif
      enddo
C .................


CC  Lanczos loop ................
C      do k=1,itmax
C        do j=1,nbasis
C           do i=1,numh(j)
C              ii=listh(listhptr(j)+i)
C              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
C           enddo
C        enddo
C        a0=0.
C        do i=1,nbasis
C           a0=a0+u0(i)*v0(i)
C        enddo
C        norm=0.
C        do i=1,nbasis
C           v0(i)=v0(i)-a0*u0(i)
C           norm=norm+v0(i)**2
C        enddo
C        b1=sqrt(norm)
C
C        do i=1,nbasis
C           hv0(i)=0.0
C        enddo
C
C        do j=1,nbasis
C           do i=1,numh(j)
C              ii=listh(listhptr(j)+i)
C              hv0(ii)=hv0(ii)+hbar(i,j)*v0(j)
C           enddo
C        enddo
C        b2=0.
C        do j=1,nbasis
C           b2=b2+u0(j)*hv0(j)/b1
C        enddo
C        c1=0.
C        do i=1,nbasis
C           c1=c1+v0(i)*hv0(i)/norm
C        enddo
C
C
CC eigenvalues and eigenvectors ...
C
C        wr(1) = 0.5d0*((a0+c1) 
C     .          + sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
C        wr(2) = 0.5d0*((a0+c1) 
C     .          - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
C
C        eivec(1)=1/sqrt(1+((a0+b1-wr(opt))/(b2+c1-wr(opt)))**2)
C        eivec(2)=-eivec(1)*(a0+b1-wr(opt))/(b2+c1-wr(opt))
C
C        ener=wr(opt)
CC ...
C
C 
C        norm=0.
C        diff=abs((wr(opt)-ecomp)/wr(opt))
C        if (diff.gt.tol) then
C          do i=1,nbasis 
C             u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
C             norm=norm+u0(i)**2
C             v0(i)=0.
C          enddo
C          do j=1,nbasis
C            u0(j)=u0(j)/sqrt(norm)
C          enddo
C          ecomp=wr(opt)
C        else 
C          goto 20
C        endif
C      enddo
CC .................

      if (Node.eq.0) then
        write(6,*) 'WARNING: lanc1 not converged after ',itmax,
     .           ' iterations'
      endif

20    continue

C  Deallocate local memory

      call de_alloc( hv0, name='hv0' )
      call de_alloc( v0, name='v0' )
      call de_alloc( u0, name='u0' )

      end   subroutine lanc1

      subroutine lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .                 ener,Node)
C *********************************************************************
C Routine to calculate the eigenvalue closest to cero for a given
C sparse Hamiltonian H, using the the Folded Spectrum Method
C (by (2nd order) the Lanczos Method).
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C real*8 hbar(maxnh,nbasis)         : Sparse Hamiltonian
C integer nhmax                   : Lower dimension of hbar
C integer numh(nhmax)             : control vector of hbar
C integer listhptr(nhmax)         : control vector of hbar
C integer listh(maxnh)            : control vector of hbar
C integer maxnh                   : Maximum number of nonzero elements of
C                                   the hamiltonian
C integer nbasis                  : number of basis orbitals
C integer Node                    : local node number
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue of H
C **********************************************************************

      use precision, only : dp

      implicit none

      integer
     .  nhmax, nbasis, Node, maxnh

      integer 
     .  listh(maxnh), numh(nbasis), listhptr(nbasis)

      real(dp)
     .  ener, hbar(nhmax,nbasis)

C  Internal variables ...

      integer
     .  itmax
      parameter (itmax=500)

      real(dp)
     .  tol
      parameter (tol=0.000001d0)

      integer
     .  i, ii, ij, j, k

      real(dp)
     .  a0, b1, b2, c1, diff, ecomp, eivec(2), 
     .  mod, norm, randomg, wr

      real(dp), dimension(:), pointer :: hv0, u0, v0, v00

      external randomg

C ...

C  Allocate local memory
      nullify( hv0 )
      call re_alloc( hv0, 1, nbasis, name='hv0', routine='lanc2' )
      nullify( u0 )
      call re_alloc( u0, 1, nbasis, name='u0', routine='lanc2' )
      nullify( v0 )
      call re_alloc( v0, 1, nbasis, name='v0', routine='lanc2' )
      nullify( v00 )
      call re_alloc( v00, 1, nbasis, name='v00', routine='lanc2' )

C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.0d0
      do i=1,nbasis
        u0(i)=2.0d0*randomg(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      v0(1:nbasis)=0.0d0
      v00(1:nbasis)=0.0d0

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v00(j)=v00(j)+hbar(i,j)*u0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v0(j)=v0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        a0=0.0d0
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.0d0
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        hv0(1:nbasis)=0.0d0
        v00(1:nbasis)=0.0d0

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v00(j)=v00(j)+hbar(i,j)*v0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              hv0(j)=hv0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        b2=0.0d0
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.0d0
        do i=1,nbasis
           c1=c1+v0(i)*hv0(i)/norm
        enddo

c  minimum eigenvalue ...
        wr = 0.5d0*((a0+c1) 
     .                 - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
c ...

c eigenvector ...
        eivec(1)=1/sqrt(1+((a0+b1-wr)/(b2+c1-wr))**2)
        eivec(2)=-eivec(1)*(a0+b1-wr)/(b2+c1-wr)
C ...

        norm=0.0d0

        do i=1,nbasis 
          u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
          norm=norm+u0(i)**2
          v0(i)=0.0d0
          v00(i)=0.0d0
        enddo
        do j=1,nbasis
          u0(j)=u0(j)/sqrt(norm)
        enddo

        diff=abs(wr-ecomp)
        ecomp=wr
        if (diff.lt.tol) goto 20
      enddo
C .................

      if (Node.eq.0) then
        write(6,*) 'WARNING: lanc2 not converged after ',itmax,
     .           ' iterations'
      endif

20    continue

C Calculate eigehvalue of H ...
      ener=0.0d0
      v0(1:nbasis)=0.0d0
      do i=1,nbasis
        do j=1,numh(i)
          ij=listh(listhptr(i)+j)
          v0(i)=v0(i)+hbar(j,i)*u0(ij)
        enddo
      enddo
      do i=1,nbasis
        ener=ener+u0(i)*v0(i)
      enddo
      
C  Deallocate local memory

      call de_alloc( hv0, name='hv0' )
      call de_alloc( u0, name='u0' )
      call de_alloc( v0, name='v0' )
      call de_alloc( v00, name='v00' )

      end subroutine lanc2

      end module m_chempot
