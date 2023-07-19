module m_tbt_kpts


implicit none

public :: get_kp_on_node, reclat
private

contains

subroutine get_kp_on_node(fname,nkpar,kpar,wkpar,kscellfdf,kdisplfdf)

use parallel, only : Node, Nodes, IONode
use sys, only : die
use m_tbt_options, only : RNode

#ifdef MPI
use mpi_siesta
use m_mpi_utils
#endif MPI

implicit    none  
character*75 fname
integer :: funit,i,j, nkpar
integer, dimension(3,3) ::  kscell, kscellfdf 
logical :: fexist
real*8 :: cutoff
real*8, dimension(3) ::  kdispl, kdisplfdf
real*8, dimension(3,3) ::  cell
real*8, dimension(:), allocatable :: wktmp
real*8, dimension(:,:), allocatable ::ktmp
real*8, pointer :: kpar(:,:),wkpar(:)
! Parallel
integer, pointer :: nkOnNodes(:),ik
real*8, pointer :: kOnNodes(:,:,:),wkOnNodes(:,:)


! To read from file ...
integer :: nua, nuotot, notot, nspin, maxnh
integer, allocatable, dimension (:) :: isa
double precision, allocatable, dimension (:,:) :: xa

#ifdef MPI
      integer :: MPIerror
#endif MPI


      if (IOnode) then
        inquire(file=fname,exist=fexist)
        if (.not. fexist) then

         write(*,*) 'File :',fname,'not found in TSiokp'
         call die('Stopping code !')

        else
  
         call io_assign(funit)
         open(funit,file=fname,form='unformatted')
       
         read(funit) nua, nuotot, notot, nspin, maxnh
         allocate(xa(3,nua))
         allocate(isa(nua)) 
         read(funit) xa
         read(funit) isa
         read(funit) cell

         deallocate(xa,isa)

         kscell=kscellfdf
         kdispl=kdisplfdf

         do i=1,3
           kscell(i,3)=0
           kscell(3,i)=0
         end do
         kscell(3,3)=1
         kdispl(3)=0d0
         cutoff=0d0

         nkpar = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) + &
                  kscell(2,1) * kscell(3,2) * kscell(1,3) +  &
                  kscell(3,1) * kscell(1,2) * kscell(2,3) -  &
                  kscell(1,1) * kscell(3,2) * kscell(2,3) -  &
                  kscell(2,1) * kscell(1,2) * kscell(3,3) -  &
                  kscell(3,1) * kscell(2,2) * kscell(1,3) )        
         
        end if 

      end if ! IOnode 
#ifdef MPI
      call broadcast(nkpar)
#endif MPI
      allocate(ktmp(3,nkpar))
      allocate(wktmp(nkpar))

#ifdef MPI       
      call broadcast(cell(:,:))
      call broadcast(kscell(:,:))
      call broadcast(kdispl(:))
      call broadcast(cutoff)
#endif MPI

      call kgrid( cell, kscell, kdispl, &
             cutoff, nkpar, ktmp, wktmp )    

      nullify(nkOnNodes,kOnNodes,wkOnNodes)
      if (Nodes > 1) then
         call dist_kpts_on_nodes(nkpar,ktmp,wktmp,nkOnNodes,kOnNodes,&
           wkOnNodes)
      end if


      if ( Nodes > 1 ) then
        nkpar=nkOnNodes(Node)
      end if

      nullify(kpar,wkpar)
      allocate(kpar(3,nkpar))
      allocate(wkpar(nkpar))
           
      if ( Nodes > 1 ) then
        kpar(1:3,1:nkpar)=kOnNodes(1:3,1:nkpar,Node)
        wkpar(1:nkpar)=wkOnNodes(1:nkpar,Node)
      else 
        kpar(:,1:nkpar)=ktmp(:,1:nkpar)
        wkpar(1:nkpar)=wktmp(1:nkpar)
      end if ! Nodes > 1

 
      end subroutine get_kp_on_node


      subroutine dist_kpts_on_nodes(nkpar,ktmp,wktmp,nkOnNodes,kOnNodes,&
                wkOnNodes)

 
      use parallel, only : Node, Nodes, IONode
      use sys, only : die
      use m_tbt_options, only : RNode      


#ifdef MPI
use mpi_siesta
use m_mpi_utils
#endif MPI


      integer, intent(in) :: nkpar
      real*8, dimension(nkpar), intent(in) :: wktmp
      real*8, dimension(3,nkpar), intent(in) ::ktmp
      integer, pointer :: nkOnNodes(:)
      real*8, pointer :: kOnNodes(:,:,:),wkOnNodes(:,:)

! Auxs ...
      integer :: iNode,MPIerror,ik,nkmod,nkmin,nnodesmin,nkparSum, &
                 nkmax
      integer :: ij


      allocate(nkOnNodes(0:Nodes-1))

      if ( Nodes <= nkpar ) then
         nkmod=mod(nkpar,Nodes)
         nkmin=(nkpar-nkmod)/Nodes
         nkOnNodes(0:Nodes-1)=nkmin
         nkmax=nkmin
         if ( nkmod > 0 ) then
           nkmax=nkmin+1
           do ik = 1,nkmod
              nkOnNodes(ik)=nkOnNodes(ik)+1
           enddo
         end if 
      else
        if (IOnode) then
           write(*,*) 'WARNING !!!!'
           write(*,*) 'Number of Nodes greater than number of k-Points !!'
           write(*,*) 'Some nodes will be innactive in k-point parallelization !!'
        end if 
        nkmax=1
        nkOnNodes(0)=0
        nkOnNodes(1:nkpar)=1
        nkOnNodes(nkpar+1:Nodes-1)=0
      end if ! Nodes < nkpar

! Verify if sum of Number of kpoints on nodes matches total number 
      nkparSum=0
      do iNode = 0,Nodes-1
         nkparSum=nkparSum+nkOnNodes(iNode)
      enddo 
      if ( nkparSum /= nkpar ) then
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI
         if (IOnode) then
            write(*,*) 'Sum of K Points over nodes does not match Total number !!'
            write(*,*) 'Sum = ',nkparSum,'Total = ',nkpar
            call die('Stopping code !!')
         end if
      end if ! nkparSum /= nkpar

! Distribute k Pts and weights Over Nodes
      allocate(kOnNodes(3,nkmax,0:Nodes-1))
      allocate(wkOnNodes(nkmax,0:Nodes-1))
      nkparSum=0
      do iNode = 0,Nodes-1
         do ik = 0,nkOnNodes(iNode)-1
            nkparSum=nkparSum+1
            kOnNodes(1:3,ik+1,iNode)=ktmp(1:3,nkparSum)
            wkOnNodes(ik+1,iNode)=wktmp(nkparSum)
         enddo
      enddo

      if (IONode) then
         write(*,*)
         write(*,'(a)') repeat('*',72)
         write(*,*) 'K-Points On Nodes:'
         write(*,*) 
         write(*,'(a6,a14,a13)') ' Node ',' Numb. K Pts. ',' Tot. K Pts. '
         do iNode = 0,Nodes-1
            write(*,'(I3,I9,I13)') iNode,nkOnNodes(iNode),nkpar
         enddo
         write(*,'(a)') repeat('-',72)
         write(*,*) 'K-Points Coordinates and Weights (Node,kx,ky,kz,wk),:'
         write(*,*)
         do iNode = 0,Nodes-1
            do ik = 0,nkOnNodes(iNode)-1
               write(*,'(I9,4F12.5)') iNode,(kOnNodes(ij,ik+1,iNode),ij=1,3), &
                                      wkOnNodes(ik+1,iNode)
            enddo
         enddo
         write(*,'(a)') repeat('*',72)
         write(*,*)
      end if


      end subroutine dist_kpts_on_nodes



      subroutine kgrid( cell, kscell, displ, &
          cutoff, nk, points, weight )

! **********************************************************************
! Finds Monkhost-Pack k-point coordinates and weights.
! This version assumes no symmetry except time reversal, i.e.
! inversion in reciprocal space.
! Refs: H.J.Monkhorst and J.D.Pack, Phys Rev B 13, 5188 (1976)
!       J.Moreno and J.M.Soler, Phys Rev B 45, 13891 (1992)
! Written by J.M.Soler. July 1997.
! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! integer kscell(3,3): Supercell reciprocal of k-grid unit cell
!                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
! real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
!                      origin(ix) = sum_j gridk(ix,j)*displ(j)
! real*8  cutoff     : Minimum k-grid cutoff required
!                      Not used unless det(kscell)=0
! integer nk         : Dimension of arrays points and weight
! ***************** OUTPUT *********************************************
! integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
!                        Only if det(kscell)=0 on input
! real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
!                        Only if det(kscell)=0 on input
! real*8  cutoff       : Actual k-grid cutoff
! integer nk           : Actual number of irreducible k-points
! real*8  points(3,nk) : K-point cartesian coordinates
!                        Only if input_nk .ge. output_nk
! real*8  weight(nk)   : K-point weights 
!                        Only if input_nk .ge. output_nk
! ***************** UNITS **********************************************
! cutoff must be in the same units as cell
! points returned in units reciprocal of cell
! ***************** BEHAVIOUR ******************************************
! - If det(kscell).ne.0, input cutoff is not used
! - If det(kscell).eq.0, kscell and displ are generated according to
!   input cutoff
! - If det(kscell).eq.0 .AND. cutoff.le.0.d0, they are readed
!   from the input fdf data file.
!   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
!   If both are present, kgrid_Monkhorst_Pack has priority. If none is
!   present, the cutoff defect is cero, producing only the gamma point.
!   Examples of fdf data specifications:
!     kgrid_cutoff  50. Bohr
!     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
!     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
!     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
!     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
!     %endblock kgrid_Monkhorst_Pack
! - If input_nk < output_nk, points and weight are not generated and
!   a warning is printed before return
! **********************************************************************
      implicit          none
      integer           kscell(3,3), nk
      double precision  cell(3,3), cutoff, displ(3), &
                        points(3,*), weight(*)
!      external          idiag, reclat
! ----------------------------------------------------------------------

! Internal variables
      integer           i, i1, i2, i3, igmax(3), igmin(3), &
                       ir, ik, iu, ix, j, &
                       kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3), &
                       ng(3), ni, nkmax, nkr(3), nktot
      double precision  d(3), defcut, dkg(3), dkx(3), dscell(3,3), &
                       gridk(3,3), gscell(3,3), huge, pi, &
                       scell(3,3), tiny, vmod, w1, wtot
      parameter (defcut = 0.d0)
      parameter (huge   = 1.d30)
      parameter (tiny   = 1.d-12)

! Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) + &
                  kscell(2,1) * kscell(3,2) * kscell(1,3) + &
                  kscell(3,1) * kscell(1,2) * kscell(2,3) - &
                  kscell(1,1) * kscell(3,2) * kscell(2,3) - &
                  kscell(2,1) * kscell(1,2) * kscell(3,3) - &
                  kscell(3,1) * kscell(2,2) * kscell(1,3) )



! Find kscell from required cutoff
      if (nktot .eq. 0) then
        nktot = 1
        do j = 1,3
          do i = 1,3
            kscell(j,i) = 0
          enddo
          vmod = sqrt( cell(1,j)**2 + cell(2,j)**2 + cell(3,j)**2 )
          kscell(j,j) = int(cutoff/(vmod/2.d0)) + 1
          if (mod(kscell(j,j),2) .eq. 0) displ(j) = 0.5d0
          nktot = nktot * kscell(j,j)
        enddo
      endif

! Find k-grid supercell
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) + &
                       cell(ix,2) * kscell(2,i) + &
                       cell(ix,3) * kscell(3,i)
        enddo
      enddo

! Find actual cutoff
      cutoff = huge
      do i = 1,3
        vmod = sqrt( scell(1,i)**2 + scell(2,i)**2 + scell(3,i)**2 )
        cutoff = min( cutoff, vmod/2.d0 )
      enddo

! Find equivalent diagonal supercell
      call idiag( 3, kscell, kdsc, ml, mr, maux )
      do i = 1,3
        do ix = 1,3
          dscell(ix,i) = scell(ix,1) * mr(1,i) + &
                        scell(ix,2) * mr(2,i) + &
                        scell(ix,3) * mr(3,i)
        enddo
      enddo

! Find k-grid unit vectors
! FDN Changed to 1 --> 2Pi
      call reclat( dscell, gridk, 1 )
! FDN

! Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) + &
                 gscell(ix,2) * displ(2) + &
                 gscell(ix,3) * displ(3)
      enddo

! Find grid origin in gridk coordinates
      pi = 4.d0 * atan(1.d0)
      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) + &
                  dkx(2) * dscell(2,i) + &
                  dkx(3) * dscell(3,i) ) / (2*pi)
      enddo

! Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
        igmin(j) = -( (ng(j)-1) / 2)
        igmax(j) = ng(j) / 2
      enddo

! Find number of points with time-reversal (inversion) symmetry,
! after reflection on each alternative plane
      do j = 1,3
        ni = ng(j)
        if (abs(dkg(j)) .lt. tiny) then
          ni = ng(j)/2 + 1
        elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
          ni = (ng(j)-1)/2 + 1
        endif
        nkr(j) = ni * nktot / kdsc(j,j)
      enddo

! Select reflection plane
      ir = 3
      if (nkr(2) .lt. nkr(ir)) ir = 2
      if (nkr(1) .lt. nkr(ir)) ir = 1
      igmin(ir) = 0
      if (abs(dkg(ir)-0.5d0) .lt. tiny) &
       igmax(ir) = (ng(ir)-1)/2

! Find k points and weights
      nkmax = nk
      nk = nkr(ir)
      if (nk .le. nkmax) then
        w1 = 1.d0 / nktot
        nk = 0
        do i3 = igmin(3),igmax(3)
        do i2 = igmin(2),igmax(2)
        do i1 = igmin(1),igmax(1)
          nk = nk + 1
          d(1) = i1 + dkg(1)
          d(2) = i2 + dkg(2)
          d(3) = i3 + dkg(3)
          if (d(1) .gt. 0.5d0*ng(1)+tiny) d(1) = d(1) - ng(1)
          if (d(2) .gt. 0.5d0*ng(2)+tiny) d(2) = d(2) - ng(2)
          if (d(3) .gt. 0.5d0*ng(3)+tiny) d(3) = d(3) - ng(3)
          do ix = 1,3
            points(ix,nk) = gridk(ix,1)*d(1) + &
                           gridk(ix,2)*d(2) + &
                           gridk(ix,3)*d(3)
          enddo
          if ( abs(d(ir))              .lt. tiny .or. &
              abs(d(ir)-0.5d0*ng(ir)) .lt. tiny) then
            weight(nk) = w1
          else
            weight(nk) = 2.d0 * w1
          endif
        enddo
        enddo
        enddo
      else
        write(6,'(/,a,i6,/)') &
        'kgrid: dimension nk too small. Must be at least', nk
      endif

! A couple of tests for debugging
      if (nk .le. nkmax) then
        if (nk .ne. nkr(ir)) &
          write(6,*) 'kgrid: ERROR: nk, nkr(ir) =', nk, nkr(ir)
        wtot = 0.d0
        do ik = 1,nk
          wtot = wtot + weight(ik)
        enddo
        if (abs(wtot-1.d0) .gt. nk*tiny) &
         write(6,*) 'kgrid: ERROR: wtot =', wtot
      endif
      return
      end subroutine kgrid


! $Id: mkqgrid.f,v 1.2 2002/12/04 12:32:20 mbr Exp $

      SUBROUTINE RECLAT (A,B,IOPT)

!  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
!  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3,3)
      DOUBLE PRECISION PI, C,CI
      INTEGER IOPT,I

      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO 20 I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
  20  CONTINUE

      return
      END SUBROUTINE RECLAT

! $Id: mkqgrid.f,v 1.2 2002/12/04 12:32:20 mbr Exp $

      SUBROUTINE IDIAG (NN,MOLD,MNEW,MLEFT,MRIGHT,MAUX)

! GIVEN A SQUARE INTEGER MATRIX MOLD, FINDS A DIAGONAL MATRIX MNEW,
! AND TWO MATRICES MLEFT AND MRIGHT OF DETERMINANT ONE, SUCH THAT
! MNEW = MLEFT * MOLD * MRIGHT
! Written by J.Moreno and J.M.Soler

      INTEGER NN,IMIN,JMIN,I,J,N,ITER,MMIN
      INTEGER MOLD(NN,NN),MNEW(NN,NN),MLEFT(NN,NN),MRIGHT(NN,NN), &
             MAUX(NN,NN,2)

      INTEGER, PARAMETER :: NITER=50,IBIG=9999999
      IMIN=-10000
      JMIN=-10000
      DO 20 J=1,NN
        DO 10 I=1,NN
          MNEW(I,J)=MOLD(I,J)
          MLEFT(I,J)=0
          MRIGHT(I,J)=0
          MAUX(I,J,1)=0
   10   CONTINUE
        MLEFT(J,J)=1
        MRIGHT(J,J)=1
        MAUX(J,J,1)=1
   20 CONTINUE
      DO 60 N=NN,2,-1
        DO 50 ITER=1,NITER
          MMIN=IBIG
          DO 30 J=1,N
          DO 30 I=1,N
            IF ((I.NE.N.AND.J.NE.N) .OR. (I.EQ.N.AND.J.EQ.N)) GOTO 30
            IF (MNEW(I,J).EQ.0 .OR. ABS(MNEW(I,J)).GE.MMIN) GOTO 30
               IMIN=I
               JMIN=J
               MMIN=ABS(MNEW(I,J))
   30     CONTINUE
          IF (MMIN.EQ.IBIG) GOTO 60
          I=MIN(IMIN,JMIN)
          MAUX(I,I,1)=0
          MAUX(N,N,1)=0
          MAUX(I,N,1)=SIGN(1,MNEW(IMIN,JMIN))
          MAUX(N,I,1)=SIGN(1,MNEW(IMIN,JMIN))
          IF (IMIN.LT.JMIN) THEN
            CALL IMXM (MAUX,MNEW,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MAUX,MLEFT,NN,NN,NN,MLEFT,MAUX(1,1,2))
          ELSE
            CALL IMXM (MNEW,MAUX,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MRIGHT,MAUX,NN,NN,NN,MRIGHT,MAUX(1,1,2))
          ENDIF
          MAUX(I,N,1)=0
          MAUX(N,I,1)=0
          MAUX(I,I,1)=1
          MAUX(N,N,1)=1
          DO 40 I=1,N-1
            MAUX(I,N,1)=-(MNEW(I,N)/MNEW(N,N))
            CALL IMXM (MAUX,MNEW,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MAUX,MLEFT,NN,NN,NN,MLEFT,MAUX(1,1,2))
            MAUX(I,N,1)=0
            MAUX(N,I,1)=-(MNEW(N,I)/MNEW(N,N))
            CALL IMXM (MNEW,MAUX,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MRIGHT,MAUX,NN,NN,NN,MRIGHT,MAUX(1,1,2))
            MAUX(N,I,1)=0
   40     CONTINUE
   50   CONTINUE
          WRITE(6,*)'IDIAG: ERROR. ITERATION HAS NOT CONVERGED. N=',N
          STOP 'IDIAG: ERROR. ITERATION HAS NOT CONVERGED.'
   60 CONTINUE
      END SUBROUTINE IDIAG



      SUBROUTINE IMXM (MA,MB,N1,N2,N3,MC,MAUX)

! MULTIPLIES TWO INTEGER MATRICES. ARGUMENTS:
!   MA(N1,N2),MB(N2,N3) : INPUT MATRICES OF DIMENSIONS AS INDICATED
!   MC(N1,N3) : OUTPUT PRODUCT MATRIX. M! MAY BE THE SAME AS MA OR MB
!              (OR BOTH) FOR 'IN PLACE' MULTIPLICATION
!   MAUX : AUXILIARY ARRAY OF MINIMUM SIZE N1*N3. IF M! IS DIFFERENT
!          FROM BOTH MA AND MB, YOU CAN MAKE MAUX=MC
! WRITTEN BY JOSE SOLER. 22/5/90

      INTEGER I1,I2,I3,N1,N2,N3,I,J
      INTEGER MA(N1,N2),MB(N2,N3),MC(N1,N3),MAUX(N1,N3)
      DO 20 I3=1,N3
      DO 20 I1=1,N1
        MAUX(I1,I3)=0
        DO 10 I2=1,N2
          MAUX(I1,I3)=MAUX(I1,I3)+MA(I1,I2)*MB(I2,I3)
   10   CONTINUE
   20 CONTINUE
      DO 30 I=N3,1,-1
      DO 30 J=N1,1,-1
        MC(J,I)=MAUX(J,I)
   30 CONTINUE
      return
      END SUBROUTINE IMXM






end module m_tbt_kpts
