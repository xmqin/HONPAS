C FDN cell, kscell and displ added as dummy
      subroutine mkqgrid(joutfile,TRS,NA1,NA2,nq,q,wq,
     . cell,kscell,displ )
C FDN

      implicit none

c     INPUT
      integer joutfile          !output file unit
      logical TRS               ! Use time-reversal symmetry or not
      integer NA1,NA2           !no. repetitions of simple unitcell
c                               !in A1,A2 directions

c     OUTPUT
      integer nq                !No. q-points <= NA1*NA2
      real*8, pointer:: q(:,:),wq(:)


c     Helpers
      logical PRINTALOT
      parameter(PRINTALOT=.false.)

      real*8 cell(3,3),displ(3),cutoff
      real*8 ,allocatable ::     ktmp(:,:),wqtmp(:) 


      integer kscell(3,3),i,j,iq
      
c      BEGIN

C FDN
C      do j=1,3
C         do i=1,3
C            cell(i,j)=0d0
C            kscell(i,j)=0
C         end do
C         cell(j,j)=1d0
C         displ(j)=0d0
C      end do
C
C
C      cutoff=0d0
C      nq=NA1*NA2                !initial value 
C      kscell(1,1)=NA1
C      kscell(2,2)=NA2
C      kscell(3,3)=1


      do i=1,3
        kscell(i,3)=0
        kscell(3,i)=0
      end do
      kscell(3,3)=1
      displ(3)=0d0
      cutoff=0d0

      nq = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )



      allocate(ktmp(3,nq))
      allocate(wqtmp(nq))


C      allocate(ktmp(3,NA1*NA2))    
C      allocate(wqtmp(NA1*NA2))    

C      do i=1,NA1*NA2
C         ktmp(1,i)=0d0
C         ktmp(2,i)=0d0
C         ktmp(3,i)=0d0
C         wqtmp(i)=0d0
C      end do
    
      ktmp=0d0
      wqtmp=0d0

C                                                                                                                                                                                 
C FDN                                                                                                                                                                             
C 

      

      if(TRS) then
c         write(joutfile,*) 'cell, kscell, displ,
c     .        cutoff, nq, ktmp, wqtmp'
c         write(joutfile,*) cell, kscell, displ,
c     .        cutoff, nq, ktmp, wqtmp

         call kgrid( cell, kscell, displ,
     .        cutoff, nq, ktmp, wqtmp )
         
      else                      ! do not use time-reversal symmetry

      iq=0
      do i=1,NA1
         do j=1,NA2
            iq=iq+1
            ktmp(1,iq)=1d0*(i-1)/NA1
            ktmp(2,iq)=1d0*(j-1)/NA2
            wqtmp(iq)=1.0/(NA1*NA2)
         end do
      end do
      nq=iq
      end if                    !TRS

C FDN q(2,nq) ==> q(3,nq)
      allocate(q(3,nq))
C FDN
      allocate(wq(nq))


      if(PRINTALOT) write(joutfile,*) 'q-points: ',nq
      do i=1,nq
         q(1,i)=ktmp(1,i)
         q(2,i)=ktmp(2,i)
C FDN
         q(3,i)=ktmp(3,i)
C FDN
         wq(i) =wqtmp(i)
C FDN
         if(PRINTALOT)  write(joutfile,'(3F10.7,2x,g10.5)') 
     &        q(1,i),q(2,i),q(3,i),wq(i)
C FDN
      end do

      deallocate(ktmp)
      deallocate(wqtmp)
c=============================================================
      return
      end
c=============================================================



      subroutine kgrid( cell, kscell, displ,
     .     cutoff, nk, points, weight )

c **********************************************************************
c Finds Monkhost-Pack k-point coordinates and weights.
c This version assumes no symmetry except time reversal, i.e.
c inversion in reciprocal space.
c Refs: H.J.Monkhorst and J.D.Pack, Phys Rev B 13, 5188 (1976)
c       J.Moreno and J.M.Soler, Phys Rev B 45, 13891 (1992)
c Written by J.M.Soler. July 1997.
c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c integer kscell(3,3): Supercell reciprocal of k-grid unit cell
c                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
c real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
c                      origin(ix) = sum_j gridk(ix,j)*displ(j)
c real*8  cutoff     : Minimum k-grid cutoff required
c                      Not used unless det(kscell)=0
c integer nk         : Dimension of arrays points and weight
c ***************** OUTPUT *********************************************
c integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
c                        Only if det(kscell)=0 on input
c real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
c                        Only if det(kscell)=0 on input
c real*8  cutoff       : Actual k-grid cutoff
c integer nk           : Actual number of irreducible k-points
c real*8  points(3,nk) : K-point cartesian coordinates
c                        Only if input_nk .ge. output_nk
c real*8  weight(nk)   : K-point weights 
c                        Only if input_nk .ge. output_nk
c ***************** UNITS **********************************************
c cutoff must be in the same units as cell
c points returned in units reciprocal of cell
c ***************** BEHAVIOUR ******************************************
c - If det(kscell).ne.0, input cutoff is not used
c - If det(kscell).eq.0, kscell and displ are generated according to
c   input cutoff
c - If det(kscell).eq.0 .AND. cutoff.le.0.d0, they are readed
c   from the input fdf data file.
c   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
c   If both are present, kgrid_Monkhorst_Pack has priority. If none is
c   present, the cutoff defect is cero, producing only the gamma point.
c   Examples of fdf data specifications:
c     kgrid_cutoff  50. Bohr
c     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
c     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
c     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
c     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
c     %endblock kgrid_Monkhorst_Pack
c - If input_nk < output_nk, points and weight are not generated and
c   a warning is printed before return
c **********************************************************************
      implicit          none
      integer           kscell(3,3), nk
      double precision  cell(3,3), cutoff, displ(3), 
     .                  points(3,*), weight(*)
      external          idiag, reclat
c ----------------------------------------------------------------------

c Internal variables
      integer           i, i1, i2, i3, igmax(3), igmin(3),
     .                  ir, ik, iu, ix, j,
     .                  kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3),
     .                  ng(3), ni, nkmax, nkr(3), nktot
      double precision  d(3), defcut, dkg(3), dkx(3), dscell(3,3),
     .                  gridk(3,3), gscell(3,3), huge, pi,  
     .                  scell(3,3), tiny, vmod, w1, wtot
      parameter (defcut = 0.d0)
      parameter (huge   = 1.d30)
      parameter (tiny   = 1.d-12)

c Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )

c Look for kscell or cutoff in input fdf file
c      if ( nktot.eq.0 .and. cutoff.lt.tiny ) then
c        if ( fdf_block('kgrid_Monkhorst_Pack',iu) ) then
c          do i = 1,3
c            read(iu,*) (kscell(j,i),j=1,3), displ(i)
c          enddo
c          nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
c     .                 kscell(2,1) * kscell(3,2) * kscell(1,3) +
c     .                 kscell(3,1) * kscell(1,2) * kscell(2,3) -
c     .                 kscell(1,1) * kscell(3,2) * kscell(2,3) -
c     .                 kscell(2,1) * kscell(1,2) * kscell(3,3) -
c     .                 kscell(3,1) * kscell(2,2) * kscell(1,3) )
c        else
c         The second argument is the default value
c          cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
c        endif
c      endif

c Find kscell from required cutoff
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

c Find k-grid supercell
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) +
     .                  cell(ix,2) * kscell(2,i) +
     .                  cell(ix,3) * kscell(3,i)
        enddo
      enddo

c Find actual cutoff
      cutoff = huge
      do i = 1,3
        vmod = sqrt( scell(1,i)**2 + scell(2,i)**2 + scell(3,i)**2 )
        cutoff = min( cutoff, vmod/2.d0 )
      enddo

c Find equivalent diagonal supercell
      call idiag( 3, kscell, kdsc, ml, mr, maux )
      do i = 1,3
        do ix = 1,3
          dscell(ix,i) = scell(ix,1) * mr(1,i) +
     .                   scell(ix,2) * mr(2,i) +
     .                   scell(ix,3) * mr(3,i)
        enddo
      enddo

c Find k-grid unit vectors
C FDN Changed to 1 --> 2Pi
      call reclat( dscell, gridk, 1 )
C FDN

c Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) +
     .            gscell(ix,2) * displ(2) +
     .            gscell(ix,3) * displ(3)
      enddo

c Find grid origin in gridk coordinates
      pi = 4.d0 * atan(1.d0)
      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) +
     .             dkx(2) * dscell(2,i) +
     .             dkx(3) * dscell(3,i) ) / (2*pi)
      enddo

c Some printout for debugging
*     write(6,'(/,a,/,(3f12.6,i6,f12.6))') 'kgrid: gridK,ng,dg =',
*    .  ((gridk(ix,i),ix=1,3),kdsc(i,i),dkg(i),i=1,3)

c Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
        igmin(j) = -( (ng(j)-1) / 2)
        igmax(j) = ng(j) / 2
      enddo

c Find number of points with time-reversal (inversion) symmetry,
c after reflection on each alternative plane
      do j = 1,3
        ni = ng(j)
        if (abs(dkg(j)) .lt. tiny) then
          ni = ng(j)/2 + 1
        elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
          ni = (ng(j)-1)/2 + 1
        endif
        nkr(j) = ni * nktot / kdsc(j,j)
      enddo

c Select reflection plane
      ir = 3
      if (nkr(2) .lt. nkr(ir)) ir = 2
      if (nkr(1) .lt. nkr(ir)) ir = 1
      igmin(ir) = 0
      if (abs(dkg(ir)-0.5d0) .lt. tiny)
     .  igmax(ir) = (ng(ir)-1)/2

c Find k points and weights
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
            points(ix,nk) = gridk(ix,1)*d(1) + 
     .                      gridk(ix,2)*d(2) +
     .                      gridk(ix,3)*d(3)
          enddo
          if ( abs(d(ir))              .lt. tiny .or.
     .         abs(d(ir)-0.5d0*ng(ir)) .lt. tiny) then
            weight(nk) = w1
          else
            weight(nk) = 2.d0 * w1
          endif
        enddo
        enddo
        enddo
      else
        write(6,'(/,a,i6,/)')
     .   'kgrid: dimension nk too small. Must be at least', nk
      endif

c A couple of tests for debugging
      if (nk .le. nkmax) then
        if (nk .ne. nkr(ir))
     .     write(6,*) 'kgrid: ERROR: nk, nkr(ir) =', nk, nkr(ir)
        wtot = 0.d0
        do ik = 1,nk
          wtot = wtot + weight(ik)
        enddo
        if (abs(wtot-1.d0) .gt. nk*tiny)
     .    write(6,*) 'kgrid: ERROR: wtot =', wtot
      endif
      return
      end


C $Id: mkqgrid.f,v 1.2 2002/12/04 12:32:20 mbr Exp $

      SUBROUTINE RECLAT (A,B,IOPT)

C  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
C  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3,3)
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
      END

C $Id: mkqgrid.f,v 1.2 2002/12/04 12:32:20 mbr Exp $

      SUBROUTINE IDIAG (NN,MOLD,MNEW,MLEFT,MRIGHT,MAUX)

C GIVEN A SQUARE INTEGER MATRIX MOLD, FINDS A DIAGONAL MATRIX MNEW,
C AND TWO MATRICES MLEFT AND MRIGHT OF DETERMINANT ONE, SUCH THAT
C MNEW = MLEFT * MOLD * MRIGHT
C Written by J.Moreno and J.M.Soler

      INTEGER MOLD(NN,NN),MNEW(NN,NN),MLEFT(NN,NN),MRIGHT(NN,NN),
     .        MAUX(NN,NN,2)
      PARAMETER (NITER=50,IBIG=9999999)
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
      END



      SUBROUTINE IMXM (MA,MB,N1,N2,N3,MC,MAUX)

C MULTIPLIES TWO INTEGER MATRICES. ARGUMENTS:
C   MA(N1,N2),MB(N2,N3) : INPUT MATRICES OF DIMENSIONS AS INDICATED
C   MC(N1,N3) : OUTPUT PRODUCT MATRIX. MC MAY BE THE SAME AS MA OR MB
C              (OR BOTH) FOR 'IN PLACE' MULTIPLICATION
C   MAUX : AUXILIARY ARRAY OF MINIMUM SIZE N1*N3. IF MC IS DIFFERENT
C          FROM BOTH MA AND MB, YOU CAN MAKE MAUX=MC
C WRITTEN BY JOSE SOLER. 22/5/90

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
      END


