C ##################################################################
C ##            Transmission Eigenchannel Analysis                ##
C ##                           By                                 ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

      subroutine teigchana(nc,tt,neigch,Teig)

      implicit none


C ======================================================================
c     Note: Indices are starting at 0 in this routine
C ======================================================================
C INPUT:

      integer nc                !nc: no. states in contact region
      complex*16 tt(0:nc*nc-1)  !Eqv. to t^dagger.t
      integer ispin,nspin
C ======================================================================
C OUTPUT:


      real*8 Teig(0:nc-1)     ! The MAXEIG highest contrib. channels
c      complex*16 UR(0:nc*nc-1)  ! Right eigenvectors 
c      complex*16 UL(0:nc*nc-1)  ! Left eigenvectors 

      

C ======================================================================
c Helpers, tempos ...
c      complex*16 eig1(0:maxuo-1) ! eigenvalues 
c      integer eigord(0:maxuo-1)     ! final ordering      
      integer i,neigch

C --------------------------------------------------------------------
C ESSL and LAPACK:
c      logical DUMSELCT(maxuo)          ! used for ESSL-diagonalization 
c      complex*16  aux(3*maxuo)  !         - - 
C LAPACK diagononalization:
c      COMPLEX*16 WORKLAP(4*maxuo) ! used for LAPACK matrix diag.
c      REAL*8     RWORK(4*maxuo)   !         - - 
      INTEGER DIAGINFO             !         - -       

      complex*16, dimension (:), allocatable :: worklap,UR,UL, eig1
      real*8, allocatable ::  rwork(:)
      integer, allocatable ::  eigord(:)

c=======================================================================
C BEGIN

C *** Eigenchannel calculation:

      do i=0,nc-1
         tt(i+i*nc)=tt(i+i*nc)+dcmplx(1d0,0d0)
      end do

      allocate(worklap(4*nc))
      allocate(rwork(4*nc))
      allocate(ur(0:nc*nc-1))
      allocate(ul(0:nc*nc-1))
      allocate(eig1(0:nc-1))
      allocate(eigord(0:nc-1))



            
C **  Diagonalizing t^dagger.t 
C **  eigenvectors in UR, eigenvalues in eig1
c  ESSL BEGIN
c      call zgeev(1,tt,nc,eig1,UR,nc,DUMSELCT,nc,aux,3*maxuo)
c     make sure these eigenvectors are normalized - not certain with ESSL:
c      call normeigv(nc,UR) 
c NOTE: ESSL destroys tt in general !!
c ESSL END

C LAPACK BEGIN
      CALL ZGEEV('N','V',nc,tt,nc,eig1,UL,nc,UR,nc,WORKLAP,
     &     4*nc,RWORK,DIAGINFO)
      IF(DIAGINFO.LT. 0) THEN
         DIAGINFO= -DIAGINFO
         write(*,*) 'ERROR: DIAG. FAILED: element no. ',DIAGINFO, ' had
     &        illegal value. '
         STOP
      END IF
      IF(DIAGINFO.GT. 0) THEN
         write(*,*) 'ERROR: DIAG. FAILED: only the ',DIAGINFO,'+1:',nc,' 
     &        elements converged.'
      END IF
c LAPACK END

      do i=0,nc-1
         eig1(i)=eig1(i)-dcmplx(1d0,0d0)
      end do


c ==============================================================
C     find orderingvector for eigenchannels
c ==============================================================
c eigord1: After size of |T| 
      CALL eigsort(eig1,eigord,nc,nc) 
     
c ==============================================================

      do i=0,nc-1
         Teig(i)=CDABS(eig1(eigord(i)))
      end do
      
 1234 continue

      deallocate(worklap)
      deallocate(rwork)
      deallocate(ur)
      deallocate(ul)
      deallocate(eig1)
      deallocate(eigord)



C ====================================================================
      RETURN
      END
C ====================================================================
c
c
c
c
c
C ====================================================================
      SUBROUTINE eigsort(eig,eigord,N,MAXN)

c  sorting eigenvectors and eigenvalues so the MAXN eigenvectors with
c  largest eig comes first. 
C ====================================================================
      IMPLICIT NONE
      INTEGER N,eigord(0:N-1),maxN
      COMPLEX*16 eig(0:N-1), wrk(0:N-1), ctmp

      INTEGER i,j,k,P(0:N-1),itmp

      DO i=0,N-1 
         P(i)=i
         wrk(i)=eig(i)
         eigord(i)=0
      END DO
      
      DO 2000 j=0,N-2
         DO 2001 i=N-1,j+1,-1
            IF((CDABS(wrk(i))).LT.(CDABS(wrk(j)))) THEN
               ctmp=wrk(j)
               wrk(j)=wrk(i)
               wrk(i)=ctmp
               itmp=P(j)
               P(j)=P(i)
               P(i)=itmp
            END IF
 2001    CONTINUE
 2000 CONTINUE
      
C *** we want largest first:
c       DO i=0,N-1
c          eigord(i)=P(N-1-i)
c       END DO

C *** we want the MAXN largest first keeping the sequence among them:
      k=0
      DO i=0,N-1          
         DO j=0,MAXN-1
            if(P(N-1-j).EQ.i) then 
               eigord(k)=P(N-1-j)
               k=k+1
            end if
         END DO
      END DO
      if(k.NE.MAXN) write(*,*) 'sorting trouble k=',k



      RETURN
C ====================================================================
      END ! of eigsort
C ====================================================================
c
c
c
c
c
C ====================================================================
      SUBROUTINE normeigv(N,ev)

c  normalizing eigenvectors 
C ====================================================================
      IMPLICIT NONE
      INTEGER N,i,j
      COMPLEX*16 ev(0:N*N-1), norm
c     the eigenvector are the columns of ev_rs=ev(r+N*s)
c     FORTRAN is column-major i.e. column runs slowest
c     (major = slow runner - like in the army)      

      do j=0,N-1
         norm=dcmplx(0d0,0d0)
         do i=0,N-1
            norm=norm + dconjg(ev(i+j*N))*ev(i+j*N)
         end do
         norm = sqrt(norm) 
         
         do i=0,N-1
            ev(i+j*N)=ev(i+j*N)/norm
         end do
      end do                    ! eigenvector j

      RETURN
C ====================================================================
      END                       ! of normeig
C ====================================================================






