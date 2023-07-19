C ##################################################################
C ##         The Tight-Binding Transmission Calculator            ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

      
      subroutine transmission(USEBULK,joutfile,no,
     &     Hbig,Sbig,
     &  nc,NGL,HLL,NGR,HRR,ZEnergy,GF,GFRGF,tottrans,tt)
     
      implicit none


C ======================================================================
c     Note: Indices are starting at 0 in this routine
C ======================================================================
C INPUT:
      logical USEBULK           ! if true self-energy only is input else
c                                 z*S-H-Sigma for bulk is in sfe

      integer joutfile          !outfile unit
      integer no                !no : total no. states
      integer nc                !no. : states in contact region
      complex*16 Hbig(0:no*no-1),Sbig(0:no*no-1)
      complex*16 ZEnergy 
      integer NGL,NGR
      complex*16 HLL(NGL,NGL)   ! Left electrode GF
      complex*16 HRR(NGR,NGR)   ! Right electrode GF
      
C ======================================================================
C OUTPUT:

      complex*16, dimension(:), allocatable ::
     . SigmaL,SigmaR
  
       integer,dimension(:), allocatable :: ipvt


      complex*16  tt(0:nc*nc-1)
      complex*16 GF(0:nc*nc-1)  ! 1/(Sc*Energy - Hc)
      complex*16 GFRGF(0:nc*nc-1)  ! GF^dagger.GammaR.GF

      real*8 tottrans

C ======================================================================
c Helpers, tempos ...
      complex*16, dimension(:), allocatable ::
     . U,iU,Hc,tiU 

      real*8 Energy !Real part of ZEnergy

      complex*16 a,b
      
      complex*16 csum
      integer nu
      integer i,j,ii,jj,il,jl,ir,jr

C --------------------------------------------------------------------
c      INTEGER ipvt(no)      ! pivoting vector for matrix inv.
c----------------------------------------------------------------------


c=======================================================================
C BEGIN
c=======================================================================






      if(nc .NE.  no - (NGL+NGR)) then
         write(joutfile,*) 'ERROR : transmission : nc = ',nc,'
     &        , expected:', no - (NGL+NGR)
         stop 'ERROR : transmission : nc not OK'
      end if


      Energy = DREAL(ZEnergy)

      allocate(ipvt(no))

c
c     SigmaL:
c

      allocate(U(0:no*no-1))
      allocate(iU(0:no*no-1))


      if(USEBULK) then

         do j = 0, NGL-1
            jj=j+1
            do i = 0, NGL-1
               ii=i+1
               U(i + j*NGL) = HLL(ii,jj)
            end do              !i
         end do
      else

         do j = 0, NGL-1
            jj=j+1
            do i = 0, NGL-1
               ii=i+1
               U(i + j*NGL) = Energy*Sbig(i+j*no) - Hbig(i+j*no)
     .              -  HLL(ii,jj) 
            end do              !i
         end do                 !j
      endif

    
      do j = 0, nc-1
         do i = 0, NGL-1
            jj=j+NGL
             iU(i+j*NGL) = 
     .         Energy*Sbig(i+jj*no)-Hbig(i+jj*no)
         end do                 !i
      end do                    !j





       call csolveg(NGL,nc,U,iU,ipvt,joutfile)

      do j = 0, NGL-1
         do i = 0, nc-1
            ii=i+NGL
            U(i + j*nc) =
     .         Energy*Sbig(ii+j*no)-Hbig(ii+j*no)
         end do                 !i
      end do                    !j





      a=(0.d0,0.d0)
      b=(1.d0,0.d0)

      allocate(sigmaL(0:nc*nc-1))
 
      call zgemm('N','N',nc,nc,ngl,b,U,nc,iU,
     &  ngl,a,sigmaL,nc)      





c
c     SigmaR:
c

      nu = NGL + nc



      if(USEBULK) then

      do j = 0, NGR-1
         jj=j+1
         do i = 0, NGR-1
            ii=i+1
            U(i + j*NGR) = HRR(ii,jj)
         end do                 !i
      end do                    !j

      else

      do j = 0, NGR-1
         jj=j+1
         jr = j + nu
         do i = 0, NGR-1
            ii=i+1
            ir = i + nu
            U(i + j*NGL) = Energy*Sbig(ir+jr*no) - Hbig(ir+jr*no)
     .         -  HRR(ii,jj)
         end do                 !i
      end do                    !j

      endif



      do j = 0, nc-1
         do i = 0, NGR-1
            jj=j+NGL
            ii = i + nu
             iU(i+j*NGR) =
     .          Energy*Sbig(ii+jj*no)-Hbig(ii+jj*no)
         end do                 !i
      end do                    !j

      call csolveg(NGR,nc,U,iU,ipvt,joutfile)


      do j = 0, NGR-1
         do i = 0, nc-1
            ii=i+NGL
            jj=j+nu
            U(i + j*nc) =
     .          Energy*Sbig(ii+jj*no)-Hbig(ii+jj*no)
         end do                 !i
      end do                    !j
      allocate(sigmaR(0:nc*nc-1))

      call zgemm('N','N',nc,nc,ngr,b,U,nc,iU,
     &  ngr,a,sigmar,nc)     




c
c     Total Hc:
c
      allocate(Hc(0:nc*nc-1))


      do j = 0, nc-1
         jj=j+NGL
         do i = 0, nc-1
            ii=i+NGL
            Hc(i + j*nc) =
     .           Energy*Sbig(ii+jj*no)-Hbig(ii+jj*no)
     .         - SigmaR(i+nc*j) - SigmaL(i+nc*j)
            GF(i+nc*j) = dcmplx(0d0,0d0)
         end do                 !i
         GF(j+nc*j) = dcmplx(1d0,0d0)
      end do                    !j





      call csolve(nc,Hc,GF,ipvt,joutfile)






c -----------       
c GammaL -> [U], GammaR -> [iU]

      do j = 0, nc-1
         do i = 0,j-1 
            U(i+nc*j)  = 
     .          SigmaL(i+nc*j) - dconjg(SigmaL(j+nc*i))
            U(j+nc*i)  = -dconjg( U(i+nc*j) )
            iU(i+nc*j) = 
     .          SigmaR(i+nc*j) - dconjg(SigmaR(j+nc*i))
            iU(j+nc*i) = -dconjg( iU(i+nc*j) )
        end do
            U(i+nc*i)  = 
     .          SigmaL(i+nc*i) -dconjg(SigmaL(i+nc*i))
            iU(i+nc*i) = 
     .          SigmaR(i+nc*i) -dconjg(SigmaR(i+nc*i))
      end do



 



c ...alternatively..
c      do j = 0, nc-1
c         do i = 0,nc-1
c            U(i+nc*j)  = 
c     .          SigmaL(i+nc*j) - dconjg(SigmaL(j+nc*i))
c            iU(i+nc*j) = 
c     .          SigmaR(i+nc*j) - dconjg(SigmaR(j+nc*i))
c         end do
c      end do



       deallocate(sigmaL)
       deallocate(sigmaR)

c =============================================================
c *** TRANSMISSION
c =============================================================
c tt  =  GammaL.(GF*.GammaR.GF):

c     -GammaR.GF == -[iU].GF -> [Hc]
      CALL zgemm('N','N',nc,nc,nc,dcmplx(-1d0,0d0),
     &     iU,nc,GF,nc,dcmplx(0d0,0d0),Hc,nc) 
      
c     GF*.GammaR.GF == GF*.[Hc] -> [iU]
      CALL zgemm('C','N',nc,nc,nc,dcmplx(1d0,0d0),
     &     GF,nc,Hc,nc,dcmplx(0d0,0d0),iU,nc) 





       deallocate(Hc)


c     GammaL.GF*.GammaR.GF == [U].[iU] -> tt




       CALL zgemm('N','N',nc,nc,nc,dcmplx(1d0,0d0),
     &      U,nc,iU,nc,dcmplx(0d0,0d0),tt,nc)




      do i=0,nc*nc - 1
         GFRGF(i) = dcmplx(0.5d0,0d0)*iU(i)
      end do




      deallocate(U)
      deallocate(iU)

c ==============================================================

c Trace of tt:
      csum=dcmplx(0d0,0d0)
      do i=0,nc-1
         csum = csum + tt(i+nc*i)
       end do
       
       tottrans = DREAL(csum)
c       open(56,position='append')
c       write(56,'(F9.5,e16.8)') Energy,tottrans
c       close(56)


 1234  continue ! inversion failed
C ====================================================================
      deallocate(ipvt)
      RETURN
      END
C ====================================================================
