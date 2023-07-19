C ##################################################################
C ## Subroutine which read-in GF for the 1x1 surface cell         ##
C ## and returns the GF for the NA1xNA2 surface cell              ##
C ## by summing over the appropriate q-points                     ## 
C ##                                                              ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

      subroutine getsfe(USEBULK,joutfile,jgfu,
     &     HAA,SAA,ZEnergy,nq1,
     &     NG2,nua,lasto,NA1,NA2,
     &     NG,SFE,GAMMA)


      use m_tbt_gf

c-----------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------
      real*8 EPS
      parameter(EPS=1d-7)

      real*8 Pi
      parameter(Pi=3.14159265358979323846264338327950288419717d0)

 
c ==================================================================
c INPUT
      logical USEBULK           ! if true return self-energy only else
c                                 return z*S-H-Sigma for bulk
      integer joutfile          ! fileunit for OUT

      integer jgfu              !fileunit from which we read GF

      integer NA1,NA2           ! repetitions of surface UC.
      integer NG                ! size of GF in NA1xNA2 supercell

      integer nua               ! no. uc (1x1) atoms
      integer lasto(0:nua)      ! index of last orbitals of uc atoms.

      complex*16 ZEnergy 

      integer nq1                !no. q1-points <= NA1*NA2

c ----------- 
c 1x1 surface cell
c
      integer NG2               ! size of used part of GAA in UC.
c     GF at fixed energy: ZEnergy 
c      complex*16  GAA(NG2,NG2,nq1)

c Note: HAA,SAA is 
c  read-in first time and *must* be saved in calling routine!
c     Ham,Ovl:
      complex*16 HAA(NG2,NG2,nq1),SAA(NG2,NG2,nq1)
c -----------

c======================================================================
c OUTPUT: 
      complex*16 SFE(NG,NG) !Electrode atoms self-energy for USEBULK true
c                           !else z*S - H - Sigma for ideal surface
      complex*16  GAMMA(NG,NG)  ! 0.5*(Sigma - Sigma^dagger) 

c=====================================================================
c Helpers, tempos etc...
      
      complex*16 ZEi, wgfi
      integer iEni
      integer iq1
      integer io,ia1,ia2,iau,iuo
      integer jo,ja1,ja2,jau,juo
      complex*16 cphase,ccphase,ctmp

      complex*16, allocatable :: WRK1(:,:),WRK2(:,:),GAA(:,:,:)
      integer,dimension(:), allocatable :: ipvt

      integer i,j

c=======================================================================
c     BEGIN:
c=======================================================================



     
c
c     Init:
c
       allocate(WRK1(NG,NG))
       allocate(WRK2(NG,NG))
       allocate(ipvt(NG))

      do io=1,NG
         do jo=1,NG
            WRK1(jo,io) = dcmplx(0d0,0d0)
            WRK2(jo,io) = dcmplx(0d0,0d0)
            SFE(jo,io)  = dcmplx(0d0,0d0)
         end do                 !jo
         SFE(io,io) = dcmplx(1d0,0d0) !to be used in inv.
      end do                    !io


        allocate(GAA(NG2,NG2,nq1)) 
        if (LJob) then

           do j=1,NG2
              do i=1,NG2
                 GAA(i,j,nq1)=LGAAq(i+(j-1)*NG2,1) 
                 HAA(i,j,nq1)=LHAA(i+(j-1)*NG2,1)
                 SAA(i,j,nq1)=LSAA(i+(j-1)*NG2,1)
              enddo
           enddo

        else if (RJob) then
       
           do j=1,NG2
              do i=1,NG2
                 GAA(i,j,nq1)=RGAAq(i+(j-1)*NG2,1) 
                 HAA(i,j,nq1)=RHAA(i+(j-1)*NG2,1)
                 SAA(i,j,nq1)=RSAA(i+(j-1)*NG2,1)
              enddo
           enddo
           
        else

    
        end if    

      ctmp=dcmplx(0d0,0d0)
   
C FDN
       WRK1(:,:)=GAA(:,:,1)
       WRK2(:,:)=ZEnergy*SAA(:,:,1)-HAA(:,:,1)
C FDN



c      if(iEni.eq.1) write(*,*) 'hola getsfe: H(1,28)=',ctmp

         
c     Find inverse of GAA --> SFE
      call csolve(NG,WRK1,SFE,ipvt,joutfile)


c     SFE = Z*S - HAA - Sigma ---> WRK2 - SFE = Sigma -> WRK1:
      do io=1,NG
         do jo=1,NG
            WRK1(jo,io) = WRK2(jo,io)-SFE(jo,io)
         end do                 !jo
      end do                    !io

c
c Gamma 
c
      do io=1,NG
         do jo=1,NG
            GAMMA(jo,io) = 
     .           dcmplx(0d0,0.5d0)*(WRK1(jo,io) - 
     .           dconjg(WRK1(io,jo)))
         end do                 !jo
      end do                    !io
      
      if(.not.USEBULK) then
c      
c    Sigma=(Z*S - H) - inv[GF] ---> SFE
c
         do io=1,NG
            do jo=1,NG
               SFE(jo,io) = WRK1(jo,io)
            end do              !jo
         end do                 !io
      end if                    !USEBULK


       deallocate(WRK1)
       deallocate(WRK2)
       deallocate(GAA)
       deallocate(ipvt)

c======================================================================
      return
      end
c======================================================================







