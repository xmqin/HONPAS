       subroutine calc_green(nv,zenergy,h00,s00,h01,s01,
     &  gs,zdos,joutfile,tleft,tdos)

       use m_tbt_gf, only : WGFFiles

       implicit none
       integer nv

       integer nv2,iter
       integer i,j,ic,ic2
       integer joutfile
       logical tleft,tdos

       complex*16 a,b,zdos
       double precision accur,ro
c       parameter (accur=5.d-17)
       parameter (accur=1.d-15)
       complex*16 ZEnergy 
       complex*16 h00(0:nv*nv-1),s00(0:nv*nv-1)
       complex*16 h01(0:nv*nv-1),s01(0:nv*nv-1)
       complex*16 gs(0:nv*nv-1)

       integer, dimension (:), allocatable:: ipvt
       complex*16, dimension (:), allocatable:: 
     .    rh,rh1,rh3,alpha,beta,ab,ba,gb,gs2

       real*8 Pi
       parameter(Pi=3.14159265358979323846264338327950288419717d0)

       allocate(ipvt(nv))
       allocate(rh(0:2*nv*nv))
       allocate(rh1(0:2*nv*nv))
       allocate(rh3(0:4*nv*nv))
       allocate(alpha(0:nv*nv-1))
       allocate(beta(0:nv*nv-1))
       allocate(ba(0:nv*nv-1))
       allocate(ab(0:nv*nv-1))
       allocate(gb(0:nv*nv-1))
       allocate(gs2(0:nv*nv-1))

       a=(1.d0,0.d0)
       b=(0.d0,0.d0)
       nv2 =2*nv

        
      do i=0,nv*nv-1
        gb(i) = zenergy*s00(i)-h00(i)
        alpha(i) = h01(i)-zenergy*s01(i)
      end do

      do i=0,nv*nv-1
        gs(i) = gb(i)
        gs2(i) = gb(i)
      end do                

      do j=0,nv-1
       do i=0,nv-1
        ic = i + nv*j
        ic2 = j + nv*i
        beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
      end do


        iter=0
1000    continue
        iter=iter+1


      do j=0,nv2-1
         do i=0,nv-1
              ic =i + j*nv
              ic2=i + (j - nv)*nv
              if(j.lt.nv)then
               rh(ic) = alpha(ic)
              else
               rh(ic) = beta(ic2)
              endif                
         end do                 
      end do                  

      do i=0,nv*nv-1
        rh3(i) = gb(i)
      end do                


      call csolveg(nv,nv2,rh3,rh,ipvt,joutfile) 

      do j=0,nv-1
         do i=0,nv2-1
             ic =i + j*nv
             ic2 =i-nv + j*nv
              if(i.lt.nv)then
               rh1(i + nv2*j) = alpha(ic)
              else
               rh1(i + nv2*j) = beta(ic2)
             end if
         end do                
      end do                   
      call zgemm('N','N',nv2,nv2,nv,a,rh1,nv2,rh,nv,b,rh3,nv2)

      do j=0,nv-1
         do i=0,nv2-1
             ic =i + j*nv
             ic2 =i-nv + j*nv
              if(i.lt.nv)then
                alpha(ic) = rh3(i + nv2*j) 
              else
                ba(ic2) = -rh3(i + nv2*j) 
             end if
         end do                 
      end do                   
      do j=nv,nv2-1
         do i=0,nv2-1
              ic=i + (j - nv)*nv
              ic2=i - nv + (j - nv)*nv
              if(i.lt.nv)then
                ab(ic)= -rh3(i + nv2*j) 
              else
                beta(ic2)= rh3(i + nv2*j) 
             end if
         end do                 
      end do                   
 
      do i=0,nv*nv-1
       gb(i) =  gb(i) + ba(i) + ab(i)
       gs(i) =  gs(i) + ab(i) 
       gs2(i) =  gs2(i) + ba(i) 
      end do                    

      ro =-1.0
      do j =0,nv*nv-1
        ro =max(ro,dreal(ab(j))**2+dimag(ab(j))**2)
      end do                   
      ro =dsqrt(ro)

      if(ro.gt.accur) go to 1000

      do i=0,nv*nv-1
        rh3(i) = gs(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,joutfile)

      do i=0,nv*nv-1
        gs(i) = rh(i)
      end do

      do i=0,nv*nv-1
        rh3(i) = gs2(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,joutfile)

      do i=0,nv*nv-1
        gs2(i) = rh(i)
      end do


      do i=0,nv*nv-1
        rh3(i) = gb(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,joutfile)

      do i=0,nv*nv-1
        gb(i) = rh(i)
      end do


c      ----      DOS     -----

      do i=0,nv*nv-1
        alpha(i) = h01(i)-zenergy*s01(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,gs2,nv,alpha,nv,b,ab,nv)
      do i=0,nv*nv-1
        alpha(i) = ab(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,alpha,nv,gb,nv,b,ab,nv)


      do j=0,nv-1
       do i=0,nv-1
        ic = i + nv*j
        ic2 = j + nv*i
        beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
      end do
      call zgemm('N','N',nv,nv,nv,a,gs,nv,beta,nv,b,ba,nv)
      do i=0,nv*nv-1
        beta(i) = ba(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,beta,nv,gb,nv,b,ba,nv)


      do i=0,nv*nv-1
        rh3(i) = 0.0
      end do

      call zgemm('N','N',nv,nv,nv,a,gb,nv,s00,nv,b,rh3,nv)
      call zgemm('N','C',nv,nv,nv,a,ab,nv,s01,nv,a,rh3,nv)
      call zgemm('N','N',nv,nv,nv,a,ba,nv,s01,nv,a,rh3,nv)


      zdos =0.0

      do j=0,nv-1
        zdos = zdos + (rh3(j*(nv+1)))
      end do



      if(tleft) then

       do i=0,nv*nv-1
        gs(i) =  gs2(i) 
       end do

      

      endif

       deallocate(ipvt)
       deallocate(rh)
       deallocate(rh1)
       deallocate(rh3)
       deallocate(alpha)
       deallocate(beta)
       deallocate(ba)
       deallocate(ab)
       deallocate(gb)
       deallocate(gs2)

      return
      end
