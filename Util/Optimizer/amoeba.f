      module m_amoeba
      use precision, only: dp
      use vars_module, only: constrained

      implicit none

      public :: amoeba
      private

      CONTAINS

      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,req_itmax)
      INTEGER iter,mp,ndim,np
      REAL(dp)  ftol,p(mp,np),y(mp)
      integer, intent(in), optional :: req_itmax

      integer, parameter :: ITMAX_default = 2000

      interface
         function funk(x) result(value)
         use precision, only: dp
         real(dp), dimension(:), intent(in) :: x
         real(dp)                           :: value
      end function funk
      end interface

      integer :: itmax

      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL(dp)  rtol,sum,swap,ysave,ytry
      real(dp), parameter :: tiny = 1.0e-12_dp

      REAL(dp)  :: psum(ndim)  ! Automatic

      if (present(req_itmax)) then
         itmax = req_itmax
      else
         itmax = itmax_default
      endif

      iter=0
1     continue
      do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)) + tiny)
      print *, "Fractional dispersion: ", rtol
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      if (iter.ge.ITMAX) RETURN

      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0_dp)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0_dp)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5_dp)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum(1:ndim))
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END subroutine amoeba
!------------------------------------------------------------------------------
C
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np
      REAL(dp)  amotry,fac,p(mp,np),psum(np),y(mp)

      INTEGER j
      REAL(dp)  fac1,fac2,ytry
      
      REAL(dp)  :: ptry(size(psum))  ! Automatic

      interface
         function funk(x) result(val)
         use precision, only: dp
         real(dp), dimension(:), intent(in) :: x
         real(dp)                           :: val
         end function funk
      end interface


      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=constrained(psum(j)*fac1-p(ihi,j)*fac2,j)
11    continue

      ytry=funk(ptry(1:ndim))
      print *, "New point: ", ptry(1:ndim), " --- ",  ytry

      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END function amotry

C  (C) Copr. 1986-92 Numerical Recipes Software P#5,t$D#)K#)?.

      end module m_amoeba
