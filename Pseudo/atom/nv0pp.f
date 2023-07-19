c
c $Id: nv0pp.f,v 1.1 1991/12/14 00:34:49 alberto Exp $
c
c $Log: nv0pp.f,v $
c Revision 1.1  1991/12/14 00:34:49  alberto
c Initial revision
c
      double precision function nv0pp(xgamma)
c
      implicit none
c
      include 'coeffs.h'
c
      double precision xgamma
c
      double precision accuracy
      parameter (accuracy = 1.d-12)
c
      double precision x1, x2
      logical bracketed
c
      double precision chg_mism, brent
      external chg_mism, brent
c
c     Propagate gamma to the common block.
c
      gamma = xgamma
c
c     Bracket delta
c
      x1 = -1.d0
      x2 =  1.d0
      call brac(chg_mism,x1,x2,bracketed)
      if (.not. bracketed) then
        write(6,9000)
 9000   format(//1x,'Error in nv0pp - delta not found')
        call ext(860+ang_moment+1)
      endif
c
c     Get the delta which is consistent with gamma and the {alpha}
c
      delta = brent(chg_mism,x1,x2,accuracy)
c
c     Compute V''(0)
c
      nv0pp = 8*( (2*ang_moment +5)*alpha + xgamma**2 )
c
      return
c
      end
