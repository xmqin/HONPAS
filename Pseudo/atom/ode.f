c
      subroutine ode(x,y,dxdy)
c
      implicit none
c
      include 'radial.h'
      include 'ode_blk.h'
c
c     Returns the derivatives (right hand sides) for the
c     system of differential equations:
c
c     dy1/dx = y2
c
c     dy2/dx = ( V(x) - energy ) y1
c
c     appropriate for the solution of the Schroedinger equation.
c
      double precision x
      double precision y(2), dxdy(2)
      double precision potx, delta_pot
c
c     Interpolate (three points used)
c
      call hunt(r,nr,x,jint)
      call polint(r(jint),v(jint),3,x,potx,delta_pot)
c
      dxdy(1) = y(2)
      dxdy(2) = ( potx - energ ) * y(1)
c
      return
c
      end
