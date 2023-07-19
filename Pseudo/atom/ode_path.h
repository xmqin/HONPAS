c-----
c  Ode_path.h
c
c     See Numerical Recipes (p.559)
c
      integer nmax      
      parameter (nmax = 10)
      integer kmax, kount
      double precision dxsav
      double precision xp(200), yp(nmax,200)
c
      common /path/ kmax, kount, dxsav, xp, yp
      save /path/
c-----
