c-----
c Ode_blk.h
c
c     Provides the variables needed for the solution of
c     the Schroedinger equation by an adaptive stepsize method.
c
      integer jint
      double precision energ
      double precision v(nrmax)
c
      common /ode_blk/ energ, v, jint
      save /ode_blk/
c----
