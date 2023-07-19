c------
c     General operating parameters
c
      character ispp*1, icorr*2, nameat*2
      character irel*3, nicore*4
      double precision rsh, zel, znuc, zsh
c
      integer indd(5), indu(5)
      double precision rc(5), rc_input(5), cfac, rcfac,
     &                 logder_radius
      integer scheme, ncore, job, ifcore
      logical normal, polarized, relativistic, is_gga
c
      common /param/ scheme, job, ifcore, ncore,
     &               rsh, zel, znuc, zsh, rc, rc_input, cfac, rcfac,
     &               logder_radius, indu, indd
      common /par_char/ nameat, ispp, icorr, irel, nicore
      common /par_log/  normal, polarized, relativistic, is_gga
      save /param/, /par_char/, /par_log/
c------
