! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
c Common arrays for plrho routines
      integer maxx, maxy
      parameter ( maxx = 1240 )
      parameter ( maxy = 1024 )
      integer ixmin, ixmax, iymin, iymax,
     .        pixmap(0:maxx-1,0:maxy-1)
      real     zbuff(0:maxx-1,0:maxy-1)
      common /compix/ ixmin, ixmax, iymin, iymax, pixmap, zbuff

