! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
c Common arrays for plrho routines
      integer maxx, maxy
      parameter ( maxx = 1240 )
      parameter ( maxy = 1024 )
      integer ixmin, ixmax, iymin, iymax,
     .        pixmap(0:maxx-1,0:maxy-1)
      real     zbuff(0:maxx-1,0:maxy-1)
      common /compix/ ixmin, ixmax, iymin, iymax, pixmap, zbuff

