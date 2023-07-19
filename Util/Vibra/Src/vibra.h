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
c  maxa    =  maximum number of atoms in the unit cell
c  maxx    =  maximum number of shells repeating the first lattice vector
c                   in the supercell
c  maxy    =  maximum number of shells repeating the second lattice vector
c                   in the supercell
c  maxz    =  maximum number of shells repeating the third lattice vector
c                   in the supercell
c  maxasc  =  maximum number of atoms in the supercell
c  maxd    =  maximum size of the dynamical matrix

      integer maxa,maxasc,maxx,maxy,maxz,maxd
      parameter (maxa = 8)
      parameter (maxx = 2)
      parameter (maxy = 2)
      parameter (maxz = 2)
      parameter (maxasc = maxa * (2*maxx+1) * (2*maxy+1) * (2*maxz+1))
      parameter (maxd = 3*maxa)
