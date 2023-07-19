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
c $Id: constr.f,v 1.6 2003/06/23 09:46:16 ordejon Exp $

      subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
c *****************************************************************
c User-written routine to implement specific geometric constraints,
c by orthogonalizing the forces and stress to undesired changes.
c Arguments:
c real*8  cell(3,3)    : input lattice vectors (Bohr)
c integer na           : input number of atoms
c integer isa(na)      : input species indexes
c real*8  amass(na)    : input atomic masses
c real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
c real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
c real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
c integer ntcon        : total number of positions constr. imposed
c *****************************************************************
      implicit         none
      integer          na, isa(na), ntcon
      double precision amass(na), cell(3,3), fa(3,na),
     .                 stress(3,3), xa(3,na)

c Write here your problem-specific code.

      end

