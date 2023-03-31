! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
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

