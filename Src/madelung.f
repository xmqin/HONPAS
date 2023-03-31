! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine madelung(cell, shape, charnet, Emad)

c *******************************************************************
c Finds Madelung correction for charged systems
c (Makov & Payne, PRB 51, 4014 (1995))
c (Values of Madelung constant: 
C  Coldwell-Horsfall & Maradudin, J. Math. Phys. 1, 395 (1960))
c
c Written by P. Ordejon, March 1999.
c ********* INPUT ***************************************************
c real*8       cell(3,3): Lattice vectors (in Bohr)
c character*10 shape    : Cell shape (atom,molecule,chain,slab,bulk)
c real*8       charnet  : Net charge of the system
c ********* OUTPUT **************************************************
c real*8       Emad     : Madelung correction energy (in Ry)
c *******************************************************************
C
C  Modules
C
      use precision
      use parallel,   only : Node
      use sys, only: die

      implicit none

      real(dp)
     .  cell(3,3), charnet, Emad

      character
     .  shape*10

C Local variables.........

      real(dp)
     .  alpha, lv
 
      character
     .  ctype*4

C If system is charged, find if energy correction terms can be applied
C (cluster or molecule, with SC, FCC or BCC cell) .....................
      if (shape .ne. 'molecule' .and. shape .ne. 'atom') then
        if (Node.eq.0) then
          write(6,'(2a)')
     .      'madelung: WARNING: Charged system, but not'
     .      ,' an atom or molecule.'
          write(6,'(2a)')
     .      'madelung: WARNING: Energy correction terms'
     .      ,' can not be applied. '
          write(6,'(2a)')
     .      'madelung: WARNING: Continue only if you really know'
     .      ,' what you are doing.'
        endif
        Emad = 0.0_dp
        return
      else
        call typecell(cell,ctype,lv)
        if (ctype .eq. 'none') then
          if (Node.eq.0) then
            write(6,'(2a)')
     .        'madelung: WARNING: Charged system, but not'
     .        ,' SC, FCC or BCC cell.'
            write(6,'(2a)')
     .        'madelung: WARNING: Energy correction terms'
     .        ,' can not be applied. '
            write(6,'(2a)')
     .        'madelung: WARNING: Continue only if you really know'
     .        ,' what you are doing.'
          endif
          Emad = 0.0_dp
          return
        else if (ctype .eq. 'sc') then
          alpha = 5.6745947_dp
          if (Node.eq.0) then
            write(6,'(a)')
     .        'madelung: Charged system with SC cell'
          endif
        else if (ctype .eq. 'bcc') then
          alpha = 7.278473_dp
          if (Node.eq.0) then
            write(6,'(a)')
     .        'madelung: Charged system with BCC cell'
          endif
        else if (ctype .eq. 'fcc') then
          alpha = 9.1697536_dp
          if (Node.eq.0) then
            write(6,'(a)')
     .        'madelung: Charged system with FCC cell'
          endif
        else
          call die('madelung: ERROR: Wrong type of cell')
        endif
      endif

      Emad = charnet**2 * alpha / (2.0_dp * lv)
C ...................
      return
      end

