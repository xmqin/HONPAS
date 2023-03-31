! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine redcel( alat, ucell, scell, mscell )

c *********************************************************************
c Reads the unit cell (and supercell, if present) vectors from the
c data file.
c Written by J.M.Soler
c ********* Output ****************************************************
c real*8 alat       : Lattice constant
c real*8 ucell(3,3) : Unit cell vectors ucell(ixyz,ivec)
c real*8 scell(3,3) : Supercell vectors (equal to ucell if no supercell)
c integer mscell(3,3): Supercell vectors in units of ucell, i.e.
c                      scell(ix,i) = Sum_j( ucell(ix,j)*mscell(j,i) )
c *********************************************************************
C
C  Modules
C
      use precision, only : dp
      use parallel,  only : Node
      use sys,       only : die
      use fdf
      implicit          none

      integer, intent(out)  ::      mscell(3,3)
      real(dp), intent(out) ::      alat, scell(3,3), ucell(3,3)

C Internal variables
      integer           i, ix, j
      real(dp)          alp, alplp, betlp, blp, cell(3,3), clp,
     .                  gamlp, pi, xxx

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Lattice constant
      alat = fdf_physical('LatticeConstant',0.0_dp,'Bohr')

C Lattice vectors

      if ( fdf_isblock('LatticeParameters') .and.
     .    fdf_isblock('LatticeVectors') ) then
         call die('redcel: ERROR: Lattice vectors doubly ' //
     $       'specified: by LatticeVectors and by LatticeParameters.')
      endif

      if ( fdf_block('LatticeParameters',bfdf) ) then
        if (.not. fdf_bline(bfdf, pline))
     .    call die('redcel: ERROR in LatticeParameters block')
        alp = fdf_bvalues(pline,1)
        blp = fdf_bvalues(pline,2)
        clp = fdf_bvalues(pline,3)
        alplp = fdf_bvalues(pline,4)
        betlp = fdf_bvalues(pline,5)
        gamlp = fdf_bvalues(pline,6)
*       write(6,'(a)')
*    .   'redcel: Lattice Parameters (units of Lattice Constant) ='
*       write(6,'(a,3f10.5,3f9.3)')
*    .   'redcel: ',alp,blp,clp,alplp,betlp,gamlp
        pi = acos(-1.d0)
        alplp = alplp * pi/180.0_dp
        betlp = betlp * pi/180.0_dp
        gamlp = gamlp * pi/180.0_dp
        cell(1,1) = alp
        cell(2,1) = 0.0_dp
        cell(3,1) = 0.0_dp
        cell(1,2) = blp * cos(gamlp)
        cell(2,2) = blp * sin(gamlp)
        cell(3,2) = 0.0_dp
        cell(1,3) = clp * cos(betlp)
        xxx = (cos(alplp) - cos(betlp)*cos(gamlp))/sin(gamlp)
        cell(2,3) = clp * xxx
        cell(3,3) = clp * sqrt(sin(betlp)*sin(betlp) - xxx*xxx)

        ! Make sure that we do not get very small numbers
        where (abs(cell) < 1.0e-14) cell = 0.0_dp

      elseif ( fdf_block('LatticeVectors',bfdf) ) then

        do i = 1,3
          if (.not. fdf_bline(bfdf, pline))
     .      call die('redcel: ERROR in LatticeVectors block')
          cell(1,i) = fdf_bvalues(pline,1)
          cell(2,i) = fdf_bvalues(pline,2)
          cell(3,i) = fdf_bvalues(pline,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.0_dp
          enddo
          cell(i,i) = 1.0_dp
        enddo
      endif
      call fdf_bclose(bfdf)

      ! Note that cell will be zero if alat is not specified
      ! This is not very intuitive

      cell = alat * cell

C     Super cell 

      if ( fdf_block('SuperCell',bfdf) ) then
        if ( alat .eq. 0.0_dp ) then
          call die('redcel: ERROR: LatticeConstant required' //
     .             ' to define SuperCell')
        endif
        do i = 1,3
          if (.not. fdf_bline(bfdf, pline))
     .      call die('redcel: ERROR in SuperCell block')
          mscell(1,i) = fdf_bvalues(pline,1)
          mscell(2,i) = fdf_bvalues(pline,2)
          mscell(3,i) = fdf_bvalues(pline,3)
        enddo
        call fdf_bclose(bfdf)
      else
        do i = 1,3
          do j = 1,3
            mscell(j,i) = 0
          enddo
          mscell(i,i) = 1
        enddo
      endif

      ucell = cell

      do i = 1,3
        do ix = 1,3
          scell(ix,i) = ucell(ix,1) * mscell(1,i) +
     .                  ucell(ix,2) * mscell(2,i) +
     .                  ucell(ix,3) * mscell(3,i)
        enddo
      enddo

      end
