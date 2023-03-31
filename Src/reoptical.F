! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine reoptical(DoOptical,nmeshk,displ,wmin,wmax,smooth,
     .                     efield,nbands,option,escissor)

C *******************************************************************
C Reading k_space grids, frequency range, vector defining 
C the electric field (or propagation direction of the ligth)
C for the calculation of the optical absorption, and dielectric 
C function. 
C
C Reads fdf block. 
C    First the frequency range is specified, and the width
C    of the 'delta' function ( with its units)
C    secondly the number of bands to be used and 
C    the option: 'polarized', 'unpolarized', 'polycrystal'
C    the k-grid for Brillouin-Zone integration, ('yes' if the grid
C    is displaced, this implies that Gamma is not included in 
C    samplings with an odd number of points)
C    and finally a vector, only with options  'polarized' and 
C    'unpolarized', which is interpreted as the electric field
C    for the 'polarized' option, or as the propagation direction
C    for the 'unpolarized' option.  
C    For example: 
C        %block  Optical
C            0.0  5.0 0.1 eV 
C            5  polarized
C            5  5  5  yes
C            1.000 0.000 0.500
C        %enblock Optical
C
C      
C Written by D. Sanchez-Portal. August 1999
C Restyled for f90 version by JDG. June 2004
C ********* INPUT/OUTPUT ********************************************
C integer   nbands        : On input, the number of bands in the 
C                           material and on output equals the number 
C                           of bands to be used.
C ********* OUTPUT **************************************************
C logical   DoOptical     : if .true. then do optical calculation
C integer   nmeshk(3)     : Grids for the integration in reciprocal 
C                           space
C real*8    displ         : Displacement of the bidimensional 
C                           integration grid from the origin
C real*8    wmin          : Minimum frequency
C real*8    wmax          : Maximum frequency
C real*8    smooth        : width of the 'delta' function 
C real*8    efield(3)     : Vector specifying the direction of 
C                           the electric field
C character*12  option    : Polarization option
C *******************************************************************

      use precision
      use parallel,    only : Node, Nodes
      use fdf
      use sys,         only : die
      use m_spin, only: NonCol, SpOrb
      
      implicit          none

C Passed variables 
      integer           nmeshk(3), nbands
      logical           DoOptical
      real(dp)          displ, wmin, wmax, efield(3), smooth
      real(dp)          escissor
      character*12      option

C Internal variables 
      integer           ni, nn, nr, ix
      logical           offset
      logical           present

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Initialise outputs to defaults
      do ix = 1,3
        nmeshk(ix) = 0
        efield(ix) = 0.0d0 
      enddo 
      wmin = 0.0d0
      wmax = 0.0d0
      displ = 0.0d0
      smooth = 0.0d0
      option = 'polycrystal'

C Read flag for optical calculation
      DoOptical = fdf_boolean('OpticalCalculation',.false.)
      if ( NonCol .or. SpOrb ) then
         if ( DoOptical ) then
            if (Node.eq.0) then
               write(6,'(/2a)') 'reoptical: Not possible ',
     + 'with non-collinear spin configuration.'
            end if
         end if
         DoOptical = .false.
      end if

      if (DoOptical) then
C Read energy range
        wmin = fdf_physical('Optical.EnergyMinimum', 0.0d0,'Ry')
        wmax = fdf_physical('Optical.EnergyMaximum',10.0d0,'Ry')
        wmin = abs(wmin)
        wmax = abs(wmax)

C Read Gaussian smearing width
        smooth = fdf_physical('Optical.Broaden',0.0d0,'Ry')
        smooth = abs(smooth)

C Read scissor operator
        escissor = fdf_physical('Optical.Scissor',0.0d0,'Ry')
        escissor = abs(escissor)

C Read number of bands
        nn = nbands
        nbands = fdf_integer('Optical.NumberOfBands',nn)
        nbands = min(nbands,nn)
        nbands = max(nbands,1)

C Read whether to offset the mesh from gamma
        offset = fdf_boolean('Optical.OffsetMesh',.false.)
        if (offset) displ = 0.5d0

C Read BZ sampling grid
        present = fdf_block('Optical.Mesh',bfdf)
        if (present) then
          if (.not. fdf_bline(bfdf,pline))
     +      call die('reoptical: ERROR in Optical.Mesh block')
          ni = fdf_bnintegers(pline)
          if (ni .lt. 3) then 
            call die('reoptical: Check syntax in Optical.Mesh block')
          endif  
          nmeshk(1) = fdf_bintegers(pline,1)
          nmeshk(2) = fdf_bintegers(pline,2)
          nmeshk(3) = fdf_bintegers(pline,3)
          call fdf_bclose(bfdf)
        endif

C Read whether this is a polarised, unpolarised or polycrystal calc
        option = fdf_string('Optical.PolarizationType','polycrystal')
        if ((option.ne.'polycrystal') .and. (option.ne.'unpolarized')
     +      .and. (option.ne.'polarized')) then
          if (Node.eq.0) then
            write(6,'(/a)') 'reoptical: WARNING '
            write(6,'(a,a)')
     +        'reoptical: read Optical.PolarizationType: ', option
            write(6,'(a,/a,/a)')
     +      'reoptical: input Optical.PolarizationType must be one of ',
     +      'reoptical: the followings: ', 
     +      'reoptical: polarized, unpolarized, or polycrystal' 
            write(6,'(a)')
     +      'reoptical: Assuming Optical.PolarizationType = polycrystal'
            option = 'polycrystal'
          endif
        endif

C Read datablock for electric field
        if (option .ne. 'polycrystal') then
          present = fdf_block('Optical.Vector',bfdf)
          if (present) then
            if (.not. fdf_bline(bfdf,pline))
     +        call die('reoptical: ERROR in Optical.Vector block')
            nr = fdf_bnreals(pline)
            if (nr .lt. 3) then   
              call die('reoptical: Check syntax in ' //
     +                 'Optical.Vector block')
            endif
            efield(1) = fdf_breals(pline,1)
            efield(2) = fdf_breals(pline,2)
            efield(3) = fdf_breals(pline,3)
          endif
        endif

      endif

      return
      end
