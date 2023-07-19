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
      use parsing,     only : parse
      use sys,         only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

C Passed variables 
      integer           nmeshk(3), nbands
      logical           DoOptical
      real(dp)          displ, wmin, wmax, efield(3), smooth
      real(dp)          escissor
      character*12      option

C Internal variables 
      character         line*130, names*80
      integer           ni, nn, nr, nv, iu, ix
      integer           integs(4), lastc, lc(0:3)
      logical           offset
      logical           present
      real(dp)          reals(4), values(4)
#ifdef MPI
      integer           MPIerror
#endif

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

      if (Node.eq.0) then
C Read flag for optical calculation
        DoOptical = fdf_boolean('OpticalCalculation',.false.)

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
          present = fdf_block('Optical.Mesh',iu)
          if (present) then
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (.not.(ni.ge.3)) then 
              call die('reoptical: Check syntax in Optical.Mesh block')
            endif  
            nmeshk(1) = integs(1)
            nmeshk(2) = integs(2)
            nmeshk(3) = integs(3)
          endif

C Read whether this is a polarised, unpolarised or polycrystal calc
          option = fdf_string('Optical.PolarizationType','polycrystal')
          if(option.ne.'polycrystal' .and. option.ne.'unpolarized' 
     +         .and. option.ne.'polarized' ) then 
               write(6,'(/a)')
     + 'reoptical: WARNING '
              write(6,'(a,a)')
     + 'reoptical: read Optical.PolarizationType: ',option
              write(6,'(a,/a,/a)')
     + 'reoptical: input Optical.PolarizationType must be one of ',
     + 'reoptical: the followings: ', 
     + 'reoptical: polarized, unpolarized, or polycrystal' 
              write(6,'(a)')
     + 'reoptical: Assuming Optical.PolarizationType = polycrystal'
              option='polycrystal'
          endif

C Read datablock for electric field
          if (option.ne.'polycrystal') then
            present = fdf_block('Optical.Vector',iu)
            if (present) then
              read(iu,'(a)', end=50) line
              lastc = index(line,'#') - 1
              if (lastc .le. 0) lastc = len(line)
              call parse( line(1:lastc), nn, lc, names, nv, values,
     .                    ni, integs, nr, reals )
              if (.not.(nr.ge.3)) then   
                call die('reoptical: Check Optical.Vector block')
              endif  
              efield(1) = reals(1)
              efield(2) = reals(2)
              efield(3) = reals(3)
            endif
          endif

        endif

      endif

  50  continue 

#ifdef MPI
C Distribute information to all nodes
      if (Nodes.gt.1) then
        call MPI_Bcast(DoOptical,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
        if (DoOptical) then
          call MPI_Bcast(nmeshk,3,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(nbands,1,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(displ,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(wmin,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(wmax,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(smooth,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(escissor,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(efield,3,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(option,12,MPI_character,0,
     .      MPI_Comm_World,MPIerror)
        endif
      endif
#endif

      return
      end
