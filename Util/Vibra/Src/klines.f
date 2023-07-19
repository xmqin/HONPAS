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
      subroutine klines( maxk, nk, nlines, lastk, label, kpoint )
C *********************************************************************
C Finds k-points at given lines to compute band structures
C Reads the information in FDF format, and passes it to calling routine
C
C Extracted from bands.f (written by J.Soler, August 1997),
C with some modifications in behaviour.
C P. Ordejon, August 1998
C **************************** INPUT **********************************
C integer maxk                : Last dimension of kpoint and ek
C                               (maximum number of k points)
C *************************** OUTPUT **********************************
C integer nk                  : Number of band k-points
C integer nlines              : Number of band k-lines
C integer lastk(maxlin)       : Index of last k-point in each k-line
C character*8 label(maxlin)   : Label of each k-line
C real*8  kpoint(3,maxk)      : k point vectors
C **************************** UNITS **********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C ***************** BEHAVIOUR *****************************************
C - k-points are read from labels BandLines and 
C   BandLinesScale of the input fdf data file. If these labels are 
C   not present, it returns with nk=0
C - Allowed values for BandLinesScale are ReciprocalLatticeVectors and
C   pi/a (default). If another value is given, it returns with nk=0
C   after printing a warning.
C - If nk>maxk, a warning is printed before return.
C ***************** USAGE *********************************************
C Example of fdf band lines specification for an FCC lattice.
C Last column is an optional LaTex label (for plot)
C     BandLinesScale  pi/a
C     %block BandLines                  # These are comments
C      1  0.000  0.000  0.000  \Gamma   # Begin at Gamma
C     25  2.000  0.000  0.000     X     # 25 points from Gamma to X
C     10  2.000  1.000  0.000     W     # 10 points from X to W
C     15  1.000  1.000  1.000     L     # 15 points from W to L
C     20  0.000  0.000  0.000  \Gamma   # 20 points from L to Gamma
C     25  1.500  1.500  1.500     K     # 25 points from Gamma to K
C     %endblock BandLines
C
C Example for BCC:
C     BandLinesScale  pi/a
C     %block BandLines
C      1  0.000  0.000  0.000  \Gamma
C     20  2.000  0.000  0.000     H
C     15  1.000  1.000  0.000     N
C     15  0.000  0.000  0.000  \Gamma
C     20  1.000  1.000  1.000     P
C     10  1.000  1.000  0.000     N
C     10  1.000  1.000  1.000     P
C     20  2.000  2.000  2.000     H
C     %endblock BandLines
C
C Example for HCP (an angle of 120 deg is assumed between reciprocal
C lattice vectors, what implies an angle of 60 deg between the first 
C two vectors of cell argument):
C     BandLinesScale  ReciprocalLatticeVectors
C     %block BandLines
C      1  0.000000000  0.000000000  0.000000000  \Gamma
C     20  0.666666667  0.333333333  0.000000000     K 
C     10  0.500000000  0.000000000  0.000000000     M
C     20  0.000000000  0.000000000  0.000000000  \Gamma
C     15  0.000000000  0.000000000  0.500000000     A
C     20  0.666666667  0.333333333  0.500000000     H
C     10  0.500000000  0.000000000  0.500000000     L
C     %endblock BandLines
C
C If only given points (not lines) are desired, simply specify 1 as 
C the number of points along the line.
C *********************************************************************
      use fdf

      implicit          none

      integer maxlin
      parameter (maxlin = 1000)

      integer           lastk(maxlin), maxk, nk, nlines
      double precision  kpoint(3,maxk)
      character         label(maxlin)*8

      external          parse_vibra
C *********************************************************************

C  Internal variables 

      integer
     .  i, ik, il, integs(4), iu, ix, 
     .  lastc, lc(0:3), 
     .  ni, nkl, nn, nr, nv

      double precision
     .  alat, 
     .  pi, rcell(3,3), reals(4),
     .  ucell(3,3), values(4)

      character 
     .  line*130, names*80, scale*30

      logical
     .  found, frstme

      save frstme
      data frstme /.true./

C     Find if there are band-lines data
      if ( fdf_block('BandLines',iu) ) then

C         Find lattice constant
          alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
          if (alat .eq. 0.d0) then
            write(6,'(a)') 'klines: ERROR: Lattice constant required'
            goto 999
          endif

C         Find scale used in k point data
          scale = fdf_string( 'BandLinesScale', 'pi/a' )
          if (scale .eq. 'pi/a') then
            pi = 4.d0 * atan(1.d0)
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            if ( fdf_block('LatticeVectors',iu) ) then
              do i = 1,3
                read(iu,*) (ucell(ix,i),ix=1,3)
                do ix = 1,3
                  ucell(ix,i) = alat * ucell(ix,i)
                enddo
              enddo
            else
              do i = 1,3
                do ix = 1,3
                  ucell(ix,i) = 0.d0
                enddo
                ucell(i,i) = alat
              enddo
            endif
            call reclat( ucell, rcell, 1 )
          else
            write(6,'(a,/,2a,/,a)')
     .        'klines: WARNING: Invalid value for BandLinesScale',
     .        'klines: Allowed values are pi/a and',
     .              ' ReciprocalLatticeVectors',
     .        'klines: No band calculation performed'
          endif

C         Loop on data lines
          nk = 0
          do il = 1,maxlin

C           Read and parse data line
            if (il .eq. 1) found = fdf_block('BandLines',iu)
            read(iu,'(a)',end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse_vibra( line(1:lastc), nn, lc,
     $           names, nv, values, ni, integs, nr, reals )

C           Check if data are already finished
            if (nv .ge. 3) then

C             Add to total number of k points
              nkl = integs(1)
              nk = nk + nkl

C             If there is room to store k points
              if (nk .le. maxk) then

C               Find last point in line
                if (scale .eq. 'pi/a') then
                  kpoint(1,nk) = values(2) * pi / alat
                  kpoint(2,nk) = values(3) * pi / alat
                  kpoint(3,nk) = values(4) * pi / alat
                elseif (scale .eq. 'ReciprocalLatticeVectors') then
                  do ix = 1,3
                    kpoint(ix,nk) = rcell(ix,1) * values(2) +
     .                              rcell(ix,2) * values(3) +
     .                              rcell(ix,3) * values(4)
                  enddo
                endif

C               Find points along the line
                do ik = 1,nkl-1
                  do ix = 1,3
                    kpoint(ix,nk-nkl+ik) =
     .                  kpoint(ix,nk-nkl) * dble(nkl-ik) / dble(nkl) + 
     .                  kpoint(ix,nk)     * dble(ik)     / dble(nkl)
                  enddo
                enddo

C               Find point label
                if (nn .gt. 0) then
                  label(il) = names(1:lc(1))
                else
                  label(il) = ' '
                endif
                lastk(il) = nk

              endif
            else
C             No more lines to read => Exit do loop
              goto 50
            endif
          enddo
            write(6,'(a)') 'klines: ERROR. Parameter maxlin too small'
   50     continue
          nlines = il - 1
      else
C No k-point data available => return with nk=0
        write(6,'(a)') 'klines: WARNING: No k-points specified!!!'
        goto 999
      endif

C Check parameter maxk 
      if (nk .gt. maxk) then
        write(6,'(/,a)')
     .    'klines: WARNING: parameter maxk too small'
        goto 999
      endif

C This is the only exit point 
  999 continue
      frstme = .false.
      end


