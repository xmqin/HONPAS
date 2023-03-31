! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module band
      use precision, only: dp
      use sys,       only: die
      use alloc,     only: re_alloc, de_alloc
      implicit none
      private

      integer, public :: nbk, maxbk

      integer, parameter :: maxlin = 1000

      real(dp), pointer, save, public :: bk(:,:)
      logical, public , save          :: gamma_bands

      character, save :: label(maxlin)*8
      integer,   save :: lastk(maxlin), nlines

      public  :: setup_bands, bands, reset_bands
      
      contains

      subroutine setup_bands( )
        implicit none
        ! Find number of band k-points
        maxbk = 1
        nullify(bk)
        call re_alloc( bk, 1, 3, 1, maxbk, 'bk', 'band',
     .                copy=.false. )

        call initbands( maxbk, nbk, bk )

        if (nbk .gt. maxbk) then
          ! If there wasn't enough space to store bands on first call 
          ! correct the dimensions and repeat the initialisation
          maxbk = max(nbk,1)
          call re_alloc( bk, 1, 3, 1, maxbk, 'bk', 'band',
     .                  copy=.false. )
          call initbands( maxbk, nbk, bk )
        endif
        
        gamma_bands = (nbk == 0)   ! Flag for use of auxiliary cell

      end subroutine setup_bands

      subroutine reset_bands( )
      implicit none
      if (associated(bk)) call de_alloc( bk, 'bk', 'band' )
      end subroutine reset_bands

      subroutine initbands( maxk, nk, kpoint )
C *********************************************************************
C Finds band k-points
C Based on initialisation part of subroutine bands from original code.
C Written by J.Soler, August 1997 and August 1998.
C BandPoints option by EA, November 2006.
C **************************** INPUT **********************************
C integer maxk           : Last dimension of kpoint
C *************************** OUTPUT **********************************
C integer nk             : Number of band k points
C real*8  kpoint(3,maxk) : k point vectors
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C ***************** BEHAVIOUR *****************************************
C - If nk=0 on input, k-points are read from labels BandLines or
C   BandPoints, and BandLinesScale of the input fdf data file. 
C   If these labels are not present, it returns with nk=0.
C - BandLines for band-like representation, BandPoints for arbitrary
C   k points. BandPoints supersede BandLines.
C - Allowed values for BandLinesScale are ReciprocalLatticeVectors and
C   pi/a (default). If another value is given, it returns with nk=0
C   after printing a warning.
C - If nk>maxk, k points and bands are not calculated and no warning
C   is printed before return
C - If BandPoints, nlines is given a negative value, for routine bands
C   to do the appropriate writing
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
C
C  Modules
C
      use precision
      use parallel,     only : Node
      use fdf
      use m_get_kpoints_scale, only: get_kpoints_scale
      
      implicit          none

      integer, intent(in) :: maxk
      integer, intent(out) :: nk
      real(dp), pointer ::  kpoint(:,:)

C *********************************************************************

      logical BandLinesPresent, BandPointsPresent
      logical :: Nlines_overflow

      integer
     .  ik, il, ix, ierr,
     .  ni, nkl, nn, nr, nv

      real(dp) caux(3,3), rcell(3,3)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Start time counter 
      call timer( 'bands', 1 )

C Initialise the number of band lines
      nlines = 0

C Find if there are band-lines data
        BandLinesPresent = fdf_defined('BandLines')
        BandPointsPresent = fdf_defined('BandPoints')

        if ( BandLinesPresent .or. BandPointsPresent ) then

           call get_kpoints_scale('BandLinesScale',rcell,ierr)

           if (ierr /= 0) then
              nk = 0
              goto 999
           endif
        endif
!
        if ( BandLinesPresent .and. BandPointsPresent ) then
           if (Node .eq. 0) then
            write(6,'(a)')
     .        'bands: WARNING: BandPoints block supersedes BandLines'
            endif
            BandLinesPresent = .false.
         endif

C For BandLines: Loop on data lines
CEA
        if ( BandLinesPresent) then

          il = 1
          nk = 0
          Nlines_overflow = .FALSE.
          BandLinesPresent = fdf_block('BandLines',bfdf)

          do while(fdf_bline(bfdf, pline))

C Read and parse data line
            nn = fdf_bnnames(pline)
            nv = fdf_bnvalues(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

C Add to total number of k points
            nkl = fdf_bintegers(pline,1)
            nk = nk + nkl

! If there is room, store k points, otherwise skip this section
! and wait for a second call after proper allocation

            if (nk .le. maxk) then

C Find last point in line
               do ix = 1,3
                  kpoint(ix,nk) = rcell(ix,1) * fdf_bvalues(pline,2) +
     .                            rcell(ix,2) * fdf_bvalues(pline,3) +
     .                            rcell(ix,3) * fdf_bvalues(pline,4)
               enddo

C Find points along the line
              do ik = 1,nkl-1
                do ix = 1,3
                  kpoint(ix,nk-nkl+ik) =
     .              kpoint(ix,nk-nkl) * dble(nkl-ik) / dble(nkl) + 
     .              kpoint(ix,nk)     * dble(ik)     / dble(nkl)
                enddo
              enddo

C Find point label
              if (nn .gt. 0) then
                label(il) = fdf_bnames(pline,1)
              else
                label(il) = ' '
              endif
              lastk(il) = nk

            endif

            il = il + 1
            if (il .gt. maxlin) then
              ! No room to store bands --- hard coded limit...
              Nlines_overflow = .TRUE.
            endif
          enddo

          if ((Nlines_overflow) .and. (Node .eq. 0)) then
            write(6,'(a)') 'bands: ERROR. Parameter maxlin too small'
          endif
          nlines = il - 1
        endif

C For BandPoints: Loop on data points
CEA
        if ( BandPointsPresent ) then

          il = 0
          nk = 0
          Nlines_overflow = .FALSE.
          BandPointsPresent = fdf_block('BandPoints',bfdf)
          do while(fdf_bline(bfdf, pline))

C Read and parse data line
            nn = fdf_bnnames(pline)
            nv = fdf_bnvalues(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

C Add to total number of k points
            nk = nk + 1

! If there is room, store k points, otherwise skip this section
! and wait for a second call after proper allocation
            if (nk .le. maxk) then

C Accumulate new point
                do ix = 1,3
                  kpoint(ix,nk) = rcell(ix,1) * fdf_bvalues(pline,1) +
     .                            rcell(ix,2) * fdf_bvalues(pline,2) +
     .                            rcell(ix,3) * fdf_bvalues(pline,3)
                enddo
            endif

            il = il + 1
            if (il .gt. maxlin) then
            !  No room to store bands
              Nlines_overflow = .TRUE.
            endif
          enddo

          if ((Nlines_overflow) .and. (Node .eq. 0)) then
            write(6,'(a)') 'bands: ERROR. Parameter maxlin too small'
          endif
          nlines = - 3
        endif

C This is the only exit point 
  999 continue
      call timer( 'bands', 2 )

      end subroutine initbands


      subroutine bands( no_s, spin,
     .                  no_u, no_l,maxnh, maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  writeb, nk, kpoint, ek, occtol, getPSI )

C Finds band energies at selected k-points.
C Optionally write WFS to file
C Written by J.Soler, August 1997 and August 1998.
C Initialisation moved into a separate routine, JDG Jan 2000.
C WFS options by A. Garcia, April 2012
C Cleaned for tSpin by N. Papior, August 2019

C **************************** INPUT **********************************
C integer no_s                : Number of basis orbitals in supercell
C type(tSpin) spin            : Containing all spin-information
C integer maxnh               : Maximum number of orbitals interacting
C                               with any orbital
C integer maxk                : Last dimension of kpoint and ek
C integer numh(no_l)          : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(no_l)      : Pointer to start of each row of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row
C real*8  H(maxnh,spin%H)     : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  ef                  : Fermi energy
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer no_u                : First dimension of ek
C integer no_l                : Second dimension of H and S
C integer indxuo(no_s)        : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.no_l, with
C                               no_l the number of orbitals in unit cell
C real*8  ef                  : Fermi energy
C logical writeb              : This routine must write bands?
C integer no_u                : Total number of orbitals in unit cell
C integer nk                  : Number of band k points
C real*8  kpoint(3,maxk)      : k point vectors
C real*8  occtol              : Occupancy threshold for DM build
C *************************** OUTPUT **********************************
C real*8  ek(no_u,spin%spinor,maxk) : Eigenvalues
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C ***************** BEHAVIOUR *****************************************
C - When writeb=true, bands are saved in file sys_name.bands, where
C   sys_name is the value of fdf label SystemLabel, or 'siesta'
C   by default.
C - When nlines is negative, eigenvalues are written for the
C   preselected k points (as introduced through BandPoints).
C   Otherwise written as lines for band-like representation (BandLines)
C *********************************************************************
C
C  Modules
C
      use precision,    only : dp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs
      use densematrix,  only : allocDenseMatrix, resetDenseMatrix
      use densematrix,  only : Haux, Saux, psi
      use alloc,        only : re_alloc, de_alloc
      use files,        only : slabel, label_length
      use atomlist,     only : iaorb, iphorb
      use siesta_geom,  only : isa
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use writewave,    only : wfs_filename

      use t_spin, only: tSpin

      use m_diag_option, only: ParallelOverK, Serial
      use m_diag, only: diag_init


      implicit          none

      type(tSpin), intent(in) :: spin
      integer  ::   maxk, maxnh, no_u, no_l, nk, no_s,
     .              indxuo(no_s), listh(maxnh), numh(no_l),
     .              listhptr(*)
      logical  ::   writeb
      real(dp) ::   ef, ek(no_u,spin%spinor,maxk),
     .              H(maxnh,spin%H), kpoint(3,maxk),
     .              S(maxnh), xij(3,maxnh), occtol
      logical, intent(in) :: getPSI

      external ::       io_assign, io_close, memory


      character(len=label_length+6) :: fname
      character(len=10)             :: string

      logical, parameter :: getD = .false.
      logical, parameter :: fixspin = .false.

      integer :: ik, il, io, ispin, iu, iu_wfs, iuo, naux, nhs, j
      logical :: SaveParallelOverK
      
      real(dp)
     .  Dnew, qs(2), e1, e2, efs(2), emax, emin, Enew, eV, qk, qtot,
     .  path, temp, wk, Entropy

C Dynamic arrays
      real(dp), dimension(:), pointer :: aux => null()

      parameter ( eV = 1.d0 / 13.60580d0 )

      save Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

C Start time counter
      call timer( 'bands', 1 )

C Check parameter maxk
      if (nk .gt. maxk) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'bands: WARNING: parameter maxk too small',
     .       'bands: No bands calculation performed'
        endif
        goto 999
      endif

C     Allocate local arrays - only aux is relevant here

      if ( spin%Grid == 4 ) then
         nhs = 2 * (2*no_u) * (2*no_l)
      else
         nhs = 2 * no_u*no_l
      endif
      naux = 2*no_u*5
      call allocDenseMatrix(nhs, nhs, nhs)
      call re_alloc( aux, 1, naux, 'aux', 'bands' )

      if (getPSI .AND. (Node.eq.0)) then

        call io_assign( iu_wfs )
        open(iu_wfs,file=wfs_filename,
     $              form="unformatted",status='unknown')

        rewind (iu_wfs)

        write(iu_wfs) nk, .false. ! nk, Gamma, same file-format in WFS as for Gamma-point
        write(iu_wfs) spin%Grid
        write(iu_wfs) no_u
        write(iu_wfs) (iaorb(j),labelfis(isa(iaorb(j))),
     .            iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
     .            symfio(isa(iaorb(j)),iphorb(j)), j=1,no_u)

        call io_close(iu_wfs)
      endif
      
C Handle parallel over K points option which is not allowed for here
      SaveParallelOverK = ParallelOverK
      if ( ParallelOverK ) then
        if ( Node .eq. 0 ) then
          write(*,"(a)")
     $        "*** Note: ParallelOverK option not used for Bands"
        end if
        ParallelOverK = .false.
        Serial = Nodes == 1
        call diag_init()
      end if

C Find the band energies
      if ( spin%none .or. spin%Col ) then
C fixspin and qs are not used in diagk, since getD=.false. ...
        qs(1) = 0.0_dp
        qs(2) = 0.0_dp

        call diagk( spin%H, no_l, no_s, spin%spinor, 
     .              maxnh, maxnh, 
     .              no_u, numh, listhptr, listh, numh, listhptr,
     .              listh, H, S, getD, getPSI, fixspin, qtot, qs, temp,
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              ek, qk, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux,
     .              no_u, occtol, 1, no_u )

      elseif ( spin%NCol ) then
        call diag2k(no_l, no_s, maxnh, maxnh, no_u,
     .              numh, listhptr, listh, numh, listhptr,
     .              listh, H, S, getD, getPSI, qtot, temp, e1, e2,
     .              xij, indxuo, nk, kpoint, wk,
     .              ek, qk, Dnew, Enew, ef, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux,
     .              no_u, occtol, 1, no_u )

      elseif ( spin%SO ) then
        call diag3k(no_l, no_s, maxnh, maxnh, no_u, numh,
     .              listhptr, listh, numh, listhptr, listh,
     .              H, S, getD, getPSI, qtot, temp, e1, e2, xij,
     .              indxuo, nk, kpoint, wk, ek, qk, Dnew,
     .              Enew, ef, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux,
     .              no_u, occtol, 1, no_u )
      else
        call die( 'bands: ERROR: incorrect value of spin%H')
      endif

      ! Restore parallel over k
      ParallelOverK = SaveParallelOverK
      if ( ParallelOverK ) Serial = .true.

C Write bands
      if (writeb.and.Node.eq.0) then

C Find name of output file and open it
        fname = trim(slabel) // '.bands'
        call io_assign(iu)
        open( iu, file=fname, status='unknown')

C Write Fermi energy
        write(iu,*) Ef/eV

C Find and write the ranges of k and ek
        path = 0.d0
        emax = ek(1,1,1)
        emin = ek(1,1,1)
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          do ispin = 1, spin%spinor
            do io = 1, no_u
              emax = max( emax, ek(io,ispin,ik) )
              emin = min( emin, ek(io,ispin,ik) )
            enddo
          enddo
        enddo
C Write path for BandLines but not for BandPoints
        if (nlines .ge. 0) write(iu,*) 0.d0, path
        write(iu,*) emin/eV, emax/eV

C Write eigenvalues
        if ( spin%Grid == 4 ) then
          write(iu,*) 2*no_u, 1, nk     ! A single spin channel, double number of bands
        else 
          write(iu,*) no_u, spin%spinor, nk
        endif

        ! For non-collinear calculations, the ek array will contain the first no_u
        ! bands in ek(:,1,ik), and the next no_u bands in ek(:,2,ik)   
        path = 0.d0
        do ik = 1,nk
          if (nlines .ge. 0) then
            if (ik .gt. 1)
     .        path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                            (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                            (kpoint(3,ik)-kpoint(3,ik-1))**2 )
            write(iu,'(f10.6,10f12.4,/,(10x,10f12.4))')
     .         path, ((ek(io,ispin,ik)/eV,io=1,no_u),
     .         ispin=1,spin%spinor)
          else
            write(iu,'(3f9.5,8f12.4,/,(27x,8f12.4))')
     .        kpoint(1,ik), kpoint(2,ik), kpoint(3,ik),
     .        ((ek(io,ispin,ik)/eV,io=1,no_u),ispin=1,spin%spinor)
          endif
        enddo

C Write abscisas of line ends and their labels
        if ( nlines .ge. 0 ) then
          write(iu,*) nlines
          il = 1
          path = 0.d0
          do ik = 1,nk
            if (ik .gt. 1)
     .        path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                            (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                            (kpoint(3,ik)-kpoint(3,ik-1))**2 )
            if (ik .eq. lastk(il)) then
C Put label between quotes
              if (label(il) .eq. ' ') then
                string = "' '"
              else
                string = "'"//trim(label(il))//"'"
              endif
              write(iu,'(f12.6,3x,a)') path, string
              il = il + 1
            endif
          enddo
        endif

C Close output file
        call io_close(iu)
      endif

C Free local arrays
      call de_alloc( aux, 'aux', 'bands' )
      call resetDenseMatrix()
      
C This is the only exit point
  999 continue
      call timer( 'bands', 2 )

      end subroutine bands

      end module band
