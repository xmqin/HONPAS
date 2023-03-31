! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module writewave

      USE precision, only : dp
      use alloc,     only : re_alloc, de_alloc
      use sys, only: die, message
      implicit none

      private
      character(len=192), save, public :: wfs_filename

      integer,    save :: nwf = 1
      logical,    save :: wwf = .false.
      logical,    save :: scf_set = .false.
      logical,    save :: debug = .false.


      integer, save, dimension(:), pointer   :: nwflist => null()
      integer, save, dimension(:,:), pointer :: iwf => null()

      integer, public, save           :: nwk
      real(dp), public, pointer, save :: wfk(:,:) => null()
      logical, public, save           :: gamma_wavefunctions

      real(dp),   save :: wfs_energy_min = -huge(1.0_dp)
      real(dp),   save :: wfs_energy_max =  huge(1.0_dp)
      logical,    save :: wfs_energy_window  = .true.

      public :: setup_wfs_list, setup_wf_kpoints, wwave, writew 

      CONTAINS

      subroutine setup_wf_kpoints()
      USE alloc, only: re_alloc
      USE sys, only: die
      USE atomlist, only: no_u
      USE m_spin, only: spin

      implicit none

      integer :: max_n_states
      
      ! Find number of k-points for wavefunction printout

      nullify(wfk)

      nwk = 0
      if ( spin%Grid == 4 ) then
        max_n_states = 2 * no_u
      else
        max_n_states = no_u
      end if
      call initwave( max_n_states, nwk, wfk )

      if (nwk .eq. 0) then
         gamma_wavefunctions = .true.
      else if (nwk .eq. 1) then
         if (dot_product(wfk(:,1),wfk(:,1)) < 1.0e-20_dp) then
            gamma_wavefunctions = .true.
         else
            gamma_wavefunctions = .false.
         endif
      else  
         gamma_wavefunctions = .false.
      endif

      end subroutine setup_wf_kpoints

      subroutine setup_wfs_list(nk,norb,wfs_band_min,wfs_band_max,
     $                          use_scf_weights,use_energy_window)
      use fdf
      !
      ! Initialize structures if the number of bands
      ! to output is the same for all k-points
      ! (e.g., for the SCF k-point set, or for the new 'bands' option

      integer, intent(in) :: nk, norb
      integer, intent(in) :: wfs_band_min, wfs_band_max  ! Band range
      logical, intent(in) :: use_scf_weights            
      logical, intent(in) :: use_energy_window           

      integer :: j, nwfs

      call re_alloc(nwflist, 1, nk, name="nwflist",
     $     routine="setup_wfs_list")
      call re_alloc(iwf, 1, nk, 1, norb, name="iwf",
     $     routine="setup_wfs_list")

      scf_set = use_scf_weights
      wfs_energy_window = use_energy_window
!
      debug = fdf_get("WriteWaveDebug", .false.)

!     Here we request output of all the eigenstates in a k-point
!     except if the user selects a band or energy range.
!
      if (wfs_band_min < 1) then
         call die("WFS.BandMin < 1")
      endif
      if (wfs_band_max > norb) then
         call die("WFS.BandMax > maximum number of states")
      endif
      if (wfs_band_max < wfs_band_min) then
         call message("WARNING","WFS.BandMax < WFS.BandMin." //
     $        " Nothing written to file")
      endif

      nwfs = max(0,wfs_band_max - wfs_band_min + 1)
      nwflist(1:nk) = nwfs
      do j=1, nwfs
         iwf(1:nk,j) = wfs_band_min - 1 + j
      enddo

!     The energy window, if allowed, will take effect inside the routine
!     that writes to file, as the energies are not yet known

      wfs_energy_min = fdf_get("WFS.EnergyMin",
     $                                -huge(1.0_dp),"Ry")
      wfs_energy_max = fdf_get("WFS.EnergyMax",
     $                                 huge(1.0_dp),"Ry")

      end subroutine setup_wfs_list


      subroutine initwave( norb, nk, kpoint )
C *********************************************************************
C Finds k-points for wavefunction printout
C Based on initband routine by J.Soler
C Written by P. Ordejon, June 2003
C **************************** INPUT **********************************
C integer norb           : Number of orbitals (actually max number of bands)
C *************************** OUTPUT **********************************
C integer nk             : Number k points to compute wavefunctions
C real*8  kpoint(3,maxk) : k point vectors
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C ***************** BEHAVIOUR *****************************************
C - If nk=0 on input, k-points are read from labels WaveFuncKPoints and
C   WaveFuncKPointsScale from the input fdf data file. If these labels
C   are not present, it returns with nk=0.
C - Allowed values for WaveFuncKPointsScale are ReciprocalLatticeVectors
C   and pi/a (default). If another value is given, it returns with nk=0
C   after printing a warning.
C - If nk>maxk, k-points and wavefunctions are not calculated and no
C   warning is printed before return
C ***************** USAGE *********************************************
C Example of fdf wavefunction k-points specification for an FCC lattice.
C
C     WaveFuncKPointsScale  pi/a
C     %block WaveFuncKPoints              # These are comments
C     0.000  0.000  0.000  from 1 to 10   # eigenstates 1-10 of Gamma
C     2.000  0.000  0.000  1 3 5          # eigenstates 1,3,5 of X
C     1.500  1.500  1.500                 # all eigenstates of K
C     %endblock WaveFuncKPoints
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
      use sys,          only : die
      use alloc,        only : re_alloc
      use m_get_kpoints_scale, only: get_kpoints_scale
      
      implicit          none

      integer           :: norb, nk
      real(dp), pointer :: kpoint(:,:)
      external          :: memory
C *********************************************************************

      character
     .  name1*10, name2*10, name3*10,
     .  msg*120

      logical       :: outlng
      logical       :: WaveFuncPresent

      integer
     .  i, ix, iw, iw1, iw2, iw3, ierr,
     .  ni, nn, nr, nv

      real(dp)
     .  caux(3,3), 
     .  rcell(3,3)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      wfs_energy_window = .true.
C Find if there are k-points data
      outlng = fdf_boolean('LongOutput', .false.)
      wwf = fdf_boolean('WriteWaveFunctions',outlng)

C Find if there are k-points data
      WaveFuncPresent = fdf_block('WaveFuncKPoints',bfdf)
      if ( WaveFuncPresent ) then

         call get_kpoints_scale('WaveFuncKPointsScale',rcell,ierr)

         if (ierr /= 0) then
            nk = 0
            RETURN
         endif

C Find max number of k points and max number of bands per k-point
C Loop on data lines
        nk = 0
        nwf = 0
        do while(fdf_bline(bfdf,pline))
C Read and parse data line
          nn = fdf_bnnames(pline)
          nv = fdf_bnvalues(pline)
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)

C Check if syntax line is correct
          if (nv .ge. 3) then

C Check syntax
            if (nr .ne. 3) then
              write(msg,'(a,/,a)')
     .          'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .          trim(fdf_getline(bfdf%mark))
              call die(msg)
            endif

C Add this point to total number of k points
            nk = nk + 1

C Find which eigenvectors should be printed
C Store indexes of wave functions to printout

C #X.XXX #Y.YYY #Z.ZZZ from #start to #end [step] #step
            if (nn .ge. 1) then

C Check that line contains 'from', 'to' and maybe 'step'
              if ((nn .ne. 2) .and. (nn .ne. 3)) then
                write(msg,'(a,/,a)')
     .            'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .            trim(fdf_getline(bfdf%mark))
                call die(msg)
              endif
              name1 = fdf_bnames(pline,1)
              name2 = fdf_bnames(pline,2)
              if (nn .eq. 3) name3 = fdf_bnames(pline,3)
              if (name1 .ne. 'from' .or. name2 .ne. 'to') then
                write(msg,'(a,/,a)')
     .            'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .            trim(fdf_getline(bfdf%mark))
                call die(msg)
              endif
              if (nn .eq. 3) then
                if (name3 .ne. 'step') then
                  write(msg,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .             trim(fdf_getline(bfdf%mark))
                  call die(msg)
                endif
              endif

              iw1 = fdf_bintegers(pline,1)
              iw2 = fdf_bintegers(pline,2)
              if (iw1 .lt. 0) iw1 = norb + iw1 + 1
              if (iw2 .lt. 0) iw2 = norb + iw2 + 1
              if (nn .eq. 3) then
                iw3 = abs(fdf_bintegers(pline,3))
              else
                iw3 = 1
              endif
              ni = 0
              do iw = min(iw1,iw2), max(iw1,iw2), iw3
                ni = ni + 1
              enddo
              nwf = max(nwf,ni)

C #X.XXX #Y.YYY #Z.ZZZ #eigen1 [#eigen2 ...]
            elseif (ni .ne. 0) then
              nwf = max(nwf,ni)

C #X.XXX #Y.YYY #Z.ZZZ
            elseif (ni .eq. 0) then
              nwf = max(nwf,norb)
            endif
          else
C Bad syntax
            write(msg,'(a,/,a)')
     .        'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .        trim(fdf_getline(bfdf%mark))
            call die(msg)
          endif
        enddo

        ! Minimize to the maximum number of eigenvalues
        nwf = min(nwf, norb)
        
C Allocate nwflist, iwf and kpoints structures according to nk and nwf
        call re_alloc( nwflist, 1, nk, 'nwflist', 'initwave' )
        call re_alloc( iwf, 1, nk, 1, nwf, 'iwf', 'initwave' )

        call re_alloc( kpoint, 1, 3, 1, nk, 'kpoint', 'initwave' )

C Fill nwflist, iwf and kpoints structures
C Loop on data lines
        nk = 0
        call fdf_brewind(bfdf)
        do while(fdf_bline(bfdf,pline))
C Read and parse data line
          nn = fdf_bnnames(pline)
          nv = fdf_bnvalues(pline)
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)

C Add this point to total number of k points
          nk = nk + 1

C Find coordinates of k point

            do ix = 1,3
              kpoint(ix,nk) = rcell(ix,1) * fdf_bvalues(pline,1) +
     .                        rcell(ix,2) * fdf_bvalues(pline,2) +
     .                        rcell(ix,3) * fdf_bvalues(pline,3)
            enddo

C Find which eigenvectors should be printed
C Store indexes of wave functions to printout

C #X.XXX #Y.YYY #Z.ZZZ from #start to #end [step] #step
          if (nn .ge. 2) then

            name1 = fdf_bnames(pline,1)
            name2 = fdf_bnames(pline,2)
            if (nn .eq. 3) name3 = fdf_bnames(pline,3)

            iw1 = fdf_bintegers(pline,1)
            iw2 = fdf_bintegers(pline,2)
            if (iw1 .lt. 0) iw1 = norb + iw1 + 1
            if (iw2 .lt. 0) iw2 = norb + iw2 + 1
            ! Reduce to maximally the number of orbitals
            iw2 = min(iw2, norb)
            if (nn .eq. 3) then
              iw3 = abs(fdf_bintegers(pline,3))
            else
              iw3 = 1
            endif
            ni = 0
            do iw = min(iw1,iw2), max(iw1,iw2), iw3
              ni = ni + 1
              if (iw .lt. 0) then
                iwf(nk,ni) = norb + iw + 1
              else
                iwf(nk,ni) = iw
              endif
            enddo

C #X.XXX #Y.YYY #Z.ZZZ #eigen1 [#eigen2 ...]
          elseif (ni .ne. 0) then
            do i= 1, ni
              iw = fdf_bintegers(pline,i)
              if (iw .lt. 0) then
                iw = norb + iw + 1
              endif
              if ( iw <= norb ) then
                iwf(nk,i) = iw
              end if
            enddo

C #X.XXX #Y.YYY #Z.ZZZ
          elseif (ni .eq. 0) then
            ni = norb
            do i= 1, ni
              iwf(nk,i) = i
            enddo
          endif

          nwflist(nk) = ni
        enddo

      endif

      end subroutine initwave

      subroutine wwave( no, spin, maxo, maxuo, maxnh, 
     .                  maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  gamma, nk, kpoint, nuotot, occtol)
C *********************************************************************
C Finds wavefunctions at selected k-points.
C Written by P. Ordejon, June 2003
C from routine 'bands' written by J.M.Soler
C Changed to use tSpin, N. Papior, August 2019
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C tSpin spin                  : Contains spin information
C integer maxo                : First dimension of ek
C integer maxuo               : Second dimension of H and S
C integer maxnh               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxk                : Last dimension of kpoint and ek
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to start of each row of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,spin%H)     : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  ef                  : Fermi energy
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C logical Gamma               : whether only the Gamma-point is sampled
C integer nk                  : Number of band k points
C real*8  kpoint(3,maxk)      : k point vectors
C integer nuotot              : Total number of orbitals in unit cell
C real*8  occtol              : Occupancy threshold for DM build
C *************************** OUTPUT **********************************
C None; output is dumped to wave functions file SystemLabel.WFSX
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C
C  Modules
C
      use precision
      use parallel,    only : Node, Nodes
      use fdf
      use densematrix, only : allocDenseMatrix, resetDenseMatrix
      use densematrix, only : Haux, Saux, psi
      use alloc
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use atomlist,     only : iaorb, iphorb
      use siesta_geom,  only : isa
      use t_spin, only: tSpin
      use units, only : eV

#ifdef MPI
      use parallelsubs, only : GetNodeOrbs
#endif
      use m_diag_option, only: Serial, ParallelOverK
      use m_diag, only: diag_init

      implicit          none

      type(tSpin), intent(in) :: spin
      logical, intent(in) :: Gamma
      integer           maxk, maxnh, maxo, maxuo, nk, no, 
     .                  nuotot, indxuo(no), listh(maxnh), numh(*), 
     .                  listhptr(*)
      real(dp)          ef, H(maxnh,spin%H), kpoint(3,maxk), 
     .                  S(maxnh), xij(3,maxnh), occtol

      external          io_assign, io_close, memory
C *********************************************************************

      logical
     .  fixspin

      integer
     .  io, iu, iuo, naux, nuo, j, nhs

      real(dp)
     .  Dnew, qs(2), e1, e2, efs(2), Enew, qk, qtot, 
     .  temp, wk, Entropy

C Dynamic arrays
      logical, parameter :: getD = .false.
      logical, parameter :: getPSI = .true.
      real(dp), allocatable :: aux(:)
      real(dp), allocatable :: ek(:,:,:)

      save Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

      logical :: SaveParallelOverK
      

C Get local number of orbitals
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

C Start time counter 
      call timer( 'writewave', 1 )

C Check parameter maxk 
      if (nk .gt. maxk) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'writewave: WARNING: parameter maxk too small',
     .       'writewave: No wavefunction calculation performed'
        endif
        goto 999
      endif

C Check spin

C Allocate local arrays - only aux is relevant here
      if ( spin%Grid == 4 ) then
         nhs = 2 * (2*nuotot) * (2*nuo)
      else
         nhs = 2 * nuotot * nuo
      endif
      call allocDenseMatrix(nhs, nhs, nhs) 

      allocate(ek(nuotot,spin%spinor,nk))
      call memory('A','D',nuotot*spin%spinor*nk,'writewave')
      naux = 2*nuotot*5
      allocate(aux(naux))
      call memory('A','D',naux,'writewave')

C Open file

      if (Node.eq.0) then

        call io_assign( iu )
        open(iu, file=wfs_filename,form="unformatted",status='unknown')

        rewind (iu)

        if (wwf) then
          write(6,*)
          write(6,'(a)') 'writewave: Wave Functions Coefficients'
          write(6,*)
          write(6,'(a26,2x,i6)') 'Number of k-points = ', nk
          write(6,'(a26,2x,i6)') 'Number of Spins = ', spin%H
          write(6,'(a26,2x,i6)') 'Number of basis orbs = ',nuotot
          write(6,*)
        endif
        write(iu) nk, gamma
        write(iu) spin%H   ! This will be 1, 2, 4 or 8 (to tell NC from SOC)
        write(iu) nuotot
        write(iu) (iaorb(j),labelfis(isa(iaorb(j))),
     .            iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
     .            symfio(isa(iaorb(j)),iphorb(j)), j=1,nuotot)

        call io_close(iu)

      endif

C Find the eigenvectors
c fixspin and qs are not used in diagk, since getD=.false. ...
      fixspin = .false.
      qs(1) = 0.0d0
      qs(2) = 0.0d0

C Handle parallel over K points option which is not allowed for here
      SaveParallelOverK = ParallelOverK
      if ( ParallelOverK ) then
        if ( Node .eq. 0 ) then
          write(*,"(a)")
     $        "*** Note: ParallelOverK option not used for "//
     &        "WF info generation"
        end if
        ParallelOverK = .false.
        Serial = (Nodes == 1)
        call diag_init()
      end if

C Call appropriate diagonalization routine

      if ( spin%none .or. spin%Col ) then
       if (gamma) then
         call diagg( spin%H, nuo, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, nuotot, occtol, 1, nuotot )
       else
         call diagk( spin%H, nuo, no, spin%spinor, maxnh, maxnh,
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, xij, indxuo, nk, kpoint, wk,
     .               ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, Haux, Saux, aux, nuotot,
     .               occtol, 1, nuotot )
       endif

      elseif ( spin%NCol) then
        if (gamma) then
          call diag2g( nuo, no, maxnh, maxnh, maxo, numh,
     .                 listhptr, listh, numh, listhptr, listh,
     .                 H, S, getD, getPSI, qtot, temp, e1, e2, ek, qk,
     .                 Dnew, Enew, ef, Entropy,
     .                 Haux, Saux, psi, aux,
     .                 nuotot, occtol, 1, nuotot )
        else
          call diag2k( nuo, no, maxnh, maxnh, maxo, numh,
     .                 listhptr, listh, numh, listhptr, listh,
     .                 H, S, getD, getPSI, qtot, temp, e1, e2, xij,
     .                 indxuo, nk, kpoint, wk, ek, qk, Dnew,
     .                 Enew, ef, Entropy,
     .                 Haux, Saux, psi, Haux, Saux, aux,
     .                 nuotot, occtol, 1, nuotot )
        endif
      elseif ( spin%SO ) then
        if (gamma) then
          call diag3g( nuo, no, maxnh, maxnh, maxo, numh,
     .                 listhptr, listh, numh, listhptr, listh,
     .                 H, S, getD, getPSI, qtot, temp, e1, e2, ek, qk,
     .                 Dnew, Enew, ef, Entropy,
     .                 Haux, Saux, psi, aux,
     .                 nuotot, occtol, 1, nuotot )
        else
          call diag3k( nuo, no, maxnh, maxnh, maxo, numh,
     .                 listhptr, listh, numh, listhptr, listh,
     .                 H, S, getD, getPsi, qtot, temp, e1, e2, xij,
     .                 indxuo, nk, kpoint, wk, ek, qk, Dnew,
     .                 Enew, ef, Entropy,
     .                 Haux, Saux, psi, Haux, Saux, aux,
     .                 nuotot, occtol, 1, nuotot )
        endif
      endif

      ! Restore ParallelOverK
      ParallelOverK = SaveParallelOverK
      if ( ParallelOverK ) Serial = .true.

      call resetDenseMatrix()
C Free local arrays 
      call memory('D','D',size(aux),'writewave')
      deallocate(aux)
      call memory('D','D',size(ek),'writewave')
      deallocate(ek)

C This is the only exit point 
  999 continue
      call timer( 'writewave', 2 )

      end subroutine wwave


      subroutine writew(nuotot,nuo,ik,k,ispin,eo,psi,
     $                  gamma,non_coll,blocksize)

      use precision
      use sys
      use parallel,     only : Node, Nodes
      use parallel,     only : Module_Blocksize=>BlockSize
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use units,        only : eV
      use kpoint_grid,  only : kweight
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

#ifdef MPI
      integer  MPIerror, mpistatus(MPI_STATUS_SIZE), tag
#endif

      integer, intent(in)  ::  nuotot, nuo, ispin, ik
      logical, intent(in)  ::  non_coll, gamma
      real(dp), intent(in) ::  eo(*), k(3)
      real(dp), intent(in) ::  psi(:)
      integer, intent(in)  ::  blocksize


C  Internal variables .............................................
      integer  BNode, iie, iw, indwf, j, ind, iu
      integer  :: blocksize_saved
      integer number_of_wfns
      integer, allocatable  :: ind_wfn(:)
      real(dp)  :: kpoint_weight
      real(SP), dimension(:,:), allocatable :: aux   !! NOTE SP

      external io_assign, io_close

C ...................

C Allocate auxiliary arrays

      if (non_coll) then
         allocate(aux(4,nuotot))
         call memory('A','D',4*nuotot,'writewave')
      else
         if (gamma) then
            allocate(aux(1,nuotot))
            call memory('A','D',nuotot,'writewave')
         else
            allocate(aux(2,nuotot))
            call memory('A','D',2*nuotot,'writewave')
         endif
      endif

C Find file name, and open for Node 0

      if (Node .eq. 0) then
        call io_assign( iu )
        open(iu,file=wfs_filename,form="unformatted",position='append',
     .        status='old')
      endif

      if (wfs_energy_window) then  
         allocate(ind_wfn(nwflist(ik)))
         number_of_wfns = 0
         do iw = 1,nwflist(ik)
            indwf = iwf(ik,iw)
            if (eo(indwf) >= wfs_energy_min .AND.
     $          eo(indwf) <= wfs_energy_max ) then
               number_of_wfns = number_of_wfns + 1
               ind_wfn(number_of_wfns) = indwf
            endif
         enddo
         ! Replace indexes with new values
         nwflist(ik) = number_of_wfns
         do iw = 1, nwflist(ik)
            iwf(ik,iw) = ind_wfn(iw)
         enddo
         deallocate(ind_wfn)
      endif

C First print the index and value of k-point

      if (Node .eq. 0) then
        if (wwf) then
          write(6,*)
          write(6,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik, k(1:3)
          if (non_coll) then
             write(6,'(a22,2x)') 'Spinors have mixed spins'
          else
             write(6,'(a22,2x,i6)') 'Spin component = ',ispin
          endif
          write(6,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist(ik)
          write(6,'(a)') 'Use readwfsx utility to print '
     $                   // "wavefunction coefficients from WFSX file"
        endif
        if (scf_set) then
           kpoint_weight = kweight(ik)
        else
           kpoint_weight = 1.0_dp
        endif
        write(iu) ik, k(1:3), kpoint_weight
        write(iu) ispin
        write(iu) nwflist(ik)
      endif

      ! Loop over wavefunctions that should be stored

      ! The effective blocksize to use is passed as an argument.
      ! Module_blocksize is the variable in the 'parallel' module that
      ! holds the value that all the helper functions (WhichNodeOrb,
      ! GetNodeOrb, etc) see.  Hence, we must temporarily set
      ! Module_blocksize to the passed blocksize.

      ! Examples: non-collinear calculations need a wf blocksize double
      ! the 'orbital' one

      blocksize_saved = Module_blocksize
      Module_blocksize = blocksize
      
      if (Node .eq. 0 .and. debug) then
         print *, " -*- Will print ", nwflist(ik), " wfns."
         print *, " -*- ", iwf(ik,1:nwflist(ik))
      endif
      do iw = 1,nwflist(ik)
         indwf = iwf(ik,iw)

        ! Determine which node handles this wavefunction
        ! Note that the distribution is block cyclic,
        ! so the same 'orbital' helper functions apply

        call WhichNodeOrb(indwf,Nodes,BNode)

        if (Node .eq. 0 .and. debug) then
           print *, " -*- About to print wfn ", indwf,
     $                    "from node ", Bnode
        endif

        if (Node.eq.BNode) then

          ! Determine the index of the wavefunction in the local node
          call GlobalToLocalOrb( indwf, BNode, Nodes, iie)

          ! Save wavefunction in aux array

          ! Despite being passed as a one-dimensional array, psi has
          ! different underlying structures for different cases, so the
          ! indexing must be handled differently

          if (.not. non_coll) then
             if (gamma) then
                ! Real functions
                do j = 1,nuotot
                   ind = j + (iie-1)*nuotot 
                   aux(1,j) = real(psi(ind),kind=sp)
                enddo
             else
                ! Complex functions: two reals per orbital coefficient
                do j = 1,nuotot
                   ind = 1+(j-1)*2+(iie-1)*2*nuotot
                   aux(1,j) = real(psi(ind),kind=sp)
                   aux(2,j) = real(psi(ind+1),kind=sp)
                enddo
             endif
          else                  ! non-collinear case, including SOC
             ! Each 'orbital' coefficient is actually four reals: complex up/down
             ! (spinor)
             do j = 1,nuotot
                ind = 1 + (j-1)*4 + (iie-1)*4*nuotot 
                aux(1,j) = real(psi(ind),kind=sp)
                aux(2,j) = real(psi(ind+1),kind=sp)
                aux(3,j) = real(psi(ind+2),kind=sp)
                aux(4,j) = real(psi(ind+3),kind=sp)
             enddo
          endif  !  non-coll or not

       endif          ! Node == BNode
       
#ifdef MPI
C Pass the wf to the master node
        if (BNode.ne.0) then
           tag = indwf
           if (Node .eq. BNode) then
              call MPI_Send(aux(1,1),size(aux),MPI_real,
     .             0,tag,MPI_Comm_World,MPIerror)
              if (debug) print *, Node, " sent msg to 0 with tag ",tag
           elseif (Node .eq. 0) then
              call MPI_Recv(aux(1,1),size(aux),MPI_real,
     .             BNode,tag,MPI_Comm_World,MPIStatus,MPIerror)
              if (debug) print *, "Got msg with tag ",tag,
     $             " from node ", BNode, ". status: ", mpistatus
           endif
           call MPI_Barrier(MPI_Comm_World, MPIError) ! Just to be safe 
        endif
#endif

C eigenvector is now stored in aux in Node 0, and can be printed

        if (Node .eq. 0) then
           write(iu) indwf
           write(iu) eo(indwf)/eV
           write(iu) (aux(1:,j), j=1,nuotot)
        endif

      enddo
      
C Close output file

      if (Node .eq. 0) then
        call io_close(iu)
      endif

      ! Restore blocksize value in module
      Module_blocksize = blocksize_saved

C Deallocate auxiliary arrays

      call memory('D','D',size(aux),'writewave')
      deallocate(aux)
        
      end subroutine writew

      end module writewave

