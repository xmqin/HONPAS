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
      module writewave

      use precision, only: dp
      use sys, only: die, message
      use alloc, only: re_alloc
      use m_fdf_global, only: fdf_global_get

      implicit none

      integer :: maxlin
      parameter (maxlin=1000)

      character,  save :: label(maxlin)*8
      character(len=200), save, public :: wfs_filename
      integer,    save :: lastk(maxlin)

      integer,    save :: nwf = 1
      logical,    save :: wwf = .false.
      logical,    save :: scf_set = .false.
      logical,    save :: debug = .false.

      integer, save, dimension(:), pointer   :: nwflist => null()
      integer, save, dimension(:,:), pointer :: iwf => null()
 
      integer, public, save               :: nwk
      integer, public, save               :: maxwk
      real(dp), public, pointer, save     :: wfk(:,:) => null()
      logical, public, save               :: gamma_wavefunctions

      real(dp),   save :: wfs_energy_min = -huge(1.0_dp)
      real(dp),   save :: wfs_energy_max =  huge(1.0_dp)
      logical,    save :: wfs_energy_window  = .true.

      private 
      public :: setup_wfs_list, setup_wf_kpoints, wwave, writew

      CONTAINS

      subroutine setup_wf_kpoints()
      USE alloc, only: re_alloc
      USE sys, only: die
      USE atomlist, only: no_u

      implicit none

      ! Find number of k-points for wavefunction printout ....................

      logical  :: overflow

      nwk = 0
      maxwk = 1
      nullify(wfk)
      call re_alloc(wfk,1,3,1,maxwk,name='wfk',
     $              routine='setup_wf_kpoints')

      call initwave( maxwk, no_u, nwk, wfk, overflow )

      if (overflow) then
! If there wasn't enough space to store bands on first call correct
! the dimensions and repeat the initialisation        
         if (nwk .gt. maxwk) then
            maxwk = max(nwk,1)
            call re_alloc(wfk,1,3,1,maxwk,name='wfk',
     $           routine='setup_wf_kpoints')
         endif
         nwk = 0
         call initwave( maxwk, no_u, nwk, wfk, overflow )
         if (overflow)
     $      call die('siesta: ERROR: Unsuccessful initialization of' //
     $     ' list of wavefunctions to print')
      endif

      if (nwk == 0) then
         gamma_wavefunctions = .true.
      else if (nwk == 1) then
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
      call fdf_global_get(debug,"WriteWaveDebug", .false.)

!     Here we request output of all the eigenstates in a k-point
!     except if the user selects a band or energy range.
!
      if (wfs_band_min < 1) then
         call die("WFS.BandMin < 1")
      endif
      if (wfs_band_max > norb) then
         call die("WFS.BandMax > number of orbitals")
      endif
      if (wfs_band_max < wfs_band_min) then
         call message("WFS.BandMax < WFS.BandMin." //
     $        " Nothing written to file")
      endif

      nwfs = max(0,wfs_band_max - wfs_band_min + 1)
      nwflist(1:nk) = nwfs
      do j=1, nwfs
         iwf(1:nk,j) = wfs_band_min - 1 + j
      enddo

!     The energy window, if allowed, will take effect inside the routine
!     that writes to file, as the energies are not yet known

      call fdf_global_get(wfs_energy_min,"WFS.EnergyMin",
     $                                -huge(1.0_dp),"Ry")
      call fdf_global_get(wfs_energy_max,"WFS.EnergyMax",
     $                                 huge(1.0_dp),"Ry")

      end subroutine setup_wfs_list

      subroutine initwave( maxk, norb, nk, kpoint, overflow)
C *********************************************************************
C Finds k-points for wavefunction printout
C Based on initband routine by J.Soler
C Written by P. Ordejon, June 2003
C **************************** INPUT **********************************
C integer maxk           : Last dimension of kpoint
C integer norb           : Number of orbitals
C *************************** OUTPUT **********************************
C integer nk             : Number k points to compute wavefunctions
C real*8  kpoint(3,maxk) : k point vectors
C logical overflow       : true if the internal memory was not sufficient
C                          in which case it should be called again
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
      use parsing
#ifdef MPI
      use mpi_siesta
#endif
      implicit          none

      integer           maxk, norb, nk
      real(dp)          kpoint(3,maxk)
      logical           overflow
      external          redcel, memory
C *********************************************************************

      character 
     .  line*130, names*80, name1*10, name2*10, name3*10,
     .  scale*30

      logical, save :: Frstime = .true.
      logical       :: outlng
      logical       :: overread
      logical       :: WaveFuncPresent


      integer
     .  i, il, iu, ix, iw, iw1, iw2, iw3,
     .  lastc, lc(0:3), mscell(3,3), ni, nn, nr, nv, nmk

#ifdef MPI
      integer
     .  MPIerror
#endif

      real(dp)
     .  alat, caux(3,3), pi,
     .  rcell(3,3), reals(4), ucell(3,3)

      real(dp), dimension(:), allocatable :: values(:)
      integer :: maxnv 
   
      integer, dimension(:), allocatable :: integs(:)
      integer :: maxni 
 
C Start time counter 
      call timer( 'writewave', 1 )

      overflow = .false.

      if (Frstime) then
         call re_alloc(nwflist, 1, maxk, name="nwflist",
     $        routine="initwave")
         call re_alloc(iwf, 1, maxk, 1, 1, name="iwf",
     $        routine="initwave")
      endif

      wfs_energy_window = .true.
      call fdf_global_get(wfs_energy_min,"WFS.EnergyMin",
     $                                -huge(1.0_dp),"Ry")
      call fdf_global_get(wfs_energy_max,"WFS.EnergyMax",
     $                                 huge(1.0_dp),"Ry")

C Find k points if they are not given in argument 
      if (nk .le. 0) then

C Find if there are k-points data
        if (Node.eq.0) then
          outlng = fdf_boolean('LongOutput', .false.)
          wwf = fdf_boolean('WriteWaveFunctions',outlng)
        endif

C Find if there are k-points data
        if (Node.eq.0) then
          WaveFuncPresent = fdf_defined('WaveFuncKPoints')
          debug = fdf_boolean('WriteWaveDebug',.false.)
        endif
#ifdef MPI
        call MPI_Bcast(WaveFuncPresent,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
        call MPI_Bcast(debug,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif
        if ( WaveFuncPresent ) then

C Find lattice constant
          if (Node.eq.0) then
            alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
          endif
#ifdef MPI
          call MPI_Bcast(alat,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
#endif
          if (alat .eq. 0.d0) then
            if (Node.eq.0) then
              write(6,'(a)') 
     .          'writewave: ERROR: Lattice constant required'
            endif
            goto 999
          endif

C Find scale used in k point data
          if (Node.eq.0) then
            scale = fdf_string( 'WaveFuncKPointsScale', 'pi/a' )
          endif
#ifdef MPI
          call MPI_Bcast(scale,30,MPI_character,0,
     .      MPI_Comm_World,MPIerror)
#endif

          if (scale .eq. 'pi/a') then
            pi = 4.d0 * atan(1.d0)
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            if (Node.eq.0) then
              call redcel( alat, caux, ucell, mscell )
            endif
#ifdef MPI
            call MPI_Bcast(alat,1,MPI_double_precision,0,
     .        MPI_Comm_World,MPIerror)
            call MPI_Bcast(ucell(1,1),9,MPI_double_precision,0,
     .        MPI_Comm_World,MPIerror)
#endif
            call reclat( ucell, rcell, 1 )
          elseif (Node.eq.0) then
            write(6,'(a,/,2a,/,a)')
     .     'writewave: WARNING: Invalid value for WaveFuncKPointsScale',
     .     'writewave: Allowed values are pi/a and',
     .                ' ReciprocalLatticeVectors',
     .     'writewave: No band calculation performed'
          endif

C Loop on data lines
          if (Node.eq.0) then
            nk = 0
            WaveFuncPresent = fdf_block('WaveFuncKPoints',iu)
            do il = 1,maxlin
C Read and parse data line
              read(iu,'(a)',end=50) line
              lastc = index(line,'#') - 1
              if (lastc .le. 0) lastc = len(line)
 45           continue
C The maximum number of integers that could be introduced 
C in the input file is the 
C number of bands per k-point = number of atomic orbitals per cell
              if( .not. allocated(integs)) then
                maxni = norb
                allocate(integs(maxni))
                integs(:) = 0
              endif
C The maximum number of values that could be introduced 
C in the input file is the number of bands plus
C the coordinates of the k-point (three real numbers).
              if( .not. allocated(values)) then
                maxnv = norb + 3
                allocate(values(maxnv))
                values(:) = 0.d0
              endif

              call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )

              overread = .false.
              if( nv .gt. maxnv ) then
                maxnv    = nv
                overread = .true.
              endif
              if( ni .gt. maxni ) then
                maxni    = ni
                overread = .true.
              endif
              if( overread ) then
                if( allocated(values) ) then
                  deallocate(values)
                  allocate(values(maxnv))
                endif
                if( allocated(integs) ) then
                  deallocate(integs)
                  allocate(integs(maxni))
                endif
                goto 45
              endif

C Check if data are already finished
              if (nv .ge. 3) then

C Check syntax
                if (nr.ne.0 .and. nr.ne.3) then
                  write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                  goto 10
                endif


C Add this point to total number of k points
                nk = nk + 1

C If there is room to store k points
                if (nk .le. maxk) then

C Find coordinates of k point
                  if (scale .eq. 'pi/a') then
                    kpoint(1,nk) = values(1) * pi / alat
                    kpoint(2,nk) = values(2) * pi / alat
                    kpoint(3,nk) = values(3) * pi / alat
                  elseif (scale .eq. 'ReciprocalLatticeVectors') then
                    do ix = 1,3
                      kpoint(ix,nk) = rcell(ix,1) * values(1) +
     .                                rcell(ix,2) * values(2) +
     .                                rcell(ix,3) * values(3)
                    enddo
                  endif
                else
                  overflow = .true.
                  if (.not. Frstime)
     $                 call die('writewave: invalid initialization')
                endif


C Do the next even if there is no room to store k points
C so that the information on number of vectors to print
C is available on next call

C Find which eigenvectors should be printed

                if (nn .ge. 1) then

C Check that line contains 'from', 'to' and maybe 'step'
                  if (nn .ne. 2 .and. nn .ne. 3) then
                    write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                    goto 10
                  endif
                  name1 = names(lc(0)+1:lc(1))
                  name2 = names(lc(1)+1:lc(2))
                  if (nn .eq. 3) name3=names(lc(2)+1:lc(3))
                  if (name1 .ne. 'from' .or. name2 .ne. 'to') then
                    write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                    goto 10
                  endif
                  if (nn .eq. 3) then
                    if (name3 .ne. 'step') then
                      write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .                line(1:lastc)
                      goto 10
                    endif
                  endif

                  iw1=integs(1)
                  iw2=integs(2)
                  if (iw1 .lt. 0) iw1 = norb + iw1 +1
                  if (iw2 .lt. 0) iw2 = norb + iw2 +1
                  if (nn .eq. 3) then
                    iw3 = abs(integs(3))
                  else
                    iw3 = 1
                  endif
                  ni = 0
                  do iw = min(iw1,iw2),max(iw1,iw2),iw3
                    ni = ni + 1
                    integs(ni) = iw
                  enddo

                endif

C Store indexes of wave functions to printout

                if (ni .eq. 0) then
                  ni = norb
                  iw1=1
                  iw2=norb
                  do i=1,ni
                    integs(i)=i
                  enddo
                endif


                if (ni .le. nwf) then

                  if (nk .le. maxk) then
                    nwflist(nk) = ni
                    do i = 1,ni
                      iw = integs(i)
                      if (iw .lt. 0) iw = norb +iw +1
                      iwf(nk,i) = iw
                    enddo

                  else 

                    overflow = .true.
                    nmk=nk
                    call re_alloc(nwflist, 1, nmk, name="nwflist",
     $                   routine="initwave")
                    call re_alloc(iwf, 1, nmk, 1, nwf, name="iwf",
     $                   routine="initwave")

                  endif

                else 

                  overflow = .true.
                  if (.not. Frstime)
     $                 call die('writewave: invalid initialization')

                  nwf = ni
                  nmk = max(nk,maxk)
                  call re_alloc(nwflist, 1, nmk, name="nwflist",
     $                 routine="initwave")
                  call re_alloc(iwf, 1, nmk, 1, nwf, name="iwf",
     $                 routine="initwave")

                endif
              else
C No more lines to read => Exit do loop
                goto 50
              endif
   10         continue
            enddo
            write(6,'(a)') 'writewave: ERROR Parameter maxlin too small'
   50       continue
          endif
        else
C No k-point data available => go to exit point
          goto 999
        endif
      endif

C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(nk,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nwf,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nmk,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(overflow,1,MPI_logical,0,MPI_Comm_World,MPIerror)

      if (.not. overflow) then
         call MPI_Bcast(kpoint(1,1),3*nk,MPI_double_precision,0,
     .     MPI_Comm_World,MPIerror)
      endif

      if (overflow) then
         if (Node .ne. 0) then
            call re_alloc(nwflist, 1, nmk, name="nwflist",
     $           routine="initwave")
            call re_alloc(iwf, 1, nmk, 1, nwf, name="iwf",
     $           routine="initwave")
         endif
       endif

      if (.not. overflow) then
        call MPI_Bcast(nwflist(1),maxk,MPI_integer,0,
     .     MPI_Comm_World,MPIerror)
        call MPI_Bcast(iwf(1,1),maxk*nwf,MPI_integer,0,
     .     MPI_Comm_World,MPIerror)
      endif
 
#endif


C This is the only exit point 
  999 continue

      if( allocated(values) ) then
         deallocate(values)
      endif
                                                                                
      if( allocated(integs) ) then
         deallocate(integs)
      endif


      Frstime = .false.
      call timer( 'writewave', 2 )

      end subroutine initwave


      subroutine wwave( no, nspin, maxspn, maxo, maxuo, maxnh, maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  nk, kpoint, nuotot, gamma, occtol)
C *********************************************************************
C Finds wavefunctions at selected k-points.
C Written by P. Ordejon, June 2003
C from routine 'bands' written by J.M.Soler
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C integer maxspn              : Second dimension of ek
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
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  ef                  : Fermi energy
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of band k points
C real*8  kpoint(3,maxk)      : k point vectors
C integer nuotot              : Total number of orbitals in unit cell
C logical gamma               : Indicates if this is a Gamma point run
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
      use parallel,    only : ParallelOverK, ResetFirstCall
      use fdf
      use densematrix, only : Haux, Saux, psi
      use alloc
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use atomlist,     only : iaorb, iphorb
      use siesta_geom,  only : isa

#ifdef MPI
      use parallelsubs, only : GetNodeOrbs
#endif

      implicit          none

      integer           maxk, maxnh, maxspn, maxo, maxuo, nk, no, nspin,
     .                  nuotot, indxuo(no), listh(maxnh), numh(*), 
     .                  listhptr(*)
      logical           gamma
      real(dp)          ef, H(maxnh,nspin), kpoint(3,maxk), 
     .                  S(maxnh), xij(3,maxnh), occtol

      external          io_assign, io_close, memory
C *********************************************************************

      logical
     .  fixspin

      integer
     .  io, iu, iuo, naux, nuo, j

      real(dp)
     .  Dnew, qs(2), e1, e2, efs(2), Enew, eV, qk, qtot, 
     .  temp, wk, Entropy

C Dynamic arrays
      integer,  allocatable,   save :: muo(:)
      logical,                 save :: getD = .false.
      logical,                 save :: getPSI = .true.
      real(dp), allocatable,   save :: aux(:)
      real(dp), allocatable,   save :: ek(:,:,:)

      parameter ( eV = 1.d0 / 13.60580d0 )
      save Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

      logical  SaveParallelOverK

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
      if (nspin .gt. 2) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'writewave: WARNING: wavefunctions printed is only '
          write(6,'(/,a,/,a)')
     .       '           implemented for nspin = 1 or 2'
        endif
        goto 999
      endif

C Allocate local arrays - only aux is relevant here
      call re_alloc(Haux,1,2*nuotot*nuo,name='Haux',routine='writewave')
      call re_alloc(Saux,1,2*nuotot*nuo,name='Saux',routine='writewave')
      call re_alloc(psi,1,2*nuotot*nuo,name='psi',routine='writewave')
      allocate(ek(nuotot,nspin,nk))
      call memory('A','D',nspin*nk*nuotot,'writewave')
      naux = 2*nuotot*5
      allocate(aux(naux))
      call memory('A','D',naux,'writewave')
      allocate(muo(nuotot))
      call memory('A','I',nuotot,'writewave')

C Check indxuo 
      do iuo = 1,nuotot
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'writewave: invalid index: io, indxuo =',
     .        io, indxuo(io)
            call die('writewave: invalid indxuo')
          else
            call die()
          endif
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuotot
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
          write(6,'(/,2a,3i6)') 'writewave: ERROR: inconsistent indxuo.'
     .     ,' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            call die('writewave: ERROR: inconsistent indxuo.')
          else
            call die()
          endif
        endif
      enddo


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
          write(6,'(a26,2x,i6)') 'Number of Spins = ', nspin
          write(6,'(a26,2x,i6)') 'Number of basis orbs = ',nuotot
          write(6,*)
        endif
        write(iu) nk, gamma
        write(iu) nspin
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

C Call appropriate diagonalization routine

      if (gamma) then
         call diagg( nspin, nuo, no, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, aux, nuotot, occtol, 1, nuotot )
      else
C Handle parallel over K points option which is not allowed for here
         SaveParallelOverK = ParallelOverK
         if ((Node .eq. 0) .and. ParallelOverK) then
            write(6,"(a)")
     $           "*** Note: ParallelOverK option " //
     $           "not used for WF info generation"
         endif
         ParallelOverK = .false.
         ResetFirstCall = .true.

            
         call diagk( nspin, nuo, no, maxspn, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, xij, indxuo, nk, kpoint, wk,
     .               ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, Haux, Saux, aux, nuotot,
     .               occtol, 1 )
         ParallelOverK = SaveParallelOverK
         ResetFirstCall = .false.

      endif

C Free local arrays 
      call memory('D','I',size(muo),'writewave')
      deallocate(muo)
      call memory('D','D',size(aux),'writewave')
      deallocate(aux)
      call memory('D','I',size(ek),'writewave')
      deallocate(ek)

C This is the only exit point 
  999 continue
      call timer( 'writewave', 2 )

      end subroutine wwave


      subroutine writew(nuotot,nuo,ik,k,ispin,eo,psi,gamma)

      use precision
      use sys
      use parallel,     only : Node, Nodes
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
      real(dp), intent(in) ::  eo(*), psi(*), k(3)

      logical gamma

C  Internal variables .............................................
      integer  BNode, ntot, iie, iw, indwf, j, ind, iu

      integer number_of_wfns
      integer, allocatable  :: ind_wfn(:)

      real(dp)  :: kpoint_weight
      real(SP), dimension(:,:), allocatable :: aux   !! NOTE SP

      external io_assign, io_close

C ...................

C Allocate auxiliary arrays

      if (gamma) then
         allocate(aux(1,nuotot))
         call memory('A','D',nuotot,'writewave')
      else
         allocate(aux(2,nuotot))
         call memory('A','D',2*nuotot,'writewave')
      endif

C Find file name, and open for Node 0

      if (Node .eq. 0) then
        call io_assign( iu )
        open(iu,file=wfs_filename,form="unformatted",position='append',
     .        status='old')
      endif

C Check that the total number of orbitals is correct

#ifdef MPI
      if (Nodes.gt.1) then
        call MPI_AllReduce(nuo,ntot,1,MPI_integer,MPI_sum,
     .    MPI_Comm_World,MPIerror)
      else
        ntot = nuo
      endif
#else
      ntot = nuo
#endif

      if (ntot .ne. nuotot) call die('Inconsistent number of orbitals')

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
          write(6,'(a72)')     ' ******'
          write(6,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik, k(1:3)
          write(6,'(a22,2x,i6)') 'Spin component = ',ispin
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

C Loop over wavefunctions that should be printed

      if (Node .eq. 0 .and. debug) then
         print *, " -*- Will print ", nwflist(ik), " wfns."
         print *, " -*- ", iwf(ik,1:nwflist(ik))
      endif
      do iw = 1,nwflist(ik)
        indwf = iwf(ik,iw)

C Determine which node handles this wavefunction
        call WhichNodeOrb(indwf,Nodes,BNode)

        if (Node .eq. 0 .and. debug) then
           print *, " -*- About to print wfn ", indwf,
     $                    "from node ", Bnode
        endif

        if (Node.eq.BNode) then

C Determine the index of the orbital in the local node

          call GlobalToLocalOrb( indwf, BNode, Nodes, iie)

C Save wavefunction in aux array

C psi has different structure in diagk and diagg, so the indexing
C must be handled differently

          if (gamma) then
            do j = 1,ntot
              ind = j + (iie-1)*nuotot 
              aux(1,j) = real(psi(ind),kind=sp)
            enddo
          else
            do j = 1,ntot
              ind = 1+(j-1)*2+(iie-1)*2*nuotot
              aux(1,j) = real(psi(ind),kind=sp)
              aux(2,j) = real(psi(ind+1),kind=sp)
            enddo
          endif

        endif

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

C eigenvector is now stored in aux in all processors, and can be printed

        if (Node .eq. 0) then
          write(iu) indwf
          write(iu) eo(indwf)/eV
          write(iu) (aux(1:,j), j=1,ntot)
       endif

      enddo

C Close output file

      if (Node .eq. 0) then
        close (iu)
        call io_close(iu)
      endif

C Deallocate auxiliary arrays

      call memory('D','D',size(aux),'writewave')
      deallocate(aux)
        
      end subroutine writew

      end module writewave

