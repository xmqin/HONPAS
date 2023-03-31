! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_dynamics

      use precision
      use parallel,    only : Node, IONode
      use m_ioxv,      only : xv_file_read
      use sys,         only : die
      use atomlist,    only : iza
      use units,       only : Ang, eV
      use m_mpi_utils, only : broadcast
      use files,       only : slabel
      use alloc,       only : re_alloc, de_alloc

      implicit none


      public :: npr, nose, verlet2, pr, anneal
      private

      real(dp), parameter  :: tol = 1.0e-12_dp
      logical, parameter   :: debug = .false.
      character(len=60)    :: restart_file
      logical              :: first_time = .true.

      CONTAINS

      subroutine npr(istep,iunit,natoms,fa,stress,tp,tt,dt,
     .               ma,mn,mpr,ntcon,va,xa,hdot,h,kin,kn,kpr,vn,vpr,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED PRESSURE AND TEMPERATURE.
C The temperature is controlled with a NOSE thermostat.
C The Pressure is controlled with the PARRINELLO-RAHMAN method.
C (See Allen-Tildesley, Computer Simulations of Liquids, pp 227-238
C (Oxford Science Publications), and references therein).
C It allows for VARIABLE CELL SHAPE to accomodate the external        
C pressure.
C The combined use of Nose and Parrinello-Rahman dynamics provide 
C trajectories which sample the isothermal-isobaric ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Nose and Parrinello-Rahman variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     Nose and PR masses in eV * fs**2
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     Nose and PR masses in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces
C real*8 stress(3,3)    : Stress tensor components
C real*8 tp             : Target pressure
C real*8 tt             : Target temperature
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses
C real*8 mn             : Mass of Nose thermostat
C real*8 mpr            : Mass of Parrinello-Rahman variables
C integer ntcon         : Total number of position constraints imposed
C ******************* INPUT AND OUTPUT ****************************************
C real*8 va(3,natoms)   : Atomic velocities
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C real*8 hdot(3,3)      : Matrix of time derivatives of
C                         the vectors defining the unit cell
C real*8 h(3,3)         : Matrix of the vectors defining the unit
C                         cell at the current time step
C                         h(i,j) is the ith component of jth basis vector
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system
C real*8 kn             : Kinetic energy of Nose variable
C real*8 kpr            : Kinetic energy of Parrinello-Rahman variable
C real*8 vn             : Potential energyy of Nose variable
C real*8 vpr            : Potential energyy of P-R variables
C real*8 temp           : Instantaneous system temperature
C real*8 pressin        : Instantaneous system pressure
C *****************************************************************************
C

      implicit none

      integer, intent(in) ::  natoms,ntcon,istep, iunit

      real(dp), intent(in) ::  dt,fa(3,natoms), ma(natoms),
     $                         mn, mpr,stress(3,3),tp, tt
      real(dp), intent(inout) :: h(3,3), hdot(3,3),
     $                           va(3,natoms),xa(3,natoms)
      real(dp), intent(out)   :: kin, kn, kpr, vpr, vn, temp, pressin

C Internal variables .............................................................

      logical  ::  found

      integer  ::  ct, i, ia, info, j, k

      real(dp) :: aux1(3,3), aux2(3,3), diff, dt2, dtby2, gi(3,3),
     &            f(3,3), fi(3,3), fovermp, g(3,3), gdot(3,3),
     &            hi(3,3), hnew(3,3), hlast(3,3), hs(3), pgas,
     &            press(3,3), tdiff, tekin, twodt, vol, xdot, xlast,
     &            xnew
      real(dp) :: xdum(3), xaold(3), suncdot(3), raux(3)

      real(dp), pointer :: s(:,:)=>null(), sdot(:,:)=>null(),
     $                     snew(:,:)=>null(), sunc(:,:)=>null()
      real(dp)          :: rr

      real(dp), pointer, save :: sold(:,:)=>null()
      real(dp), save ::  hold(3,3), x, xold

      integer :: old_natoms, dummy_iza, iacc

      real(dp) :: volcel
      external :: volcel, memory

C ......................................................................

      restart_file = trim(slabel) // '.NPR_RESTART'
      if (iunit .ne. 1 .and. iunit .ne. 2)
     $   call die('pr: Wrong iunit option;  must be 1 or 2')

C Allocate local memory
      if (first_time) then
        nullify(sold)
        call re_alloc( sold, 1, 3, 1, natoms, 'sold', 'npr' )
        first_time =.false.
      endif

      nullify( s, sdot, snew, sunc )
      call re_alloc( s, 1, 3, 1, natoms, 's', 'npr' )
      call re_alloc( sdot, 1, 3, 1, natoms, 'sdot', 'npr' )
      call re_alloc( snew, 1, 3, 1, natoms, 'snew', 'npr' )
      call re_alloc( sunc, 1, 3, 1, natoms, 'sunc', 'npr' )

      ct = 3 + ntcon            ! center of mass constraints
      if (natoms .eq. 1) ct = 0

C Define constants and conversion factors .......................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
        tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif

C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute G=HtH at current time 

!      g = matmul(transpose(h),h)
      do i= 1, 3
        do j= 1, 3
          rr = 0.0
          do k= 1, 3
            rr = rr + h(k,i)*h(k,j)
          enddo
          g(i,j) = rr
        enddo
      enddo

C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) call die('npr: INVER failed')
      call inver(g,gi,3,3,info)
      if (info .ne. 0) call die( 'npr: INVER failed')

C Calculate scaled coordinates (referred to matrix H) at current time
      do ia = 1,natoms
!         s(1:3,ia) = matmul(hi,xa(1:3,ia))
        s(1:3,ia) = 0.0
        do k= 1, 3
          rr = xa(k,ia)
          do i= 1, 3
            s(i,ia) = s(i,ia) + hi(i,k)*rr
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation

!     Will need to read hold and xaold, and x and xold
!     from file if the simulation
!     is continuing by reading an XV file.
      if (istep .eq. 1) then

         if (xv_file_read) then
           if (IOnode) inquire( file=restart_file, exist=found )
           call broadcast(found)
         else
           found=.false.
         endif

         if (.not. found) then

            x = 0.0_dp
            xold = 0.0_dp
            hold = h - dt*hdot
            sold = 0.0
            do ia = 1,natoms
               xaold = xa(1:3,ia) - dt*va(1:3,ia)
     $              + (dt2/2.0_dp) * fovermp * fa(1:3,ia) / ma(ia)
!               sold(:,ia) = matmul(hi,xaold)
               do k= 1, 3
                 rr = xaold(k)
                 do i= 1, 3
                   sold(i,ia) = sold(i,ia) + hi(i,k)*rr
                 enddo
               enddo
               if (debug .and. IOnode) print *, sold(:,ia)
            enddo

         else

            ! Need to read xaold, hold, and the Nose variables
            ! from file

            if (IOnode) then
               call io_assign(iacc)
               open(unit=iacc,file=restart_file, form='formatted',
     $              status='old', action='read', position='rewind')
               read(iacc,*) x, xold
               do i = 1,3
                  read(iacc,*) (hold(j,i),j=1,3)
               enddo
               read(iacc,*) old_natoms
               if (old_natoms .ne. natoms)
     $              call die('Wrong number of atoms in NPR_RESTART')
               sold = 0.0
               do ia = 1, natoms
                  read(iacc,*) dummy_iza, (xaold(i),i=1,3)
                  if (dummy_iza .ne. iza(ia))
     $                 call die('Wrong species number in NPR_RESTART')
!                  sold(:,ia) = matmul(hi,xaold)
                  do k= 1, 3
                    rr = xaold(k)
                    do i= 1, 3
                      sold(i,ia) = sold(i,ia) + hi(i,k)*rr
                    enddo
                  enddo

                  if (debug .and. IOnode) print *, sold(:,ia)
               enddo
               call io_close(iacc)
               write(6,*)
     $            'MD restart: Read old positions, cell',
     $            ' and Nose variables from NPR_RESTART'

            endif               ! IONODE

            call broadcast(x)
            call broadcast(xold)
            call broadcast(sold(1:3,1:natoms))
            call broadcast(hold(1:3,1:3))

         endif                  ! xv_file_read
      endif      ! istep == 1
C ..................

C Compute uncorrected next positions .....................................
      if (debug .and. IOnode)
     $          print *, 'Uncorrected new reduced coordinates'
      do ia = 1,natoms
!       xdum       =  (dt2*fovermp/ma(ia)) * matmul(hi,fa(1:3,ia))     
        xdum = 0.0
        do k= 1, 3
          rr = fa(k,ia)
          do i= 1, 3
            xdum(i) = xdum(i) + hi(i,k)*rr
          enddo
        enddo
        xdum       =  (dt2*fovermp/ma(ia)) * xdum
        sunc(:,ia) = -sold(:,ia) + 2.0_dp*s(:,ia) + xdum
        if (debug .and. IOnode) print *, sunc(:,ia)
      enddo
C ...................

C Compute initial guess for Nose and Parrinello-Rahman 
C   variables at next time step ...........................................
      xnew = 2.0_dp * x - xold
      hnew = 2.0_dp * h - hold
      if (debug .and. IOnode) print *, 'xnew:\n', xnew
      if (debug .and. IOnode) print *, 'hnew:\n', hnew

C ...................

C Start selfconsistency loop to calculate Nose and P-R variables ..........
10    continue

      xlast = xnew
      hlast = hnew
        
C xdot and hdot (time derivatives at current time), and related stuff
      xdot = (xnew - xold) / twodt
      hdot = (hnew - hold)/twodt
      if (debug .and. IOnode) print *, 'xdot:\n', xdot
      if (debug .and. IOnode) print *, 'hdot:\n', hdot

      gdot = matmul(transpose(h),hdot) +
     $       matmul(transpose(hdot),h)
      if (debug .and. IOnode) print *, 'gdot:\n', gdot

      f =  matmul(gi,gdot)
      do i = 1,3
        f(i,i) = f(i,i) + xdot
      enddo
      f = dtby2 * f
      if (debug .and. IOnode) print *, 'f:\n', f

      aux1 = f
      do i = 1,3
        aux1(i,i) = aux1(i,i) + 1.0_dp
      enddo

      call inver( aux1,fi,3,3,info )
      if (info .ne. 0) call die('npr: INVER failed')
      fi = fi/twodt
      if (debug .and. IOnode) print *, 'fi:\n', fi

C Calculate corrected velocities at current time
      if (debug .and. IOnode) print *, 'Corrected velocities:'

      do ia = 1,natoms
        sdot(:,ia) = 0.0
        raux = sunc(:,ia)-sold(:,ia)
        do i= 1, 3
          do k= 1, 3
            sdot(k,ia) = sdot(k,ia) + fi(k,i)*raux(i)
          enddo
        enddo
        if (debug .and. IOnode) print *, sdot(:,ia)
      enddo

C Calculate pressure tensor at current time and ideal gas pressure
      press = 0.0_dp
      do ia = 1,natoms
        hs = 0.0
        do i= 1, 3
          do k= 1, 3
            hs(k) = hs(k) + h(k,i)*sdot(i,ia)
          enddo
        enddo

        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      if (debug .and. IOnode) print *, 'press:\n', press

      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      press = press/vol - stress

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C  Compute Nose and Parrinello-Rahman variables for next time step 
      xnew = 2.0d0 * x - xold 
     .       + (dt2 / mn) * (3.0d0 * vol * pgas - 2.0d0 * tekin)

      aux1 = 0.0_dp
      aux2 = 0.0_dp
      do i = 1,3
        aux1(i,i) = -tp
      enddo

      aux1 = aux1 + press
      aux2 = matmul(aux1,transpose(hi))
      if (debug .and. IOnode) print *, 'aux2:\n', aux2

      hnew = ( 2.0_dp * h + (dt2 * exp(2 * x) * vol / mpr) * aux2
     .                - (1.0d0 + dtby2 * xdot) * hold )
     .                / (1.0d0 - dtby2 * xdot)
      if (debug .and. IOnode) print *, 'hnew:\n', hnew

C Check if selfconsistency has been reached
      diff = abs(xnew - xlast)
      if (xlast .eq. 0.0_dp) then
        if (diff .gt. tol)  goto 10
      else
        if (diff/abs(xlast) .gt. tol)  goto 10
      endif

      diff = sum(abs(hnew-hlast))
      tdiff = sum(abs(hlast))
      if (tdiff .eq. 0.0_dp) then
        if (diff .gt. tol) goto 10
      else
        if (diff/tdiff .gt. tol) goto 10
      endif
C ...................   This is the bottom of the effective loop

C Calculate corrected atomic coordinates at next time step ................
      if (debug .and. IOnode)
     $     print *, 'Corrected new reduced coordinates'
      do ia = 1,natoms
!        suncdot = matmul(f,sold(:,ia))
        suncdot = 0.0
        do i= 1, 3
          do k= 1, 3
            suncdot(k) = suncdot(k) + f(k,i)*sold(i,ia)
          enddo
        enddo

!        snew(:,ia) = twodt * matmul(fi,sunc(:,ia) + suncdot(:))
        raux = sunc(:,ia) + suncdot(:)
        snew(:,ia) = 0.0
        do i= 1, 3
          do k= 1, 3
            snew(k,ia) = snew(k,ia) + fi(k,i)*raux(i)
          enddo
        enddo
      enddo
      snew = twodt*snew

      if (debug .and. IOnode) then
        do ia = 1,natoms
          print *, snew(:,ia)
        enddo
      endif

!     This is the place to store the current magnitudes
!
!      do ia = 1,natoms
!         xa(:,ia) = matmul(h,s(:,ia))
!         va(:,ia) = matmul(h,sdot(:,ia))
!      enddo
!      call add_to_md_file(xa,va,cell=h,vcell=hdot,nose=x,nosedot=xdot)
C Save current atomic positions as old ones, 
C   and next positions as current ones
      sold = s
      s    = snew
      hold = h
      h    = hnew
      xold = x
      x    = xnew

C Transform back to absolute coordinates 
      xa= 0.0
      do ia = 1,natoms
!       xa(:,ia) = matmul(h,s(:,ia))
        do i= 1, 3
          do k= 1, 3
            xa(k,ia) = xa(k,ia) + h(k,i)*s(i,ia)
          enddo
        enddo
      enddo

      va = 0.0
      do ia = 1,natoms
!       va(:,ia) = matmul(h,sdot(:,ia))
        do i= 1, 3
          do k= 1, 3
            va(k,ia) = va(k,ia) + h(k,i)*sdot(i,ia)
          enddo
        enddo
      enddo

!       Save (now old) positions and cell, and
!       Nose variables,  to NPR_RESTART
!
      if (Node .eq. 0) then
        call io_assign(iacc)
        open(unit=iacc,file=restart_file, form='formatted',
     $       status='unknown', action= 'write', position='rewind')
        write(iacc,*) x, xold 
        do i = 1,3
           write(iacc,*) (hold(j,i),j=1,3)
        enddo
        write(iacc,*) natoms
        do ia = 1, natoms
!         xaold(:) = matmul(h,sold(:,ia))
          xaold = 0.0
          do i= 1, 3
            do k= 1, 3
              xaold(k) = xaold(k) + h(k,i)*sold(i,ia)
            enddo
          enddo
          write(iacc,*) iza(ia), (xaold(i),i=1,3)
        enddo
        call io_close(iacc)
      endif

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Kinetic energy of Nose variable
      kn = (1.0d0 / 2.0d0) * mn * xdot**2

C Kinetic energy of Parrinello-Rahman variables
      kpr = 0.5_dp * mpr * sum(hdot(1:3,1:3)**2) / exp(2.0_dp * x)

C Potential energy of Nose variable
      vn = 2.0d0 * tekin * xold

C Potential energy of Parrinello-Rahman variables
      vpr = tp * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call de_alloc( sunc, 'sunc', 'npr' )
      call de_alloc( snew, 'snew', 'npr' )
      call de_alloc( sdot, 'sdot', 'npr' )
      call de_alloc( s, 's', 'npr' )

      end subroutine npr
    
      subroutine pr(istep,iunit,iquench,natoms,fa,stress,tp,dt,
     .               ma,mpr,ntcon,va,xa,hdot,h,kin,kpr,vpr,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED PRESSURE.
C The Pressure is controlled with the PARRINELLO-RAHMAN method.
C (See Allen-Tildesley, Computer Simulations of Liquids, pp 227-238
C (Oxforf Science Publications), and references therein).
C It allows for VARIABLE CELL SHAPE to accomodate the external
C pressure.
C The use of Parrinello-Rahman dynamics provide 
C trajectories which sample the constant NPE ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Parrinello-Rahman variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     PR mass in eV * fs**2
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     PR mass in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:   
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity 
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 stress(3,3)    : Stress tensor components 
C real*8 tp             : Target pressure
C real*8 dt             : Length of the time step 
C real*8 ma(natoms)     : Atomic masses
C real*8 mpr            : Mass of Parrinello-Rahman variables
C integer ntcon         : Total number of position constraints imposed
C ******************* INPUT AND OUTPUT ****************************************
C real*8 va(3,natoms)   : Atomic velocities
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C real*8 hdot(3,3)      : Matrix of time derivatives of
C                         the vectors defining the unit cell
C real*8 h(3,3)         : Matrix of the vectors defining the unit
C                         cell at the current time step 
C                         h(i,j) is the ith component of jth basis vector
C                        (input: current time step; output: next time step)
C ************************* OUTPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system 
C real*8 kpr            : Kinetic energy of Parrinello-Rahman var 
C real*8 vpr            : Potential energy of P-R variables 
C real*8 temp           : Instantaneous system temperature 
C real*8 pressin        : Instantaneous system pressure 
C *****************************************************************************

      implicit none

      integer, intent(in) ::  natoms,ntcon,istep,iquench,iunit

      real(dp), intent(in) ::  dt,fa(3,natoms), ma(natoms),
     $                         mpr,stress(3,3),tp
      real(dp), intent(inout) :: h(3,3), hdot(3,3),
     $                           va(3,natoms),xa(3,natoms)
      real(dp), intent(out)   :: kin, kpr, vpr, temp, pressin

C Internal variables 

      logical :: found

      integer  :: ct, i, info, ia, j, k

      real(dp) :: a1, a2, aux1(3,3), aux2(3,3), diff, dot, dt2, dtby2,
     &            f(3,3), fi(3,3), fovermp, g(3,3), gdot(3,3), gi(3,3),
     &            hi(3,3), hnew(3,3), hlast(3,3), hs(3), pgas,
     &            press(3,3), tdiff, twodt, vol
      real(dp) :: xdum(3), xaold(3), suncdot(3), raux(3)

      real(dp), pointer :: s(:,:)=>null(), sdot(:,:)=>null(),
     $                     snew(:,:)=>null(), sunc(:,:)=>null()
      real(dp), pointer, save :: sold(:,:)=>null()
      real(dp), save ::  hold(3,3)

      integer :: old_natoms, dummy_iza, iacc

      real(dp) :: volcel
      external :: volcel, memory
C ...............................................................................
      
      restart_file = trim(slabel) // '.PR_RESTART'

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $   call die('pr: Wrong iunit option;  must be 1 or 2')
      ct = 3 + ntcon             ! center-of-mass constraint
      if (natoms .eq. 1) ct = 0 

C Allocate local memory
      call re_alloc( sold, 1, 3, 1, natoms, 'sold', 'pr' )
      call re_alloc( s, 1, 3, 1, natoms, 's', 'pr' )
      call re_alloc( sdot, 1, 3, 1, natoms, 'sdot', 'pr' )
      call re_alloc( snew, 1, 3, 1, natoms, 'snew', 'pr' )
      call re_alloc( sunc, 1, 3, 1, natoms, 'sunc', 'pr' )

C Define constants and conversion factors .......................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute G=HtH at current time 

      g = matmul(transpose(h),h)

C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) call die('pr: INVER failed')
      call inver(g,gi,3,3,info)
      if (info .ne. 0) call die( 'pr: INVER failed')

C Calculate scaled coordinates (referred to matrix H) at current time
      s = 0.0
      do ia = 1,natoms
!       s(1:3,ia) = matmul(hi,xa(1:3,ia))
        do i= 1, 3
          do k= 1, 3
            s(k,ia) = s(k,ia) + hi(k,i)*xa(i,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
!
!     Will need to read hold and xaold from file if the simulation
!     is continuing by reading an XV file.
!
      if (istep .eq. 1) then

         if (xv_file_read) then
           if (IONode) inquire( file=restart_file, exist=found )
           call broadcast(found)
         else
           found=.false.
         endif

         if (.not. found) then

            hold = h - dt * hdot
            if (debug .and. IOnode) print *, 'Old reduced coordinates'
            sold = 0.0
            do ia = 1,natoms
              xaold(1:3) = xa(1:3,ia) - dt*va(1:3,ia)
     $             + (dt2/2.0d0) * fovermp * fa(1:3,ia) / ma(ia)
!             sold(1:3,ia) = matmul(hi,xaold(1:3))
              do i= 1, 3
                do k= 1, 3
                  sold(k,ia) = sold(k,ia) + hi(k,i)*xaold(i)
                enddo
              enddo

              if (debug .and. IOnode) print *, sold(:,ia)
            enddo

         else

            ! Need to read xaold and hold from file

            if (IOnode) then
               call io_assign(iacc)
               open(unit=iacc,file=restart_file, form='formatted',
     $              status='old', action='read', position='rewind')
               do i = 1,3
                  read(iacc,*) (hold(j,i),j=1,3)
               enddo
               read(iacc,*) old_natoms
               if (old_natoms .ne. natoms)
     $           call die('Wrong number of atoms in PR_RESTART')
               sold = 0.0
               do ia = 1, natoms
                 read(iacc,*) dummy_iza, (xaold(i),i=1,3)
                 if (dummy_iza .ne. iza(ia))
     $                call die('Wrong species number in PR_RESTART')
!                sold(1:3,ia) = matmul(hi,xaold(1:3))
                 do i= 1, 3
                   do k= 1, 3
                     sold(k,ia) = sold(k,ia) + hi(k,i)*xaold(i)
                   enddo
                 enddo

                 if (debug .and. IOnode) print *, sold(:,ia)
               enddo
               call io_close(iacc)
               write(6,*)
     $                 'MD restart: Read old positions and cell',
     $                 ' from PR_RESTART'

            endif               ! IONODE

            call broadcast(sold(1:3,1:natoms))
            call broadcast(hold(1:3,1:3))

         endif                  ! xv_file_read
      endif              ! istep ==1

C Compute uncorrected next positions .....................................
!
!     Note that re-scaled forces are used.
      if (debug .and. IOnode)
     $          print *, 'Uncorrected new reduced coordinates'
      do ia = 1,natoms
!       raux = matmul(hi,fa(1:3,ia)) 
        raux = 0.0
        do i= 1, 3
          do k= 1, 3
            raux(k) = raux(k) + hi(k,i)*fa(i,ia)
          enddo
        enddo

        xdum =  (dt2*fovermp/ma(ia)) *raux    
        sunc(1:3,ia) = -sold(1:3,ia) + 2.0_dp*s(1:3,ia) + xdum
        if (debug .and. IOnode) print *, sunc(:,ia)
      enddo
C ...................

C Compute initial guess for Parrinello-Rahman 
C   variables at next time step ...........................................

      hnew = 2.0_dp * h - hold
      if (debug .and. IOnode) print *, 'hnew:\n', hnew

C Start selfconsistency loop to calculate P-R variables ..........
10    continue

      hlast = hnew
        
C hdot (time derivatives at current time), and related stuff

      hdot = (hnew - hold)/twodt
      if (debug .and. IOnode) print *, 'hdot:\n', hdot

      gdot = matmul(transpose(h),hdot) +
     $       matmul(transpose(hdot),h)
      if (debug .and. IOnode) print *, 'gdot:\n', gdot

      f(:,:) = dtby2 * matmul(gi,gdot)
      if (debug .and. IOnode) print *, 'f:\n', f

      aux1(1:3,1:3) = f(1:3,1:3)
      do i = 1,3
        aux1(i,i) = aux1(i,i) + 1.0_dp
      enddo

      call inver(aux1,fi,3,3,info)
      if (info .ne. 0) call die( 'pr: INVER failed')

      fi = fi/twodt
      if (debug .and. IOnode) print *, 'fi:\n', fi

C Calculate corrected velocities at current time

      if (debug .and. IOnode) print *, 'Corrected velocities:'
      do ia = 1,natoms
        raux = sunc(1:3,ia)-sold(1:3,ia)
!       sdot(1:3,ia) = matmul(fi,raux)
        sdot(1:3,ia) = 0.0
        do i= 1, 3
          do k= 1, 3
            sdot(k,ia) = sdot(k,ia) + fi(k,i)*raux(i)
          enddo
        enddo
        if (debug .and. IOnode) print *, sdot(:,ia)
      enddo

C Calculate pressure tensor at current time and ideal gas pressure

      press(1:3,1:3) = 0.0_dp
      do ia = 1,natoms
!        hs(1:3) = matmul(h,sdot(1:3,ia))
        hs = 0.0
        do i= 1, 3
          do k= 1, 3
            hs(k) = hs(k) + h(k,i)*sdot(i,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      if (debug .and. IOnode) print *, 'press:\n', press

      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      press(:,:) = press(:,:)/vol - stress(:,:)

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C  Compute Parrinello-Rahman variables for next time step 
      aux1 = 0.0_dp
      aux2 = 0.0_dp
      do i = 1,3
        aux1(i,i) = -tp
      enddo
      aux1(1:3,1:3) = aux1(1:3,1:3) + press(1:3,1:3)

      aux2(1:3,1:3) = matmul(aux1,transpose(hi))
      if (debug .and. IOnode) print *, 'aux2:\n', aux2

      hnew = 2.0_dp * h + (dt2 * vol / mpr) * aux2 - hold
      if (debug .and. IOnode) print *, 'hnew:\n', hnew


C Check if selfconsistency has been reached
      diff = sum(abs(hnew(1:3,1:3)-hlast(1:3,1:3)))
      tdiff = sum(abs(hlast(1:3,1:3)))
      if (tdiff .eq. 0.0d0) then
        if (diff .gt. tol) goto 10
      else
        if (diff/tdiff .gt. tol) goto 10
      endif
C ...................   This is the bottom of the effective loop

C Calculate corrected atomic coordinates at next time step ................
      if (debug .and. IOnode)
     $     print *, 'Corrected new reduced coordinates'

      snew = 0.0
      do ia = 1,natoms
!       suncdot(:) = matmul(f,sold(:,ia))
        suncdot = 0.0
        do i= 1, 3
          do k= 1, 3
            suncdot(k) = suncdot(k) + f(k,i)*sold(i,ia)
          enddo
        enddo

        raux = sunc(:,ia) + suncdot
!       snew(:,ia) = twodt * matmul(fi,raux)
        do i= 1, 3
          do k= 1, 3
            snew(k,ia) = snew(k,ia) + fi(k,i)*raux(i)
          enddo
        enddo
        snew(:,ia) = twodt*snew(:,ia)

        if (debug .and. IOnode) print *, snew(:,ia)
      enddo

C Quench option if iquench = 0 ..............................................

C Quench velocity components going uphill
      if (iquench .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            a1 = 0.0d0
            a2 = 0.0d0
            do j = 1,3
              a1 = a1 + hi(i,j) * fovermp * fa(j,ia) / ma(ia) 
              a2 = a2 - f(i,j) * sdot(j,ia) / dtby2
            enddo
            dot = a1 * sdot(i,ia)
            if (dot .lt. 0.0) then
              sdot(i,ia) = 0.0
              snew(i,ia) = s(i,ia)
            endif
          enddo
        enddo
  
        do i = 1,3
          do j = 1,3
            dot = hdot(i,j) * aux2(i,j)
            if (dot .le. 0.0) then
              hdot(i,j) = 0.0
              hnew(i,j) = h(i,j)
            endif
          enddo
        enddo
            
C Compute gas pressure again, in case quench has happened
        press = 0.0_dp
        do ia = 1,natoms
!         hs(:) = matmul(h,sdot(:,ia))
          hs = 0.0
          do i= 1, 3
            do k= 1, 3
              hs(k) = hs(k) + h(k,i)*sdot(i,ia)
            enddo
          enddo
          do j = 1,3
             do i = 1,3
                press(i,j) = press(i,j) +
     $                        ma(ia) * hs(i) * hs(j) / fovermp
             enddo
          enddo
        enddo
        pgas = 0.0d0
        do i = 1,3
          pgas = pgas + press(i,i) / vol
        enddo

      endif  ! quench
C ....................
          
!     This is the place to store the current magnitudes
!
!      do ia = 1,natoms
!         xa(:,ia) = matmul(h,s(:,ia))
!         va(:,ia) = matmul(h,sdot(:,ia))
!      enddo
!      call add_to_md_file(xa,va,cell=h,vcell=hdot)

C Save current atomic positions as old ones, 
C   and next positions as current ones

      sold = s
      s    = snew
      hold = h
      h    = hnew

C Transform back to absolute coordinates 
      xa= 0.0
      do ia = 1,natoms
!       xa(:,ia) = matmul(h,s(:,ia))
        do i= 1, 3
          do k= 1, 3
            xa(k,ia) = xa(k,ia) + h(k,i)*s(i,ia)
          enddo
        enddo
      enddo

      va = 0.0
      do ia = 1,natoms
!       va(:,ia) = matmul(h,sdot(:,ia))
        do i= 1, 3
          do k= 1, 3
            va(k,ia) = va(k,ia) + h(k,i)*sdot(i,ia)
          enddo
        enddo
      enddo

!       Save (now old) positions and cell to PR_RESTART
!
      if (Node .eq. 0) then
        call io_assign(iacc)
        open(unit=iacc,file=restart_file, form='formatted',
     $       status='unknown', action= 'write', position='rewind')
        do i = 1,3
           write(iacc,*) (hold(j,i),j=1,3)
        enddo
        write(iacc,*) natoms
        do ia = 1, natoms
!          xaold(:) = matmul(h,sold(:,ia))
          xaold = 0.0
          do i= 1, 3
            do k= 1, 3
              xaold(k) = xaold(k) + h(k,i)*sold(i,ia)
            enddo
          enddo
          write(iacc,*) iza(ia), (xaold(i),i=1,3)
        enddo
        call io_close(iacc)
      endif


C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Kinetic energy of Parrinello-Rahman variables
      kpr = 0.5_dp * mpr * sum(hdot(1:3,1:3)**2)

C Potential energy of Parrinello-Rahman variables
      vpr = tp * vol
      if (debug .and. IOnode) print *, 'kpr, vpr: ', kpr, vpr

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif
C .....................

C Deallocate local memory
      call de_alloc( sunc, 'sunc', 'pr' )
      call de_alloc( snew, 'snew', 'pr' )
      call de_alloc( sdot, 'sdot', 'pr' )
      call de_alloc( s, 's', 'pr' )

      end subroutine pr
    
      subroutine nose( istep, iunit, natoms, fa, tt, dt, ma, mn, ntcon,
     .                 va, xa, kin, kn, vn, temp )
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED TEMPERATURE.
C The temperature is controlled with a NOSE thermostat.
C The use of Nose dynamics provides trajectories which sample the 
C isothermal ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Nose variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C     Nose mass in eV * fs**2
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C     Nose mass in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces
C real*8 tt             : Target temperature
C real*8 dt             : Length of the time step 
C real*8 ma(natoms)     : Atomic masses
C real*8 mn             : Mass of Nose thermostat
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system
C real*8 kn             : Kinetic energy of Nose variable 
C real*8 vn             : Potential energyy of Nose var
C real*8 temp           : Instantaneous system temperature 
C *****************************************************************************
C

      integer 
     .   natoms,ntcon,istep,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,kn,
     .  ma(natoms),mn,tt,
     .  va(3,natoms),vn,xa(3,natoms)

      external
     .  memory
C Internal variables .........................................................

      logical
     .  found

      integer
     .  ct,i,ia

      integer  :: iacc, dummy_iza, old_natoms
      real(dp) :: old_dt

      save x,xold

      real(dp)
     .  diff,dt2,dtby2,fact,fovermp,
     .  tekin,temp,twodt,
     .  x,xdot,xlast,xnew,xold

      real(dp), pointer, save :: xanew(:,:)=>null(), xaold(:,:)=>null()
C .............................................................................

      restart_file = trim(slabel) // '.NOSE_RESTART'

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $    call die('nose: Wrong iunit option;  must be 1 or 2')
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory and initialise
      nullify(xanew)
      call re_alloc( xanew, 1, 3, 1, natoms, 'xanew', 'nose' )
      if (.not.associated(xaold)) then
        nullify(xaold)
        call re_alloc( xaold, 1, 3, 1, natoms, 'xaold', 'nose' )
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia)=0.0d0
          enddo
        enddo
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
         tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then

         if (xv_file_read) then
           if (IONode) inquire( file=restart_file, exist=found )
           call broadcast(found)
         else
           found=.false.
         endif

         if (.not. found) then

C     Compute old positions in terms of current positions and velocities
C     if the time step is the first of the simulation 
!     and we start from x(t), v(t) *at the same time*.
!     (e.g., when the velocities are constructed from
!      the Boltzmann distribution).
!     In this case the algorithm works out well.
!     Nose variables are set to zero, as there is currently no
!     better way to initialize them...

            x = 0.0d0
            xold = 0.0d0
            do ia = 1,natoms
               do i = 1,3
                  xaold(i,ia) = xa(i,ia) - dt * va(i,ia)
     .                 + (dt2/2.0d0) * fovermp * fa(i,ia) / ma(ia)
               enddo
            enddo

         else

!         For restarts, we need information about the old 
!         positions, and the Nose variables
!
           if (Node .eq. 0) then
            call io_assign(iacc)
            open(unit=iacc,file=restart_file, form='formatted',
     $           status='old', action='read', position='rewind')
            read(iacc,*) old_natoms, old_dt
            read(iacc,*) x, xold
            if (old_natoms .ne. natoms)
     $           call die('Wrong number of atoms in NOSE_RESTART')
            do ia = 1, natoms
               read(iacc,*) dummy_iza, (xaold(i,ia),i=1,3)
               if (dummy_iza .ne. iza(ia))
     $              call die('Wrong species number in NOSE_RESTART')
            enddo
            call io_close(iacc)
            write(6,*)
     $         'MD restart: Read old positions and Nose variables',
     $              ' from NOSE_RESTART'
            if (abs(old_dt - dt) .gt. 1.0d-8) then
               write(6,*) '**WARNING: Timestep has changed. Old: ',
     $              old_dt, ' New: ', dt
               write(6,*) '**WARNING: Approximating old positions.'
                  ! First order, using the positions and velocities 
                  ! at t-old_dt (positions from NOSE_RESTART, velocities
                  !              from XV file)
               xaold(1:3,1:natoms) = xaold(1:3,1:natoms) -
     $                              (dt-old_dt) * va(1:3,1:natoms)
            endif               ! dt /= old_dt

           endif                  ! IONode
          
            call broadcast(x)
            call broadcast(xold)
            call broadcast(xaold(1:3,1:natoms))

         endif     ! xv_file_read
      endif        ! istep == 1
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) =  2.0d0 * xa(i,ia) - xaold(i,ia) +
     .                   dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Compute uncorrected velocities and kinetic energy ......................
      kin = 0.d0
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
C ..................

C Compute initial guess for Nose variables at next time step .............
      xnew = 2.0_dp * x - xold
C ...................

C Start selfconsistency loop to calculate Nose variable ..................
10    continue

      xlast = xnew
        
C xdot and hdot (time derivatives at current time), and related stuff
      xdot = (xnew - xold) / twodt
      fact = (1.0_dp/(1.0_dp+xdot*dtby2))

C  Compute Nose variable for next iteration
      xnew = 2.0_dp * x - xold 
     .       + (dt2/mn) * 2.0_dp * (fact**2 * kin - tekin)

C Check if selfconsistency has been reached
      diff = abs(xnew - xlast)
      if (xlast .eq. 0.0d0) then
        if (diff .gt. tol)  goto 10
      else
        if (diff/abs(xlast) .gt. tol)  goto 10
      endif
C ...................

C Calculate corrected atomic coordinates at next time step, 
C and corrected velocities and kinetic energy at current time step .........
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) = fact * ( xanew (i,ia) +
     .                   dtby2 * xdot * xaold(i,ia))
          va(i,ia) = fact * va(i,ia)
        enddo
      enddo
      kin = kin * fact**2 
C ...................
!     Here we can save x, xa, va for MD  (experimental)
!
!     call add_to_md_file(xa,va,nose=x,nosedot=xdot)
C Save current atomic positions as old ones, 
C   and next positions as current ones

      xaold(1:3,1:natoms) = xa(1:3,1:natoms)
      xa(1:3,1:natoms) = xanew(1:3,1:natoms)

      xold = x
      x = xnew

C Calculate Kinetic and potential energies ................................
C Kinetic energy of Nose variable
      kn = (1.0d0 / 2.0d0) * mn * xdot**2

C Potential energy of Nose variable (in eV)
      vn = 2.0d0 * tekin * xold

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = kin / (0.5d0 * (3.d0 * natoms - ct) * 8.617d-5)
      else
        temp = kin / (0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * eV)
      endif

      if (Node .eq. 0) then
!
!       Save (now old) positions and nose variables to NOSE_RESTART
!
         call io_assign(iacc)
         open(unit=iacc,file=restart_file, form='formatted',
     $        status='unknown', action= 'write', position='rewind')
         write(iacc,*) natoms, dt
         write(iacc,*) x, xold
         do ia = 1, natoms
            write(iacc,*) iza(ia), (xaold(i,ia),i=1,3) 
         enddo
         call io_close(iacc)
      endif

C .....................
      end subroutine nose

      subroutine anneal(istep,iunit,ianneal,taurelax,bulkm,
     .               natoms,fa,stress,tp,tt,dt,
     .               ma,ntcon,va,xa,h,kin,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with a TARGET PRESSURE AND TEMPERATURE.
C The system is driven to a desired temperature and pressure in
C a given time, by rescaling the velocities and the cell shape and size.
C It needs an estimate of the bulk modulus of the system, to determine
C the rate of change of cell shape to accomodate to the target
C pressure in the required time. A wrong estimate will simply
C drive the system to the desired pressure, but in a different time
C than the especified in the input. (See Kittel for representative
C values of bulk moduli for materials).
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     Bulk modulus in eV/A**3
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     Bulk modulus in Ry/Bohr**3
C ************************* INPUT *********************************************
C integer istep       : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer ianneal     : Work mode option:
C                       1 = reach target temperature only
C                       2 = reach target pressure only
C                       3 = reach target temperature and pressure
C real*8 taurelax     : Relaxation time to reach desired T and P
C real*8 bulkm        : Estimate of the Bulk Modulus of the system
C integer natoms      : Number of atoms in the simulation cell
C real*8 fa(3,natoms) : Atomic forces 
C real*8 stress(3,3)  : Stress tensor components 
C real*8 tp           : Target pressure 
C real*8 tt           : Target temperature 
C real*8 dt           : Length of the time step 
C real*8 ma(natoms)   : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms) : Atomic velocities
C                       (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms) : Atomic coordinates 
C                      (input: current time step; output: next time step)
C real*8 h(3,3)       : Matrix of the vectors defining the unit
C                       cell at the current time step 
C                       h(i,j) is the ith component of jth basis vector
C                      (input: current time step; output: next time step)
C ************************* OUTPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system 
C real*8 temp           : Instantaneous system temperature 
C real*8 pressin        : Instantaneous system pressure 
C *****************************************************************************

      logical
     .   found

      integer 
     .   natoms,ntcon,istep,ianneal,iunit,iacc

      real(dp)
     .  bulkm,dt,fa(3,natoms),h(3,3),kin,
     .  ma(natoms),stress(3,3),taurelax,tp,tt,
     .  va(3,natoms),xa(3,natoms)

C Internal variables .............................................................

      integer  :: ct, i, ia, info, j, k
      real(dp) :: dt2, fovermp, hi(3,3), hs(3), pgas, press(3,3),
     &            pressin, rfac, rfac2, tekin, temp, twodt, vol, volcel

      real(dp), pointer :: s(:,:)=>null(), sdot(:,:)=>null(),
     $                     snew(:,:)=>null(), sunc(:,:)=>null()
      real(dp), pointer, save :: sold(:,:)=>null()

      real(dp), dimension(3) :: xdum, xaold
      logical :: have_t_target, have_p_target
       external
     &  volcel, memory
C ....................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $    call die('anneal: Wrong iunit option;  must be 1 or 2')

!!!      if (taurelax/dt .lt. 0.1) return

      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      call re_alloc( sold, 1, 3, 1, natoms, 'sold', 'anneal' )
      call re_alloc( s, 1, 3, 1, natoms, 's', 'anneal' )
      call re_alloc( sdot, 1, 3, 1, natoms, 'sdot', 'anneal' )
      call re_alloc( snew, 1, 3, 1, natoms, 'snew', 'anneal' )
      call re_alloc( sunc, 1, 3, 1, natoms, 'sunc', 'anneal' )

C Define constants and conversion factors .......................................
      dt2   = dt**2
      twodt = dt*2.0d0

      have_t_target = (ianneal .eq. 1 .or. ianneal .eq. 3)
      have_p_target = (ianneal .eq. 2 .or. ianneal .eq. 3)

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
        tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif

C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute Inverse of H at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) call die('anneal: INVER failed')

C Calculate scaled coordinates (referred to matrix H) at current time
      s = 0.0
      do ia = 1,natoms
!       s(1:3,ia) = matmul(hi,xa(1:3,ia))
        do i= 1, 3
          do k= 1, 3
            s(k,ia) = s(k,ia) + hi(k,i)*xa(i,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then
        if (IOnode) then
          restart_file = trim(slabel) // '.ANNEAL_RESTART'
          inquire(file=restart_file, exist=found)
        end if
        call broadcast(found)
        if (found) then
          if (IOnode) then
            call io_assign(iacc)
            open(unit=iacc, file=restart_file, form='formatted',
     $           status='old', action='read', position='rewind')
            read(iacc,*) sold
            call io_close(iacc)
          endif
          call broadcast(sold)
        else
          do ia = 1,natoms
            xaold(1:3) = xa(1:3,ia) - dt*va(1:3,ia)
     $           + (dt2/2.0_dp) * fovermp * fa(1:3,ia) / ma(ia)
!           sold(1:3,ia) = matmul(hi,xaold(1:3))
            do i= 1, 3
              do k= 1, 3
                sold(k,ia) = sold(k,ia) + hi(k,i)*xaold(i)
              enddo
            enddo
          enddo
        endif
      endif
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
!         xdum(1:3) =  (dt2*fovermp/ma(ia)) * matmul(hi,fa(1:3,ia))
        xdum = 0.0
        do i= 1, 3
          do k= 1, 3
            xdum(k) = xdum(k) + hi(k,i)*fa(i,ia)
          enddo
        enddo
        xdum = (dt2*fovermp/ma(ia))*xdum
        sunc(:,ia) = -sold(:,ia) + 2.0_dp*s(:,ia) + xdum
      enddo
C ...................

C Calculate uncorrected velocities at current time
      sdot = (sunc - sold) / twodt

C Calculate pressure tensor at current time and ideal gas pressure
      press = 0.0_dp
      do ia = 1,natoms
!       hs(1:3) = matmul(h,sdot(1:3,ia))
        hs = 0.0
        do i= 1, 3
          do k= 1, 3
            hs(k) = hs(k) + h(k,i)*sdot(i,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0

      press = press/vol -stress

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0d0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C Compute kinetic energy
      kin = (3.0d0 / 2.0d0) * pgas * vol

      if (IOnode) write(*,*) 'Anneal: Kinetic Energy= ', kin

      if (have_t_target) then
C Correct velocities to reach target termperature
      if (kin .eq. 0.0) then
        rfac2 = 1.0d0 + dt/taurelax
      else
        rfac2 = (1.0d0 + dt/taurelax * (tekin/kin -1.0d0))
      endif
      if (rfac2 .le. 0.0) call die('Wrong anneal parameter')
      rfac = sqrt(rfac2)
      if (IOnode) write(*,*) 'Anneal: Velocity scale factor = ', rfac

      sdot(1:3,1:natoms) = rfac * sdot(1:3,1:natoms)

C Compute again pressure, with corrected velocities

      press(1:3,1:3) = 0.0_dp

      do ia = 1,natoms
!        hs(1:3) = matmul(h,sdot(1:3,ia))
        hs = 0.0
        do i= 1, 3
          do k= 1, 3
            hs(k) = hs(k) + h(k,i)*sdot(i,ia)
          enddo
        enddo

        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0

      press(1:3,1:3) = press(1:3,1:3)/vol - stress(1:3,1:3)

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0
      endif


C Correct new possitions according to corrected velocities
      snew = sold + twodt * sdot

      if (have_p_target) then

C Correct cell shape to reach target pressure
      do i = 1,3
        do j = 1,3
          if (i .ne. j) then
            rfac2 = 1.0 + (dt / taurelax) * press(i,j) 
     .                    / (0.5 * bulkm)
          else
            rfac2 = 1.0 + (dt / taurelax) * (press(i,i) - tp) 
     .                    / (0.5 * bulkm)
          endif
          if (rfac2 .le. 0.0) call die('Wrong anneal parameter')
          rfac = sqrt(rfac2)
          h(i,j) = rfac * h(i,j)
        enddo
      enddo
      if (IOnode)  write(*,*) 'Anneal: Cell scale factor = ', rfac
      endif


!     This is the place to store the current magnitudes
!
!      do ia = 1,natoms
!         xa(:,ia) = matmul(h,s(:,ia))
!         va(:,ia) = matmul(h,sdot(:,ia))
!      enddo
!      call add_to_md_file(xa,va,cell=h)


C Save current atomic positions as old ones, 
C   and next positions as current ones

      sold = s
      s = snew

C Transform back to absolute coordinates 
      xa= 0.0
      do ia = 1,natoms
!       xa(:,ia) = matmul(h,s(:,ia))
        do i= 1, 3
          do k= 1, 3
            xa(k,ia) = xa(k,ia) + h(k,i)*s(i,ia)
          enddo
        enddo
      enddo

      va = 0.0
      do ia = 1,natoms
!       va(:,ia) = matmul(h,sdot(:,ia))
        do i= 1, 3
          do k= 1, 3
            va(k,ia) = va(k,ia) + h(k,i)*sdot(i,ia)
          enddo
        enddo
      enddo

C ....................

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms 
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call de_alloc( sunc, 'sunc', 'anneal' )
      call de_alloc( snew, 'snew', 'anneal' )
      call de_alloc( sdot, 'sdot', 'anneal' )
      call de_alloc( s, 's', 'anneal' )

!!!      taurelax = taurelax - dt

      if (IOnode) then
        restart_file = trim(slabel) // '.ANNEAL_RESTART'
        call io_assign(iacc)
        open(unit=iacc, file=restart_file, form='formatted',
     $       status='unknown', action='write', position='rewind')
        write(iacc,*) sold
        call io_close(iacc)
      endif

      end subroutine anneal


      subroutine verlet1(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp)
C *************************************************************************
C Subroutine for MD simulations using the Original Verlet Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 78)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:   
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity 
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities 
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature 
C *****************************************************************************
C
      integer 
     .   natoms,ntcon,istep,iquench,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia

      real(dp)
     .  dot,dt2,fovermp,temp,twodt

      real(dp), pointer, save :: xanew(:,:)=>null(),
     $                           xaold(:,:)=>null()

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $     call die('verlet1: Wrong iunit option;  must be 1 or 2')

      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      call re_alloc( xanew, 1, 3, 1, natoms, 'xanew', 'verlet1' )
      if (.not.associated(xaold)) then
        call re_alloc( xaold, 1, 3, 1, natoms, 'xaold', 'verlet1' )
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038  * Ang**2 / eV
      endif

C ........................

C Compute old coordinates if the time step is the first of the simulation .....
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia) = xa(i,ia) - dt * va(i,ia)
     .                + (dt2/2.0d0) * fovermp * fa(i,ia) / ma(ia)
          enddo
        enddo
      endif

C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) = - xaold(i,ia) + 2.0d0 * xa(i,ia)
     .                  + dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Calculate velocities at current time .....................................
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
        enddo
      enddo
C ...................

C Quench option if iquench = 0 ..............................................
      if (iquench .eq. 1) then

C Quench velocity components going uphill
        do ia = 1,natoms
          do i = 1,3
            dot = va(i,ia) * fa(i,ia)
            if (dot .lt. 0.0) then
              va(i,ia) = 0.0
              xanew(i,ia) = xa(i,ia)
            endif
          enddo
        enddo

      endif
C......................


C Save current atomic positions as old ones, 
C   and next positions as current ones .....................................
      do i = 1,3
        do ia = 1,natoms
          xaold(i,ia) = xa(i,ia)
          xa(i,ia) = xanew(i,ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call de_alloc( xanew, 'xanew', 'verlet1' )

      end subroutine verlet1
    
      subroutine verlet2(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp)
C *************************************************************************
C Subroutine for MD simulations using the velocity-Verlet Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 81)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to zero velocity
C                                 components opposite to force)
C                              2 = 'fire' quenching, as in
C                                  Bitzek et al, PRL 97, 170201 (2006)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************** OUTPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature
C *****************************************************************************

      integer 
     .   natoms,ntcon,istep,iquench,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      logical
     .  found

      integer
     .  ct,i,ia

      integer :: old_natoms, iacc, dummy_iza
      real(dp) :: old_dt

      real(dp)
     .  dot,dt2,dtby2,fovermp,temp

      real(dp), pointer, save :: accold(:,:)=>null(), vold(:,:)=>null()

C Related to FIRE quenching
      integer, parameter  :: firenmin = 5
      integer, save       :: firenpos
      real(dp), parameter :: firefinc = 1.1_dp, firefdec = 0.5_dp,
     .                       firealf0 = 0.1_dp, firefalf = 0.99_dp,
     .                       dtmax = 10.0_dp
      real(dp), save      :: firealf
      real(dp)            :: magv, magf
C ........................

      restart_file = trim(slabel) // '.VERLET_RESTART'
      if (iunit .ne. 1 .and. iunit .ne. 2)
     $     call die('verlet2: Wrong iunit option;  must be 1 or 2')

      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.associated(accold)) then
        call re_alloc( accold, 1, 3, 1, natoms, 'accold', 'verlet2' )
      endif
      if (.not.associated(vold)) then
        call re_alloc( vold, 1, 3, 1, natoms, 'vold', 'verlet2' )
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C ........................

      
      if (istep .eq. 1) then

C Initialise FIRE quench if that is the option
         if (iquench .eq. 2) then
            firealf = firealf0
            firenpos = 0
         endif

         if (xv_file_read) then
           if (IONode) inquire( file=restart_file, exist=found )
           call broadcast(found)
           if (.not. found) old_dt=dt
         else
           found=.false.
         endif

         if (.not. found) then

C     Compute old accelerations and velocities 
C     if the time step is the first of the simulation ...........................
!     and we start from x(t), v(t) *at the same time*.
!     (e.g., when the velocities are constructed from
!      the Boltzmann distribution).
!     In this case the algorithm works out well.

            do ia = 1,natoms
               do i = 1,3
                  accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
                  vold(i,ia) = va(i,ia) - dt * accold(i,ia)
               enddo
            enddo

         else

!         For restarts, we need information about the old 
!         forces, in order to match the velocities
!         correctly (the velocities in the XV file are
!         one time step behind, so they are already the
!         'old' velocities).
!
          if (Node .eq. 0) then
           call io_assign(iacc)
           open(unit=iacc,file=restart_file, form='formatted',
     $          status='old', action='read', position='rewind')
           read(iacc,*) old_natoms, old_dt
           if (old_natoms .ne. natoms)
     $          call die('Wrong number of atoms in VERLET_RESTART')
           do ia = 1, natoms
              read(iacc,*) dummy_iza, (accold(i,ia),i=1,3) ! forces
              if (dummy_iza .ne. iza(ia)) 
     $             call die('Wrong species number in VERLET_RESTART')
              accold(:,ia) = fovermp * accold(:,ia) / ma(ia)
              vold(:,ia)  = va(:,ia)
           enddo
           call io_close(iacc)
           write(6,*) 'MD restart: Read old forces from VERLET_RESTART'
           if (abs(old_dt - dt) .gt. 1.0d-8) then
              write(6,*) 'Timestep has changed. Old: ', old_dt,
     $             ' New: ', dt
           endif

          endif             ! IONode

          call broadcast(old_dt)
          call broadcast(accold(1:3,1:natoms))
          call broadcast(vold(1:3,1:natoms))

       endif      ! XV file read

      endif    ! first step

C ....................
C Compute velocities at current time step, 
! using the previous step's velocities and the previous and current forces.

      if ((istep .eq. 1) .and. xv_file_read) then
!
!        Use old time step in case it is different, only in 
!        the first step.
!
         do ia = 1,natoms
            va(:,ia) = vold(:,ia) + 0.5_dp * old_dt
     .           * (accold(:,ia) + fovermp * fa(:,ia) / ma(ia))
         enddo
         
      else
!
!        Current timestep.
!
         do ia = 1,natoms
            va(:,ia) = vold(:,ia) + dtby2 
     .           * (accold(:,ia) + fovermp * fa(:,ia) / ma(ia))
         enddo
      endif    

C Quench option if iquench = 1 ..............................................
      if (iquench .eq. 1) then

C Quench velocity components going uphill
         do ia = 1,natoms
           do i = 1,3
             dot = va(i,ia) * fa(i,ia)
             if (dot .lt. 0.0) va(i,ia) = 0.0
           enddo
         enddo

      endif
C ................

C Fire quench option if iquench = 2 .........................................
      if (iquench .eq. 2) then

C Compute power
         dot = 0.0
         do ia = 1,natoms
           do i = 1,3
             dot = dot + va(i,ia) * fa(i,ia)
           enddo
         enddo

         if (dot .le. 0.0) then
C If uphill: quench, reduce time step & go back to initial damping alpha
            va(:,:) = 0.0
            dt = firefdec * dt
            firealf = firealf0
            firenpos = 0
! FIXME:
! Need to set magv and magf to avoid writing undef values.
! Set to zero here - but does it make more sense to do the
! double summation below for both uphill and downhill cases? 
! AMW - 3 / 7 / 2008.
            magv = 0.0
            magf = 0.0
         else
C If downhill: Compute magnitudes of v and F, and ...
            magv = 0.0
            magf = 0.0
            do ia = 1,natoms
              do i = 1,3
                magv = magv + va(i,ia) * va(i,ia)
                magf = magf + fa(i,ia) * fa(i,ia)
              enddo
            enddo
            magv = sqrt( magv )
            magf = sqrt( magf )
C ... damp velocities, and ...
            do ia = 1,natoms
              do i = 1,3
                va(i,ia) = (1.0_dp - firealf)*va(i,ia) + 
     .                     firealf * magv * fa(i,ia) / magf
              enddo
            enddo
C ... increase dt and decreas firealf if Npos > Nmin
            if ( firenpos .gt. firenmin ) then
               dt = dt * firefinc
               if (dt .gt. dtmax) dt = dtmax
               firealf = firealf * firefalf
            endif
            firenpos = firenpos + 1
         endif

         if (Node .eq. 0)
     $     write(6,"('FIRE: istep, dt, firenpos, firealf, magv'/
     .           'FIRE: ', i6, f10.4, i4, f10.4, f12.4)")
     .            istep, dt, firenpos, firealf, magv

      endif
C ................



!     This is the place to store the current magnitudes
!
!      call add_to_md_file(xa,va,cell=h,vcell=hdot)


C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = xa(i,ia) + dt * va(i,ia) 
     .                  + dt2 / 2.0_dp * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Save current velocities and accelerations as old ones .....................
      do i = 1,3
        do ia = 1,natoms
          vold(i,ia) = va(i,ia)
          accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms 
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5_dp * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

!
!       Save (now old) forces to VERLET_RESTART
!
      if (Node .eq. 0) then
         call io_assign(iacc)
         open(unit=iacc,file=restart_file, form='formatted',
     $        status='unknown', action= 'write', position='rewind')
         write(iacc,*) natoms, dt
         do ia = 1, natoms
            write(iacc,*) iza(ia), (fa(i,ia),i=1,3) ! forces
         enddo
         call io_close(iacc)
      endif
C .....................

      end subroutine verlet2

      end module m_dynamics

    
