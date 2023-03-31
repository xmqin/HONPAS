! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine vmb(nat,ttemp,mass,xa,isa,va)
C *************************************************************************
C This routine creates a velocity distribution according to the
C Maxwell-Boltzman distribution with a given temperature.
C It also imposes the constraint of zero total velocity 
C (to avoid center of mass drifts).
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C Modified by P. Ordejon to take into account constraints, June 2003
C ************************** INPUT ****************************************
C integer nat              : Number of atoms
C real*8  ttemp            : Temperature desired for the distribution (K)
C real*8  mass(nat)        : Atomic masses (in amu)
C real*8  xa(3,nat)        : Atomic positions
C integer isa(nat)         : Species indexes
C ************************* OUTPUT ****************************************
C real*8  va(3,nat)        : Atomic Velocities
C *************************************************************************
C
C  Modules
C
      use precision,   only : dp
      use parallel,    only : Node
      use m_fixed,     only : fixed
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer  MPIerror
#endif

      integer  nat, isa(nat)

      real(dp)  mass(nat), ttemp, va(3,nat), xa(3,nat)

C Internal variables ..................

      integer  iseed, i, ix, ntcon

      real(dp)
     .  cell(3,3), stress(3,3), fa(3,nat), 
     .  massi, tempe, velo, pcm(3), mtot,
     .  tmpstress(3,3), tmpfa(3,nat)

      external velo
C .....................

C Only generate the velocities on one Node in case they are truely random
      if (Node.eq.0) then

        iseed=-17

C Loop over atoms to assign velocities .................
        mtot = 0.0_dp
        pcm = 0.0_dp
        do i = 1,nat
          massi = mass(i)
          mtot = mtot + massi
           
          va(1,i) = velo(iseed,ttemp,massi)
          va(2,i) = velo(iseed,ttemp,massi)
          va(3,i) = velo(iseed,ttemp,massi)

          pcm(:) = pcm(:) + massi * va(:,i)
        enddo

C Impose constraint on zero center of mass velocity ....................
        do i = 1,nat
          do ix = 1,3
            va(ix,i) = va(ix,i) - pcm(ix)/mtot
          enddo
        enddo

      endif


C Distribute the velocities over the Nodes

#ifdef MPI
      call MPI_Bcast(va(1,1),3*nat,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
#endif


C Impose other constraints for MD, as specified in fixed.F ..................
C Note that fixed.F imposes constraints on forces; here we
C need those on velocities, so must transform multiplying by atomic mass
C

C Zero dummy matrices for stress and cell (needed in
C call to constrains routine)

      cell = 0.0_dp
      stress = 0.0_dp
      do i = 1,nat
        do i x =1,3
          fa(ix,i) = va(ix,i)*mass(i)
        enddo
      enddo

      ! NAG's F95 compiler does not like reuse of an
      ! argument with intent(in) and intent(out), such
      ! as stress and fa in the call to fixed. Quick 
      ! workaround is to make the use of a temporary
      ! variable thus:
      ! The "forces" are already mass-scaled,
      ! so regardless of the dynamics the call is correct.
      call fixed( cell, stress, nat, isa, mass, xa, fa,
     .                  tmpstress, tmpfa, ntcon )
      stress = tmpstress
      fa = tmpfa
      ! But is it better to use intent(inout) in fixed? 
      ! If this is done, remember to remove the two tmp
      ! variable definitions.

      do i = 1,nat
        do ix = 1,3
          va(ix,i) = fa(ix,i)/mass(i)
        enddo
      enddo

      if (nat .le. 1) return

C Correct velocity to exact temperature .....................
      call temp(2,Nat,mass,va,ntcon,tempe)

      if (abs(tempe-ttemp) .ge. 1.e-4 .and. tempe .ge. 1.e-4) then
        do i = 1,Nat
          do ix = 1,3
            va(ix,i) = va(ix,i) * dsqrt(ttemp/tempe)
          enddo
        enddo
      endif

      return
      end



      function velo(iseed,temp,mass) 
C *************************************************************************
C This function assigns velocities to atoms according to the 
C Maxwell-Boltzman distribution.
C It generates random numbers drawn from the normal distribution,
C using as input random numbers distributed uniformly from 0 to 1,  
C which are provided by routine randomg.
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C ************************** INPUT ****************************************
C integer iseed            : Seed for random number generator
C real*8  temp             : Temperature desired for the distribution (K)
C real*8  mass             : Atomic mass (in amu)
C ************************* OUTPUT ****************************************
C real*8  velo             : Velocity component 
C *************************************************************************

      use precision

      implicit none

      integer 
     . iseed

      real(dp)
     . mass,temp

      real(dp) :: velo

C Internal variables .................

      real(dp)
     .  arg1, arg2, gauss, med, pi, randomg, var
      
      external
     .  randomg

C ...........

C  For other distributions med may be different from cero.
      med = 0.0
      pi = acos(-1.0d0)
      var = sqrt(temp/mass)
C  conversion factor to bohr/fs
      var = var * 0.00172309d0
 
      arg1 = sqrt((-2.) * log(randomg(iseed)))

      arg2 = 2.0d0 * pi * randomg(iseed)
      gauss = arg1 * cos(arg2)

      velo = med + var * gauss

      return
      end



      subroutine temp(iunit,natoms,ma,va,ntcon,tempe)
C *************************************************************************
C Subroutine to calculate instantaneous temperature
C
C Written by P.Ordejon, November'97
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Atomic masses in atomic mass units
C
C Space units depend on input option:
C
C   Option iunit = 1:
C     Distances are in Angstrom
C   Option iunit = 2:
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 ma(natoms)     : Atomic masses 
C real*8 va(3,natoms)   : Atomic velocities 
C integer               : Total number of position constraints imposed
C ************************* OUTOPUT *******************************************
C real*8 tempe          : Instantaneous system temperature 
C *****************************************************************************

      use precision, only : dp
      use sys,       only : die

      implicit none

      integer 
     .   natoms,iunit,ntcon

      real(dp)
     .  ma(natoms),tempe,va(3,natoms)

C Internal variables ..........................................................
 
      integer
     .  i,ia

      real(dp)
     .  Ang,eV,fovermp,kin

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2)
     $     call die('temp: Wrong iunit option;  must be 1 or 2')

C Define constants and conversion factors .....................................

      Ang = 1.0d0 / 0.529177d0
      eV  = 1.0d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038d0
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038d0 *Ang**2 / eV
      endif

C ........................

C Calculate kinetic energy and temperature ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        tempe = 2.0d0*kin/(3.0d0*natoms-3.0d0-ntcon)/8.617d-5
      else
        tempe = 2.0d0*kin/(3.0d0*natoms-3.0d0-ntcon)/8.617d-5/eV
      endif
C .....................

      return
      end
    
