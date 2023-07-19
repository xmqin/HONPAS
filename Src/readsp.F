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
      subroutine readsp(q, spiral)
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Writen by J. Ferrer. June 2002
C **************************** OUTPUT *********************************
C real*8 q(3)              : Pitch wave vector for spiral configuration
C logical spiral           : True = Spiral arrangement of spins
C **********************************************************************

C
C  Modules
C
      use precision
      use parallel,   only : Node
      use sys,        only : die
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      real(dp)
     .  q(3)

      logical
     .  spiral

      external
     .  redcel

C  Internal variables .................................................

      character
     .  lattice*30

      integer 
     .  iu, mscell(3,3)

#ifdef MPI
      integer
     .  MPIerror
#endif

      logical
     .  lSpinSpiral

      real(dp)
     .  pi, alat, ucell(3,3), scell(3,3), rcell(3,3)

C Read lattice constant and unit cell vectors .........................

      if (Node.eq.0) then
        call redcel( alat, ucell, scell, mscell )
        lSpinSpiral = fdf_block('SpinSpiral',iu)
      endif
#ifdef MPI
      call MPI_Bcast(alat,1,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(ucell(1,1),9,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
      call MPI_Bcast(lSpinSpiral,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C...

C Spiral arrangement for spins ........................................
      pi = 2.d0 * asin(1.d0)

      if ( lSpinSpiral ) then
         spiral = .true.
         read(iu,*) lattice
         read(iu,*) q(1), q(2), q(3)
         if ( lattice .eq. 'Cubic')  then
            q(1) = pi * q(1) / alat
            q(2) = pi * q(2) / alat
            q(3) = pi * q(3) / alat
         else if ( lattice .eq. 'ReciprocalLatticeVectors' ) then
            call reclat( ucell, rcell, 1 )
            q(1) = rcell(1,1)*q(1) + rcell(1,2)*q(2) + rcell(1,3)*q(3)
            q(2) = rcell(2,1)*q(1) + rcell(2,2)*q(2) + rcell(2,3)*q(3)
            q(3) = rcell(3,1)*q(1) + rcell(3,2)*q(2) + rcell(3,3)*q(3)
         else
            write(6,*) 'redata: ERROR: ReciprocalCoordinates must be
     .      ''Cubic'' or ''ReciprocalLatticeVectors'' '
          call die('redata: ERROR: ReciprocalCoordinates for Spiral')
         endif
      else
         spiral = .false.
         q(1) = 0.d0
         q(2) = 0.d0
         q(3) = 0.d0
         lattice = 'unpolarized or collinear spins'
      endif

c     write(6,*) 'redata: Basis for vector q:    ', lattice
c     write(6,'(a,4(f9.2,3x))')
c    .  'redata: q(3)    ', q(1), q(2), q(3)

C...

      return
      end
