! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine poison( CELL, N1, N2, N3, Mesh, RHO, U, V, STRESS, 
     &                   NSM )

C *********************************************************************
C SOLVES POISSON'S EQUATION.
C ENERGY AND POTENTIAL RETURNED IN RYDBERG UNITS.
C WRITTEN BY J.M.SOLER, JUNE 1995.
C **************** INPUT **********************************************
C REAL*8  CELL(3,3)     : UNIT CELL VECTORS
C INTEGER N1,N2,N3      : NUMBER OF MESH DIVISIONS IN EACH CELL VECTOR
C INTEGER Mesh(3)       : Number of global mesh divisions
C REAL*4  RHO(N1,N2,N3) : DENSITIY AT MESH POINTS
C **************** OUTPUT *********************************************
C REAL*8  U             : ELECTROSTATIC ENERGY (IN RY)
C REAL*4  V(N1,N2,N3)   : ELECTROSTATIC POTENTIAL (IN RY)
C                         V AND RHO MAY BE THE SAME PHYSICAL ARRAY
C REAL*8  STRESS(3,3) : Electrostatic-energy contribution to stress
C                       tensor (in Ry/Bohr**3) assuming constant density
C                       (not charge), i.e. r->r' => rho'(r') = rho(r)
C                       For plane-wave and grid (finite difference)
C                       basis sets, density rescaling gives an extra
C                       term (not included) equal to -2*U/cell_volume
C                       for the diagonal elements of stress. For other
C                       basis sets, the extra term is, in general:
C                       IntegralOf( V * d_rho/d_strain ) / cell_volume
C INTEGER NSM           : Number of sub-mesh points per mesh point
C                       : along each axis
C *********************************************************************

C     Modules
      use precision,   only : dp, grid_p
      use parallel,    only : Node, Nodes, ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc
      use m_fft,       only : fft     ! 3-D fast Fourier transform
      use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
      use cellsubs,    only : volcel  ! Finds unit cell volume
      use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh
#ifdef MPI
      use mpi_siesta
#endif
      implicit          none
C     Input/output variables
      integer               :: N1, N2, N3, Mesh(3), NSM
      real(grid_p)          :: RHO(N1*N2*N3), V(N1*N2*N3)
      real(dp)              :: CELL(3,3), STRESS(3,3), U
C     Local variables
      integer               :: I, I1, I2, I3, IX, J, J1, J2, J3, JX,
     &                         NP, NG, NG2, NG3, 
     &                         ProcessorZ, Py, Pz, J2min, J2max,
     &                         J3min, J3max, J2L, J3L, NRemY, NRemZ,
     &                         BlockSizeY, BlockSizeZ
      real(dp)              :: C, B(3,3), DU, G(3), G2, G2MAX, 
     &                         PI, VG, VOLUME, PI8
      real(grid_p), pointer :: CG(:,:)
      real(dp),   parameter :: K0(3)= (/ 0.0, 0.0, 0.0 /), TINY= 1.0e-15
#ifdef MPI
      integer               :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE POISON' )
#endif
#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 2 )
#endif

C     Start time counter
      call timer( 'POISON', 1 )

C     Allocate local memory
      nullify( CG )
      call re_alloc( CG, 1, 2, 1, n1*n2*n3, 'CG', 'poison' )

C     Find unit cell volume
      VOLUME = VOLCEL( CELL )

C     Find reciprocal lattice vectors
      call reclat(CELL, B, 1 )

C     Find maximun planewave cutoff
      NP = N1 * N2 * N3
      G2MAX = 1.0e30_dp
      call CHKGMX( K0, B, Mesh, G2MAX )

C     Copy density to complex array
!$OMP parallel do default(shared), private(I)
      do I = 1, NP
        CG(1,I) = RHO(I)
        CG(2,I) = 0.0_grid_p
      enddo
!$OMP end parallel do

C     Find fourier transform of density
      call fft( CG, Mesh, -1 )

C     Initialize stress contribution
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = 0.0_dp
        enddo
      enddo

C     Work out processor grid dimensions
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes) 
     &     call die('ERROR: ProcessorY must be a factor of the' //
     &     ' number of processors!')
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1

C     Multiply by 8*PI/G2 to get the potential
      PI  = 4.0_dp * atan(1.0_dp)
      PI8 = PI * 8._dp
      U = 0.0_dp
      NG2 = Mesh(2)
      NG3 = Mesh(3)
      BlockSizeY = ((NG2/NSM)/ProcessorY)*NSM
      BlockSizeZ = ((NG3/NSM)/ProcessorZ)*NSM
      NRemY = (NG2 - BlockSizeY*ProcessorY)/NSM
      NRemZ = (NG3 - BlockSizeZ*ProcessorZ)/NSM
      J2min = (Py-1)*BlockSizeY + NSM*min(Py-1,NRemY)
      J2max = J2min + BlockSizeY - 1
      if (Py-1.lt.NRemY) J2max = J2max + NSM
      J2max = min(J2max,NG2-1)
      J3min = (Pz-1)*BlockSizeZ + NSM*min(Pz-1,NRemZ)
      J3max = J3min + BlockSizeZ - 1
      if (Pz-1.lt.NRemZ) J3max = J3max + NSM
      J3max = min(J3max,NG3-1)

!$OMP parallel do default(shared), 
!$OMP&private(J3,J3L,I3,J2,J2L,I2,J1,I1,G,G2), 
!$OMP&private(J,VG,DU,C), reduction(+:U,STRESS)
      do J3 = J3min,J3max
        J3L = J3 - J3min 
        if (J3.gt.NG3/2) then
          I3 = J3 - NG3
        else
          I3 = J3
        endif
        do J2 = J2min,J2max
          J2L = J2 - J2min 
          if (J2.gt.NG2/2) then
            I2 = J2 - NG2
          else
            I2 = J2
          endif
          do J1 = 0,N1-1
            if (J1.gt.N1/2) then
              I1 = J1 - N1
            else
              I1 = J1
            endif
            G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
            G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
            G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3
            G2 = G(1)**2 + G(2)**2 + G(3)**2
            J = 1 + J1 + N1 * J2L + N1 * N2 * J3L
            if (G2.LT.G2MAX .AND. G2.GT.TINY) then
              VG = PI8 / G2
              DU = VG * ( CG(1,J)**2 + CG(2,J)**2 )
              U = U + DU
              C = 2.0_dp * DU / G2
              DO IX = 1,3
                DO JX = 1,3
                  STRESS(JX,IX) = STRESS(JX,IX) + C * G(IX) * G(JX)
                ENDDO
              ENDDO
              CG(1,J) = VG * CG(1,J)
              CG(2,J) = VG * CG(2,J)
            else
              CG(1,J) = 0.0_dp
              CG(2,J) = 0.0_dp
            endif
          enddo
        enddo
      enddo
!$OMP end parallel do

      NG = Mesh(1)*Mesh(2)*Mesh(3)
      U = 0.5_dp * U * VOLUME / DBLE(NG)**2
      C = 0.5_dp / DBLE(NG)**2
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = C * STRESS(JX,IX)
        enddo
        STRESS(IX,IX) = STRESS(IX,IX) + U / VOLUME
      enddo
 
C     Go back to real space
      call fft( CG, Mesh, +1 )

C     Copy potential to array V
!$OMP parallel do default(shared), private(I)
      do I = 1, NP
        V(I) = CG(1,I)
      enddo
!$OMP end parallel do
 
C     Free local memory
      call de_alloc( CG, 'CG', 'poison' )

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif

C     Stop time counter
      call timer( 'POISON', 2 )

#ifdef DEBUG
      call write_debug( '    POS POISON' )
#endif
      end subroutine poison
