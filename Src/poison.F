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
      subroutine poison( CELL, N1, N2, N3, Mesh, RHO, U, V, STRESS, 
     .                   NSM )

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

C
C  Modules
C
      use precision,   only : dp, grid_p
      use parallel,    only : Node, Nodes, ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

#ifdef MPI
      integer           MPIerror
#endif

      integer           N1, N2, N3, Mesh(3), NSM
      real(grid_p)      RHO(N1*N2*N3), V(N1*N2*N3)
      real(dp)          CELL(3,3), STRESS(3,3), U

      real(grid_p), dimension(:,:), pointer  :: CG

      integer           I, I1, I2, I3, IX, J, J1, J2, J3, JX,
     .                  NP, NG, NG2, NG3, 
     .                  ProcessorZ, Py, Pz, J2min, J2max,
     .                  J3min, J3max, J2L, J3L, NRemY, NRemZ,
     .                  BlockSizeY, BlockSizeZ
      real(dp)          C, B(3,3), DU, G(3), G2, G2MAX, 
     .                  K0(3), PI, TINY, VG, VOLUME, VOLCEL
      externaL          CHKGMX, FFT, RECLAT, VOLCEL
      save              K0, TINY

      data K0   /3*0.0_dp/
      data TINY /1.0e-15_dp/

#ifdef MPI
C Place MPI barrier here so that any waiting time for all
C nodes to reach this point isn't included in the time for
C the FFT which is otherwise misleading.
      call MPI_Barrier(MPI_Comm_World, MPIerror)
#endif

C Start timer
      call timer('POISON',1)

C Allocate local memory
      nullify( CG )
      call re_alloc( CG, 1, 2, 1, n1*n2*n3, name='CG', 
     &                 routine='poison' )

C Find unit cell volume
      VOLUME = VOLCEL( CELL )

C Find reciprocal lattice vectors
      call reclat(CELL, B, 1 )

C Find maximun planewave cutoff
      NP = N1 * N2 * N3
      G2MAX = 1.0e30_dp
      call chkgmx( K0, B, Mesh, G2MAX )

C Copy density to complex array
      do I = 1, NP
        CG(1,I) = RHO(I)
        CG(2,I) = 0.0_grid_p
      enddo

C Find fourier transform of density
      call fft( CG, Mesh, -1, NSM )
 
C Initialize stress contribution
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = 0.0_dp
        enddo
      enddo

C Work out processor grid dimensions
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1

C Multiply by 8*PI/G2 to get the potential
      PI = 4.0_dp * atan(1.0_dp)
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

      do J3 = J3min,J3max
        if (J3.gt.NG3/2) then
          I3 = J3 - NG3
        else
          I3 = J3
        endif
        do J2 = J2min,J2max
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
            J2L = J2 - J2min 
            J3L = J3 - J3min 
            J = 1 + J1 + N1 * J2L + N1 * N2 * J3L
            if (G2.LT.G2MAX .AND. G2.GT.TINY) then
              VG = 8.0_dp * PI / G2
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

      NG = Mesh(1)*Mesh(2)*Mesh(3)
      U = 0.5_dp * U * VOLUME / DBLE(NG)**2
      C = 0.5_dp / DBLE(NG)**2
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = C * STRESS(JX,IX)
        enddo
        STRESS(IX,IX) = STRESS(IX,IX) + U / VOLUME
      enddo
 
C Go back to real space
      call fft( CG, Mesh, +1, NSM )

C Copy potential to array V
      do I = 1, NP
        V(I) = CG(1,I)
      enddo
 
C Stop time counter
  999 continue
      call timer('POISON',2)

C Free local memory
      call de_alloc( CG, name='CG' )

      return
      end

      subroutine nfft( N )

      use sys, only: die
      use parallel, only: ionode

C CHANGES N INTO THE NEXT INTEGER ALLOWED BY THE FFT ROUTINE
C WRITTEN BY J.M.SOLER. MAY/95.

      integer :: np, nmax, nmin, n, nrem, ip
      parameter (NP = 3, NMAX = 1000000)
      integer :: IPRIME(NP)
      data IPRIME / 2, 3, 5 /

      NMIN = N
      DO N = NMIN, NMAX
        NREM = N
        DO IP = 1,NP
   10     CONTINUE 
          IF ( MODULO( NREM, IPRIME(IP) ) .EQ. 0 ) THEN
            NREM = NREM / IPRIME(IP)
            GOTO 10
          ENDIF
        ENDDO
        IF (NREM .EQ. 1) RETURN
      ENDDO
      if (Ionode)
     $    write(6,*) 'NFFT: NO SUITABLE INTEGER FOUND FOR N =', NMIN
      call die()
      END
