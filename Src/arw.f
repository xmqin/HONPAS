! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!***********************************************************************
!    This file arw.F contains modified versions of routines written by 
!    A.R.Williams around 1985 (possibly based on previous work).
!    They are used for the solution of the radial Schrodinger's and
!    Poisson's equations in the atomic program. It also contains some
!    routines written by J.M.Soler in collaboration with A.R.Williams
!    and based on algorithms developed by ARW.
!    J. M. Soler, Nov. 1998 and April 2003
!***********************************************************************
! The routines contained in this file are:
!     SUBROUTINE EGOFV
!     SUBROUTINE YOFE
!     SUBROUTINE NRMLZG
!     SUBROUTINE BCORGN
!     SUBROUTINE BCRMAX
!     SUBROUTINE NUMIN
!     SUBROUTINE NUMOUT
!     SUBROUTINE VHRTRE
!     SUBROUTINE QVLOFZ  ---> moved to periodic_table.f
!     SUBROUTINE LMXOFZ  ---> moved to periodic_table.f
!     SUBROUTINE CNFIG   ---> moved to periodic_table.f
!***********************************************************************

      SUBROUTINE EGOFV(H,S,NR,E,G,Y,L,Z,A,B,RMAX,NPRIN,NNODE,DGDR)
!***********************************************************************
! Finds the radial atomic wavefunction and energy for a given potential, 
! number of nodes, and logarithmic derivative at the boundary.
! Input:
!   real*8  H(NR) : Radial hamiltonian
!   real*8  S(NR) : Metric function. H and S are defined for a radial
!                   variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                   where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                   For a logarithmic mesh (see below), they are
!                    S(j) = (A*r(j))^2
!                    H(j) = S(j)*V(r(j)) + A^2/4
!                   where V(r) is the total effective radial potential
!                   in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer NR    : Number of radial points (including r(1)=0)
!   integer L     : Angular momentum
!   real*8  Z     : Atomic number
!   real*8  A,B   : Log-mesh parameters:
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,NR
!   real*8  RMAX  : Maximum radius. Must be equal to r(NR)
!   integer NPRIN : Principal quantum number
!   integer NNODE : Required number of nodes for wavefunction
!   real*8  DGDR  : Required logarithmic derivative dlogG/dr at RMAX
! Output:
!   real*8  E     : Solution energy
!   real*8  G(NR) : Radial wavefunction
! Auxiliary:
!   real*8  Y(NR) : Intermediate array used in the Numerov method
!***********************************************************************
      use parallel,  only: ionode
      use sys,       only: die
      use precision, only: dp

      IMPLICIT NONE
      INTEGER          :: NR, L, NPRIN, NNODE
      real(dp)         :: H(NR), S(NR), E, G(NR), Y(NR),
     &                    Z, A, B, RMAX, DGDR

      INTEGER,         PARAMETER :: MAXITER = 40    ! Max. iterations
      real(dp)        ,PARAMETER :: TOL     = 1.D-5 ! Convergence tol.
      INTEGER          :: I, NCOR, NITER, NN, NN1, NN2
      real(dp)         :: DE, DE12, E1, E2, T

      NCOR=NPRIN-L-1
      NN=NNODE
      NN1=NNODE
      NN2=NNODE-1
      E1=E
      E2=E
!     Initial bracketing energy range
      DE12=0.5D0
      DE=0.D0
      DO NITER = 1,MAXITER+1
!       New energy estimate from YOFE
        E=E+DE
        IF (E.GT.E1 .AND. E.LT.E2 .AND.
     &      NN.GE.NNODE-1 .AND. NN.LE.NNODE) THEN
!         New energy and number of nodes are within desired ranges
          IF (ABS(DE).LT.TOL) EXIT ! NITER loop converged
        ELSE
!         Use bisection to force energy within range
          E=0.5D0*(E1+E2)
        END IF
!       Find wavefunction Y for energy E, and new energy estimate E+DE
        CALL YOFE(E, DE, DGDR, RMAX, H, S, Y, NR, L, NCOR, NN, Z, A, B)
        IF (NN.LT.NNODE) THEN
!         Too few nodes. Increase lower energy bound
          E1=E
          NN1=NN
          IF (NN2.LT.NNODE) THEN
!           Energy not yet bracketed. Increase higher energy bound
            DE12=DE12*2.D0
            E2=E1+DE12
          END IF
        ELSE
!         Too many nodes. Decrease higher energy bound
          E2=E
          NN2=NN
          IF (NN1.GE.NNODE) THEN
!           Energy not yet bracketed. Decrease lower energy bound
            DE12=DE12*2.D0
            E1=E2-DE12
          END IF
        END IF
      END DO ! NITER

!     No-convergence error exception
      IF (NITER.GT.MAXITER) THEN
        IF (IOnode) WRITE(6,'(A,/,A,F3.0,2(A,I2),2(A,F12.5))')
     &   ' EGOFV: ERROR: Too many iterations. Stopping.',
     &   ' Z=',Z,'  L=',L,'  NNODE=',NNODE,'  E=',E,'  DE=',DE
        call die()
      END IF

!     Find true waveftn G from auxiliary function Y and normalize
      G(1) = 0.D0
      DO I=2,NR
        T=H(I)-E*S(I)
        G(I)=Y(I)/(1.D0-T/12.D0)
      END DO
      CALL NRMLZG(G,S,NR)
      END SUBROUTINE EGOFV


      SUBROUTINE YOFE( E, DE, DGDR, RMAX, H, S, Y, NR, L,
     &                 NCOR, NNODE, Z, A, B )
!***********************************************************************
! Integrates the radial Scroedinger equation for a given energy
! Input:
!   real*8  E     : Energy of the required wavefunction
!   real*8  DGDR  : Required logarithmic derivative dlogG/dr at RMAX
!   real*8  RMAX  : Maximum radius. Must be equal to r(NR)
!   real*8  H(NR) : Radial hamiltonian
!   real*8  S(NR) : Metric function. H and S are defined for a radial
!                   variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                   where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                   For a logarithmic mesh (see below), they are
!                     S(j) = (A*r(j))^2
!                     H(j) = S(j)*V(r(j)) + A^2/4
!                   where V(r) is the total effective radial potential
!                   in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer NR    : Number of radial points (including r(1)=0)
!   integer L     : Angular momentum
!   integer NCOR  : Number of core states of same L below the given one
!   integer NNODE : Required number of nodes for wavefunction
!   real*8  Z     : Atomic number
!   real*8  A,B   : Log-mesh parameters:
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,NR
! Output:
!   real*8  DE    : Estimate of energy change required to eliminate kink
!   real*8  Y(NR) : Numerov function, related to G above by
!                     Y(j)=G(j)*(12-H(j)+E*S(j))/12
!   integer NNODE : Actual number of nodes of wavefunction
!***********************************************************************
      use precision, only: dp
      IMPLICIT NONE
      INTEGER,           intent(in) :: NR, L, NCOR
      INTEGER,          intent(out) :: NNODE
      real(dp)        ,  intent(in) :: E, DGDR, RMAX, H(NR), S(NR),
     &                                 Z, A, B
      real(dp)        , intent(out) :: DE, Y(NR)

      real(dp)        ,PARAMETER:: TMAX  =1.D0 ! Max negative kin engy
      real(dp)        ,PARAMETER:: DLGMAX=1.D3 ! Max log deriv at Rmax
      INTEGER          :: KNK, NR0, NNDIN
      real(dp)         :: GIN, GOUT, GSG, GSGIN, GSGOUT, RATIO, 
     &                    T, XIN, XOUT, Y2, YN, ZDR

      ! Find where the negative kinetic energy H-ES becomes too large
      ! and make Y=0 beyond that point
      DO NR0 = NR,1,-1
        IF ( H(NR0)-E*S(NR0) .LT. TMAX ) EXIT
        Y(NR0)=0.D0
      END DO ! NR0

      ! Find boundary condition at origin
      ZDR = Z*A*B
      CALL BCORGN(E,H,S,L,ZDR,Y2)

      ! Integrate Schroedinger equation outwards
      KNK=NR0  ! A new value for the kink position KNK will be returned
      CALL NUMOUT(E, H, S, Y, NCOR, KNK, NNODE, Y2, GOUT, GSGOUT, XOUT)

      ! Find if kinetic energy is sufficiently non negative to use
      ! Numerov at Rmax. Otherwise use zero-value boundary condition
      IF (NR0.EQ.NR .AND. ABS(DGDR).LE.DLGMAX) THEN
        CALL BCRMAX(E,DGDR,RMAX,H,S,NR0,YN,A,B)
      ELSE
        YN=0.D0
      END IF

      ! Integrate Schroedinger equation inwards
      CALL NUMIN(E,H,S,Y,NR0,NNDIN,YN,GIN,GSGIN,XIN,KNK)

      ! Make wavefunction continuous at R(KNK)
      RATIO = GOUT/GIN
      Y(KNK:NR0) = Y(KNK:NR0)*RATIO
      XIN = XIN*RATIO
      GSGIN = GSGIN*RATIO**2

      ! Estimate the energy change required to eliminate kink
      GSG=GSGOUT+GSGIN
      T=H(KNK)-E*S(KNK)
      DE=GOUT*(XOUT+XIN+T*GOUT)/GSG

      ! Find the number of nodes
      NNODE=NNODE+NNDIN
      IF (DE.LT.0.D0) NNODE=NNODE+1
      END SUBROUTINE YOFE



      SUBROUTINE NRMLZG(G,S,N)
!***********************************************************************
! Normalizes the radial wavefunction, using Simpson's rule for the norm
! Input:
!   real*8  G(N) : Wavefunction to be normalized. G is related to the
!                  true radial wavefunction psi by
!                     G(j) = (dr/dj)^(3/2) * r(j) * psi(r(j))
!   real*8  S(N) : Metric function defined for a logarithmic mesh
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,N
!                  as  S(j) = (dr/dj)^2 = (A*r(j))^2
!   integer N    : Number of radial points (including r(1)=0)
! Output:
!   real*8  G(N) : Normalized wavefunction
!***********************************************************************
      use alloc,     only: re_alloc, de_alloc
      use precision, only: dp
      IMPLICIT NONE
      INTEGER           :: N
      real(dp)          :: G(N), S(N)

      INTEGER           :: I
      real(dp)          :: NORM, SRNRM
      real(dp), POINTER :: gaux(:)

      nullify(gaux)
      call re_alloc( gaux, 1, N, name='gaux', routine='NRMLZG' )
      gaux = g*g

      call integrator(gaux,s,n,norm)

      call de_alloc( gaux, name='gaux', routine='NRMLZG' )

      ! Normalize wavefunction
      SRNRM = SQRT(NORM)
      DO I=1,N
         G(I) = G(I)/SRNRM
      END DO

      END SUBROUTINE NRMLZG


      SUBROUTINE BCORGN(E,H,S,L,ZDR,Y2)
!************************************************************************
! Finds the boundary condition at the origin (actually the wavefunction 
! at the second point) to integrate the radial Schroedinger equation.
! Input:
!   real*8  E    : Energy of the required wavefunction
!   real*8  H(3) : Radial hamiltonian (only the first 3 points are used)
!   real*8  S(3) : Metric function. H and S are defined for a radial
!                   variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                   where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                   For a logarithmic mesh (see below), they are
!                     S(j) = (A*r(j))^2
!                     H(j) = S(j)*V(r(j)) + A^2/4
!                   where V(r) is the total effective radial potential
!                   in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer L    : Angular momentum
!   real*8  ZDR  : Atomic charge times (dr/dj) at r=0
! Output:
!   real*8  Y2   : Value at j=2 of the Numerov function, related to G
!                  above by   Y(j)=G(j)*(12-H(j)+E*S(j))/12
! Notes:
! - Local variable D(j) is the inverse of the diagonal of the 
!   tri-diagonal Numerov matrix
! - For L=0, G(j) vanishes at the origin (j=1), but its first and second 
!   derivatives are finite, making Y(1) finite
! - For L=1, G and G' vanish at the origin, but G'' and Y are finite
! - For L>1, G, G', G'', and Y all vanish at the origin
!************************************************************************
      use precision, only: dp
      IMPLICIT NONE
      INTEGER          :: L
      real(dp)         :: E, H(3), S(3), ZDR, Y2

      real(dp)         :: C0, C1, C2, D2, T2, T3

      T2=H(2)-E*S(2)
      D2=-((24.D0+10.D0*T2)/(12.D0-T2))
      IF (L<2) THEN
        IF (L==0) THEN
          C0=ZDR/6.D0
          C0=C0/(1.D0-0.75D0*ZDR)
        ELSE
          C0=1.D0/12.D0
          C0=(-C0)*8.D0/3.D0
        END IF
        C1=C0*(12.D0+13.D0*T2)/(12.D0-T2)
        T3=H(3)-E*S(3)
        C2=(-5.D-1)*C0*(24.D0-T3)/(12.D0-T3)
        D2=(D2-C1)/(1.D0-C2)
      END IF
      Y2=(-1.D0)/D2

      END SUBROUTINE BCORGN



      SUBROUTINE BCRMAX(E,DGDR,RMAX,H,S,N,YN,A,B)
!************************************************************************
! Finds the boundary condition at Rmax (actually the wavefunction at the
! last point) to integrate the radial Schroedinger equation inwards.
! Input:
!   real*8  E     : Energy of the required wavefunction
!   real*8  DGDR  : Required logarithmic derivative dlogG/dr at RMAX
!   real*8  RMAX  : Maximum radius. Must be equal to r(N)
!   real*8  H(N+1): Radial hamiltonian (only the last 3 points are used)
!   real*8  S(N+1): Metric function. H and S are defined for a radial
!                   variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                   where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                   For a logarithmic mesh (see below), they are
!                     S(j) = (A*r(j))^2
!                     H(j) = S(j)*V(r(j)) + A^2/4
!                   where V(r) is the total effective radial potential
!                   in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer N     : Number of radial points (including r(1)=0)
!   real*8  A,B   : Log-mesh parameters:
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,NR
! Output:
!   real*8  YN   : Value at j=N of the Numerov function, related to G
!                  above by   Y(j)=G(j)*(12-H(j)+E*S(j))/12
!************************************************************************
      use precision, only: dp
      IMPLICIT NONE
      INTEGER          :: N
      real(dp)         :: E, DGDR, RMAX, H(N+1), S(N+1), YN, A, B

      real(dp)         :: BETA, C1, C2, C3, DG, DN, TN, TNM1, TNP1

      TNM1=H(N-1)-E*S(N-1)
      TN  =H(N  )-E*S(N  )
      TNP1=H(N+1)-E*S(N+1)
      BETA=1.D0+B/RMAX
      DG=A*BETA*(DGDR+1.D0-0.5D0/BETA)

      C2=24.D0*DG/(12.D0-TN)
      DN=-(24.D0+10.D0*TN)/(12.D0-TN)

      C1= (12.D0-2.D0*TNM1)/(12.D0-TNM1)
      C3=-(12.D0-2.D0*TNP1)/(12.D0-TNP1)
      YN=-(1.D0-C1/C3)/(DN-C2/C3)

      END SUBROUTINE BCRMAX



      SUBROUTINE NUMIN(E,H,S,Y,NR,NNODE,YN,G,GSG,DY,KNK)
!***********************************************************************
! Integrates Schroedinger's equation inwards, using Numerov's method
! Input:
!   real*8  E     : Energy of the required wavefunction
!   real*8  H(NR) : Radial hamiltonian
!   real*8  S(NR) : Metric function. H and S are defined for a radial
!                   variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                   where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                   For a logarithmic mesh (see below), they are
!                     S(j) = (A*r(j))^2
!                     H(j) = S(j)*V(r(j)) + A^2/4
!                   where V(r) is the total effective radial potential
!                   in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer NR    : Number of radial points (including r(1)=0)
!   real*8  YN    : Value at j=N of the Numerov function, defined below
!   integer KNK   : Fixes the 'kink point' r(KNK) up to which the 
!                   inward integration is required
! Output:
!   real*8  Y(NR) : Numerov function related to G(j) by
!                    Y(j)=(1-T(j)/12)*G(j), with T(j)=H(j)-E*S(j)
!   integer NNODE : Number of nodes of Y(j) between KNK and NR
!   real*8  G     : Value of G, defined above, at r(KNK)
!   real*8  GSG   : Probability G*S*G=(4*pi*r**2)*psi**2 at r(KNK)
!   real*8  DY    : Value of Y(j)-Y(j-1) at j=KNKout
!   integer KNK   : Made equal to NR-2 if the input value is larger
! Algorithm: see routine NUMOUT
!***********************************************************************
      use precision, only: dp
      IMPLICIT NONE
      INTEGER          :: NR, NNODE, KNK
      real(dp)         :: E, H(NR), S(NR), Y(NR), YN, G, GSG, DY

      INTEGER          :: I
      real(dp)         :: T

      Y(NR)=YN
      T=H(NR)-E*S(NR)
      G=Y(NR)/(1.D0-T/12.D0)
      GSG=G*S(NR)*G
      I=NR-1
      Y(I)=1.D0
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      DY=Y(I)-Y(NR)
      NNODE=0
      KNK = MIN(KNK,NR-2)
      DO I = NR-2,KNK,-1
        DY=DY+T*G
        Y(I)=Y(I+1)+DY
        IF( Y(I)*Y(I+1) .LT. 0.D0) NNODE=NNODE+1
        T=H(I)-E*S(I)
        G=Y(I)/(1.D0-T/12.D0)
        GSG=GSG+G*S(I)*G
      END DO

      END SUBROUTINE NUMIN


      SUBROUTINE NUMOUT( E, H, S, Y, NCOR, KNK, NNODE, Y2, G, GSG, DY )
!***********************************************************************
! Integrates Schroedinger's equation outwards, using Numerov's method,
! up to the first maximum of psi after NCOR nodes
! Input:
!   real*8  E      : Energy of the required wavefunction
!   real*8  H(KNK) : Radial hamiltonian
!   real*8  S(KNK) : Metric function. H and S are defined for a radial
!                    variable mesh so that d2G/dj2=(H(j)-E*S(j))*G(j)
!                    where G(j)=(dr/dj)^(3/2)*r(j)*psi(r(j))
!                    For a logarithmic mesh (see below), they are
!                      S(j) = (A*r(j))^2
!                      H(j) = S(j)*V(r(j)) + A^2/4
!                    where V(r) is the total effective radial potential
!                    in Rydbergs, including the kinetic term l*(l+1)/r^2
!   integer NCOR   : Number of core states of same L below the given one
!   integer KNK    : Number of radial points (including r(1)=0)
!   real*8  Y2     : Value at j=N of the Numerov function, defined below
! Output:
!   real*8  Y(KNK) : Numerov function related to G(j) by
!                      Y(j)=(1-T(j)/12)*G(j), with T(j)=H(j)-E*S(j)
!   integer NNODE  : One plus the number of nodes up to output KNK
!   real*8  G      : Value of G, defined above, at r(KNKout)
!   real*8  GSG    : Probability G*S*G=(4*pi*r**2)*psi**2 at r(KNKout)
!   real*8  DY     : Value of Y(j)-Y(j-1) at j=KNKout
!   integer KNK    : Position of first maximum of Y(j) after NCOR nodes,
!                    or KNK=KNKin-4 (whatever is lower)
! Algorithm:
! The Numerov equation is 
!   (g(j+1) - 2*g(j) + g(j-1)) / dx2 = 
!        (t(j+1)*g(j+1) + 10*t(j)*g(j) + t(j-1)*g(j-1)) / 12
! or, for dx=1,
!   (1-t(j+1)/12)*g(j+1) = (2+10*t(j)/12)*g(j) - ((1-t(j-1)/12)*g(j-1)
! Then, defining y(j)=(1-t(j)/12)*g(j) and after some simple algebra
!   dy(j) = y(j+1)-y(j) = t(j)*g(j) + y(j)-y(j-1) = t(j)*g(j) + dy(j-1)
!***********************************************************************
      use precision, only: dp
      IMPLICIT NONE
!     Input/Output variables
      INTEGER,           intent(in) :: NCOR
      INTEGER,        intent(inout) :: KNK
      INTEGER,          intent(out) :: NNODE
      real(dp)        ,  intent(in) :: E, H(KNK), S(KNK), Y2
      real(dp)        , intent(out) :: Y(KNK), G, GSG, DY
!     Local variables
      INTEGER                       :: I
      real(dp)                      :: DYL, T

      Y(1)  = 0.D0
      Y(2)  = Y2
      T     = H(2)-E*S(2)
      G     = Y(2)/(1.D0-T/12.D0)
      GSG   = G*S(2)*G
      Y(3)  = 1.D0
      T     = H(3)-E*S(3)
      G     = Y(3)/(1.D0-T/12.D0)
      GSG   = GSG+G*S(3)*G
      DY    = Y(3)-Y(2)
      NNODE = 0
      DO I= 4, KNK-4
        DYL  = DY
        DY   = DY+T*G
        Y(I) = Y(I-1)+DY
        IF( Y(I)*Y(I-1) .LT. 0.D0) NNODE = NNODE + 1
        T   = H(I)-E*S(I)
        G   = Y(I)/(1.D0-T/12.D0)
        GSG = GSG+G*S(I)*G
        ! End loop if |Y(j)| has a maximum after NCOR nodes
        IF (NNODE.GE.NCOR .AND. DYL*DY.LE.0.D0) EXIT
      END DO
      KNK = MIN(I,KNK-4)

      END SUBROUTINE NUMOUT



      SUBROUTINE VHRTRE(R2RHO,V,R,DRDI,SRDRDI,NR,A)
!***********************************************************************
! Finds the Hartree potential created by a radial electron density, 
! using Numerov's method to integrate the radial Poisson equation:
!   d2(r*V)/dr2 = -4*pi*rho*r = -(4*pi*r2*rho)/r
! Input:
!   real*8  R2RHO(NR) : 4*pi*r**2*rho, with rho the electron density
!   real*8  R(NR)     : Logarithmic radial mesh R(i)=B*(exp(A*(i-1)-1)
!   real*8  DRDI(NR)  : dr/di at the mesh points
!   real*8  SRDRDI(NR): sqrt(dr/di) at the mesh points
!   integer NR        : Number of radial mesh points, including r(1)=0
!   real*8  A         : The parameter A in r(i)=B*(exp(A*(i-1)-1)
! Output:
!   real*8  V(NR)     : Electrostatic potential created by rho, in Ryd
!                       The constants of integration are fixed so that
!                       V=finite at the origin and V(NR)=Q/R(NR),
!                       where Q is the integral of rho up to R(NR)
! Algorithm: see routine NUMOUT
!***********************************************************************
      use precision, only: dp
      use alloc,     only: re_alloc, de_alloc
      IMPLICIT NONE
      INTEGER           :: NR
      real(dp)          :: R2RHO(NR),V(NR),R(NR),DRDI(NR),SRDRDI(NR),A

      INTEGER           :: IR
      real(dp)          :: A2BY4, BETA, DV, DY, DZ, Q, QBYY, QPARTC, QT,
     .                     T, V0, Y, YBYQ
      real(dp), POINTER :: gaux(:)


      ! Find some constants
      A2BY4=A*A/4.D0
      YBYQ=1.D0-A*A/48.D0
      QBYY=1.D0/YBYQ

      ! Use Simpson's rule to find the total charge QT, and the 
      ! potential at the origin V0:
      ! QT = Int(4*pi*r**2*rho*dr) = Int((4*pi*r**2*rho)*(dr/di)*di)
      ! V0 = Int(4*pi*r*rho*dr) =  Int((4*pi*r**2*rho)/r*(dr/di)*di)

      call integrator(r2rho(2),drdi(2),nr-1,QT)

      nullify(gaux)
      call re_alloc( gaux, 2, NR, name='gaux', routine='VHRTRE' )

      gaux = r2rho(2:nr)/r(2:nr)
      call integrator( gaux, drdi(2), nr-1, V0 )

      call de_alloc( gaux, name='gaux', routine='VHRTRE' )

      ! Fix V(1) and V(2) to start Numerov integration. To find a 
      ! particular solution of the inhomog. eqn, V(2) is fixed by 
      ! setting rV(2)=0. Notice that V=finite => rV=0 at r=0
      V(1)=2.D0*V0    ! Factor 2 because we use Rydbergs
      T=SRDRDI(2)/R(2)
      BETA=DRDI(2)*T*R2RHO(2)
      DY=0.D0
      Y=0.D0
      Q=(Y-BETA/12.D0)*QBYY
      V(2)=2.D0*T*Q

      ! Integrate Poisson's equation outwards, using Numerov's method
      DO IR = 3,NR
        DY=DY+A2BY4*Q-BETA
        Y=Y+DY
        T=SRDRDI(IR)/R(IR)
        BETA=T*DRDI(IR)*R2RHO(IR)
        Q=(Y-BETA/12.D0)*QBYY
        V(IR)=2.D0*T*Q
      END DO

      ! Add a solution (finite at r=0) of the homogeneous equation 
      ! d2(r*V)/dr2=0 => rV=const*r => V=const, to ensure that 
      ! V(NR)=Q/R(NR). Notice that V(1) is set independently
      QPARTC = R(NR)*V(NR)/2.D0
      DZ=QT-QPARTC
      DV=2.D0*DZ/R(NR)
      DO IR=2,NR
         V(IR)=V(IR)+DV
      ENDDO

      END SUBROUTINE VHRTRE

      SUBROUTINE INTEGRATOR(F,S,NP,VAL)
!***********************************************************************
! Integrates a radial function tabulated on a logarithmic grid,
! using a generalized Simpson's rule valid for both even and odd
! number of points. Note that the "h" is 1 as the reparametrization
! involves a mapping of integers to reals.
!
! Alberto Garcia, Dec. 2006, based on code by Art Williams.
!
! Input:
!   real*8  F(NP) : Function to be integrated.
!   real*8  S(NP) : Metric function defined for a logarithmic mesh
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,NP
!                  as  S(j) = (dr/dj)^2 = (A*r(j))^2
!   integer NP    : Number of radial points (including r(1)=0)
! Output:
!   real*8  VAL   : Value of the integral
!***********************************************************************
      use precision, only: dp
      IMPLICIT NONE
      INTEGER          :: NP
      real(dp)         :: F(NP), S(NP)

      INTEGER          :: I, N
      real(dp)         :: VAL

      IF (MOD(NP,2).EQ.1) THEN        
         N = NP               ! ODD
      ELSE
         IF (NP .EQ. 2) THEN
          ! Special case of trapezoidal rule
            VAL = 0.5D0 * (F(1)*S(1) + F(2)*S(2))
            RETURN
         ENDIF
         N = NP - 3           ! EVEN: TAKE A FINAL FOUR-POINT INTERVAL
      ENDIF
!
!     STANDARD EXTENDED SIMPSON RULE WITH ALTERNATING 4/3 AND 2/3 FACTORS
!     FOR THE SECTION MADE UP OF PAIRS OF INTERVALS (THREE-POINT SEGMENTS)
!
      VAL = 0.D0
      DO I = 2,N-1,2
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL * 2.D0
      DO I = 3,N-2,2
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL * 2.D0
      DO I = 1,N,N-1                ! first and last points at 1/3
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL/3.D0

      IF (MOD(NP,2).EQ.0) THEN           ! EVEN
!        ADD THE CONTRIBUTION OF THE 
!        FINAL FOUR-POINT SEGMENT USING SIMPSON'S 3/8 RULE
!        (SEE NUMERICAL RECIPES (FORTRAN), P. 105)
         I = NP - 3
         VAL = VAL +
     $      (3.D0/8.D0) * ( (F(I)*S(I) + F(I+3)*S(I+3)) +
     $         3.D0 * (F(I+1)*S(I+1) + F(I+2)*S(I+2)) )
      ENDIF
      END SUBROUTINE INTEGRATOR
