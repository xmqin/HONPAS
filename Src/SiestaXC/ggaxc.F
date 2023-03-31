!!@LICENSE
!
!******************************************************************************
! MODULE m_ggaxc
! Provides routines for GGA XC functional evaluation
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! ggaxc,       ! General subroutine for all coded GGA XC functionals
! blypxc,      ! Becke-Lee-Yang-Parr (see subroutine blypxc)
! pbexc,       ! Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
! pbesolxc,    ! Perdew et al, PRL, 100, 136406 (2008)
! pw91xc,      ! Perdew & Wang, JCP, 100, 1290 (1994)
! revpbexc,    ! GGA Zhang & Yang, PRL 80,890(1998)
! rpbexc,      ! Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
! am05xc,      ! Mattsson & Armiento, PRB, 79, 155101 (2009)
! wcxc,        ! Wu-Cohen (see subroutine wcxc)
! pbeJsJrLOxc  ! Reparametrizations of the PBE functional by
! pbeJsJrHEGxc !   L.S.Pedroza et al, PRB 79, 201106 (2009) and
! pbeGcGxLOxc  !   M.M.Odashima et al, J. Chem. Theory Comp. 5, 798 (2009)
! pbeGcGxHEGxc ! using 4 different combinations of criteria
! pw86x,       ! Perdew & Wang, PRB 33, 8800 (1986) (exchange only)
! pw86rx,      ! pw86x refitted by Murray, Lee & Langreth, JCTC 5, 2754 (2009)
! b88x,        ! Becke, PRA 38, 3098 (1988) (exchange only)
! b88kbmx,     ! Klimes et al, JPCM 22, 022201 (2009) (exchange only)
! c09x,        ! Cooper, PRB 81, 161104 (2010) (exchange only)
! bhx          ! Berland & Hyldgaard, PRB 89, 035412 (2014) (exchange only)
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,     only: die     ! Termination routine
! use m_ldaxc, only: exchng  ! Local exchange
! use m_ldaxc, only: pw92c   ! Perdew & Wang, PRB, 45, 13244 (1992) (Corr. only)
!
!   USED module parameters:
! use precision,  only: dp   ! Real double precision type
!
!   EXTERNAL procedures used:
! none
!
!******************************************************************************

      MODULE m_ggaxc

      ! Used module procedures
      use sys,     only: die     ! Termination routine
      USE m_ldaxc, only: exchng  ! Local exchange
      USE m_ldaxc, only: pw92c   ! Perdew & Wang, PRB, 45, 13244 (1992) correl


      ! Used module parameters
      use precision, only : dp   ! Double precision real kind

      implicit none

      PUBLIC::  
     .  ggaxc,        ! General subroutine for all coded GGA XC functionals
     .  blypxc,       ! Becke-Lee-Yang-Parr (see subroutine blypxc)
     .  pbexc,        ! Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
     .  pbesolxc,     ! Perdew et al, PRL, 100, 136406 (2008)
     .  pw91xc,       ! Perdew & Wang, JCP, 100, 1290 (1994)
     .  revpbexc,     ! GGA Zhang & Yang, PRL 80,890(1998)
     .  rpbexc,       ! Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
     .  am05xc,       ! Mattsson & Armiento, PRB, 79, 155101 (2009)
     .  wcxc,         ! Wu-Cohen (see subroutine wcxc)
     .  pbeJsJrLOxc,  ! Reparametrizations of the PBE functional by
     .  pbeJsJrHEGxc, !   L.S.Pedroza et al, PRB 79, 201106 (2009) and
     .  pbeGcGxLOxc,  !   M.M.Odashima et al, JCTC 5, 798 (2009)
     .  pbeGcGxHEGxc, ! using 4 different combinations of criteria
     .  pw86x,        ! Perdew & Wang, PRB 33, 8800 (1986) (exchange only)
     .  pw86rx,       ! pw86x refit: Murray, Lee & Langreth, JCTC 5, 2754 (2009)
     .  b88x,         ! Becke, PRA 38, 3098 (1988) (exchange only)
     .  b88kbmx,      ! Klimes et al, JPCM 22, 022201 (2009) (exchange only)
     .  c09x,         ! Cooper, PRB 81, 161104 (2010) (exchange only)
     .  bhx           ! Berland & Hyldgaard, PRB 89, 035412 (2014) (exch. only)

      PRIVATE  ! Nothing is declared public beyond this point

      CONTAINS

      SUBROUTINE GGAXC( AUTHOR, IREL, nSpin, D, GD,
     .                  EPSX, EPSC, dEXdD, dECdD, dEXdGD, dECdGD )

C Finds the exchange and correlation energies at a point, and their
C derivatives with respect to density and density gradient, in the
C Generalized Gradient Correction approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec. 96.
C Modified by V.M.Garcia-Suarez to include non-collinear spin. June 2002
C Non collinear part rewritten by J.M.Soler. Sept. 2009

      implicit none

      ! Input
      character(len=*),intent(in):: AUTHOR  ! GGA flavour (author initials)
      integer, intent(in) :: IREL           ! Relativistic exchange? 0=no, 1=yes
      integer, intent(in) :: nSpin          ! Number of spin components
      real(dp),intent(in) :: D(nSpin)       ! Local electron (spin) density
      real(dp),intent(in) :: GD(3,nSpin)    ! Gradient of electron density

      ! Output
      real(dp),intent(out):: EPSX           ! Exchange energy per electron
      real(dp),intent(out):: EPSC           ! Correlation energy per electron
      real(dp),intent(out):: dEXdD(nSpin)   ! dEx/dDens, Ex=Int(dens*epsX)
      real(dp),intent(out):: dECdD(nSpin)   ! dEc/dDens
      real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)
      real(dp),intent(out):: dECdGD(3,nSpin) ! dEc/dGrad(Dens)

      ! Internal variables and arrays
      integer    :: NS, is, ix
      real(dp)   :: DD(2), dECdDD(2), dEXdDD(2),
     .              dDDdD(2,4), dDTOTdD(4), dDPOLdD(4),
     .              dECdGDD(3,2), dEXdGDD(3,2),
     .              dGDDdD(3,2,4), dGDDdGD(2,4),
     .              dGDTOTdD(3,4), dGDPOLdD(3,4),
     .              dGDTOTdGD(4), dGDPOLdGD(4),
     .              DPOL, DTOT, GDD(3,2), GDTOT(3), GDPOL(3),
     .              VPOL
      real(dp), parameter :: DENMIN = 1.e-12_dp

      real(dp):: theta, phi, c2, s2, st, ct, cp, sp, dpolz, dpolxy

      logical , parameter :: old_scheme = .true.

      ! Handle non-collinear spin case
      if (nSpin==4) then
        NS = 2             ! Diagonal spin components

        if ( old_scheme ) then
         DTOT = D(1) + D(2)
         dpolz= D(1) - D(2)
         dpolxy= 2.0d0*sqrt(d(3)**2+d(4)**2)
         dpol  = sqrt( dpolz**2 + dpolxy**2 )
         if ( dpol.gt.1.0d-12 ) then
          theta = atan2(dpolxy,dpolz)
         else
          theta = 0.0_dp
         endif
         C2 = COS(THETA/2.0_dp)
         S2 = SIN(THETA/2.0_dp)
         ST = SIN(THETA)
         CT = COS(THETA)
         PHI = ATAN2(-D(4),D(3))
         CP = COS(PHI)
         SP = SIN(PHI)

         DD(1) = 0.5D0 * ( DTOT + DPOL )
         DD(2) = 0.5D0 * ( DTOT - DPOL )
         DO IX = 1,3
          GDD(IX,1) = GD(IX,1)*C2**2 + GD(IX,2)*S2**2 +
     .                2.d0*C2*S2*(GD(IX,3)*CP - GD(IX,4)*SP)
          GDD(IX,2) = GD(IX,1)*S2**2 + GD(IX,2)*C2**2 -
     .                2.d0*C2*S2*(GD(IX,3)*CP - GD(IX,4)*SP)
         ENDDO

        else

        ! Find eigenvalues of density matrix Dij (diagonal densities DD, i.e.
        ! up and down densities along the spin direction). Note convention: 
        ! D(1)=D11, D(2)=D22, D(3)=Re(D12)=Re(D21), D(4)=Im(D12)=-Im(D21)
        DTOT = D(1) + D(2)                           ! DensTot (DensUp+DensDn)
        DPOL = SQRT( (D(1)-D(2))**2                  ! DensPol (DensUp-DensDn)
     .              + 4*(D(3)**2+D(4)**2) )
        DD(1) = ( DTOT + DPOL ) / 2                  ! DensUp
        DD(2) = ( DTOT - DPOL ) / 2                  ! DensDn
        DPOL = max( DPOL , DENMIN )                  ! Avoid division by zero

        ! Find gradients of up and down densities
        GDTOT(:) = GD(:,1) + GD(:,2)                 ! Grad(DensTot)
        GDPOL(:) = ( (D(1)-D(2))*(GD(:,1)-GD(:,2))   ! Grad(DensPol)
     .             + 4*(D(3)*GD(:,3)+D(4)*GD(:,4)) )
     .             / DPOL
        GDD(:,1) = ( GDTOT(:) + GDPOL(:) ) / 2       ! Grad(DensUp)
        GDD(:,2) = ( GDTOT(:) - GDPOL(:) ) / 2       ! Grad(DensDn)

        ! Derivatives of Dup and Ddn with respect to input density matrix
        dDTOTdD(1:2) = 1                             ! dDensTot/dD(i)
        dDTOTdD(3:4) = 0
        dDPOLdD(1) = +( D(1) - D(2) ) / DPOL         ! dDensPol/dD(i)
        dDPOLdD(2) = -( D(1) - D(2) ) / DPOL
        dDPOLdD(3) = 4 * D(3) / DPOL
        dDPOLdD(4) = 4 * D(4) / DPOL
        dDDdD(1,:) = ( dDTOTdD(:) + dDPOLdD(:) ) / 2 ! dDensUp/dD(i)
        dDDdD(2,:) = ( dDTOTdD(:) - dDPOLdD(:) ) / 2 ! dDensDn/dD(i)

        ! Derivatives of grad(Dup) and grad(Ddn) with respect to D(i)
        dGDTOTdD(1:3,1:4) = 0                        ! dGradDensTot/dD(i)
        dGDPOLdD(:,1) = + (GD(:,1)-GD(:,2))/DPOL     ! dGradDensPol/dD(i)
     .                  - GDPOL(:) * dDPOLdD(1)/DPOL
        dGDPOLdD(:,2) = - (GD(:,1)-GD(:,2))/DPOL
     .                  - GDPOL(:) * dDPOLdD(2)/DPOL
        dGDPOLdD(:,3) = 4*GD(:,3)/DPOL
     .                  - GDPOL(:) * dDPOLdD(3)/DPOL
        dGDPOLdD(:,4) = 4*GD(:,4)/DPOL
     .                  - GDPOL(:) * dDPOLdD(4)/DPOL
        dGDDdD(:,1,:) = ( dGDTOTdD(:,:)              ! dGradDensUp/dD(i)
     .                  + dGDPOLdD(:,:) ) / 2
        dGDDdD(:,2,:) = ( dGDTOTdD(:,:)              ! dGradDensDn/dD(i)
     .                  - dGDPOLdD(:,:) ) / 2

        ! Derivatives of grad(Dup) and grad(Ddn) with respect to grad(D(i))
        dGDTOTdGD(1:2) = 1                           ! dGradDensTot/dGradD(i)
        dGDTOTdGD(3:4) = 0
        dGDPOLdGD(1) = +( D(1) - D(2) ) / DPOL       ! dGradDensPol/dGradD(i)
        dGDPOLdGD(2) = -( D(1) - D(2) ) / DPOL
        dGDPOLdGD(3) = 4 * D(3) / DPOL
        dGDPOLdGD(4) = 4 * D(4) / DPOL
        dGDDdGD(1,:) = ( dGDTOTdGD(:)                ! dGradDensUp/dGradD(i)
     .                 + dGDPOLdGD(:) ) / 2
        dGDDdGD(2,:) = ( dGDTOTdGD(:)                ! dGradDensDn/dGradD(i)
     .                 - dGDPOLdGD(:) ) / 2
       endif
        
      else if (nSpin==1 .or. nSpin==2) then ! Normal (collinear) spin
        NS = nSpin
        DD(1:NS) = max( D(1:NS), 0.0_dp ) ! ag: Avoid negative densities
        GDD(1:3,1:NS) = GD(1:3,1:NS)
      else
        call die('ggaxc: ERROR: invalid value of nSpin')
      end if ! (nSpin==4)

      ! Select functional to find energy density and its derivatives
      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NS, DD, GDD,                  ! JMS
     .              EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSE IF (AUTHOR.EQ.'RPBE'.OR.AUTHOR.EQ.'rpbe') THEN
        CALL RPBEXC( IREL, NS, DD, GDD,                 ! MVFS
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSE IF (AUTHOR.EQ.'WC'.OR.AUTHOR.EQ.'wc') THEN
        CALL WCXC( IREL, NS, DD, GDD,                   ! MVFS
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSE IF (AUTHOR.EQ.'REVPBE'.OR.AUTHOR.EQ.'revpbe'
     .                           .OR.AUTHOR.EQ.'revPBE') THEN
        CALL REVPBEXC( IREL, NS, DD, GDD,               ! EA
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSE IF (AUTHOR.EQ.'BLYP'.OR.AUTHOR.EQ.'blyp'.OR.
     .         AUTHOR.EQ.'LYP'.OR.AUTHOR.EQ.'lyp') THEN
        CALL BLYPXC( NS, DD, GDD,                       ! AG
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PW91' .OR. AUTHOR.EQ.'pw91') THEN
        CALL PW91XC( IREL, NS, DD, GDD,                 ! AG
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBESOL' .OR. AUTHOR.EQ.'pbesol' .OR.
     .        AUTHOR.EQ.'PBEsol') THEN
        CALL PBESOLXC( IREL, NS, DD, GDD,
     .                 EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'AM05' .OR. AUTHOR.EQ.'am05') THEN
        CALL AM05XC( IREL, NS, DD, GDD,
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBEJsJrLO' .OR.
     .        AUTHOR.EQ.'pbejsjrlo' .OR.
     .        AUTHOR.EQ.'PBEJSJRLO') THEN
        CALL PBEJsJrLOxc( IREL, NS, DD, GDD,
     .                  EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBEJsJrHEG' .OR.
     .        AUTHOR.EQ.'pbejsjrheg' .OR.
     .        AUTHOR.EQ.'PBEJSJRHEG') THEN
        CALL PBEJsJrHEGxc( IREL, NS, DD, GDD,
     .                  EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBEGcGxLO' .OR.
     .        AUTHOR.EQ.'pbegcgxlo' .OR.
     .        AUTHOR.EQ.'PBEGCGXLO') THEN
        CALL PBEGcGxLOxc( IREL, NS, DD, GDD,
     .                  EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBEGcGxHEG' .OR.
     .        AUTHOR.EQ.'pbegcgxheg' .OR.
     .        AUTHOR.EQ.'PBEGCGXHEG') THEN
        CALL PBEGcGxHEGxc( IREL, NS, DD, GDD,
     .                  EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PW86' .OR. AUTHOR.EQ.'pw86') THEN
        CALL PW86X( IREL, NS, DD, GDD,
     .              EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PW86R' .OR. AUTHOR.EQ.'pw86r') THEN
        CALL PW86RX( IREL, NS, DD, GDD,
     .               EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'B88' .OR. AUTHOR.EQ.'b88') THEN
        CALL B88X( IREL, NS, DD, GDD,
     .             EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'B88KBM' .OR. AUTHOR.EQ.'b88kbm') THEN
        CALL B88KBMX( IREL, NS, DD, GDD,
     .             EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'C09' .OR. AUTHOR.EQ.'c09') THEN
        CALL C09X( IREL, NS, DD, GDD,
     .             EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'BH' .OR. AUTHOR.EQ.'bh') THEN
        CALL BHX( IREL, NS, DD, GDD,
     .            EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD ) 

#ifdef HAVE_LIBINT
      ELSEIF (AUTHOR.EQ.'HSE06' .OR. AUTHOR.EQ.'hse06') THEN
        CALL HSEXC( IREL, NS, DD, GDD,
     .              EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )

      ELSEIF (AUTHOR.EQ.'PBE0' .OR. AUTHOR.EQ.'pbe0') THEN
        CALL PBE0XC( IREL, NS, DD, GDD,
     .              EPSX, EPSC, dEXdDD, dECdDD, dEXdGDD, dECdGDD )
#endif

      ELSE
        call die('GGAXC: Unknown author ' // trim(AUTHOR))
      ENDIF

      ! Find dE/dD(i) and dE/dGradD(i). Note convention:
      ! DEDD(1)=dE/dD11, DEDD(2)=dE/dD22, DEDD(3)=Re(dE/dD12)=Re(dE/dD21),
      ! DEDD(4)=Im(dE/dD12)=-Im(dE/D21)
      if (nSpin==4) then  ! Non colinear spin
        ! dE/dD(i) = dE/dDup * dDup/dD(i) + dE/dDdn * dDdn/dD(i)
        !          + dE/dGDup * dGDup/dD(i) + dE/dGDdn * dGDdn/dD(i)
        ! dE/dGradD(i) = dE/dGDup * dGDup/dGD(i) + dE/dGDdn * dGDdn/dGD(i)

        if ( old_scheme ) then
         VPOL  = (dExdDD(1)-dExdDD(2)) * ct
         dExdD(1) = 0.5D0 * ( dExdDD(1) + dExdDD(2) + VPOL )
         dExdD(2) = 0.5D0 * ( dExdDD(1) + dExdDD(2) - VPOL )
         dExdD(3) = 0.5d0 * (dExdDD(1)-dExdDD(2)) * st * cp
         dExdD(4) =-0.5d0 * (dExdDD(1)-dExdDD(2)) * st * sp
         VPOL  = (dEcdDD(1)-dEcdDD(2)) * ct
         dEcdD(1) = 0.5D0 * ( dEcdDD(1) + dEcdDD(2) + VPOL )
         dEcdD(2) = 0.5D0 * ( dEcdDD(1) + dEcdDD(2) - VPOL )
         dEcdD(3) = 0.5d0 * (dEcdDD(1)-dEcdDD(2)) * st * cp
         dEcdD(4) =-0.5d0 * (dEcdDD(1)-dEcdDD(2)) * st * sp
C Gradient terms
         DO IX = 1,3
          dExdGD(IX,1) = dExdGDD(IX,1)*C2**2 + dExdGDD(IX,2)*S2**2
          dExdGD(IX,2) = dExdGDD(IX,1)*S2**2 + dExdGDD(IX,2)*C2**2
          dExdGD(IX,3) = 0.5D0*(dExdGDD(IX,1) - dExdGDD(IX,2))*ST*CP
          dExdGD(IX,4) =-0.5D0*(dExdGDD(IX,1) - dExdGDD(IX,2))*ST*SP
          dEcdGD(IX,1) = dEcdGDD(IX,1)*C2**2 + dEcdGDD(IX,2)*S2**2
          dEcdGD(IX,2) = dEcdGDD(IX,1)*S2**2 + dEcdGDD(IX,2)*C2**2
          dEcdGD(IX,3) = 0.5D0*(dEcdGDD(IX,1) - dEcdGDD(IX,2))*ST*CP
          dEcdGD(IX,4) =-0.5D0*(dEcdGDD(IX,1) - dEcdGDD(IX,2))*ST*SP
         enddo

        else

        do is = 1,4
          dEXdD(is) = sum( dEXdDD(:) * dDDdD(:,is) )
     .              + sum( dEXdGDD(:,:) * dGDDdD(:,:,is) )
          dECdD(is) = sum( dECdDD(:) * dDDdD(:,is) )
     .              + sum( dECdGDD(:,:) * dGDDdD(:,:,is) )
          do ix = 1,3
            dEXdGD(ix,is) = sum( dEXdGDD(ix,:) * dGDDdGD(:,is) )
            dECdGD(ix,is) = sum( dECdGDD(ix,:) * dGDDdGD(:,is) )
          end do
        end do
        ! Divide by two the non-diagonal derivatives. This is necessary
        ! because DEDD(3:4) intend to be Re(dE/dD12)=Re(dE/dD21) and 
        ! Im(dE/dD12)=-Im(dE/dD21), respectively. However, both D12 and D21
        ! depend on D(3) and D(4), and we have derived directly dE/dD(3) and
        ! dE/dD(4). Although less trivially, the same applies to dE/dGD(3:4).
        dEXdD(3:4) = dEXdD(3:4) / 2
        dECdD(3:4) = dECdD(3:4) / 2
        dEXdGD(:,3:4) = dEXdGD(:,3:4) / 2
        dECdGD(:,3:4) = dECdGD(:,3:4) / 2

       endif
      else   ! Collinear spin => just copy derivatives to output arrays
        dEXdD(1:nSpin) = dEXdDD(1:nSpin)
        dECdD(1:nSpin) = dECdDD(1:nSpin)
        dEXdGD(:,1:nSpin) = dEXdGDD(:,1:nSpin)
        dECdGD(:,1:nSpin) = dECdGDD(:,1:nSpin)
      end if ! (nSpin==4)

      END SUBROUTINE GGAXC



      SUBROUTINE PBEformXC( beta, mu, kappa, iRel, nSpin, Dens, GDens,
     .                      EX, EC, dEXdD, dECdD, dEXdGD, dECdGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation
C functional form, but with variable values for parameters beta, mu, and
C kappa. Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C Written by L.C.Balbas and J.M.Soler. December 1996. 
C ******** INPUT ******************************************************
C REAL*8  BETA           : Parameter beta of the PBE functional
C REAL*8  MU             : Parameter mu of the PBE functional
C REAL*8  KAPPA          : Parameter kappa of the PBE functional
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit none

C Input
      real(dp),intent(in) :: beta           ! Parameter of PBE functional
      real(dp),intent(in) :: mu             ! Parameter of PBE functional
      real(dp),intent(in) :: kappa          ! Parameter of PBE functional
      integer, intent(in) :: iRel           ! Relativistic exchange? 0=no, 1=yes
      integer, intent(in) :: nSpin          ! Number of spin components
      real(dp),intent(in) :: Dens(nSpin)    ! Local electron (spin) density
      real(dp),intent(in) :: GDens(3,nSpin) ! Gradient of electron density

C Output
      real(dp),intent(out):: EX             ! Exchange energy per electron
      real(dp),intent(out):: EC             ! Correlation energy per electron
      real(dp),intent(out):: dEXdD(nSpin)   ! dEx/dDens, Ex=Int(dens*epsX)
      real(dp),intent(out):: dECdD(nSpin)   ! dEc/dDens
      real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)
      real(dp),intent(out):: dECdGD(3,nSpin) ! dEc/dGrad(Dens)

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  A, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  ECUNIF, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KF, KFS, KS, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      GAMMA = (1 - LOG(TWO)) / PI**2

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + MU * S**2 / KAPPA
        F = 1 + KAPPA - KAPPA / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DD = 2 * (F1-1) * DSDD / S
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = 2 * MU * S * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE PBEformXC



      SUBROUTINE PBEXC( IREL, nspin, Dens, GDens,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
C Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C Modified to call PBEformXC by J.M.Soler. December 2009.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.066725_dp       ! From grad. exp. for correl. rs->0
      MU = BETA * PI**2 / 3    ! From Jell. response for x+c
      KAPPA = 0.804_dp         ! From general Lieb-Oxford bound

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBEXC



      SUBROUTINE REVPBEXC( IREL, nspin, Dens, GDens,
     .                     EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements revPBE: revised Perdew-Burke-Ernzerhof GGA.
C Ref: Y. Zhang & W. Yang, Phys. Rev. Lett. 80, 890 (1998).
C Same interface as PBEXC.
C revPBE parameters introduced by E. Artacho in January 2006
C Modified to call PBEformXC by J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.066725_dp       ! From grad. exp. for correl. rs->0
      MU = BETA * PI**2 / 3    ! From Jell. response for x+c
      KAPPA = 1.245_dp         ! From fit of molecular energies

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE REVPBEXC



      SUBROUTINE PBESOLXC( IREL, NSPIN, DENS, GDENS,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
C with the revised parameters for solids (PBEsol).
C Ref: J.P.Perdew et al, PRL 100, 136406 (2008)
C Same interface as PBEXC.
C Modified by J.D. Gale for PBEsol. May 2009.
C Modified to call PBEformXC by J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.046_dp          ! From fit of Jell. surf. energies
      MU = 10.0_dp / 81.0_dp   ! From grad. exp. for exchange.
      KAPPA = 0.804_dp         ! From general Lieb-Oxford bound

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBESOLXC



      SUBROUTINE PBEJsJrLOxc( IREL, NSPIN, DENS, GDENS,
     .                      EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation
C functional form, with revised parameters of Capelle et al:
C   Js refers to Jellium surface energies, that fix parameter beta
C   Jr refers to Jellium response, that fixes parameter mu
C   LO refers to Lieb-Oxford bound, that fixes parameter kappa
C Refs: L.S.Pedroza et al, PRB 79, 201106 (2009)
C       M.M.Odashima et al, J. Chem. Theory Comp. 5, 798 (2009)
C Same interface as PBEXC. J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.046_dp          ! From fit of Jell. surf. energies
      MU = BETA * PI**2 / 3    ! From Jell. response for x+c
      KAPPA = 0.804_dp         ! From general Lieb-Oxford bound

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBEJsJrLOxc



      SUBROUTINE PBEJsJrHEGxc( IREL, NSPIN, DENS, GDENS,
     .                      EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation
C functional form, with revised parameters of Capelle et al:
C   Js refers to Jellium surface energies, that fix parameter beta
C   Jr refers to Jellium response, that fixes parameter mu
C   HGE refers to the Lieb-Oxford bound for the low-density limit of
C       the homogeneous electron gas, that fixes parameter kappa
C Refs: L.S.Pedroza et al, PRB 79, 201106 (2009)
C       M.M.Odashima et al, J. Chem. Theory Comp. 5, 798 (2009)
C Same interface as PBEXC. J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.046_dp          ! From fit of Jell. surf. energies
      MU = BETA * PI**2 / 3    ! From Jell. response for x+c
      KAPPA = 0.552_dp         ! From Lieb-Oxford bound for HEG

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBEJsJrHEGxc



      SUBROUTINE PBEGcGxLOxc( IREL, NSPIN, DENS, GDENS,
     .                      EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation
C functional form, with revised parameters of Capelle et al:
C   Gc refers to gradient exp. for correl., that fixes parameter beta
C   Gx refers to grad. expansion for exchange, that fixes parameter mu
C   LO refers to Lieb-Oxford bound, that fixes parameter kappa
C Refs: L.S.Pedroza et al, PRB 79, 201106 (2009)
C       M.M.Odashima et al, J. Chem. Theory Comp. 5, 798 (2009)
C Same interface as PBEXC. J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.066725_dp       ! From grad. exp. for correl. rs->0
      MU = 10._dp / 81._dp     ! From grad. exp. for exchange.
      KAPPA = 0.804_dp         ! From general Lieb-Oxford bound

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBEGcGxLOxc



      SUBROUTINE PBEGcGxHEGxc( IREL, NSPIN, DENS, GDENS,
     .                      EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation
C functional form, with revised parameters of Capelle et al:
C   Gc refers to gradient exp. for correl., that fixes parameter beta
C   Gx refers to grad. expansion for exchange, that fixes parameter mu
C   HGE refers to the Lieb-Oxford bound for the low-density limit of
C       the homogeneous electron gas, that fixes parameter kappa
C Refs: L.S.Pedroza et al, PRB 79, 201106 (2009)
C       M.M.Odashima et al, J. Chem. Theory Comp. 5, 798 (2009)
C Same interface as PBEXC. J.M.Soler. December 2009.
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      real(dp):: BETA, KAPPA, MU, PI

C Fix values for PBE functional parameters
      PI = 4 * ATAN(1._dp)
      BETA = 0.066725_dp       ! From grad. exp. for correl. rs->0
      MU = 10._dp / 81._dp     ! From grad. exp. for exchange.
      KAPPA = 0.552_dp         ! From Lieb-Oxford bound for HEG

C Call PBE routine with appropriate values for beta, mu, and kappa.
      CALL PBEformXC( BETA, MU, KAPPA, IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

      END SUBROUTINE PBEGcGxHEGxc



      SUBROUTINE PW91XC( IREL, nspin, Dens, GDens,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Wang91 Generalized-Gradient-Approximation.
C Ref: JCP 100, 1290 (1994)
C Written by J.L. Martins  August 2000
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit          none
      INTEGER           IREL, nspin
      real(dp)          Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KF, KFS, KS, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA
     
      real(dp)          F5, F6, F7, F8, ASINHS
      real(dp)          DF5DD,DF6DD,DF7DD,DF8DD
      real(dp)          DF1DS, DF2DS, DF3DS, DFDS, DF7DGD

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4.0_dp * ATAN(1.0_dp)
      BETA = 15.75592_dp * 0.004235_dp
      GAMMA = BETA**2 / (2.0_dp * 0.09_dp)

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      S = GDMT / (2 * KF * DT)
      T = GDMT / (2 * KS * DT)
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      F5 = 0.002568D0 + 0.023266D0*RS + 7.389D-6*RS**2
      F6 = 1.0D0 + 8.723D0*RS + 0.472D0*RS**2 + 0.07389D0*RS**3
      F7 = EXP(-100.0D0 * S**2 * PHI**4)
      F8 =  15.75592D0*(0.001667212D0 + F5/F6 -0.004235D0 + 
     .          3.0D0*0.001667212D0/7.0D0)
      H = GAMMA * PHI**3 * LOG( 1 + F4 ) + F8 * T**2 * F7
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - THD * RS / DT
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - 1 / DT - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DSDD = -S*DKFDD/KF - S/DT    ! JMS: corrected, May.2014
        DTDD = -T*DKSDD/KS - T/DT
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = - F2 * DF1DD
        DADD = - A * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DF5DD = (0.023266D0 + 2.0D0*7.389D-6*RS)*DRSDD
        DF6DD = (8.723D0 + 2.0D0*0.472D0*RS
     .            + 3.0D0*0.07389D0*RS**2)*DRSDD
        DF7DD = -200.0D0 * S * PHI**4 * DSDD * F7
     .         -100.0D0 * S**2 * 4.0D0* PHI**3 * DPDD * F7
        DF8DD = 15.75592D0 * DF5DD/F6 - 15.75592D0*F5*DF6DD / F6**2
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DHDD = DHDD + DF8DD * T**2 * F7
        DHDD = DHDD + F8 * 2*T*DTDD *F7
        DHDD = DHDD + F8 * T**2 * DF7DD
        
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD
        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DSDGD = (S / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DF7DGD = -200.0D0 * S * PHI**4 * DSDGD * F7
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DHDGD = DHDGD + F8 * 2*T*DTDGD *F7 + F8 * T**2 *DF7DGD
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(1) = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(1))**THD
        S = GDMS / (2 * KFS * DS(1))
        F4 = SQRT(1.0D0 + (7.7956D0*S)**2)
        ASINHS = LOG(7.7956D0*S + F4)
        F1 = 1.0D0 + 0.19645D0 * S * ASINHS
        F2 = 0.2743D0 - 0.15084D0*EXP(-100.0D0*S*S)
        F3 = 1.0D0 / (F1 + 0.004D0 * S*S*S*S)
        F = (F1 + F2 * S*S ) * F3
     .       
        CALL EXCHNG( IREL, 1, DS, EXUNIF, VXUNIF )
        FX = FX + DS(1) * EXUNIF * F

        DKFDD = THD * KFS / DS(1)
        DSDD = S * ( -DKFDD/KFS - 1/DS(1) )
        DF1DS = 0.19645D0 * ASINHS +
     .    0.19645D0 * S * 7.7956D0 / F4
        DF2DS = 0.15084D0*200.0D0*S*EXP(-100.0D0*S*S)
        DF3DS = - F3*F3 * (DF1DS + 4.0D0*0.004D0 * S*S*S)
        DFDS =  DF1DS * F3 + DF2DS * S*S * F3 + 2.0D0 * S * F2 * F3
     .            + (F1 + F2 * S*S ) * DF3DS   
        DFXDD(IS) = VXUNIF(1) * F + DS(1) * EXUNIF * DFDS * DSDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(1) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE PW91XC



       subroutine blypxc(nspin,dens,gdens,EX,EC,
     .                   dEXdd,dECdd,dEXdgd,dECdgd) 
c ***************************************************************
c Implements Becke gradient exchange functional (A.D. 
c Becke, Phys. Rev. A 38, 3098 (1988)) and Lee, Yang, Parr
c correlation functional (C. Lee, W. Yang, R.G. Parr, Phys. Rev. B
c 37, 785 (1988)), as modificated by Miehlich,Savin,Stoll and Preuss,
c Chem. Phys. Lett. 157,200 (1989). See also Johnson, Gill and Pople,
c J. Chem. Phys. 98, 5612 (1993). Some errors were detected in this
c last paper, so not all of the expressions correspond exactly to those
c implemented here.
c Written by Maider Machado. July 1998.
c **************** INPUT ******************************************** 
c integer nspin          : Number of spin polarizations (1 or 2)
c real*8  dens(nspin)    : Total electron density (if nspin=1) or
c                           spin electron density (if nspin=2)
c real*8  gdens(3,nspin) : Total or spin density gradient
c ******** OUTPUT *****************************************************
c real*8  ex             : Exchange energy density
c real*8  ec             : Correlation energy density
c real*8  dexdd(nspin)   : Partial derivative
c                           d(DensTot*Ex)/dDens(ispin),
c                           where DensTot = Sum_ispin( Dens(ispin) )
c                          For a constant density, this is the
c                          exchange potential
c real*8  decdd(nspin)   : Partial derivative
c                           d(DensTot*Ec)/dDens(ispin),
c                           where DensTot = Sum_ispin( Dens(ispin) )
c                          For a constant density, this is the
c                          correlation potential
c real*8  dexdgd(3,nspin): Partial derivative
c                           d(DensTot*Ex)/d(GradDens(i,ispin))
c real*8  decdgd(3,nspin): Partial derivative
c                           d(DensTot*Ec)/d(GradDens(i,ispin))
c ********* UNITS ****************************************************
c Lengths in Bohr
c Densities in electrons per Bohr**3
c Energies in Hartrees
c Gradient vectors in cartesian coordinates
c ********************************************************************
 
      implicit none

      integer nspin
      real(dp)   dens(nspin), gdens(3,nspin), EX, EC,
     .           dEXdd(nspin), dECdd(nspin), dEXdgd(3,nspin),
     .           dECdgd(3,nspin)

c Internal variables
      integer is,ix
      real(dp)   pi, beta, thd, tthd, thrhlf, half, fothd,
     .           d(2),gd(3,2),dmin, ash,gdm(2),denmin,dt, 
     .           g(2),x(2),a,b,c,dd,onzthd,gdmin,     
     .           ga, gb, gc,becke,dbecgd(3,2),
     .           dgdx(2), dgdxa, dgdxb, dgdxc,dgdxd,dbecdd(2),
     .           den,omega, domega, delta, ddelta,cf,
     .           gam11, gam12, gam22, LYPa, LYPb1,
     .           LYPb2,dLYP11,dLYP12,dLYP22,LYP,
     .           dd1g11,dd1g12,dd1g22,dd2g12,dd2g11,dd2g22,
     .           dLYPdd(2),dg11dd(3,2),dg22dd(3,2),
     .           dg12dd(3,2),dLYPgd(3,2)
  
c Lower bounds of density and its gradient to avoid divisions by zero
      parameter ( denmin=1.0e-8_dp )
      parameter ( gdmin=1.0e-8_dp )
      parameter ( dmin=1.0e-5_dp ) 

c Fix some numerical parameters 
      parameter ( thd = 1.d0/3.d0, tthd=2.d0/3.d0 )
      parameter ( thrhlf=1.5d0, half=0.5d0,
     .            fothd=4.d0/3.d0, onzthd=11.d0/3.d0)

c Empirical parameter for Becke exchange functional (a.u.)
      parameter(beta= 0.0042d0) 

c Constants for LYP functional (a.u.) 
      parameter(a=0.04918d0, b=0.132d0, c=0.2533d0, dd=0.349d0)

      pi= 4*atan(1.d0)

c Translate density and its gradient to new variables
      if (nspin .eq. 1) then
        d(1) = half * dens(1)
        d(1) = max(denmin,d(1))
        d(2) = d(1)
        dt = max( denmin, dens(1) )
        do ix = 1,3
          gd(ix,1) = half * gdens(ix,1)    
          gd(ix,2) = gd(ix,1)
        enddo 
      else
        d(1) = dens(1)
        d(2) = dens(2)
        do is=1,2
         d(is) = max (denmin,d(is))
        enddo
        dt = max( denmin, dens(1)+dens(2) )  
        do ix = 1,3
          gd(ix,1) = gdens(ix,1)
          gd(ix,2) = gdens(ix,2)
        enddo
      endif

      gdm(1) = sqrt( gd(1,1)**2 + gd(2,1)**2 + gd(3,1)**2 )
      gdm(2) = sqrt( gd(1,2)**2 + gd(2,2)**2 + gd(3,2)**2 )
 
      do is=1,2
      gdm(is)= max(gdm(is),gdmin)
      enddo

c Find Becke exchange energy
       ga = -thrhlf*(3.d0/4.d0/pi)**thd
      do is=1,2
       if(d(is).lt.dmin) then
        g(is)=ga
       else
        x(is) = gdm(is)/d(is)**fothd
        gb = beta*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1)) 
        gc = 1+6*beta*x(is)*ash        
        g(is) = ga-gb/gc
       endif
      enddo

c   Density of energy 
      becke=(g(1)*d(1)**fothd+g(2)*d(2)**fothd)/dt

c Exchange energy derivatives
       do is=1,2
        if(d(is).lt.dmin)then
         dbecdd(is)=0.
         do ix=1,3
          dbecgd(ix,is)=0.
         enddo
        else
        dgdxa=6*beta**2*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1))
        dgdxb=x(is)/sqrt(x(is)**2+1)-ash
        dgdxc=-2*beta*x(is)
        dgdxd=(1+6*beta*x(is)*ash)**2
        dgdx(is)=(dgdxa*dgdxb+dgdxc)/dgdxd
        dbecdd(is)=fothd*d(is)**thd*(g(is)-x(is)*dgdx(is))
        do ix=1,3
         dbecgd(ix,is)=d(is)**(-fothd)*dgdx(is)*gd(ix,is)/x(is)
        enddo 
        endif
       enddo

c  Lee-Yang-Parr correlation energy
      den=1+dd*dt**(-thd)
      omega=dt**(-onzthd)*exp(-c*dt**(-thd))/den
      delta=c*dt**(-thd)+dd*dt**(-thd)/den
      cf=3.*(3*pi**2)**tthd/10.
      gam11=gdm(1)**2
      gam12=gd(1,1)*gd(1,2)+gd(2,1)*gd(2,2)+gd(3,1)*gd(3,2)
      gam22=gdm(2)**2
      LYPa=-4*a*d(1)*d(2)/(den*dt)
      LYPb1=2**onzthd*cf*a*b*omega*d(1)*d(2)
      LYPb2=d(1)**(8./3.)+d(2)**(8./3.)
      dLYP11=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)
     .*d(1)/dt)-d(2)**2)
      dLYP12=-a*b*omega*(d(1)*d(2)/9.*(47.-7.*delta)
     .-fothd*dt**2)
      dLYP22=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)*
     .d(2)/dt)-d(1)**2)

c    Density of energy
      LYP=(LYPa-LYPb1*LYPb2+dLYP11*gam11+dLYP12*gam12
     .+dLYP22*gam22)/dt

c   Correlation energy derivatives
       domega=-thd*dt**(-fothd)*omega*(11.*dt**thd-c-dd/den)
       ddelta=thd*(dd**2*dt**(-5./3.)/den**2-delta/dt)

c   Second derivatives with respect to the density
       dd1g11=domega/omega*dLYP11-a*b*omega*(d(2)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2))

       dd1g12=domega/omega*dLYP12-a*b*omega*(d(2)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)

      dd1g22=domega/omega*dLYP22-a*b*omega*(1./9.*d(2)
     . *(1.-3.*delta-(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3.+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2)-2*d(1))

      dd2g22=domega/omega*dLYP22-a*b*omega*(d(1)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2))
      
      dd2g12=domega/omega*dLYP12-a*b*omega*(d(1)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)
      
      dd2g11=domega/omega*dLYP11-a*b*omega*(1./9.*d(1)
     . *(1.-3.*delta-(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2)-2*d(2))

        dLYPdd(1)=-4*a/den*d(1)*d(2)/dt*
     . (thd*dd*dt**(-fothd)/den
     . +1./d(1)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(2)*(onzthd*
     . d(1)**(8./3.)+d(2)**(8./3.)))+dd1g11*gam11+
     . dd1g12*gam12+dd1g22*gam22

       dLYPdd(2)=-4*a/den*d(1)*d(2)/dt*(thd*dd*dt**(-fothd)/den
     . +1./d(2)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(1)*(onzthd*
     . d(2)**(8./3.)+d(1)**(8./3.)))+dd2g22*gam22+
     . dd2g12*gam12+dd2g11*gam11

c second derivatives with respect to the density gradient
        do is=1,2
          do ix=1,3
           dg11dd(ix,is)=2*gd(ix,is)
           dg22dd(ix,is)=2*gd(ix,is)
          enddo
        enddo
        do ix=1,3
          dLYPgd(ix,1)=dLYP11*dg11dd(ix,1)+dLYP12*gd(ix,2)
          dLYPgd(ix,2)=dLYP22*dg22dd(ix,2)+dLYP12*gd(ix,1)
        enddo

       EX=becke
       EC=LYP
       do is=1,nspin
        dEXdd(is)=dbecdd(is)
        dECdd(is)=dLYPdd(is)
        do ix=1,3
         dEXdgd(ix,is)=dbecgd(ix,is)
         dECdgd(ix,is)=dLYPgd(ix,is)
        enddo
       enddo
       end subroutine blypxc



      SUBROUTINE RPBEXC( IREL, nspin, Dens, GDens,
     .                   EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Hammer's RPBE Generalized-Gradient-Approximation (GGA).
C A revision of PBE (Perdew-Burke-Ernzerhof) 
C Ref: Hammer, Hansen & Norskov, PRB 59, 7413 (1999) and
C J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C
C Written by M.V. Fernandez-Serra. March 2004. On the PBE routine of
C L.C.Balbas and J.M.Soler. December 1996.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit          none
      INTEGER           IREL, nspin
      real(dp)          Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      INTEGER
     .  IS, IX

      real(dp)
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
cea Hammer's RPBE (Hammer, Hansen & Norskov PRB 59 7413 (99)
cea     F1 = DEXP( - MU * S**2 / KAPPA)
cea     F = 1 + KAPPA * (1 - F1)
cea Following is standard PBE
cea     F1 = 1 + MU * S**2 / KAPPA
cea     F = 1 + KAPPA - KAPPA / F1
cea (If revPBE Zhang & Yang, PRL 80,890(1998),change PBE's KAPPA to 1.245)
        F1 = DEXP( - MU * S**2 / KAPPA)
        F = 1 + KAPPA * (1 - F1)
 
c       Note nspin=1 in call to exchng...
 
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

cMVFS   The derivatives of F  also need to be changed for Hammer's RPBE.
cMVFS   DF1DD = 2 * F1 * DSDD  * ( - MU * S / KAPPA)
cMVFS   DF1DGD= 2 * F1 * DSDGD * ( - MU * S / KAPPA)
cMVFS   DFDD  = -1 * KAPPA * DF1DD
cMVFS   DFDGD = -1 * KAPPA * DFDGD

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
c       DF1DD = 2 * (F1-1) * DSDD / S
c       DFDD = KAPPA * DF1DD / F1**2
        DF1DD = 2* F1 * DSDD * ( - MU * S / KAPPA)
        DFDD = -1 * KAPPA * DF1DD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
c         DF1DGD = 2 * MU * S * DSDGD / KAPPA
c         DFDGD = KAPPA * DF1DGD / F1**2
          DF1DGD =2*F1 * DSDGD * ( - MU * S / KAPPA)
          DFDGD = -1 * KAPPA * DF1DGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE RPBEXC



      SUBROUTINE WCXC( IREL, nspin, Dens, GDens,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Wu-Cohen Generalized-Gradient-Approximation.
C Ref: Z. Wu and R. E. Cohen PRB 73, 235116 (2006)
C Written by Marivi Fernandez-Serra, with contributions by
C Julian Gale and Alberto Garcia,
C over the PBEXC subroutine of L.C.Balbas and J.M.Soler.
C September, 2006.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit          none
      INTEGER           IREL, nspin
      real(dp)          Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      INTEGER
     .  IS, IX

      real(dp)
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  XWC, DXWCDS, CWC,
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  TEN81, 
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )
      PARAMETER ( TEN81 = 10.0d0/81.0d0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0
      CWC = 0.0079325D0

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
c
c For PBE: 
c
c       x = MU * S**2
c       dxds = 2*MU*S
c
c Wu-Cohen form:
c
        XWC= TEN81 * s**2 + (MU- TEN81) * 
     .       S**2 * exp(-S**2) + log(1+ CWC * S**4)
        DXWCDS = 2 * TEN81 * S + (MU - TEN81) * exp(-S**2) *
     .           2*S * (1 - S*S) + 4 * CWC * S**3 / (1 + CWC * S**4) 
c-------------------

        F1 = 1 +  XWC / KAPPA
        F = 1 + KAPPA - KAPPA / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DD = DXWCDS * DSDD / KAPPA 
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = DXWCDS * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE WCXC



      SUBROUTINE AM05XC( IREL, nspin, Dens, GDens,
     .                   EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements the Armiento Mattsson AM05 GGA.
C Ref: R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005)
C Written by L.C.Balbas and J.M.Soler originally for PBE. December 1996. 
C Modified by J.D. Gale for AM05. May 2009.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C am05wbs
C ********************************************************************

      use precision, only : dp
      use am05,      only : am05wbs

      implicit          none
      integer           irel, nspin
      real(dp)          Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      integer
     .  is, ix

      real(dp)
     .  D(2), DENMIN, DFXDD(2), DFCDD(2), DFCDGD(3,2), 
     .  DFXDGD(3,2), DFXDG(2), DFCDG(2),
     .  DS(2), DT, EC, EX, FX, FC,
     .  GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3)

C Lower bounds of density and its gradient to avoid divisions by zero
      parameter ( DENMIN = 1.D-12 )
      parameter ( GDMIN  = 1.D-12 )

C Translate density and its gradient to new variables
      if (nspin .eq. 1) then
        D(1) = 0.5_dp*Dens(1)
        D(2) = D(1)
        DT = max( DENMIN, Dens(1) )
        do ix = 1,3
          GD(ix,1) = 0.5_dp*GDens(ix,1)
          GD(ix,2) = GD(ix,1)
          GDT(ix) = GDens(ix,1)
        enddo
      else
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = max( DENMIN, Dens(1)+Dens(2) )
        do ix = 1,3
          GD(ix,1) = GDens(ix,1)
          GD(ix,2) = GDens(ix,2)
          GDT(ix) = GDens(ix,1) + GDens(ix,2)
        enddo
      endif
      GDM(1) = sqrt( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = sqrt( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = sqrt( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = max( GDMIN, GDMT )

      D(1) = max(D(1),denmin)
      D(2) = max(D(2),denmin)

C Call AM05 subroutine
      call am05wbs(D(1), D(2), GDM(1), GDM(2), FX, FC,
     .             DFXDD(1), DFXDD(2), DFCDD(1), DFCDD(2), 
     .             DFXDG(1), DFXDG(2), DFCDG(1), DFCDG(2))

C Convert gradient terms into vectors
      do is = 1,nspin
        do ix = 1,3
          DFXDGD(ix,is) = DFXDG(is)*GD(ix,is)
          DFCDGD(ix,is) = DFCDG(is)*GD(ix,is)
        enddo
      enddo

C Convert FX to form required by SIESTA - note factor of 1/2
C is already applied in am05 code.
      FX = FX / DT

C Set output arguments
      EX = FX
      EC = FC
      do is = 1,nspin
        DEXDD(is) = DFXDD(is)
        DECDD(is) = DFCDD(is)
        do ix = 1,3
          DEXDGD(ix,is) = DFXDGD(ix,is)
          DECDGD(ix,is) = DFCDGD(ix,is)
        enddo
      enddo

      END SUBROUTINE AM05XC



      SUBROUTINE PW86formX( a, b, c, iRel, nSpin, Dens, GDens,
     .                      EX, dEXdD, dEXdGD )

C *********************************************************************
C Implements Perdew-Wang-86 Generalized-Gradient-Approximation exchange-only
C functional form, eps_x=eps_x_LDA*(1+15*a*s**2+b*s**4+c*s**6)**(1/15),
C but with variable values for parameters a, b, and c
C Refs: J.P.Perdew & Y.Wang, PRB 33, 8800 (1986)
C       E.D.Murray, K.Lee & D.C.Langreth, JCTC 5, 2754 (2009)
C Written by J.M.Soler. March 2010.
C ******** INPUT ******************************************************
C REAL*8  a, b, c        : Parameter a, b, and c of the PW86 functional
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG
C ********************************************************************

      implicit none

C Input
      real(dp),intent(in) :: a              ! Parameter of PW86 functional
      real(dp),intent(in) :: b              ! Parameter of PW86 functional
      real(dp),intent(in) :: c              ! Parameter of PW86 functional
      integer, intent(in) :: iRel           ! Relativistic exchange? 0=no, 1=yes
      integer, intent(in) :: nSpin          ! Number of spin components
      real(dp),intent(in) :: Dens(nSpin)    ! Local electron (spin) density
      real(dp),intent(in) :: GDens(3,nSpin) ! Gradient of electron density

C Output
      real(dp),intent(out):: EX             ! Exchange energy per electron
      real(dp),intent(out):: dEXdD(nSpin)   ! dEx/dDens, Ex=Int(dens*epsX)
      real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  D(2), DENMIN, DF1DS, DFDD, DFDGD, DFDS, DFXDD(2), DFXDGD(3,2), 
     .  DKFDD, DS(2), DSDD, DSDGD, DT, ECUNIF, EXUNIF, F, F1, FX,
     .  GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  KFS, PI, S, VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = Dens(1) / 2
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = GDens(IX,1) / 2
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**(1._dp/3)
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + 15*a*S**2 + b*S**4 + c*S**6
        F = F1**(1._dp/15)
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = KFS / DS(IS) / 3
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DS = 30*a*S + 4*b*S**3 + 6*c*S**5
        DFDS = F/F1/15 * DF1DS
        DFDD = DFDS * DSDD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = FX / DT / 2

C Set output arguments
      EX = FX
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE PW86formX



      SUBROUTINE PW86X( IREL, nspin, Dens, GDens,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Wang-86 Generalized-Gradient-Approximation 
C exchange-only functional. Correlation energy returns as zero.
C Ref: J.P.Perdew & Y.Wang, PRB 33, 8800 (1986)
C Written by J.M.Soler. March 2010.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C PW86formX
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal parameters of the PW86 exchange functional
      real(dp),parameter:: a = 0.0864_dp
      real(dp),parameter:: b = 14.0_dp
      real(dp),parameter:: c = 0.2_dp

C Call PW86 routine with appropriate values for a, b, and c.
      CALL PW86formX( a, b, c, IREL, nspin, Dens, GDens,
     .                EX, DEXDD, DEXDGD )

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE PW86X



      SUBROUTINE PW86RX( IREL, nspin, Dens, GDens,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Wang-86 Generalized-Gradient-Approximation 
C exchange-only functional with the refitted parameters of
C Murray, Lee, and Langreth. Correlation energy returns as zero.
C Refs: J.P.Perdew & Y.Wang, PRB 33, 8800 (1986)
C       E.D.Murray, K.Lee & D.C.Langreth, JCTC 5, 2754 (2009)
C Written by J.M.Soler. March 2010.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C PW86formX
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal parameters of the refitted PW86 exchange functional
      real(dp),parameter:: a = 0.1234_dp
      real(dp),parameter:: b = 17.33_dp
      real(dp),parameter:: c = 0.163_dp

C Call PW86 routine with appropriate values for a, b, and c.
      CALL PW86formX( a, b, c, IREL, nspin, Dens, GDens,
     .                EX, DEXDD, DEXDGD )

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE PW86RX


      SUBROUTINE BHX( IREL, nspin, Dens, GDens,
     .                EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements the combination by Berland and Hyldgaard of Perdew-Wang-86 
C Generalized-Gradient-Approximation exchange-only functional with the 
C refitted parameters of Murray, Lee, and Langreth and Langreth-Vosko 
C screened exchange. Correlation energy returns as zero.
C Refs: J.P.Perdew & Y.Wang, PRB 33, 8800 (1986)
C       E.D.Murray, K.Lee & D.C.Langreth, JCTC 5, 2754 (2009)
C       K.Berland & P.Hyldgaard, PRB 89, 035412 (2014)
C Written by Michelle Fritz Feb. 2014.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN          : Number of spin polarizations (1 or 2)
C REAL*8  DENS(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  D(2), DENMIN, DF1DS, DFDD, DFDGD, DFDS, DFLVDS, DFPW86RDS, 
     .  DFXDD(2), DFXDGD(3,2), DKFDD, DS(2), DSDD, DSDGD, DT, ECUNIF, 
     .  EXUNIF, FPW86R, FLV, F, F1, FX, GD(3,2), GDM(2), GDMIN, GDMS, 
     .  GDMT, GDS, GDT(3), KFS, PI, S, VXUNIF(2), ZETA, MLV

C Internal parameters of the refitted PW86 exchange functional
      real(dp),parameter:: a = 0.1234_dp
      real(dp),parameter:: b = 17.33_dp
      real(dp),parameter:: c = 0.163_dp
      real(dp),parameter:: zab = -0.8491_dp
      real(dp),parameter:: alpha = 0.02178_dp
      real(dp),parameter:: beta = 1.15_dp

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)

C Translate density and its gradient to new variables
      IF (NSPIN .EQ. 1) THEN
        D(1) = DENS(1) / 2
        D(2) = D(1)
        DT = MAX( DENMIN, DENS(1) )
        DO 10 IX = 1,3
          GD(IX,1) = GDENS(IX,1) / 2
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDENS(IX,1)
   10   CONTINUE
      ELSE
        D(1) = DENS(1)
        D(2) = DENS(2)
        DT = MAX( DENMIN, DENS(1)+DENS(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDENS(IX,1)
          GD(IX,2) = GDENS(IX,2)
          GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**(1._dp/3)
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + 15*a*S**2 + b*S**4 + c*S**6
        FPW86R = F1**(1._dp/15)
        MLV = -zab/9._dp
        FLV = 1 + MLV*S**2
        F = (1 / (1 + alpha*S**6)) * FLV
        F = F + ((alpha*S**6) / (beta + alpha*S**6)) * FPW86R
        
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = KFS / DS(IS) / 3
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DS = 30*a*S + 4*b*S**3 + 6*c*S**5
        DFLVDS = 2*MLV*S
        DFPW86RDS = FPW86R/F1/15 * DF1DS
        DFDS = -((6*alpha*S**5) / (1 + alpha*S**6)**2) * FLV
        DFDS = DFDS + (1 / (1 + alpha*S**6)) * DFLVDS
        DFDS = DFDS + ((6*alpha*S**5) / (beta + alpha*S**6)) * FPW86R
        DFDS =DFDS-((6*alpha**2*S**11)/(beta + alpha*S**6)**2)*FPW86R
        DFDS = DFDS + ((alpha*S**6) / (beta + alpha*S**6)) * DFPW86RDS
        DFDD = DFDS * DSDD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = FX / DT / 2

      EX = FX
      DO 90 IS = 1,NSPIN
        DEXDD(IS) = DFXDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE BHX


      SUBROUTINE B88formX( beta, mu, c, iRel, nSpin, Dens, GDens,
     .                     EX, dEXdD, dEXdGD )

C *********************************************************************
C Implements Becke-88 Generalized-Gradient-Approximation exchange-only
C functional form, eps_x=eps_x_LDA*(1+mu*s**2/(1+beta*s*asinh(c*s))),
C but with variable values for parameters beta, mu, and c
C Refs: A.D.Becke, PRA 38, 3098 (1988)
C       J.Klimes, D.R.Bowler, and A.Michaelides, JPCM 22, 022201 (2009)
C Written by J.M.Soler. April 2010.
C ******** INPUT ******************************************************
C REAL*8  beta, mu, c    : Parameter beta mu, and c of the B88 functional,
C                          as formulated by Klimes et al
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG
C ********************************************************************

      implicit none

C Input
      real(dp),intent(in) :: beta           ! Parameter of B88 functional
      real(dp),intent(in) :: mu             ! Parameter of B88 functional
      real(dp),intent(in) :: c              ! Parameter of B88 functional
      integer, intent(in) :: iRel           ! Relativistic exchange? 0=no, 1=yes
      integer, intent(in) :: nSpin          ! Number of spin components
      real(dp),intent(in) :: Dens(nSpin)    ! Local electron (spin) density
      real(dp),intent(in) :: GDens(3,nSpin) ! Gradient of electron density

C Output
      real(dp),intent(out):: EX             ! Exchange energy per electron
      real(dp),intent(out):: dEXdD(nSpin)   ! dEx/dDens, Ex=Int(dens*epsX)
      real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  ASINHCS, CS, D(2), DASINHDS, 
     .  DENMIN, DF1DS, DFDD, DFDGD, DFDS, DFXDD(2), DFXDGD(3,2), 
     .  DKFDD, DS(2), DSDD, DSDGD, DT, ECUNIF, EXUNIF, F, F1, FX,
     .  GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  KFS, PI, S, VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = Dens(1) / 2
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = GDens(IX,1) / 2
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**(1._dp/3)
        S = GDMS / (2 * KFS * DS(IS))
        CS = C * S
        ! Next line is ASINH(CS). See Numerical Recipes
        ASINHCS = log(CS+sqrt(CS**2+1))
        F1 = 1 + beta * S * ASINHCS
        F = 1 + mu * S**2 / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = KFS / DS(IS) / 3
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DASINHDS = C * (1 + CS/sqrt(CS**2+1)) / (CS+sqrt(CS**2+1))
        DF1DS = beta * ASINHCS + beta * S * DASINHDS
        DFDS = 2*mu*S/F1 - mu*S**2/F1**2 * DF1DS
        DFDD = DFDS * DSDD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = FX / DT / 2

C Set output arguments
      EX = FX
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END SUBROUTINE B88formX



      SUBROUTINE B88X( IREL, nspin, Dens, GDens,
     .                 EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Becke-88 Generalized-Gradient-Approximation 
C exchange-only functional. Correlation energy returns as zero.
C Refs: A.D.Becke, PRA 38, 3098 (1988)
C       J.Klimes, D.R.Bowler, and A.Michaelides, JPCM 22, 022201 (2009)
C Written by J.M.Soler. April 2010.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C PW86formX
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables and arrays
      real(dp):: beta, c, mu, pi

C Internal parameters of the B88 exchange functional, written as
C by Klimes et al
      pi = acos(-1._dp)
      c = 2*(6*pi**2)**(1._dp/3)
      mu = 0.2743_dp
      beta = 9 * mu * (6/pi)**(1._dp/3) / (2*c)

C Call PW86 routine with appropriate values for a, b, and c.
      CALL B88formX( beta, mu, c, IREL, nspin, Dens, GDens,
     .               EX, DEXDD, DEXDGD )

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE B88X



      SUBROUTINE B88KBMX( IREL, nspin, Dens, GDens,
     .                    EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Becke-88 GGA exchange functional, reparametrized by Klimes
C et al for its combined use with vdW-DF nonlocal correlation.
C Correlation energy returns as zero.
C Refs: A.D.Becke, PRA 38, 3098 (1988)
C       J.Klimes, D.R.Bowler, and A.Michaelides, JPCM 22, 022201 (2009)
C Written by J.M.Soler. April 2010.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C PW86formX
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables and arrays
      real(dp):: beta, c, mu, pi

C Klimes et al parameters for the B88 exchange functional
      pi = acos(-1._dp)
      c = 2*(6*pi**2)**(1._dp/3)
      mu = 0.22_dp
      beta = mu / 1.2_dp

C Call PW86 routine with appropriate values for a, b, and c.
      CALL B88formX( beta, mu, c, IREL, nspin, Dens, GDens,
     .               EX, DEXDD, DEXDGD )

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE B88KBMX


      SUBROUTINE C09X( IREL, nspin, Dens, GDens,
     .                 EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Cooper-09 exchange-only GGA (for use with vdW-DF1 correlation). 
C Correlation energy returns as zero.
C Ref: V.R.Cooper, PRB 81, 161104(R) (2010)
C Written by J.M.Soler. April 2010.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C none
C ********************************************************************

      IMPLICIT NONE

C Passed arguments
      integer, intent(in) :: IREL, NSPIN
      real(dp),intent(in) :: DENS(NSPIN), GDENS(3,NSPIN)
      real(dp),intent(out):: EX, EC, DECDD(NSPIN), DECDGD(3,NSPIN),
     .                       DEXDD(NSPIN), DEXDGD(3,NSPIN)

C Internal variables and arrays (V.R.Cooper, PRB 81, 161104(R) (2010))
      real(dp),parameter:: alpha = 0.0483_dp
      real(dp),parameter:: kappa = 1.245_dp
      real(dp),parameter:: mu    = 0.0617

C Internal variables
      INTEGER
     .  IS, IX
      real(dp)
     .  D(2), DENMIN, DES2DS, DF1DS, DFDD, DFDGD, DFDS,  
     .  DFXDD(2), DFXDGD(3,2), DKFDD, DS(2), DSDD, DSDGD, DT, 
     .  ECUNIF, ES2, EXUNIF, F, FX,
     .  GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  KFS, PI, S, S2, VXUNIF(2)

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = Dens(1) / 2
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = GDens(IX,1) / 2
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS) = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**(1._dp/3)
        S = GDMS / (2 * KFS * DS(IS))
        S2 = S**2
        
        ! Next lines are the core of the C09 functional
        ES2 = exp(-alpha*S2/2)
        F = 1 + mu*S2*ES2**2 + kappa*(1-ES2)
        DES2DS = -alpha*S*ES2
        DFDS = 2*mu*S*ES2**2 + 2*mu*S2*ES2*DES2DS - kappa*DES2DS

        ! Note nspin=1 in call to exchng...
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = KFS / DS(IS) / 3
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DFDD = DFDS * DSDD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = FX / DT / 2

C Set output arguments
      EX = FX
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

C Set correlation energy and derivatives to zero
      EC = 0
      DECDD(:) = 0
      DECDGD(:,:) = 0

      END SUBROUTINE C09X

      END MODULE m_ggaxc
