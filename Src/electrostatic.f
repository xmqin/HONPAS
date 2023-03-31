! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module electrostatic
!
!     Electrostatic correction energy
!
      use atmparams, only: ntbmax
      use precision, only: dp
      use atm_types, only: nspecies, species, elec_corr, npairs
      use radial
      use atmfuncs, only: floating, zvalfis, psch
      use interpolation, only: polint  ! polynomial interpolation
      use m_radfft
      use sys, only: die
      use m_bessph, only: bessph   ! Spherical Bessel functions
!-----------------------------------------------
      
      public :: elec_corr_setup
!-----------------------------------------------

      CONTAINS  !=================================

      subroutine elec_corr_setup()

      integer is, is2, i
      real(dp) rchloc, rchloc2

      type(rad_func), pointer :: func

      npairs = ((nspecies+1)*nspecies)/2
      allocate(elec_corr(npairs))
      do is=1,npairs
        func=>elec_corr(is)
        call rad_zero(func)
      enddo

      do is=1,nspecies
         rchloc = species(is)%chlocal%cutoff
         do is2=is,1,-1
            rchloc2 = species(is2)%chlocal%cutoff
            i = ((is-1)*is)/2+is2
            func => elec_corr(i)
            if (floating(is) .or. floating(is2)) then
               call rad_zero(func)
            else
               call rad_alloc(func,NTBMAX)
               func%cutoff = rchloc + rchloc2 + 0.2_dp
               func%delta =  func%cutoff / (NTBMAX - 1)
               call ch_overlap(is,is2,func%cutoff,func%f)
               call rad_setup_d2(func,yp1=0.0_dp,ypn=huge(1.0_dp))
            endif
         enddo
      enddo

      end subroutine elec_corr_setup

!======================================================================
!
      SUBROUTINE CH_OVERLAP(IS1,IS2,RMX,CORR)

      integer, intent(in)   :: is1, is2
      real(dp), intent(in )    :: rmx
      real(dp), intent(out)    :: corr(:)

C     Returns a table with the difference between the electrostatic energy 
C     of two spherical charge-densities and two point charges with the 
C     same total charge as a function of the distance between the centers 
C     of these charge densities. 
C     Written by D.Sanchez-Portal. March, 1997.(from routine MATEL, written 
C     by Jose M. Soler)

C     INTEGER IS1,IS2             :  Species indexes.
C     RMX                         :  Maximum range of the correction.
C     CORR(NTBMAX)                :  Electrostatic correction energy.

C     Distances in Bohr. Energy in Rydbergs.

C     Internal precision parameters  ------------------------------------
C     NQ is the number of radial points in reciprocal space.
C     Npoint , 2npoint+1 is the number of points used by POLINT in the 
C     interpolation.
C     Q2CUT is the required planewave cutoff for the expansion of
C     the 'local-pseudopotential atomic charge density'
C     (in Ry if lengths are in Bohr).
C     CHERR is a small number to check the precision of the charge density
C     integration.

      integer nq, npoint
      real(dp)            :: q2cut, cherr
      PARAMETER ( NQ     =  512  )
      PARAMETER ( NPOINT =  4     ) 
      PARAMETER ( Q2CUT  =  2.5e3_dp )
      PARAMETER ( CHERR   =  5.e-2_dp )

      real(dp)
     .     CH(0:NQ,2),VTB(NTBMAX,2),
     .     V(0:NQ,2),
     .     GRCH(3),RX(3),RAUX(2*NPOINT+1)


      REAL(DP) cons, qmax, rmax, delt, c, dlt, z1, z2, ch1, ch2, pi
      REAL(DP) r, vd, vv1, vv2, energ1, energ2
      integer itb, nr, nmin, nmax, nn, iq, ir
      real(dp) zval1, zval2

      REAL(DP) QTMP             

      PI= 4._DP * ATAN(1._DP)       
      CONS= 1.0_dp/(2.0_dp*PI)**1.5_DP
C     
C***  CUT-OFF IN REAL AND RECIPROCAL SPACE**
C     
      QMAX =  SQRT( Q2CUT )
      RMAX = PI * NQ / QMAX
      IF(RMX.GT.RMAX) THEN  
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION',
     .           ' POINTS MUST BE INCREASED'
            write(6,'(a,2f15.6)') 'ch_overlap: rmx,rmax =', rmx, rmax
         call die
      ENDIF 
      DELT=PI/QMAX
      C=4.0_DP*PI*DELT
      DLT=RMX/(NTBMAX-1)

      ZVAL1=ZVALFIS(IS1)
      ZVAL2=ZVALFIS(IS2)

      Z1=0.0_DP
      Z2=0.0_DP

      RX(2)=0.0_DP
      RX(3)=0.0_DP 

      DO IR=0,NQ
         R=IR*DELT
         RX(1)=R
             
         CALL PSCH(IS1,RX,CH1,GRCH)
         CALL PSCH(IS2,RX,CH2,GRCH)

         CH(IR,1)=-CH1
         CH(IR,2)=-CH2

         Z1=Z1-C*CH1*R*R    
         Z2=Z2-C*CH2*R*R

      ENDDO
           
      IF((ABS(Z1-ZVAL1).GT.CHERR).OR.
     .     (ABS(Z2-ZVAL2).GT.CHERR)) THEN 
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION',
     .           ' POINTS MUST BE INCREASED'
            WRITE(6,*) 'CH_OVERLAP: Z1=',Z1,' ZVAL1=',ZVAL1
            WRITE(6,*) 'CH_OVERLAP: Z2=',Z2,' ZVAL2=',ZVAL2
         call die
      ENDIF

      DO IR=0,NQ
         CH(IR,1)=real(ZVAL1,dp)*CH(IR,1)/Z1
         CH(IR,2)=real(ZVAL2,dp)*CH(IR,2)/Z2
      ENDDO 
C
C     REAL SPACE INTEGRATION OF POISSON EQUATION
C          
          
      CALL NUMEROV(NQ,DELT,CH(0,1),V(0,1))
      CALL NUMEROV(NQ,DELT,CH(0,2),V(0,2))
           
      DO ITB=1,NTBMAX
         R=DLT*(ITB-1)
         NR=NINT(R/DELT)
         NMIN=MAX(0,NR-NPOINT)
         NMAX=MIN(NQ,NR+NPOINT)
         NN=NMAX-NMIN+1
         DO IR=1,NN
            RAUX(IR)=DELT*(NMIN+IR-1) 
         ENDDO 
         CALL POLINT(RAUX,V(NMIN,1),NN,R,VV1,VD)
         CALL POLINT(RAUX,V(NMIN,2),NN,R,VV2,VD)
 
         VTB(ITB,1)=VV1
         VTB(ITB,2)=VV2
      ENDDO 
         
C****FOURIER-TRANSFORM OF RADIAL CHARGE DENSITY****
C
      CALL RADFFT( 0, NQ, RMAX, CH(0:NQ,1), CH(0:NQ,1) )
      CALL RADFFT( 0, NQ, RMAX, CH(0:NQ,2), CH(0:NQ,2) )
      CALL RESET_RADFFT( )
C

CNEUTRALIZE CHARGE DENSITY FOR FOURIER-SPACE CALCULATION
C
      DO IQ=0,NQ
         R=IQ*QMAX/NQ
         CH1 = (CH(IQ,1)-ZVAL1*CONS)*CH(IQ,2)
         CH2=  (CH(IQ,2)-ZVAL2*CONS)*CH(IQ,1)
         CH(IQ,1) = CH1
         CH(IQ,2) = CH2
      ENDDO
C
C     THE ELECTROSTATIC ENERGY CORRECTION IS STORED IN 'CORR'
C  
      DO IR=1,NTBMAX

         R=DLT*(IR-1)
         ENERG1=0.0_dp
         ENERG2=0.0_dp


         DO IQ=0,NQ
            QTMP=IQ*QMAX/NQ
            QTMP=QTMP*R 
            ENERG1=ENERG1+BESSPH(0,QTMP)*CH(IQ,1)
            ENERG2=ENERG2+BESSPH(0,QTMP)*CH(IQ,2)
         ENDDO 

         ENERG1=ENERG1*QMAX/NQ
         ENERG2=ENERG2*QMAX/NQ
   
         ENERG2=ENERG2*4.0_DP*(2.0_dp*PI)**2
         ENERG1=ENERG1*4.0_DP*(2.0_dp*PI)**2
              
         ENERG1=-(ENERG1*R)-(ZVAL2*(VTB(IR,1)*R-ZVAL1))
         ENERG2=-(ENERG2*R)-(ZVAL1*(VTB(IR,2)*R-ZVAL2))
  
         CORR(IR)=0.5_DP*(ENERG1+ENERG2)

      ENDDO 

      END subroutine ch_overlap


      SUBROUTINE NUMEROV(NR,DELT,Q,V)
      integer, intent(in)  :: nr
      REAL(DP), intent(in)   :: delt
      REAL(DP), intent(in)   :: q(0:nr)
      REAL(DP), intent(out)  :: v(0:nr)

C   Being Q(r) a spherical charge density in a homogeneus radial mesh
C   with distance DELT between consecutive points, this routine returns
C   the electrostatic potential generated by this charge distribution.
C   Written by D. Sanchez-Portal, March 1997.

C   INTEGER NR      :    Number of radial points.
C   REAL(DP)  DELT    :    Distance between consecutive points.
C   REAL(DP)  Q(0:NR) :    Spherical charge density.
C   REAL(DP)  V(0:NR) :    Electrostatic potential at mesh points.

C   Qtot/r asimptotic behaviour is imposed.


      integer ir
      REAL(DP) pi, fourpi, qtot, r, cons

      PI=4.0_DP*ATAN(1.0_DP)
      FOURPI=4.0_DP*PI

C     NUMEROV ALGORITHM* 
C
      V(0)=0.0_DP
      V(1)=1.0_DP

      DO IR=2,NR
         V(IR)=2.0_DP*V(IR-1)-V(IR-2) - FOURPI*DELT**3*
     .        ( Q(IR)*IR+10.0_DP*Q(IR-1)*(IR-1)+Q(IR-2)*(IR-2) )/12.0_DP
      ENDDO 

C***CALCULATE TOTAL CHARGE***
   
      QTOT=0.0_DP
      DO IR=1,NR
         R=IR*DELT
         QTOT=QTOT+R*R*Q(IR)
      ENDDO
      QTOT=4.0_DP*PI*QTOT*DELT

C** FIXING QTOT/R ASYMPTOTIC BEHAVIOUR*

      CONS=(QTOT-V(NR))/(NR*DELT)
             
      DO IR=1,NR
         R=IR*DELT
         V(IR)=V(IR)/(IR*DELT)+CONS
      ENDDO 
      V(0)=(4.0_DP*V(1)-V(2))/3.0_DP

      END subroutine numerov

      end module electrostatic

