! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module periodic_table

      use precision
      use sys

      implicit none

      public :: qvlofz, lmxofz, cnfig, symbol, atmass

      private

      CONTAINS


      SUBROUTINE QVLOFZ( Z, QVAL )  
C Returns the atomic LSD ground state configurations, as an array with 
C valence population for each angular momentum L. The valence orbitals 
C for each L are those returned by routine LMXOFZ
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,  intent(in)  :: Z        ! Atomic number
      real(dp), intent(out) :: QVAL(0:) ! Valence charge for each L

      integer, parameter :: LMAX=3, ZMAX=98, NCHNG=15
      integer :: I, ICHNG, L, LMXATM, LMXCHM, LCHNG(NCHNG),
     .           N(0:LMAX), NVAL(0:5), Q(0:LMAX,0:ZMAX), ZCHNG(NCHNG)

      ! Notice: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
      DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/

      ! Notice that s valence charge must be consitent with ZCHNG
!
!     Notice: For Pd, the ground state configuration is set as
!                     4d9-5s1.
!
      DATA (Q(0,I),I= 0, 2) /0,1,2/
      DATA (Q(0,I),I= 3,10) /1,2,                      2,2,2,2,0,0/
      DATA (Q(0,I),I=11,18) /1,2,                      2,2,2,2,0,0/
      DATA (Q(0,I),I=19,36) /1,2, 2,2,2,1,2,2,2,2,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=37,54) /1,2, 2,2,1,1,2,1,1,1,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=55,71) /1,2, 2,                               14*2/
      DATA (Q(0,I),I=72,86) /       2,2,2,2,2,2,1,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=87,98) /1,2, 10*2/

      DATA (Q(1,I),I= 0, 2) / 3*0/                ! H-He
      DATA (Q(1,I),I= 3,10) / 2*0, 1,2,3,4,5,6/   ! Li-Ne
      DATA (Q(1,I),I=11,18) / 2*0, 1,2,3,4,5,6/   ! Na-Ar
      DATA (Q(1,I),I=19,36) /12*0, 1,2,3,4,5,6/   ! K-Kr
      DATA (Q(1,I),I=37,54) /12*0, 1,2,3,4,5,6/   ! Rb-Xe
      DATA (Q(1,I),I=55,86) /26*0, 1,2,3,4,5,6/   ! Cs-Rn
      DATA (Q(1,I),I=87,98) /12*0/                ! Fr-Cf

      DATA (Q(2,I),I= 0,36) /21*0, 1,2,3,5,5,6,7, 8,10,10, 6*0/
      DATA (Q(2,I),I=37,54) / 2*0, 1,2,4,5,5,7,8, 9,10,10, 6*0/
      DATA (Q(2,I),I=55,71) / 2*0, 1,                      6*0,1,6*0,1/
      DATA (Q(2,I),I=72,86) /        2,3,4,5,6,7,9,10,10, 6*0/
      DATA (Q(2,I),I=87,98) / 2*0, 1,2,1,1,3*0,1,2,1/

      DATA (Q(3,I),I= 0,71) /58*0, 2,3,4,5,6,7,7,9,10,11,12,13,14,14/
      DATA (Q(3,I),I=72,98) /18*0, 0,2,3,5,6,7,7,7,9/
 
      IF (Z.GT.ZMAX) call die('QVLOFZ: ERROR: Z out of range')

      ! Find the principal quantum numbers assigned by my own data
      DO L=0,LMAX
         N(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         N(L)=N(L)+1
      END DO

      ! Find the valence principal quantum numbers assigned by GNFIG
      CALL LMXOFZ (Z,LMXCHM,LMXATM)
      CALL CNFIG (Z,NVAL)

      ! Check size of QVAL and initialize it (for L>LMXCHM)
      IF (UBOUND(QVAL,1).LT.LMXCHM) 
     .  call die('QVLOFZ: ERROR: Size of QVAL too small')
      QVAL(:) = 0

      ! Make valence occupation consistent with CNFIG valence assignment
      DO L=0,LMXCHM
        IF (NVAL(L).GT.N(L)) THEN
          QVAL(L)=0
        ELSE
          QVAL(L) = Q(L,Z) + 2*(2*L+1)*(N(L)-NVAL(L))
        ENDIF
      END DO

      END subroutine qvlofz


      SUBROUTINE LMXOFZ( Z, LMXVAL, LMXATM ) 
C Given the atomic number Z, returns the maximum angular mometum L
C which is populated in the atomic ground state configuration:
C For the valence orbitals => LMXVAL. For any orbitals => LMXATM
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,intent(in) :: Z      ! Atomic number
      integer,intent(out):: LMXVAL ! Max. L for valence states
      integer,intent(out):: LMXATM ! Max. L for valence and core states

      integer, parameter :: NCHNG=17
      integer :: ICHNG, LCHNG(NCHNG), ZCHNG(NCHNG)

      !            B,Na,Al, K,Sc,Ga,Rb, Y,In,Cs,La,Ce,Hf,Tl,Fr,Ac,Th
      DATA ZCHNG / 5,11,13,19,21,31,37,39,49,55,57,58,72,81,87,89,90/
      DATA LCHNG / 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 3/

      LMXVAL=0
      LMXATM=0
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         LMXVAL=LCHNG(ICHNG)
         LMXATM=MAX(LMXATM,LMXVAL)
      ENDDO

      END subroutine lmxofz


      SUBROUTINE CNFIG( Z, NVAL ) 
C Returns the valence configuration for atomic ground state, i.e.
C the principal quantum number NVAL of the valence orbilas for each L
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,intent(in) :: Z        ! Atomic number
      integer,intent(out):: NVAL(0:) ! Valence electrons for each L

      integer, parameter :: LMAX=3, NCHNG=15
      integer :: ICHNG, L, LCHNG(NCHNG), ZCHNG(NCHNG)

      ! Originally: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
*     DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/

      ! Changed to: s valence orbital switched for full p occupation
      !           Li,Na,Na, K, K,Ga,Rb,Rb,In,Cs,Cs,Hf,Tl,Fr,Fr
      DATA ZCHNG / 3,11,11,19,19,31,37,37,49,55,55,72,81,87,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
      DO L=0,LMAX
         NVAL(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         NVAL(L)=NVAL(L)+1
      END DO

      END subroutine cnfig


      FUNCTION SYMBOL( Z )

!! ** This function should not be called from within an I/O statement

C Given the atomic number, returns the atomic symbol (e.g. 'Na')
C Written by J. Soler

      character(len=2)    :: SYMBOL  ! Atomic symbol
      integer, intent(in) :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) = 
     .         (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     .           'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     .           'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     .           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .           'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     .           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .           'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     .           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .           'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     .           'Md','No','Lr'/)

      IF (Z.EQ.0 .OR. Z.EQ.-100) THEN
         SYMBOL = 'BS'
      ELSE IF (ABS(Z).LE.NZ) THEN
         SYMBOL = NAME(ABS(Z))
      ELSE IF (Z.GT.200) THEN
         write(SYMBOL,'(a1,i1)') 'S', mod(Z-200,10)
      ELSE IF (Z.LT.-200) THEN
         write(SYMBOL,'(a1,i1)') 'G', mod(abs(Z)-200,10)
      ELSE
         WRITE(6,*) 'SYMBOL: ERROR: No data for Z =', Z
         SYMBOL = ' '
      ENDIF

      END function symbol


      FUNCTION ATMASS( Z )
C Returns the average atomic mass from the atomic number Z.
C Dta taken from VCH periodic table.
C Written by J.M.Soler. April'97.

      real(dp)             :: ATMASS ! Average atomic mass, in amu
      integer, intent(in)  :: Z      ! Atomic number

      integer, PARAMETER  :: NZ=94
      character(len=50) message

      real(dp), parameter ::  AMASS(0:NZ) =                (/0.00_dp,
     .     1.01_dp,  4.00_dp,  6.94_dp,  9.01_dp, 10.81_dp, 12.01_dp,
     .    14.01_dp, 16.00_dp, 19.00_dp, 20.18_dp, 22.99_dp, 24.31_dp,
     .    26.98_dp, 28.09_dp, 30.97_dp, 32.07_dp, 35.45_dp, 39.95_dp,
     .    39.10_dp, 40.08_dp, 44.96_dp, 47.88_dp, 50.94_dp, 52.00_dp,
     .    54.94_dp, 55.85_dp, 58.93_dp, 58.69_dp, 63.55_dp, 65.39_dp,
     .    69.72_dp, 72.61_dp, 74.92_dp, 78.96_dp, 79.90_dp, 83.80_dp, 
     .    85.47_dp, 87.62_dp, 88.91_dp, 91.22_dp, 92.91_dp, 95.94_dp, 
     .    98.91_dp,101.07_dp,102.91_dp,106.42_dp,107.87_dp,112.41_dp,
     .   114.82_dp,118.71_dp,121.75_dp,127.60_dp,126.90_dp,131.29_dp,
     .   132.91_dp,137.33_dp,138.91_dp,140.12_dp,140.91_dp,144.24_dp,
     .   146.92_dp,150.36_dp,151.97_dp,157.25_dp,158.93_dp,162.50_dp,
     .   164.93_dp,167.26_dp,168.93_dp,173.04_dp,174.97_dp,178.49_dp,
     .   180.95_dp,183.85_dp,186.21_dp,190.20_dp,192.22_dp,195.08_dp,
     .   196.97_dp,200.59_dp,204.38_dp,207.20_dp,208.98_dp,208.98_dp,
     .   209.99_dp,222.02_dp,223.02_dp,226.03_dp,227.03_dp,232.04_dp,
     .   231.04_dp,238.03_dp,237.05_dp,244.06_dp/)

      IF (Z.LT.0 .OR. Z.GT.NZ) THEN
         write(message,'(a,i4)') 'ATMASS: ERROR: No data for Z =',Z
         call die(message)
      ELSE
         ATMASS=AMASS(Z)
      ENDIF

      END function atmass

      end module periodic_table







