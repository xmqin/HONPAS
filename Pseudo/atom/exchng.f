      SUBROUTINE EXCHNG( IREL, NSP, DS, EX, VX )

C *****************************************************************
C  Finds local exchange energy density and potential
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input ******************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C **** Units ******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      integer irel, nsp
      real*8 ex
      real*8 DS(NSP), VX(NSP)

      real*8 zero, one, pfive, opf, c014
      PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,C014=0.014D0)

      real*8 pi, trd, ftrd, tftm, a0, alp

      real*8 d, z, fz, fzp, rs, vxp, exp_var, beta, sb, alb, vxf, exf

       PI=4*ATAN(ONE)
       TRD = ONE/3
       FTRD = 4*TRD
       TFTM = 2**FTRD-2
       A0 = (4/(9*PI))**TRD

C      X-alpha parameter:       
       ALP = 2 * TRD

       IF (NSP .EQ. 2) THEN
         D = DS(1) + DS(2)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           RETURN
         ENDIF
         Z = (DS(1) - DS(2)) / D
         FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
         FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = (3 / (4*PI*D) )**TRD
       VXP = -3*ALP/(2*PI*A0*RS)
       EXP_VAR = 3*VXP/4
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP_VAR = EXP_VAR * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = 2**TRD*VXP
       EXF = 2**TRD*EXP_VAR
       IF (NSP .EQ. 2) THEN
         VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP_VAR)
         VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP_VAR)
         EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
       ELSE
         VX(1) = VXP
         EX    = EXP_VAR
       ENDIF
      END
