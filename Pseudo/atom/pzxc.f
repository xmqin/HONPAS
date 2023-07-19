      SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC )

C *****************************************************************
C  Perdew-Zunger parameterization of Ceperley-Alder exchange and 
C  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input *****************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  EC      : correlation energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C REAL*8  VC(NSP) : (spin-dependent) correlation potential
C **** Units *******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      integer nsp, irel
      real*8 DS(NSP), VX(NSP), VC(NSP)

      real*8 zero, one, pfive, opf, pnn
      PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,PNN=.99D0)
      real*8 pthree, psevf, c0504
      PARAMETER (PTHREE=0.3D0,PSEVF=0.75D0,C0504=0.0504D0) 
      real*8 c0254, c014, c0406
      PARAMETER (C0254=0.0254D0,C014=0.014D0,C0406=0.0406D0)
      real*8 c15p9, c0666, c11p4
      PARAMETER (C15P9=15.9D0,C0666=0.0666D0,C11P4=11.4D0)
      real*8 c045, c7p8, c88, c20p59
      PARAMETER (C045=0.045D0,C7P8=7.8D0,C88=0.88D0,C20P59=20.592D0)
      real*8 c3p52, c0311, c0014
      PARAMETER (C3P52=3.52D0,C0311=0.0311D0,C0014=0.0014D0)
      real*8 c0538, c0096, c096
      PARAMETER (C0538=0.0538D0,C0096=0.0096D0,C096=0.096D0)
      real*8 c0622, c004, c0232
      PARAMETER (C0622=0.0622D0,C004=0.004D0,C0232=0.0232D0)
      real*8 c1686, c1p398, c2611
      PARAMETER (C1686=0.1686D0,C1P398=1.3981D0,C2611=0.2611D0)
      real*8 c2846, c1p053, c3334
      PARAMETER (C2846=0.2846D0,C1P053=1.0529D0,C3334=0.3334D0)

C    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.

      real*8 con1, con2, con3, con4, con5, con6, con7, con8, con9,
     $       con10, con11
      
       PARAMETER (CON1=1.D0/6, CON2=0.008D0/3, CON3=0.3502D0/3) 
       PARAMETER (CON4=0.0504D0/3, CON5=0.0028D0/3, CON6=0.1925D0/3)
       PARAMETER (CON7=0.0206D0/3, CON8=9.7867D0/6, CON9=1.0444D0/3)
       PARAMETER (CON10=7.3703D0/6, CON11=1.3336D0/3)

C      X-alpha parameter:

       real*8 alp
       PARAMETER ( ALP = 2.D0 / 3.D0 )

       real*8 ONEZ
       parameter ( ONEZ = 1.0d0 + 1.d-12)

C      Other variables converted into parameters by J.M.Soler

       real*8 pi, half, trd, ftrd, tftm, a0, crs, cxp, cxf

       PARAMETER ( PI   = 3.14159265358979312D0 )
       PARAMETER ( HALF = 0.5D0 ) 
       PARAMETER ( TRD  = 1.D0 / 3.D0 ) 
       PARAMETER ( FTRD = 4.D0 / 3.D0 )
       PARAMETER ( TFTM = 0.51984209978974638D0 )
       PARAMETER ( A0   = 0.52106176119784808D0 )
       PARAMETER ( CRS  = 0.620350490899400087D0 )
       PARAMETER ( CXP  = - 3.D0 * ALP / (PI*A0) )
       PARAMETER ( CXF  = 1.25992104989487319D0 )

cag    Do not attempt to do anything with very small densities,
cag    as vx can blow up with some compilers (pgi, old versions of g77)
cag
       real*8 threshold
       parameter (threshold = 1.0d-30)
cag
       integer isp
       real*8 d, ex, ec, z, fz, fzp, rs, vxp, exp_var, beta, sb,
     $        alb, vxf, exf, sqrs, te, be, ecp, vcp, ecf, vcf,
     $        rslog

c      Find density and polarization

       IF (NSP .EQ. 2) THEN
         D = DS(1) + DS(2)
         IF (D .LE. THRESHOLD) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
         Z = (DS(1) - DS(2)) / D

C
         FZ = ((ONEZ+Z)**FTRD+(ONEZ-Z)**FTRD-2)/TFTM
         FZP = FTRD*((ONEZ+Z)**TRD-(ONEZ-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. THRESHOLD) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VC(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = CRS / D**TRD

C      Exchange
       VXP = CXP / RS
       EXP_VAR = 0.75D0 * VXP
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP_VAR = EXP_VAR *(ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = CXF * VXP
       EXF = CXF * EXP_VAR

C      Correlation 
       IF (RS .GT. ONE) THEN  
         SQRS=SQRT(RS)
         TE = ONE+CON10*SQRS+CON11*RS
         BE = ONE+C1P053*SQRS+C3334*RS
         ECP = -C2846/BE
         VCP = ECP*TE/BE
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -C1686/BE
         VCF = ECF*TE/BE
       ELSE
         RSLOG=LOG(RS)
         ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
         VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
         ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
         VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
       ENDIF

C      Find up and down potentials
       IF (NSP .EQ. 2) THEN
         EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
         EC    = ECP + FZ*(ECF-ECP)
         VX(1) = VXP + FZ*(VXF-VXP) + (ONEZ-Z)*FZP*(EXF-EXP_VAR)
         VX(2) = VXP + FZ*(VXF-VXP) - (ONEZ+Z)*FZP*(EXF-EXP_VAR)
         VC(1) = VCP + FZ*(VCF-VCP) + (ONEZ-Z)*FZP*(ECF-ECP)
         VC(2) = VCP + FZ*(VCF-VCP) - (ONEZ+Z)*FZP*(ECF-ECP)
       ELSE
         EX    = EXP_VAR
         EC    = ECP
         VX(1) = VXP
         VC(1) = VCP
       ENDIF

C      Change from Rydbergs to Hartrees
       EX = HALF * EX
       EC = HALF * EC
       DO 10 ISP = 1,NSP
         VX(ISP) = HALF * VX(ISP)
         VC(ISP) = HALF * VC(ISP)
   10  CONTINUE
      END









