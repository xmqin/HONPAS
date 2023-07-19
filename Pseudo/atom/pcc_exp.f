      SUBROUTINE PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)

c mmga
c mmga   M.M.G. Alemany, January 2000  
c mmga
c
c     constructs the partial core correction. Core charge,
c     first and second derivatives are conserved at the point
c     icore. Interpolation is done with Lagrange formula.
c     The core function is exp(ac+bc*r**2+cc*r**4).

      IMPLICIT none


      integer icore, nr
      double precision ac, bc, cc
      double precision R(NR),CDC(NR)

      integer nn
      PARAMETER (NN = 5)
      double precision DR_K(-NN:NN),DDR_K(-NN:NN)
      double precision CDC_SCA(-NN:NN)

      integer i, in1, in2, in, jn, kn, ln
      double precision f1, f2, dr, cdcp, ddr, ddcdc, cdcpp

      IF ( (ICORE - NN).LT.1 .OR. (ICORE + NN).GT.NR ) THEN
        WRITE(6,2017)
        CALL EXT(830)
      ENDIF

      DO 2000   I = ICORE-NN,ICORE+NN
        CDC_SCA(I-ICORE)= LOG(CDC(I)) - 2.D0 * LOG(R(I))
 2000 CONTINUE

      IN1 =-NN
      IN2 = NN

      DO 2011 IN = IN1,IN2
        IF (IN.EQ.0) THEN
          DR_K(IN) = 0.D0
          DO 2012 JN = IN1,IN2
            IF (JN.NE.0) DR_K(IN) = DR_K(IN) + 1.D0/(0 - JN)
 2012     CONTINUE
        ELSE
          F1 = 1.D0
          F2 = 1.D0
          DO 2014 JN = IN1,IN2
            IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
            IF (JN.NE.IN)               F2 = F2 * (IN - JN)
 2014     CONTINUE
          DR_K(IN) = F1 / F2
        ENDIF
 2011 CONTINUE

      DR = 0.D0
      DO 2016 IN = IN1,IN2
        DR = DR + DR_K(IN) * R(ICORE+IN)
 2016 CONTINUE

      CDCP = 0.D0
      DO 2001 IN = IN1,IN2
        CDCP = CDCP + DR_K(IN) * CDC_SCA(IN)
 2001 CONTINUE
      CDCP = CDCP / DR

      DO 2002 IN = IN1,IN2
        DDR_K(IN)= 0.D0
        IF (IN.EQ.0) THEN
          DO 2003 JN = IN1,IN2
            IF (JN.NE.0) THEN
              DO 2004 KN = JN+1,IN2
                IF (KN.NE.0) DDR_K(IN) = DDR_K(IN) + 
     .                                   2.D0/((0 - JN)*(0 - KN))
 2004         CONTINUE
            ENDIF
 2003     CONTINUE
        ELSE
          F2 = 1.D0
          DO 2005 JN = IN1,IN2
            IF (JN.NE.IN) THEN
              F2 = F2 * (IN - JN)
              DO 2006 KN = JN+1,IN2
                IF (KN.NE.IN) THEN
                  F1 = 2.D0
                  DO 2007 LN = IN1,IN2
                    IF (LN.NE.IN .AND. LN.NE.JN .AND. LN.NE.KN)
     .                  F1 = F1 * (0 - LN)
 2007             CONTINUE
                  DDR_K(IN) = DDR_K(IN) + F1                    
                ENDIF
 2006         CONTINUE
            ENDIF
 2005     CONTINUE
          DDR_K(IN) = DDR_K(IN) / F2
        ENDIF
 2002 CONTINUE

      DDR   = 0.D0
      DDCDC = 0.D0
      DO 2010 IN = IN1,IN2
        DDR   = DDR   + DDR_K(IN) * R(ICORE+IN)
        DDCDC = DDCDC + DDR_K(IN) * CDC_SCA(IN)
 2010 CONTINUE

      CDCPP = (DDCDC - DDR * CDCP) / DR**2



      CC = R(ICORE) * CDCPP - CDCP
      CC = CC / (8.D0 * R(ICORE)**3)

      BC = CDCP - 4.D0 * CC * R(ICORE)**3
      BC = BC / (2.D0 * R(ICORE))

      AC = CDC_SCA(0) - BC*R(ICORE)**2 - CC*R(ICORE)**4


      DO 2009 I = 1,ICORE
          CDC(I)= R(I)*R(I) * 
     .            EXP( (CC*R(I)*R(I) + BC) * R(I)*R(I) + AC ) 
 2009 CONTINUE


 2017 FORMAT(//,' error in pcc_exp - ',/,
     . ' derivatives at r(icore) not calculated')


      RETURN

      END
