      SUBROUTINE GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )

C Finds the exchange and correlation energies at a point, and their
C derivatives with respect to density and density gradient, in the
C Generalized Gradient Correction approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
cc      DOUBLE PRECISION  D(NSPIN), DECDD, DECDGD, DEXDD, DEXDGD,
cc     .                  EPSC, EPSX, GD(3,NSPIN)
      DOUBLE PRECISION  D(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN), 
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN),
     .                  EPSC, EPSX, GD(3,NSPIN)

C after JLM
      INTEGER I

      DO I=1,NSPIN
        IF(D(I) .LT. 0.0D0) then
           write(0,*) 'density <0 in ggaxc. dens=', d(i)
           D(I) = 0.0D0
        endif
      ENDDO
C after JLM

      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE IF (AUTHOR.EQ.'RPBE' .OR. AUTHOR.EQ.'rpbe') THEN
        CALL RPBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE IF (AUTHOR.EQ.'REVPBE' .OR. AUTHOR.EQ.'revpbe') THEN
        CALL REVPBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE IF (AUTHOR.EQ.'LYP'.OR.AUTHOR.EQ.'lyp') THEN
        CALL BLYPXC(NSPIN,D,GD,EPSX,EPSC,dEXdd,dECdd,dEXdgd,dECdgd)
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END






