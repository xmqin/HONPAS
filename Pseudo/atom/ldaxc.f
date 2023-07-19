      SUBROUTINE LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, VX, VC )

C Finds the exchange and correlation energies and potentials, in the
C Local (spin) Density Approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), EPSC, EPSX, VX(NSPIN), VC(NSPIN)

C after JLM
      INTEGER I

      DO I=1,NSPIN
        IF(D(I) .LT. 0.0D0) then
           write(0,*) 'density <0 in ldaxc. dens=', d(i)
           D(I) = 0.0D0
        endif
      ENDDO
C after JLM

      IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.
     .     AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        CALL PZXC( IREL, NSPIN, D, EPSX, EPSC, VX, VC )
      ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        CALL PW92XC( IREL, NSPIN, D, EPSX, EPSC, VX, VC )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END
