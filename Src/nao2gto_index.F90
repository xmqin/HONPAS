      MODULE extended_index

      private 
      public :: extended_index_init, indexsc 

      integer, pointer, save :: INDSC1(:), INDSC2(:)
      logical, save          :: nullified_pointers = .false.

      CONTAINS

      SUBROUTINE extended_index_init( NSC, NUO )
      use alloc, only: re_alloc

      IMPLICIT NONE

      integer, intent(in)   :: nsc(3)
      integer, intent(in)   :: nuo

      INTEGER I1, I2, I3, IC, J1, J2, J3, JC, &
              KUO, LASTIO, LASTJO, NCELLS, NO
      EXTERNAL MEMORY
      
      NCELLS = NSC(1) * NSC(2) * NSC(3)
      NO = NUO * NCELLS

      if (.not. nullified_pointers) then
         nullify(indsc1, indsc2)
         nullified_pointers = .true.
      endif

      call re_alloc(INDSC1,1,8*NO,name="indsc1",routine="extended_index_init")
      call re_alloc(INDSC2,1,NO,name="indsc2",routine="extended_index_init")

      DO J3 = 0,2*NSC(3)-1
      DO J2 = 0,2*NSC(2)-1
      DO J1 = 0,2*NSC(1)-1

        I1 = MOD(J1,NSC(1))
        I2 = MOD(J2,NSC(2))
        I3 = MOD(J3,NSC(3))

        IC = I1 +   NSC(1)*I2 +   NSC(1)*NSC(2)*I3
        JC = J1 + 2*NSC(1)*J2 + 4*NSC(1)*NSC(2)*J3

        LASTIO = IC * NUO
        LASTJO = JC * NUO

        DO KUO = 1,NUO
          INDSC1(LASTJO+KUO) = LASTIO + KUO
        ENDDO

        IF (I1+NSC(1).EQ.J1 .AND. I2+NSC(2).EQ.J2 .AND. I3+NSC(3).EQ.J3) THEN
          DO KUO = 1,NUO
            INDSC2(LASTIO+KUO) = LASTJO + KUO
          ENDDO
        ENDIF

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE extended_index_init

      FUNCTION INDEXSC( IO, IUO, JO )

      IMPLICIT NONE

      integer INDEXSC
      INTEGER, intent(in) ::  IO, IUO, JO

      INDEXSC = INDSC1( INDSC2(JO) - INDSC2(IO) + INDSC2(IUO) )

      END FUNCTION INDEXSC

      END MODULE extended_index
