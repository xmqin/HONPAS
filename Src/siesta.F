      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor

      implicit none

      integer :: istep
      logical :: relaxd

!------------------------------------------------------------------------- BEGIN
      call siesta_init()

      ! Begin of coordinate relaxation iteration
      relaxd = .false.
      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))

        call siesta_forces( istep )

        call siesta_move( istep, relaxd )
        if (.not. relaxd) then
          istep = istep + 1
        endif

      ENDDO
      ! End of coordinate-relaxation loop 

      call siesta_analysis( relaxd )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
