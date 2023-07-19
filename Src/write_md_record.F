      subroutine write_md_record( istep )
        Use siesta_options
        use siesta_geom
        use atomlist,          only: iza
        use md_out,            only: md_v_format
        use parallel,          only: IOnode
        use m_energies,        only: Etot, Eharrs1, Ekinion
        use m_steps
        use m_kinetic
        use m_stress, only: tstress
#ifdef CDF
        use md_out,     only: md_netcdf
#endif

        implicit none

        integer, intent(in) ::  istep

        real(dp) :: Etot_output
        real(dp):: getot  ! Total energy including thermostat dynamic variables

        real(dp) :: Psol    ! Total pressure (static + kinetic)

        logical:: lastst  


        ! Save atomic positions and velocities accumulatively
        if (writmd.and.IOnode) then
          if (.not. harrisfun) then
            Etot_output = Etot
          else
            Etot_output = Eharrs1
          endif
          getot = Etot_output + Ekinion + kn + kpr + vn + vpr
          Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
          call iomd( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, istep, inicoor, 
     .        fincoor, tempion, Etot_output, getot,
     .        volume_of_some_cell/Ang**3, Psol/kbar)
          call md_v_format(na_u,isa,xa,ucell)
#ifdef CDF
            call md_netcdf( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, 
     .        tempion, Etot_output, getot,
     .        volume_of_some_cell/Ang**3, Psol/kbar)
#endif

        endif

        ! Accumulate coor in Xmol file for animation 
        lastst = fincoor .le. istep
        if (writpx.and.IOnode) then
          call pixmol(iza, xa, na_u, lastst)
        endif
      END
