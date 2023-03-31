! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      subroutine write_md_record( istep )
        Use siesta_options
        use siesta_geom
        use atomlist,          only: iza
        use md_out,            only: md_v_format
        use parallel,          only: IOnode
        use m_energies,        only: Etot, Eharrs, Ekinion
        use m_steps
        use m_kinetic
        use m_stress, only: tstress
        use units, only: Ang, kbar
#ifdef CDF
        use md_out,     only: md_netcdf
#endif

        implicit none

        integer, intent(in) ::  istep

        real(dp) :: Etot_output
        real(dp):: getot  ! Total energy including thermostat dynamic

        real(dp) :: Psol    ! Total pressure (static + kinetic)


        ! Save atomic positions and velocities accumulatively
        if (writmd.and.IOnode) then
          if (.not. harrisfun) then
            Etot_output = Etot
          else
            Etot_output = Eharrs
          endif
          getot = Etot_output + Ekinion + kn + kpr + vn + vpr
          Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
          call iomd( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, istep, inicoor, 
     .        tempion, Etot_output, getot,
     .        volume_of_some_cell/Ang**3, Psol/kbar)
          call md_v_format(na_u,isa,xa,ucell)
#ifdef CDF
          call md_netcdf( na_u, isa, iza, 
     .         xa, va, ucell, vcell, varcel, 
     .         tempion, Etot_output, getot,
     .         volume_of_some_cell/Ang**3, Psol/kbar)
#endif
        endif

        ! Accumulate coor in Xmol file for animation 
        if (writpx.and.IOnode) call pixmol(iza, xa, na_u )
      END
