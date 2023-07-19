      subroutine get_target_stress(tp,tstres)
      use precision, only: dp
      use parallel, only: ionode
      use fdf, only: fdf_block
      use m_mpi_utils, only: broadcast
      use units, only: kBar

      implicit none

      real(dp), intent(in) :: tp
      real(dp), intent(out) :: tstres(3,3)

      integer  :: i, j, iu
      real(dp) :: sxx, syy, szz, sxy, sxz, syz
      logical  :: tarstr

C Look for target stress and read it if found, otherwise generate it --------

          if (ionode) then
            tarstr = fdf_block('MD.TargetStress',iu)

            if (tarstr) then
               write(6,'(/a,a)')
     $              'Reading %block MD.TargetStress',
     .                        ' (units of MD.TargetPressure).'
               read(iu,*, end=50) sxx, syy, szz, sxy, sxz, syz
               tstres(1,1) = - sxx * tp
               tstres(2,2) = - syy * tp
               tstres(3,3) = - szz * tp
               tstres(1,2) = - sxy * tp
               tstres(2,1) = - sxy * tp
               tstres(1,3) = - sxz * tp
               tstres(3,1) = - sxz * tp
               tstres(2,3) = - syz * tp
               tstres(3,2) = - syz * tp
   50          continue
            else
              write(6,'(/a,a)')
     $              'No target stress found, ',
     .                'assuming hydrostatic MD.TargetPressure.'
              do i = 1, 3
                do j = 1, 3
                  tstres(i,j) = 0._dp
                enddo
                tstres(i,i) = - tp
              enddo
            endif

C Write target stress down --------------------------------------------------

            write(6,"(/a)") 'Target stress (kBar)'
            do i = 1, 3
               write(6,"(3f12.3)") 
     .            tstres(i,1)/kBar, tstres(i,2)/kBar, 
     .            tstres(i,3)/kBar
            enddo
          endif  ! node0

          call broadcast(tstres(1:3,1:3))

          end subroutine get_target_stress
