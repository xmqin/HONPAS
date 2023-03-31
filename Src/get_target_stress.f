! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      subroutine get_target_stress(tp,tstres)
      use precision, only: dp
      use parallel, only: ionode
      use fdf
      use units, only: kBar
      use sys,   only: die

      implicit none

      real(dp), intent(in) :: tp
      real(dp), intent(out) :: tstres(3,3)

      integer  :: i, j
      real(dp) :: sxx, syy, szz, sxy, sxz, syz
      logical  :: tarstr

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Look for target stress and read it if found, otherwise generate it --------
          tarstr = fdf_block('MD.TargetStress',bfdf)

          if (tarstr) then
            if (ionode) then
              write(6,'(/a,a)')
     $             'Broyden_optim: Reading %block MD.TargetStress',
     .             ' (units of MD.TargetPressure).'
            endif
            if (.not. fdf_bline(bfdf,pline))
     $        call die('get_target_stress: ERROR in ' //
     .                 'MD.TargetStress block')
            sxx = fdf_bvalues(pline,1)
            syy = fdf_bvalues(pline,2)
            szz = fdf_bvalues(pline,3)
            sxy = fdf_bvalues(pline,4)
            sxz = fdf_bvalues(pline,5)
            syz = fdf_bvalues(pline,6)
            call fdf_bclose(bfdf)
            
            tstres(1,1) = - sxx * tp
            tstres(2,2) = - syy * tp
            tstres(3,3) = - szz * tp
            tstres(1,2) = - sxy * tp
            tstres(2,1) = - sxy * tp
            tstres(1,3) = - sxz * tp
            tstres(3,1) = - sxz * tp
            tstres(2,3) = - syz * tp
            tstres(3,2) = - syz * tp
          else
            if (ionode) then
              write(6,'(/a,a)')
     $              'Broyden_optim: No target stress found, ',
     .              'assuming hydrostatic MD.TargetPressure.'
            endif
            do i= 1, 3
              do j= 1, 3
                tstres(i,j) = 0._dp
              enddo
              tstres(i,i) = - tp
            enddo
          endif

C Write target stress down --------------------------------------------------
          if (ionode) then
            write(6,"(/a)") 'Broyden_optim: Target stress (kBar)'
            do i = 1, 3
              write(6,"(a,2x,3f12.3)") 
     .            'get_target_stress:',
     .              tstres(i,1)/kBar, tstres(i,2)/kBar, 
     .            tstres(i,3)/kBar
            enddo
          endif

          end subroutine get_target_stress
