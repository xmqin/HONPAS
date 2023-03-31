! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_save_density_matrix
  
  implicit none
  public :: save_density_matrix

contains
  
  subroutine save_density_matrix(SCFconverged, when, file)
    
    ! Stores DM and EDM on files
    ! This version uses module variables.
    ! Eventually it can be a cleaner routine
    
    use precision,            only:  dp
    use class_dSpData2D,      only:  val
    use sparse_matrices,      only:  maxnh
#ifdef TIMING_IO
    use sparse_matrices,      only:  numh
    use sparse_matrices,      only:  listh, listhptr
#endif
    use files, only : slabel
    use m_iodm,               only:  write_dm
    use m_matio,              only:  write_mat
    use m_spin,               only:  spin
    use atomlist,             only:  no_l
    use siesta_geom,          only:  nsc
    use siesta_options,       only:  writedm, writedm_cdf
    use siesta_options,       only:  writedm_cdf_history
    use siesta_options,       only:  idyn, nmove, fixspin
    use m_steps,              only:  istp
    use sparse_matrices,      only:  Dscf, DM_2D
    use fdf,                  only:  fdf_get
#ifdef CDF
    use iodm_netcdf,          only:  write_dm_netcdf
#ifdef NCDF_4
    use dictionary
    use siesta_options, only:  write_cdf
    use m_ncdf_siesta, only : cdf_save_state
#endif
#endif
    use sparse_matrices,      only:  EDM_2D
    use m_ts_iodm,            only:  write_ts_dm
    use m_ts_global_vars,     only:  TSrun
    use m_ts_options,         only:  TS_DE_save
    use m_energies,           only:  Ef, Efs

    implicit none
    
    logical, intent(in) :: SCFconverged
    logical, intent(in), optional  :: when
    character(len=*), intent(in), optional   :: file

    real(dp), pointer :: EDM(:, :), DM(:,:)
    
    logical :: do_write
    
#ifdef CDF
#ifdef NCDF_4
    type(dictionary_t) :: dic_save
#endif
#endif
#ifdef TIMING_IO
    integer :: i
#endif

    ! Control the condition to write the DM externally
    ! This gives more flexibility
    ! We retain the "writedm" default condition for compatibility
    
    if ( present(when) ) then
      do_write = when
    else
      do_write = writedm
    end if

#ifdef CDF
    ! Save density matrix on disk, after mixing, to
    ! be used for re-starting the SCF cycle.
    if ( writedm_cdf_history ) then
      call write_dm_netcdf( no_l, maxnh, spin%H, Dscf, overwrite=.false. )
    else if ( writedm_cdf ) then
      call write_dm_netcdf( no_l, maxnh, spin%H, Dscf, overwrite=.true. )
    end if
    
#ifdef NCDF_4
    if ( write_cdf .and. do_write ) then
      
      if ( (idyn == 0 .and. nmove == 0) .or. &
          (idyn == 6 .and. istp == 1 ) ) then
        
        dic_save = ('Ef'.kv.1)//('DM'.kv.1)//('EDM'.kv.1)
        call cdf_save_state(trim(slabel)//'.nc',dic_save)
        
        call delete(dic_save)
      end if
      
    end if
#endif
#endif

#ifdef TIMING_IO
    if ( do_write ) then
      call timer('IO-W-DM',1)
      do i = 1 , 100
        if (fdf_get("Use.Blocked.WriteMat",.false.)) then
          call write_mat (maxnh, no_l, spin%H, &
              numh, listhptr, listh, Dscf, &
              userfile=trim(file)//'.blocked',compatible=.false.)
        else
          ! The filename falls back to the standard one     
          call write_dm(trim(slabel)//'.DM', nsc, DM_2D)
        end if
      end do
      call timer('IO-W-DM',2)
      call timer('IO-W-DM',3)
    end if
#endif
    
    ! Check wheter we should write
    select case ( idyn )
    case ( 6, 7, 9 )
      if ( istp /= 1 ) do_write = .false.
    end select

    ! Figure out when to write TSDE
    if ( TSrun ) then
      ! We should only store based on the input of the super routine.
      ! I.e. if the user does not request storing the mixed DM, then
      ! we need not save the TSDE since it is already stored from the
      ! previous iteration.
    else
      if ( do_write ) then
        ! If we are running the Siesta initialization we should
        ! store the DM file, but only if do_write allows it ;)
        call write_dm(trim(slabel)//'.DM', nsc, DM_2D)
      end if
      
      ! This is not a transiesta run, possibly the siesta initialization run
      ! Here we should only store the TSDE file if it has converged *AND*
      ! if TS_DE is requested to be saved
      do_write = SCFconverged .and. TS_DE_save
    end if

    if ( do_write ) then
      if ( fixspin ) then
        ! For fixed spin calculations we shift the energy-density according
        ! to the first spin
        DM => val(DM_2D)
        EDM => val(EDM_2D)
        call daxpy(size(DM,1),Efs(1)-Efs(2),DM(1,2),1,EDM(1,2),1)
        Ef = Efs(1)
      end if
      call write_ts_dm(trim(slabel)//'.TSDE', nsc, DM_2D, EDM_2D, Ef)
      if ( fixspin ) then
        ! Shift back
        call daxpy(size(DM,1),Efs(2)-Efs(1),DM(1,2),1,EDM(1,2),1)
      end if
    end if
    
  end subroutine save_density_matrix
  
end module m_save_density_matrix
