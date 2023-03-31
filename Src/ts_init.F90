!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_ts_init

  implicit none
  private
  public :: ts_init

contains

  subroutine ts_init(nspin, ucell, na_u, xa, lasto, no_u, inicoor, fincoor )
  ! Routine for initializing everything related to the Transiesta package.
  ! This is to comply with the SIESTA way of initializing at the beginning
  ! and make the code more managable.

  ! This routine does the following in that order:
  !  1. Read in ts_options
  !  2. Create Transiesta k-points
  !  3. Setup the contour path
  !  4. Create the GF related to the integration scheme
  
! Used modules
    use parallel, only : IONode

    use m_os, only : file_exist
    use fdf, only: fdf_get

    use m_ts_gf,        only : do_Green
    use m_ts_electrode, only : init_Electrode_HS
    
    use m_ts_kpoints, only : setup_ts_kpoint_grid
    use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
    use m_ts_kpoints, only : ts_kscell, ts_kdispl, ts_gamma
    use m_ts_cctype
    use m_ts_electype
    use m_ts_options ! Just everything (easier)
    use m_ts_method

    use m_ts_global_vars, only : TSmode, TSinit, onlyS
    use siesta_options, only : isolve, SOLVE_TRANSI, Nmove

    use m_fixed, only : is_fixed, is_constr

    implicit none
! *********************
! * INPUT variables   *
! *********************
    integer, intent(in)  :: nspin
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: lasto(0:na_u)
    integer, intent(in)  :: no_u
    integer, intent(in)  :: inicoor, fincoor
    
! *********************
! * LOCAL variables   *
! *********************
    logical :: neglect_conn
    integer :: i, ia
    integer :: nC, nTS

    if ( isolve .eq. SOLVE_TRANSI ) then
       TSmode = .true.
       ! If in TSmode default to initalization
       ! In case of 'DM.UseSaveDM TRUE' TSinit will be set accordingly
       TSinit = .true.
    end if

    ! initialize regions of the electrodes and device
    ! the number of LCAO orbitals on each atom will not change
    call ts_init_regions('TS', na_u, lasto)

    ! Read generic transiesta options
    call read_ts_generic( ucell )

    ! read the chemical potentials
    call read_ts_chem_pot( )

    ! read in the electrodes
    call read_ts_elec( ucell, na_u, xa, lasto )

    ! Read in the k-points
    call setup_ts_kpoint_grid( ucell )

    ! Read after electrode stuff
    call read_ts_after_Elec( ucell, nspin, na_u, xa, lasto, &
         ts_kscell, ts_kdispl)

    ! Print the options
    call print_ts_options( ucell )

    ! Print all warnings
    call print_ts_warnings( ts_Gamma, ucell, na_u, xa, Nmove )

    ! If we actually have a transiesta run we need to process accordingly!
    if ( .not. TSmode ) return

    ! If onlyS we do not need to do anything about the electrodes
    if ( onlyS ) return

    ! Print out the contour blocks etc. for transiesta
    call print_ts_blocks( na_u, xa )

    ! Check that an eventual CGrun will fix all electrodes and 
    ! buffer atoms
    if ( fincoor - inicoor > 0 ) then
       ! check fix
       do ia = 1 , na_u
          if ( .not. a_isBuffer(ia) ) cycle
          if ( .not. is_fixed(ia) ) then
             call die('All buffer atoms *MUST* be &
                  &fixed while doing transiesta geometry optimizations. &
                  &Please correct Geometry.Constraints for buffer atoms.')
          end if
       end do
       do i = 1 , N_Elec
          do ia = Elecs(i)%idx_a , Elecs(i)%idx_a + TotUsedAtoms(Elecs(i)) - 1
             if ( .not. ( is_constr(ia,'rigid') .or. is_constr(ia,'rigid-dir') &
                  .or. is_constr(ia,'rigid-max') .or. is_constr(ia,'rigid-max-dir') &
                  .or. is_fixed(ia) ) ) then
                call die('All electrode atoms *MUST* be &
                     &fixed while doing transiesta geometry optimizations. &
                     &Please correct Geometry.Constraints for electrode atoms.')
             end if
          end do
       end do
    end if

    ! here we will check for all the size requirements of Transiesta
    ! I.e. whether some of the optimizations can be erroneous or not

    ! Do a crude check of the sizes
    ! if the transiesta region is equal of size to or smaller 
    ! than the size of the combined electrodes, then the system
    ! is VERY WRONG...
    ! First calculate L/C/R sizes (we remember to subtract the buffer
    ! orbitals)
    nTS = no_u - no_Buf
    nC  = nTS  - sum(TotUsedOrbs(Elecs))
    if ( nC < 1 ) &
         call die("The contact region size is &
         &smaller than the electrode size. &
         &What have you done? Please correct this insanity...")
    
    if ( minval(TotUsedOrbs(Elecs)) < 2 ) &
         call die('We cannot perform sparse pattern on the electrode &
         &system.')
    
    ! Show every region of the Transiesta run
    call ts_show_regions(ucell,na_u,xa,N_Elec,Elecs)

    if ( .not. TS_Analyze ) then

       ! GF generation:
       do i = 1 , N_Elec

          ! initialize the electrode for Green function calculation
          call init_Electrode_HS(Elecs(i))

          ! Whether we can do with a single k-point (say a chain)
          if ( Elecs(i)%is_gamma ) then
            call do_Green(Elecs(i), &
                ucell,1,(/(/0._dp, 0._dp, 0._dp/)/),(/1._dp/), &
                Elecs_xa_Eps, .false. )

          else
            call do_Green(Elecs(i), &
                ucell,ts_nkpnt,ts_kpoint,ts_kweight, &
                Elecs_xa_Eps, .false. )

          end if
          
          ! clean-up
          call delete(Elecs(i))
          
       end do
    else

       neglect_conn = .false.

       do i = 1 , N_Elec
          call delete(Elecs(i)) ! ensure clean electrode
          call read_Elec(Elecs(i),Bcast=.true.)
          
          ! print out the precision of the electrode (whether it extends
          ! beyond first principal layer)
          if ( .not. check_Connectivity(Elecs(i)) ) then

             neglect_conn = .true.
             
          end if

          call delete(Elecs(i))
       end do

       if ( neglect_conn ) then
          neglect_conn = fdf_get('TS.Elecs.Neglect.Principal', .false.)

          if ( .not. neglect_conn ) then

             call die('Electrode connectivity is not perfect, refer &
                  &to the manual for achieving a perfect electrode.')

          end if

       end if

    end if

  end subroutine ts_init

end module m_ts_init

