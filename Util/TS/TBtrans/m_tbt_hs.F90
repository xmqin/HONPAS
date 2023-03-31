! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! Module to contain all TSHS files
! We use this module to interpolate between two bias
! calculations, and hence withhold all information related to the 
! calculation.

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication
module m_tbt_hs

  use precision, only : dp
  use class_OrbitalDistribution
  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  implicit none

  ! Generic attributes
  public
  save

  type :: tTSHS
     ! The system that we are currently investigating
     type(OrbitalDistribution) :: dit
     type(Sparsity) :: sp
     integer :: nspin
     type(dSpData1D) :: S_1D
     type(dSpData2D) :: H_2D
     ! We do not need the xij array
     integer :: nsc(3)
     integer, pointer :: isc_off(:,:) => null()
     real(dp), pointer :: sc_off(:,:) => null()
     ! The system coordinate information
     integer :: na_u, no_u
     real(dp), pointer :: xa(:,:) => null()
     integer, pointer :: lasto(:) => null()

     ! The voltage associated with this Hamiltonian object.
     real(dp) :: Volt

     ! The unit-cell
     real(dp) :: cell(3,3)

  end type tTSHS

  type(tTSHS) :: TSHS

  ! A list of hamiltonian files and their biases
  type :: tHSfile
     character(len=256) :: HSfile ! The Hamiltonian file
     real(dp) :: Volt ! The bias
  end type tHSfile

  ! the HSfiles
  integer :: N_HS = 0
  type(tHSfile), allocatable :: tHS(:)

  ! The Hamiltonians are associated with 
  ! a spin-component.
  ! However, we can easily imagine that we split the two calculations
  ! by doing : 1st node == ispin == 1
  !            2nd node == ispin == 2
  integer :: spin_idx

  ! The voltage we are working at
  ! It is given here for reasons questionable.
  ! However, it is related to the Hamiltonian, so it makes "some" sense
  real(dp) :: Volt
  ! whether we will use the bias-contour
  logical :: IsVolt = .false.

contains

  subroutine tbt_init_HSfile( )

    use fdf
    use files, only : slabel
    use units, only : eV
    use m_interpolate
    use m_spin, only: init_spin

    use m_ts_io, only: ts_read_TSHS_opt
    use m_ts_io_ctype, only : ts_c_bphysical, ts_c_bisphysical
    
    integer :: iHS, nspin
    ! For reading in the TSHS file block
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    type(tHSfile), allocatable :: tmpHS(:)

    ! Create a pivot table for the Hamiltonians
    integer, allocatable :: ipvt(:)
    character(len=50) :: chars
    logical :: btmp
    
    ! The bias, you can say, is associated
    ! with the Hamiltonian
    Volt = fdf_get('TS.Voltage',0._dp,'Ry') 
    Volt = fdf_get('TBT.Voltage',Volt,'Ry')
    ! Voltage situation is above 0.01 mV
    IsVolt = abs(Volt) > 0.00001_dp * eV

    ! Read in device region via the new block
    btmp = fdf_block('TBT.HS.Files',bfdf)

    if ( btmp ) then
       
       ! First we read number of entries
       N_HS = 0
       do while ( fdf_bline(bfdf,pline) ) 
          if ( fdf_bnnames(pline) == 0 ) cycle
          N_HS = N_HS + 1
       end do
       if ( N_HS == 1 ) then
          call die('You cannot ask for interpolation of &
               &one TSHS file. What were you thinking?')
       end if

       ! backspace to the block
       if ( .not. fdf_block('TBT.HS.Files',bfdf) ) then
          call die('Error on second reading of block')
       end if

       allocate(tmpHS(N_HS))
       iHS = 0
       do while ( fdf_bline(bfdf,pline) ) 
          if ( fdf_bnnames(pline) == 0 ) cycle
          
          iHS = iHS + 1
          ! we read the first name as the TSHS file
          tmpHS(iHS)%HSfile = fdf_bnames(pline,1)
          ! The first value on the line is the 
          !  bias in units
          tmpHS(iHS)%Volt = ts_c_bphysical(pline,1,'Ry')

       end do

       allocate(ipvt(N_HS))
       call crt_pivot(N_HS,tmpHS(:)%Volt,ipvt)

       chars = fdf_get('TBT.HS.Interp','spline')

       if ( leqi(chars,'linear') ) then

          ! We only do a linear interpolation...
          ! Remove all other than the two biases
          ! closest
          allocate(tHS(2))

          if ( Volt <= tmpHS(ipvt(1))%Volt ) then
             ! extrapolation to lower V
             tHS(1) = tmpHS(ipvt(1))
             tHS(2) = tmpHS(ipvt(2))
          else if ( tmpHS(ipvt(N_HS))%Volt <= Volt ) then
             ! extrapolation to HIGHER V
             tHS(1) = tmpHS(ipvt(N_HS-1))
             tHS(2) = tmpHS(ipvt(N_HS))
          else

             tHS(1)%HSfile = ' '
             ! We need to find the closest value
             do iHS = 2 , N_HS
                if ( Volt < tmpHS(ipvt(iHS))%Volt ) then
                   tHS(1) = tmpHS(ipvt(iHS-1))
                   tHS(2) = tmpHS(ipvt(iHS))
                   exit
                end if
             end do

             if ( len_trim(tHS(1)%HSfile) == 0 ) then
                call die('Error in linear sorting algorithm, &
                     &could not figure out the bias placement...')
             end if

          end if

          N_HS = 2

       else if ( leqi(chars,'spline') ) then
          
          allocate(tHS(N_HS))
          do iHS = 1 , N_HS
             tHS(iHS) = tmpHS(ipvt(iHS))
          end do
          
       else

          call die('Unknown interpolation scheme')

       end if

       deallocate(tmpHS,ipvt)

    else

       ! We expect a single file 
       N_HS = 1
       allocate(tHS(1))

       tHS(1)%HSfile = fdf_get('TBT.HS',trim(slabel)//'.TSHS')
       tHS(1)%Volt = Volt

     end if

     ! Before we can read the spin information we have to
     ! read the default from the HS file.
     ! This is important for TB calculations where
     ! one need not specify "Spin polarized"
     call ts_read_TSHS_opt(tHS(1)%HSfile, nspin=nspin)

     ! Now we can read the spin-configuration using fdf-flags
     call init_spin( default_nspin=nspin )

    ! We read in the first file and 
    ! initialize the sparsity pattern.

    ! We determine which spin of the files we should use
    if ( nspin == 1 ) then
      ! If there is only one spin, then
      ! we do not read in the option.
      spin_idx = 0
    else if ( nspin > 2 ) then
      call die("TBtrans is currently not implemented for non-collinear or spin-orbit")
    else
      spin_idx = fdf_get('TBT.Spin',0)
      if ( spin_idx > nspin ) then
        call die('You have asked for a spin index not existing')
      else if ( spin_idx <= 0 ) then
        spin_idx = 0 ! all spins are used
      end if
    end if

    ! Start by creating the Hamiltonian!
    ! just prepare the next ispin..
    if ( spin_idx == 0 ) then
      call prep_next_HS(1, Volt)
    else
      call prep_next_HS(spin_idx, Volt)
    end if

  end subroutine tbt_init_HSfile
                                
  subroutine prep_next_HS(ispin, Volt)

    use parallel, only : IONode
    use units, only : eV
    use m_ts_io
    use m_sparsity_handling
    use m_handle_sparse, only : reduce_spin_size

    ! The two indices corresponding to the indices in the 
    ! tHS array
    integer, intent(in) :: ispin
    real(dp), intent(in) :: Volt

    ! The files (we read them all in...)
    type(tTSHS), allocatable :: files(:)

    type(OrbitalDistribution), pointer :: dit
    integer :: iHS

    ! Clean up the intrinsic stuff..
    call clean_HS()

    ! In case we are only having one TSHS file...
    if ( N_HS == 1 ) then

       call read_HS(ispin, tHS(1), TSHS)
       dit => dist(TSHS%S_1D)
       TSHS%dit = dit

       allocate(TSHS%sc_off(3,size(TSHS%isc_off,dim=2)))
       TSHS%sc_off = matmul(TSHS%cell,TSHS%isc_off)

       TSHS%Volt = Volt
       return

    end if

    ! We are requested to read everything...
    allocate(files(N_HS))

    if ( IONode ) then
       write(*,'(a,f8.4,a)') 'tbt: Interpolation of Hamiltonian to ',Volt/eV,' V'
    end if

    do iHS = 1 , N_HS

       ! Read in the Hamiltonian
       call read_HS(ispin,tHS(iHS),files(iHS))

       if ( IONode ) then
          write(*,'(2a,f8.4,a)') '  Hamiltonian ',trim(tHS(iHS)%HSfile),tHS(iHS)%Volt / eV,' V'
       end if

       if ( iHS > 1 ) then
          if ( .not. equivalent(files(iHS-1)%sp,files(iHS)%sp) ) then
             call die('Currently we do not support interpolating between two &
                  &different sparsity structures.')
          end if

          ! Do a couple of simple checks...
          if ( any(files(iHS-1)%lasto /= files(iHS)%lasto) ) then
             call die('The orbitals per atom are not consistent. ' // &
                  'This is pretty grave. Do not mix different simulations.')
          end if

          ! We know the sparsity pattern is the same...
          ! So align them
          call set_vars_HS(files(1),files(iHS))

       else
          
          ! Align the distributions (S and H have same distribution
          dit => dist(files(iHS)%S_1D)
          files(iHS)%dit = dit

       end if

    end do

    ! Debug prints, delete when we now the above "trick" works.
!    call print_type(files(1)%sp)
!    call print_type(files(1)%H_2D)
!    call print_type(files(N_HS)%sp)
!    call print_type(files(N_HS)%S_1D)
!    call print_type(files(N_HS)%H_2D)

    ! Now we interpolate the 
    call SpData_interp(N_HS,files(:)%H_2D,tHS(:)%Volt,Volt)

    ! Now files(1) contains the interpolated values
    ! copy files(1) to the original one...
    TSHS%nspin = files(1)%nspin
    TSHS%dit = files(1)%dit
    TSHS%sp = files(1)%sp
    TSHS%S_1D = files(1)%S_1D
    TSHS%H_2D = files(1)%H_2D
    TSHS%nsc(:) = files(1)%nsc(:)
    TSHS%isc_off => files(1)%isc_off
    TSHS%na_u = files(1)%na_u
    TSHS%no_u = files(1)%no_u
    TSHS%xa => files(1)%xa
    TSHS%lasto => files(1)%lasto
    TSHS%cell = files(1)%cell

    do iHS = 1 , N_HS
       call delete(files(iHS)%H_2D)
       call delete(files(iHS)%S_1D)
       call delete(files(iHS)%sp)
       call delete(files(iHS)%dit)
       ! do not delete the xa, lasto, isc_off arrays
       ! of the first instance, we point to them!
       if ( iHS > 1 ) then
          deallocate(files(iHS)%isc_off)
          deallocate(files(iHS)%xa)
          deallocate(files(iHS)%lasto)
       end if
       nullify(files(iHS)%isc_off)
       nullify(files(iHS)%xa)
       nullify(files(iHS)%lasto)
    end do

    deallocate(files)

    ! We should now have cleaned everything...

    ! Save the bias
    TSHS%Volt = Volt

    ! Create the sc_off array
    allocate(TSHS%sc_off(3,size(TSHS%isc_off,dim=2)))
    TSHS%sc_off = matmul(TSHS%cell,TSHS%isc_off)

  contains

    subroutine read_HS(ispin, tHS, file)
      integer, intent(in) :: ispin
      type(tHSfile), intent(in) :: tHS
      type(tTSHS), intent(inout) :: file

      integer :: kscell(3,3), istep, ia1
      real(dp) :: kdispl(3), Qtot, Temp, Ef
      logical :: onlyS, Gamma, TSGamma
      
      ! The easy case is when they have the same number
      call ts_read_TSHS(tHS%HSfile, onlyS, Gamma, TSGamma, &
           file%cell, file%nsc, file%na_u, file%no_u, file%nspin, &
           kscell, kdispl, &
           file%xa, file%lasto, &
           file%sp, file%H_2D, file%S_1D, file%isc_off, &
           Ef, Qtot, Temp, &
           istep, ia1, Bcast = .true. )

      call reduce_spin_size(ispin,file%H_2D,file%S_1D,Ef)

    end subroutine read_HS

    subroutine set_vars_HS(def,other)
      type(tTSHS), intent(inout) :: def, other

      type(Sparsity), pointer :: sp
      type(OrbitalDistribution), pointer :: dit

      sp => spar(other%S_1D)
      sp = def%sp
      sp => spar(other%H_2D)
      sp = def%sp

      dit => dist(other%S_1D)
      dit = def%dit
      dit => dist(other%H_2D)
      dit = def%dit

      other%dit = def%dit
      other%sp = def%sp

    end subroutine set_vars_HS
    
  end subroutine prep_next_HS
  
  subroutine clean_HS()
    call delete(TSHS%dit)
    call delete(TSHS%sp)
    call delete(TSHS%S_1D)
    call delete(TSHS%H_2D)
    if ( associated(TSHS%isc_off) ) then
       ! Everything must have been allocated
       deallocate(TSHS%isc_off,TSHS%sc_off)
       deallocate(TSHS%xa,TSHS%lasto)
       nullify(TSHS%isc_off, TSHS%sc_off)
       nullify(TSHS%xa,TSHS%lasto)
    end if
  end subroutine clean_HS

end module m_tbt_hs

  
