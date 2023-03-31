! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
module m_tbt_GF

  use precision, only : dp
  use m_ts_gf, only : read_Green, check_Green, read_next_GS
  use m_ts_gf, only : reread_Gamma_Green

  implicit none

  public :: do_Green
  public :: read_Green
  public :: check_Green
  public :: reread_Gamma_Green
  public :: read_next_GS

  private

contains

  subroutine do_Green(El, &
       ucell,nkpnt,kpoint,kweight, &
       xa_EPS, CalcDOS)
    
    use parallel  , only : IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast, MPI_Integer, MPI_Logical
#endif
    use m_os, only : file_exist
    use m_ts_cctype
    use m_ts_electype
    use m_ts_electrode, only : create_Green

    use m_tbt_contour

    implicit none
    
! ***********************
! * INPUT variables     *
! ***********************
    type(Elec), intent(inout)     :: El
    integer, intent(in)           :: nkpnt ! Number of k-points
    real(dp), intent(in)          :: kpoint(3,nkpnt) ! k-points
    real(dp), intent(in)          :: kweight(nkpnt) ! weights of kpoints
    real(dp), intent(in)          :: xa_Eps ! coordinate precision check
    real(dp), dimension(3,3)      :: ucell ! The unit cell of the CONTACT
    logical, intent(in)           :: CalcDOS

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: uGF, i, iE, NEn
    logical :: errorGF, exist, cReUseGF
    complex(dp), allocatable :: ce(:)
    type(ts_c_idx) :: c
#ifdef MPI
    integer :: MPIerror
#endif

    ! fast exit if the Gf-file should not be created
    ! i.e. this means the calculation of the self-energy is
    ! performed in every iteration
    if ( .not. El%out_of_core ) return

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE do_Green' )
#endif
    
! check the file for existance
    exist = file_exist(trim(El%GFfile), Bcast = .true. )
    
    cReUseGF = El%ReUseGf
! If it does not find the file, calculate the GF
    if ( exist ) then
       if (IONode ) then
          write(*,*) "Electrode Green's function file: '"//&
               trim(El%GFfile)//"' already exist."
          if ( .not. cReUseGF ) then
             write(*,*)"Green's function file '"//&
                  trim(El%GFfile)//"' is requested overwritten."
          end if
       end if
    else
       cReUseGF = .false.
    end if

    errorGF = .false.

    ! we need to create all the contours
    NEn = N_TBT_E()
    allocate(ce(NEn))
    iE = 0
    if ( El%Eta > 0._dp ) then
      do i = 1 , NEn
        c = tbt_E(i)
        ! We ensure to add the complex imaginary value
#ifdef TBT_PHONON
        ce(i) = cmplx(real(c%e,dp)**2,El%Eta, dp)
#else
        ce(i) = cmplx(real(c%e,dp),El%Eta, dp)
#endif
      end do
    else
      do i = 1 , NEn
        ! We ensure to add the complex imaginary value
        c = tbt_E(i)
#ifdef TBT_PHONON
        ce(i) = cmplx(real(c%e,dp) ** 2, aimag(c%e)**2, dp)
#else
        ce(i) = c%e
#endif
      end do
    end if
    
    ! We return if we should not calculate it
    if ( .not. cReUseGF ) then

      call create_Green(El, &
          ucell,nkpnt,kpoint,kweight, &
          NEn,ce)

    else
      
      ! Check that the Green's functions are correct!
      ! This is needed as create_Green returns if user requests not to
      ! overwrite an already existing file.
      ! This check will read in the number of orbitals and atoms in the
      ! electrode surface Green's function.
      ! Check the GF file
      if ( IONode ) then
        call io_assign(uGF)
        open(file=El%GFfile,unit=uGF,form='UNFORMATTED')

        call check_Green(uGF,El, &
            ucell,nkpnt,kpoint,kweight, &
            NEn, ce, &
            xa_Eps, errorGF)

        write(*,'(/,4a,/)') "Using GF-file '",trim(El%GFfile),"'"

        call io_close(uGF)
      endif

    end if

!     Check the error in the GF file
#ifdef MPI
    call MPI_Bcast(errorGF,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif
    if ( errorGF ) &
         call die("Error in GFfile: "//trim(El%GFfile)//". Please move or delete")

    deallocate(ce)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS do_Green' )
#endif

  end subroutine do_Green

end module m_tbt_GF
