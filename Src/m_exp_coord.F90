module m_exp_coord
! 
! This file is part of the SIESTA package.
!
! Copyright (c) Nick Papior Andersen 2012.
! Email: nickpapior@gmail.com
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!
! Module for handling the 'MD.TypeRun ExplicitCoordinates' option.
! 
! The purpose of this module is to run through specific coordinates of interest.
! The main drive for this extension is the creation of a mean Hamiltonian and/or density matrix.
! This might not always be realistic in terms of preserving charge in the system, however, it 
! gives a new tool to utilize the mean functional of such a quantity.
!
! Furthermore this could easily be used for Inelastica etc. in more ingenious ways.
!
! At the moment it relies heavily on the NetCDF interface, and as such it can only be used with this
! extension.

  use parallel, only : IONode
  use precision, only : dp
  
  implicit none

  ! The file we will save the explicit coordinates to
  character(len=250), save :: expCoordFile

  private

#ifdef NCDF_4

  public :: exp_coord_init
  public :: exp_coord_next
  public :: exp_coord_weight
  public :: read_exp_coord_options

  ! There are only routines present if the NCDF_4 has been compiled
  ! Hence the contains *MUST* not be there in that case.

contains

  subroutine exp_coord_init(slabel,na_u,inicoor,fincoor)
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta
#endif
    character(len=*), intent(in) :: slabel
    integer, intent(in)          :: na_u
    integer, intent(inout)       :: inicoor, fincoor
    ! Local variables...
    logical :: exist
    character(len=10) :: unit
    type(hNCDF) :: ncdf
    integer :: xa_size(3), i ! We need to check that the size is correct
#ifdef MPI
    integer :: MPIerror
#endif

    ! Read in the number of coordinates in the explicit coordinates...
    call ncdf_open(ncdf,trim(expcoordFile))
    call ncdf_inq_dim(ncdf,'n',len=fincoor)
    
    ! Next we read in the start of the simulation coordinate
    ! Initialize the initial coord to 0 
    inicoor = 0
    
    ! Retrieve the current reached explicit coordinate step...
    if ( IONode ) then
       inquire(file=trim(slabel)//'.EXPCOORD',exist=exist)
       if ( exist ) then
          call io_assign(i)
          open(unit=i,file=trim(slabel)//'.EXPCOORD',status='old', &
               form='formatted')
          read(i,*) inicoor
          call io_close(i)
       end if
    end if

    ! Make sure that the coordinate is the correct size...
    call ncdf_inq_var(ncdf,'xa',size=xa_size(:))
    if ( IONode .and. (xa_size(1) /= 3 .or. xa_size(2) /= na_u) ) then
       call die('The size of the coordinates in: '//trim(expcoordFile)// &
            ' does not match the size of the system.')
    end if

    call ncdf_close(ncdf)

#ifdef MPI
    call MPI_Bcast(inicoor,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(fincoor,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif
    
    ! this will also catch the obscure case where the user has no coordinates in the file
    if ( inicoor > fincoor .or. fincoor == 0 ) then
       call die("The explicit coordinate run is completed. No need to start.")
    end if

  end subroutine exp_coord_init


! ***********************************************************
! * Created by Nick Papior Andersen                         *
! *                                                         *
! *   Routine for reading the next step of coordinates      *
! ***********************************************************
  subroutine exp_coord_next(istep,na_u,xa)
    use fdf, only : fdf_convfac
    use files, only : slabel
    use netcdf_ncdf
    use units, only : Ang
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Double_Precision
#endif
! **************************
! * INPUT variables        *
! **************************
    integer, intent(in)   :: istep, na_u
! **************************
! * OUTPUT variables       *
! **************************
    real(dp), intent(out) :: xa(3,na_u)

! **************************
! * Local variables        *
! **************************
    integer :: xa_coord_size(3), i
    character(len=10) :: unit
    logical :: exist
    type(hNCDF) :: ncdf
#ifdef MPI
    integer :: MPIerror
#endif

    ! First check for the existance of the coordinate file
    call ncdf_open(ncdf,trim(expCoordFile),mode=NF90_WRITE)
    call ncdf_get_var(ncdf,'xa',xa,start=(/1,1,istep/))
    ! Check the unit of the variable...
    
    call ncdf_close(ncdf)

    ! Save the current reached explicit coordinate step...
    if ( IONode ) then
       call io_assign(i)
       open(unit=i,file=trim(slabel)//'.EXPCOORD',status='replace', &
            form='formatted')
       write(i,*) istep
       call io_close(i)
    end if

#ifdef MPI
    call MPI_Bcast(xa,3*na_u,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#endif
    
  end subroutine exp_coord_next

! ***********************************************************
! * Created by Nick Papior Andersen                         *
! *                                                         *
! *   Routine for reading the next step of coordinates      *
! ***********************************************************
  real(dp) function exp_coord_weight(istep)
    use netcdf_ncdf
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Double_Precision
#endif
! **************************
! * INPUT variables        *
! **************************
    integer, intent(in) :: istep

! **************************
! * Local variables        *
! **************************
    type(hNCDF) :: ncdf
#ifdef MPI
    integer :: MPIerror
#endif
    
    ! First check for the existance of the coordinate file
    call ncdf_open(ncdf,trim(expCoordFile))
    call ncdf_get_var(ncdf,'w',exp_coord_weight,start=(/istep/))
    call ncdf_close(ncdf)

#ifdef MPI
    call MPI_Bcast(exp_coord_weight,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#endif

  end function exp_coord_weight


! ***************************************
! *     Routine for reading options     *
! ***************************************
  subroutine read_exp_coord_options()
    use fdf
    use m_os, only : file_exist
    ! Get the Name of the NetCDF file which holds the coordinates
    expCoordFile = fdf_get('MD.ExpCoord.File','COORDS.nc')

    ! Ensure that the coordinate file does in fact exist...
    if ( IONode .and. .not. file_exist(expCoordFile) ) then
       call die('File: '//trim(expCoordFile)//' does not exist. Please create a &
            &compliant coordinate file.')
    end if
  end subroutine read_exp_coord_options

#endif
  
end module m_exp_coord
