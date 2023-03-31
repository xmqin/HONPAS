! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been implemented by:
!   Nick Papior, 2014.
!
! Handling all the initializations
! This opens MPI channels, ensures the options are read in etc.
!
subroutine tbt_init()

  use fdf
  use sys, only : die
  use precision, only : dp
  use parallel, only : parallel_init, Node, Nodes, IONode
  use m_timer, only : timer_report
  use alloc, only   : alloc_report
  use files, only   : slabel
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock
#ifdef NCDF_4
  use netcdf_ncdf, only : ncdf_IOnode
#endif
#ifdef MPI
  use mpi_siesta, only : MPI_Barrier, MPI_Comm_World
#ifdef _OPENMP
  use mpi_siesta, only : Mpi_Thread_Funneled
#endif
#endif

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  use dictionary

  use m_ts_electrode, only : init_Electrode_HS
  use m_ts_electype

  use m_tbt_kpoint
  use m_tbt_regions
  use m_tbt_kregions, only : tbt_init_kregions, tbt_print_kRegions
  use m_tbt_hs
  use m_tbt_options
  use m_tbt_contour
  use m_tbt_gf
  use m_tbt_save
  use m_tbt_proj

  use m_sparsity_handling

#ifdef _OPENMP
  use omp_lib, only : omp_get_num_threads
  use omp_lib, only : omp_get_schedule, omp_set_schedule
  use omp_lib, only : omp_get_proc_bind
  use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
  use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
  use omp_lib, only : OMP_PROC_BIND_FALSE, OMP_PROC_BIND_TRUE
  use omp_lib, only : OMP_PROC_BIND_MASTER
  use omp_lib, only : OMP_PROC_BIND_CLOSE, OMP_PROC_BIND_SPREAD
#else
!$ use omp_lib, only : omp_get_num_threads
!$ use omp_lib, only : omp_get_schedule, omp_set_schedule
!$ use omp_lib, only : omp_get_proc_bind
!$ use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
!$ use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
!$ use omp_lib, only : OMP_PROC_BIND_FALSE, OMP_PROC_BIND_TRUE
!$ use omp_lib, only : OMP_PROC_BIND_MASTER
!$ use omp_lib, only : OMP_PROC_BIND_CLOSE, OMP_PROC_BIND_SPREAD
#endif

  implicit none

  integer :: level
  real(dp) :: threshold
#ifdef MPI
  integer :: MPIerror
#endif

  integer :: iEl, itmp, it
  type(Sparsity) :: tmp_sp
  type(dSpData1D) :: tmp_1D
  type(dSpData2D) :: tmp_2D
  character(len=300) :: sname, fname

  ! Initialise MPI and set processor number
#ifdef MPI
#ifdef _OPENMP
  call MPI_Init_Thread(MPI_Thread_Funneled, it, MPIerror)
  if ( MPI_Thread_Funneled /= it ) then
     ! the requested threading level cannot be asserted
     ! Notify the user
     write(0,'(a)') '!!! Could not assert funneled threads'
  end if
#else
  call MPI_Init( MPIerror )
#endif
#endif
  
  call parallel_init()

  ! Initialize the output
  call init_output(Node == 0)

#ifdef MPI
  if (.not. fdf_parallel()) then
     call die('tbt_init: ERROR: FDF module doesn''t have parallel support')
  endif
#endif

! Print version information ...........................................
  if (IOnode) then
     
     call prversion()

#ifdef MPI
     if (Nodes > 1) then
        write(*,'(/,a,i0,tr1,a)') '* Running on ', Nodes, &
             'nodes in parallel'
     else
        write(*,'(/,a)') '* Running in serial mode with MPI'
     endif
#else
     write(*,'(/,a)') '* Running in serial mode'
#endif
!$OMP parallel default(shared)
!$OMP master
!$    it = omp_get_num_threads()
!$    write(*,'(a,i0,a)') '* Running ',it,' OpenMP threads.'
!$    write(*,'(a,i0,a)') '* Running ',Nodes*it,' processes.'
!$    it = omp_get_proc_bind()
!$    select case ( it )
!$    case ( OMP_PROC_BIND_FALSE ) 
!$    write(*,'(a)') '* OpenMP threads NOT bound (please bind threads!)'
!$    case ( OMP_PROC_BIND_TRUE ) 
!$    write(*,'(a)') '* OpenMP threads bound'
!$    case ( OMP_PROC_BIND_MASTER ) 
!$    write(*,'(a)') '* OpenMP threads bound (master)'
!$    case ( OMP_PROC_BIND_CLOSE ) 
!$    write(*,'(a)') '* OpenMP threads bound (close)'
!$    case ( OMP_PROC_BIND_SPREAD ) 
!$    write(*,'(a)') '* OpenMP threads bound (spread)'
!$    case default
!$    write(*,'(a)') '* OpenMP threads bound (unknown)'
!$    end select
     
!$    call omp_get_schedule(it,itmp)
!$    select case ( it )
!$    case ( OMP_SCHED_STATIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule STATIC, chunks ',itmp
!$    case ( OMP_SCHED_DYNAMIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule DYNAMIC, chunks ',itmp
!$    if ( itmp == 1 ) then
!$     ! this is the default scheduling, probably the user
!$     ! have not set the value, predefine it to 32
!$     itmp = 32
!$     write(*,'(a,i0)')'** OpenMP runtime schedule DYNAMIC, chunks ',itmp
!$    end if
!$    case ( OMP_SCHED_GUIDED ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule GUIDED, chunks ',itmp
!$    case ( OMP_SCHED_AUTO ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule AUTO, chunks ',itmp
!$    case default
!$    write(*,'(a,i0)') '* OpenMP runtime schedule UNKNOWN, chunks ',itmp
!$    end select
!$OMP end master
!$OMP end parallel
!$    call omp_set_schedule(it,itmp)
     call timestamp('Start of run')
     call wallclock('Start of run')
  endif

  ! Start timer .....................................................
  call timer('tbtrans', 0)
  call timer('tbtrans', 1)
  
  ! Initialise read .................................................
  call tbt_reinit( sname , slabel )

  ! Initialize save-options so that we can figure
  ! out the saving directory
  call init_save_options( )

  ! Set timer report file and threshold .............................
  threshold = fdf_get('timer_report_threshold', 0._dp)
  ! Get file name
  call name_save(1,1,fname,end='times')
  call timer_report( file=trim(fname), threshold=threshold )

  ! Set allocation report level .........................................
  ! variables level and threshold imported from module siesta_options
  level = fdf_get('alloc_report_level', 0)
  threshold = fdf_get('alloc_report_threshold', 0._dp)
  call name_save(1,1,fname,end='alloc')
  call alloc_report( level=level, file=trim(fname), &
       threshold=threshold, printNow=.false. )

#ifdef NCDF_4
  ! In case the user wants to utilize the ncdf library
  call ncdf_IOnode(IONode)
#endif

  ! Initialization now complete. Flush stdout.
  if ( IOnode ) call pxfflush( 6 )

  ! Initialize the HSfiles
  ! This will read in the HSfile and determine whether we should
  ! do interpolation due to bias not matching any TSHS files
  ! passed to the program.
  ! This will also read in the required information about the system
  call tbt_init_HSfile( )

  ! Read in generic options
  call read_tbt_generic(TSHS%na_u, TSHS%lasto)

  ! Read chemical potential
  call read_tbt_chem_pot( )

  ! Read electrodes
  call read_tbt_elec(TSHS%cell, TSHS%na_u, TSHS%xa, TSHS%lasto)

  ! Read k-points
  call setup_kpoint_grid( TSHS%cell )
  if ( sum(TSHS%nsc) == 3 .and. .not. Gamma ) then
     write(*,'(a)')'Please see flag: ForceAuxCell'
     call die('Transiesta calculation was a Gamma calculation &
          &while you request transmission k-points.')
  end if

  ! Read remaining options
  call read_tbt_after_Elec(TSHS%nspin, TSHS%cell, TSHS%na_u, TSHS%lasto, &
       TSHS%xa, TSHS%no_u, kscell, kdispl)

  call read_proj_options( save_DATA )

  ! Print options
  call print_tbt_options( TSHS%nspin )

  ! Print warnings
  call print_tbt_warnings( Gamma )

  ! Print information regarding the device
  if ( IONode ) then
     write(*,'(a)') 'Device information (full):'
     call print_type(TSHS%sp)
     write(*,*) ! newline
  end if
  
  if ( IONode ) write(*,'(a)') 'Electrode information:'

  ! We have the contour now, so we can create the GF files
  do iEl = 1 , N_Elec

     if ( IONode ) write(*,*) ! newline

     ! initialize the electrode for Green's function calculation
     call init_Electrode_HS(Elecs(iEl))

     if ( Elecs(iEl)%is_gamma ) then
       call do_Green(Elecs(iEl), &
           TSHS%cell, 1, (/(/0._dp, 0._dp, 0._dp/)/), (/1._dp/), &
           Elecs_xa_Eps, .false. )
     else
       call do_Green(Elecs(iEl), &
           TSHS%cell,nkpnt,kpoint,kweight, &
           Elecs_xa_Eps, .false. )
     end if
     
     ! clean-up
     call delete(Elecs(iEl))
     
  end do

  if ( IONode ) write(*,*) ! newline

  if ( stop_after_GS ) then
     if ( IONode ) then
        write(*,'(a)')'tbt: Stopping program per user request.'
        write(*,'(a)')'tbt: Done creating all GF files.'
     end if
#ifdef MPI
     call MPI_Barrier(MPI_Comm_World,iEl)
#endif

     call tbt_end()

  end if

  ! Initialize tbtrans regions
  ! We pass a copy of the sparsity pattern as the sparsity pattern
  ! returned is the sparsity pattern minus the buffer atoms!
  ! Hence, in order to change the sparsity patterns of the data
  ! we need to retain both!
  tmp_sp = TSHS%sp
  call tbt_init_regions(N_Elec,Elecs,TSHS%cell, &
       TSHS%na_u,TSHS%xa,TSHS%lasto, &
       TSHS%dit,tmp_sp, &
       product(TSHS%nsc),TSHS%isc_off)
  call tbt_init_kregions(r_aBuf,N_Elec,Elecs,TSHS%cell, &
       TSHS%dit,tmp_sp,TSHS%na_u,TSHS%xa,TSHS%lasto, &
       TSHS%nsc,TSHS%isc_off)

  if ( Node == 0 ) then
     itmp = nnzs(TSHS%sp) - nnzs(tmp_sp)
     write(*,'(/,a,i0,/)')'tbt: Reducing matrix (H, S) &
          &sparsity patterns by: ', itmp
  end if

  ! Change the data 
  call SpData_to_Sp(TSHS%S_1D,tmp_sp,tmp_1D)
  TSHS%S_1D = tmp_1D
  call delete(tmp_1D)
  call SpData_to_Sp(TSHS%H_2D,tmp_sp,tmp_2D)
  TSHS%H_2D = tmp_2D
  call delete(tmp_2D)
  TSHS%sp = tmp_sp
  call delete(tmp_sp)

  ! Create the device region sparsity pattern
  call tbt_region_options( TSHS%sp, save_DATA )

  call tbt_print_regions( TSHS%na_u, TSHS%lasto, N_Elec, Elecs)

  call tbt_print_kRegions( TSHS%cell )

#ifdef NCDF_4

  ! Initialize the projections here.
  call init_proj( TSHS%na_u , TSHS%lasto , r_aDev , r_oDev, save_DATA )
  call init_proj_T( N_Elec, Elecs , save_DATA )

  call proj_print( N_Elec, Elecs )

#endif

  if ( N_eigen /= 0 ) then
     ! if eigen-value calculation, reduce eigen-values calculated
     ! to a sensible number
     if ( N_eigen > 0 ) then
        itmp = N_eigen
     else
        itmp = huge(1)
     end if

     ! Reduce to the minimum size
     do iEl = 1 , N_Elec
        itmp = min(itmp,Elecs(iEl)%o_inD%n)
     end do
#ifdef NCDF_4
     if ( N_proj_ME > 0 ) then
        do it = 1 , N_proj_ME
           itmp = min(itmp,proj_ME(it)%mol%orb%n)
        end do
     end if
#endif
     if ( IONode ) then
        if ( N_eigen > 0 ) then
           if ( itmp /= N_eigen ) then
              write(*,'(/,a)')'tbt: *** Correcting number of T eigenvalues...'
           end if
        else
           write(*,'(/,a,i0)')'tbt: *** Maximizing number of T eigenvalues to ',itmp
        end if
     end if
     N_eigen = itmp
  end if
  ! Update N_eigen
  if ( N_eigen /= 0 ) then
     save_DATA = save_DATA // ('T-eig'.kv.N_eigen)
  end if

  ! Now we have the sparsity patterns in the correct sparsity
  ! information and we have deleted all un-needed information.

end subroutine tbt_init
                       
