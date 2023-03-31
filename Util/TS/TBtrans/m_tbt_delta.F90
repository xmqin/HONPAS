! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2017
!
! Please attribute the original author in case of dublication.

! Enable the feature of adding \delta terms to the Hamiltonian.
! This implementation is the generic implementation of the
! dH and dSE routines.

! There are 4 levels of usage:
!   1. constant delta
!      This will enter in a common equation
!   2. k-dependent delta
!   3. energy-dependent delta
!   4. k and E dependent delta

! Note that the highest level has precedence above the others,
! so specifying both 1. and 4. will only be the equivalent of
! using level 4.

! This is particularly handy for TB calculations
! but can prove just as useful for regular DFT
! calculations as one can use an already existing
! SIESTA.nc file, and then do several "case" studies
! on such a file, thus actually having FULL control
! over EVERYTHING.

! Regardless of the delta (dH or dSE) this module also
! implements the additive routines to the trimat classes.

module m_tbt_delta

  use precision, only : dp

  use class_zSpData1D
  use m_tbt_save, only : tNodeE

  implicit none
 
  private

  type :: tDelta

     ! The file-name
     character(len=256) :: fname
     
     ! Designator of the current level of the 
     ! quantities in the dH designator
     integer :: lvl = -1

     !> The spin component that should be read in
     integer :: ispin = 1
     
     ! As this is a reciprocal cell, k-point,
     ! we initialize it immediately
     ! This should, for level 2 calculations 
     ! speed things up as we do not need to read the 
     ! sparse matrices always.
     real(dp) :: bkpt(3) = 2.12345_dp

     ! To ease the handling of different elements
     ! we only keep the complex quantity
     ! YES, we could limit this to only the real
     ! part for purely real designators.
     ! However, I do not suspect that users
     ! will use this to change extreme amounts
     ! of elements. (200.000 elements ~450 dense matrix
     ! take up 3 MB)
     type(zSpData1D) :: d

     ! Information about the contained data
     !   lvls(<>) == 0, the level does not exist
     !   lvls(<>) == 1, the level exists and is REAL
     !   lvls(<>) == 2, the level exists and is COMPLEX
     integer :: lvls(4) = 0

     ! Information for each of the levels
     real(dp), pointer :: bkpt2(:,:) => null()
     real(dp), pointer :: E3(:) => null()
     real(dp), pointer :: bkpt4(:,:) => null()
     real(dp), pointer :: E4(:) => null()

  end type tDelta

  ! Logical to determine whether the
  ! file should be read independently of the nodes
  ! or whether, there should be an IO-node.
  logical, save :: cdf_r_parallel = .true.

  public :: tDelta

  public :: init_delta_options

  public :: delta_has_level
#ifdef NCDF_4
  public :: read_delta_Sp
  public :: read_delta_next, clean_delta, delete_delta
#endif
  public :: add_zdelta_TriMat, add_zdelta_Mat

contains

  subroutine init_delta_options( opt, delta)

    use parallel, only : Node, Nodes
    use fdf
    use m_os, only : file_exist

    use dictionary

#ifdef NCDF_4
    use netcdf_ncdf, ncdf_parallel => parallel
#endif
    
#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Comm_World
    use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
#endif

    ! Option name used (dH or dSE)
    character(len=*), intent(in) :: opt
    type(tDelta), intent(inout) :: delta
    
#ifdef NCDF_4
    type(hNCDF) :: ndelta, grp
#endif
    logical :: exists, is_real

    integer :: n_k, n_E
#ifdef MPI
    integer :: MPIerror
#endif

    ! Set all levels to non-existing
    delta%lvls(:) = 0

#ifdef NCDF_4

    ! just set a file-name that should never be created by any user :).
    delta%fname = fdf_get('TBT.'//opt,'NONE 1234567890')

    ! If the file exists, use it
    if ( .not. file_exist(delta%fname, Bcast = .true.) ) then
       
       delta%fname = ' '
       return
       
    end if

    ! Ok, the file exists, lets see if all can see the file...
    if ( file_exist(delta%fname, all = .true.) ) then
       
       ! We can check whether we should read in parallel
       ! We default it to read parallelly
       cdf_r_parallel = fdf_get('TBT.'//opt//'.Parallel',.true.) .and. cdf_r_parallel

    else

       cdf_r_parallel = .false.

    end if

    ! The user cannot decide if only one core
    if ( Nodes == 1 ) cdf_r_parallel = .true.

    ! Read in options
    if ( cdf_r_parallel ) then
       call ncdf_open(ndelta,delta%fname, mode = NF90_SHARE , parallel = .true. )
    else
       call ncdf_open(ndelta,delta%fname, mode = NF90_NOWRITE )
    end if

    ! Check if level 4 exists
    exists = .false. ! required for cdf_r_parallel == .false.
    call ncdf_inq_grp(ndelta,'LEVEL-4', exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndelta,'LEVEL-4',grp)

       call ncdf_inq_var(grp,'delta', exist=is_real)
       
       call ncdf_inq_dim(grp,'ne', len = n_E)
       call ncdf_inq_dim(grp,'nkpt', len = n_k)
       
       if ( n_E == 0 .or. n_k == 0 ) then

          ! do nothing, we just skip dH for level 4
          if ( Node == 0 ) then
             write(*,'(4a)')'tbt-',opt, ': Level 4 exists, but has no E or k points. ', &
                  'It has been discarded.'
          end if

       else

          ! Write out that we are using level 4
          if ( Node == 0 ) then
             write(*,'(3a,2(i0,a))')'tbt-',opt, ': Level 4 will be used, and has ', &
                  n_k,' k-points and ',n_E,' energy points.'
          end if

          if ( is_real ) then
             delta%lvls(4) = 1
          else
             delta%lvls(4) = 2
          end if

          ! Allocate and read
          allocate(delta%bkpt4(3,n_k))
          allocate(delta%E4(n_E))
          call ncdf_get_var(grp, 'kpt', delta%bkpt4)
          call ncdf_get_var(grp, 'E', delta%E4)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) &
         call MPI_Bcast(delta%lvls(4), 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    if ( delta%lvls(4) > 0 .and. .not. cdf_r_parallel ) then
       call MPI_Bcast(n_k, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
       call MPI_Bcast(n_E, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
       if ( Node /= 0 ) then
          allocate(delta%bkpt4(3,n_k))
          allocate(delta%E4(n_E))
       end if
       call MPI_Bcast(delta%bkpt4(1,1),3*n_k,MPI_Double_Precision, 0, &
            MPI_Comm_World, MPIerror)
       call MPI_Bcast(delta%E4(1),n_E,MPI_Double_Precision, 0, &
            MPI_Comm_World, MPIerror)
    end if
#endif

    
    ! Check if level 3 exists
    exists = .false.
    call ncdf_inq_grp(ndelta,'LEVEL-3',exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndelta,'LEVEL-3',grp)

       call ncdf_inq_var(grp, 'delta', exist=is_real)

       ! Read in quantities for level 3
       call ncdf_inq_dim(grp, 'ne' , len = n_E)
       
       if ( n_E == 0 ) then

          ! do nothing, we just skip dH for level 3
          if ( Node == 0 ) then
             write(*,'(3a)')'tbt-', opt, ': Level 3 exists, but has no E points. &
                  &It has been discarded.'
          end if

       else

          ! Write out that we are using level 3
          if ( Node == 0 ) then
             write(*,'(3a,i0,a)')'tbt-',opt, ': Level 3 will be used, and has ', &
                  n_E,' energy points.'
          end if

          if ( is_real ) then
             delta%lvls(3) = 1
          else
             delta%lvls(3) = 2
          end if

          ! Allocate and read
          allocate(delta%E3(n_E))
          call ncdf_get_var(grp, 'E', delta%E3)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) &
         call MPI_Bcast(delta%lvls(3), 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    if ( delta%lvls(3) > 0 .and. .not. cdf_r_parallel ) then
       call MPI_Bcast(n_E, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
       if ( Node /= 0 ) then
          allocate(delta%E3(n_E))
       end if
       call MPI_Bcast(delta%E3(1),n_E,MPI_Double_Precision, 0, &
            MPI_Comm_World, MPIerror)
    end if
#endif

    
    ! Check if level 2 exists
    exists = .false.
    call ncdf_inq_grp(ndelta, 'LEVEL-2', exist = exists)
    if ( exists ) then

       call ncdf_open_grp(ndelta, 'LEVEL-2', grp)

       call ncdf_inq_var(grp, 'delta', exist=is_real)

       ! Read in quantities for level 2
       call ncdf_inq_dim(grp, 'nkpt', len = n_k)
       
       if ( n_k == 0 ) then

          ! do nothing, we just skip dH for level 2
          if ( Node == 0 ) then
             write(*,'(3a)')'tbt-', opt, ': Level 2 exists, but has no k-points. &
                  &It has been discarded.'
          end if

       else

          ! Write out that we are using level 2
          if ( Node == 0 ) then
             write(*,'(3a,i0,a)')'tbt-', opt, ': Level 2 will be used, and has ', &
                  n_k,' k-points.'
          end if

          if ( is_real ) then
             delta%lvls(2) = 1
          else
             delta%lvls(2) = 2
          end if

          ! Allocate and read
          allocate(delta%bkpt2(3,n_k))
          call ncdf_get_var(grp, 'kpt', delta%bkpt2)

       end if
       
    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) &
         call MPI_Bcast(delta%lvls(2), 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    if ( delta%lvls(2) > 0 .and. .not. cdf_r_parallel ) then
       call MPI_Bcast(n_k, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
       if ( Node /= 0 ) then
          allocate(delta%bkpt2(3,n_k))
       end if
       call MPI_Bcast(delta%bkpt2(1,1),3*n_k,MPI_Double_Precision, 0, &
            MPI_Comm_World, MPIerror)
    end if
#endif

    
    ! Check if level 1 exists
    exists = .false.
    call ncdf_inq_grp(ndelta, 'LEVEL-1', exist = exists)
    if ( exists ) then
       
       call ncdf_open_grp(ndelta, 'LEVEL-1', grp)

       call ncdf_inq_var(grp, 'delta', exist=is_real)

       ! Write out that we are using level 1
       if ( Node == 0 ) then
          write(*,'(3a)')'tbt-', opt, ': Level 1 will be used.'
       end if

       if ( is_real ) then
          delta%lvls(1) = 1
       else
          delta%lvls(1) = 2
       end if

    end if

#ifdef MPI
    if ( .not. cdf_r_parallel ) &
         call MPI_Bcast(delta%lvls(1), 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
#endif

    call ncdf_close(ndelta)

#endif

  end subroutine init_delta_options


  function delta_has_level(delta, lvl) result(has)
    type(tDelta), intent(in) :: delta
    integer, intent(in) :: lvl
    logical :: has
    has = delta%lvls(lvl) > 0
  end function delta_has_level

#ifdef NCDF_4

  subroutine read_delta_Sp(delta, no_u, sp)

    use class_Sparsity
    use class_OrbitalDistribution

    use m_sparsity_handling, only : Sp_union
    use netcdf_ncdf, ncdf_parallel => parallel

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use m_ncdf_io, only : cdf_r_Sp

    type(tDelta), intent(in) :: delta

    ! Read in the sparsity pattern
    integer, intent(in) :: no_u
    type(Sparsity), intent(inout) :: sp

    ! Temporary sparsity pattern
    type(Sparsity) :: sp_tmp
    type(OrbitalDistribution) :: fdist

    type(hNCDF) :: ndelta, grp

    character(len=7) :: igrp
    integer :: i 

    call delete(sp)

    ! Quick escape if there are no levels to read.
    if ( all(delta%lvls == 0) ) return

#ifdef MPI
    call newDistribution(no_u,MPI_Comm_Self,fdist,name='TBT-fake dist')
#else
    call newDistribution(no_u,-1           ,fdist,name='TBT-fake dist')
#endif

    if ( cdf_r_parallel ) then
       call ncdf_open(ndelta, delta%fname, mode = IOR(NF90_SHARE,NF90_NOWRITE) , &
            parallel = .true. )
    else
       call ncdf_open(ndelta, delta%fname, mode = NF90_NOWRITE )
    end if

    do i = 1 , 4
       if ( .not. delta_has_level(delta, i) ) cycle
       
       write(igrp,'(a,i0)') 'LEVEL-',i

       call ncdf_open_grp(ndelta,igrp,grp)
       call cdf_r_Sp(grp,no_u,sp_tmp, 'sp-delta', Bcast = .not. cdf_r_parallel )
       call Sp_union(fdist,sp_tmp,sp,sp)

    end do

    call delete(fdist)
    call delete(sp_tmp)

    call ncdf_close(ndelta)
    
  end subroutine read_delta_Sp

  ! Updates the d type to the current energy-point
  subroutine read_delta_next(opt, delta, no_u, bkpt, nE)

    use parallel, only : IONode, Node, Nodes

#ifdef MPI
    use mpi_siesta, only : MPI_Gather
    use mpi_siesta, only : MPI_Comm_World, MPI_Integer
#endif

    use units, only: eV
    use m_verbosity, only: verbosity

    ! Input variables
    character(len=*), intent(in) :: opt
    type(tDelta), intent(inout) :: delta
    integer, intent(in) :: no_u
    real(dp), intent(in) :: bkpt(3)
    type(tNodeE), intent(in) :: nE

    integer :: ik, iE, iN
    integer, allocatable :: nlvl(:)

#ifdef MPI
    integer :: MPIerror
#endif

    character(len=*), parameter :: f1 = '(a,i0,a,tr1,a)'
    character(len=*), parameter :: f2 = '(a,i0,a,tr1,2a,"[",2(" ",f7.4,", "),f7.4,"]")'
    character(len=*), parameter :: f3 = '(a,i0,a,tr1,2a,f8.4,tr1,a)'
    character(len=*), parameter :: f4 = '(a,i0,a,tr1,2a,"[",2(" ",f7.4,", "),f7.4,"]",a,f8.4,tr1,a)'

#ifdef TBTRANS_TIMING
    call timer('read-delta',1)
#endif

    ! Everybody figures out which level they correspond to
    ik = 0
    iE = 0
    allocate(nlvl(0:Nodes-1))
    nlvl(Node) = 0
    if ( delta_has_level(delta, 4) ) then
       ik = idx_k(bkpt,delta%bkpt4)
       if ( ik > 0 ) then
          iE = idx_E(nE%E(Node),delta%E4)
          if ( iE == 0 ) ik = 0
       end if
       if ( ik + iE /= 0 ) nlvl(Node) = 4
    end if
    if ( delta_has_level(delta, 3) .and. nlvl(Node) == 0 ) then
       iE = idx_E(nE%E(Node),delta%E3)
       if ( iE /= 0 ) nlvl(Node) = 3
    end if
    if ( delta_has_level(delta, 2) .and. nlvl(Node) == 0 ) then
       ik = idx_k(bkpt,delta%bkpt2)
       if ( ik /= 0 ) nlvl(Node) = 2
    end if
    if ( delta_has_level(delta, 1) .and. nlvl(Node) == 0 ) then
       nlvl(Node) = 1
    end if
    
    !print *,Node,nlvl(node),ik,bkpt,iE,nE%E(Node) * 13.60580_dp

    ! In case there is only one IO node
    if ( .not. cdf_r_parallel ) then
#ifdef MPI
       ! Gather levels on IO
       call MPI_Gather(nlvl(Node),1,MPI_Integer, &
            nlvl(0),1,MPI_Integer,0,MPI_Comm_World, MPIerror)
#endif

       ! We can only all read the same level
       ! This is because a change in sparsity pattern
       ! might prohibit the Bcast mechanism for the sparsity
       ! patterns.
       if ( IONode .and. .not. all(nlvl == nlvl(0)) ) then
          write(*,'(a,1000(tr1,i2))')'Node levels: ',nlvl
          write(*,'(3a)')'Error in using ', opt, ' functionality'
          write(*,'(a)')'When using non-parallel reading of a delta file you must &
               &ensure that at each iteration each core will use the same level.'
          write(*,'(a)')'For (easy) full functionality please see if you can place the delta-file &
               &so that all MPI-cores can see it.'
          call die('Differing level designation and non-MPI IO, please see output...')
       end if

    end if

    if ( nlvl(Node) /= delta%lvl ) then

       ! Ensure it is emptied if the level is changed.
       call clean_delta( delta )

    end if

    select case ( nlvl(Node) )
    case ( 1 )

       if ( delta%lvl /= 1 ) then

          ! Read the delta term
          call sub_read_delta(delta,1,delta%lvls(1) == 1,0,0)
          if ( verbosity > 7 ) then
             write(*,f1) 'Level 1 (',Node,')',opt
          end if

       end if

    case ( 2 )
       
       if ( sum(abs(delta%bkpt(:) - bkpt(:))) > 0.00001_dp ) then
          
          ! The k-point has changed, read the new delta term
          delta%bkpt(:) = bkpt(:)
          call sub_read_delta(delta,2,delta%lvls(2) == 1,ik,0)

          if ( verbosity > 7 ) then
             write(*,f2) 'Level 2 (',Node,')',opt, ', kpt = ',bkpt
          end if

       end if

    case ( 3 )
       
       ! We should not have two same energy-points
       ! consecutively
       call sub_read_delta(delta,3,delta%lvls(3) == 1,0,iE)

       if ( verbosity > 7 ) then
          write(*,f3) 'Level 3 (',Node,')',opt,', E = ', nE%E(Node) / eV, 'eV'
       end if
       
    case ( 4 )

       call sub_read_delta(delta,4,delta%lvls(4) == 1,ik,iE)

       if ( verbosity > 7 ) then
          write(*,f4) 'Level 4 (',Node,')',opt,', kpt = ',bkpt,', E = ',nE%E(Node) / eV, 'eV'
       end if

    end select

    ! Inform the delta type to the current level
    delta%lvl = nlvl(Node)

    deallocate(nlvl)

#ifdef TBTRANS_TIMING
    call timer('read-delta',2)
#endif

  contains
    
    subroutine sub_read_delta(delta, lvl, is_real, ik, iE)

      use parallel, only : Node

      use class_Sparsity
      use class_OrbitalDistribution
      use netcdf_ncdf, ncdf_parallel => parallel
      use m_ncdf_io, only : cdf_r_Sp

#ifdef MPI
      use mpi_siesta, only : MPI_Send, MPI_Recv
      use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self, MPI_Status_Size
      use mpi_siesta, only : MPI_Integer
      use mpi_siesta, only : MPI_Double_Precision, MPI_Double_Complex
#endif

      type(tDelta), intent(inout) :: delta
      integer, intent(in) :: lvl
      logical, intent(in) :: is_real
      integer, intent(in) :: ik, iE

      real(dp), allocatable :: rM(:)
      complex(dp), pointer :: zM(:)

      type(hNCDF) :: grp
      character(len=7) :: igrp

      type(Sparsity) :: sp
      type(OrbitalDistribution) :: fdist
      
      integer :: nnz, oE
      integer :: start(4)

#ifdef MPI
      ! We must figure out which level each node
      ! lives on
      integer :: MPIerror, status(Mpi_status_size)
#endif

      ! At this point we know which levels we should read
      ! on each node
      write(igrp,'(a,i0)') 'LEVEL-',lvl
      
      !print *,Node,lvl,ik,iE
      
      if ( cdf_r_parallel ) then
         call ncdf_open(grp,delta%fname, mode = IOR(NF90_SHARE,NF90_NOWRITE) , &
              group = igrp , parallel = .true. )
      else
         call ncdf_open(grp,delta%fname, mode = NF90_NOWRITE , group = igrp )
      end if
      
      ! First read in sparsity pattern (the user can have
      ! different sparsity patterns for each level)
      if ( delta%lvl /= lvl ) then
         
         call cdf_r_Sp(grp,no_u,sp, 'sp-delta', Bcast = .not. cdf_r_parallel )
         
#ifdef MPI
         call newDistribution(no_u,MPI_Comm_Self,fdist,name='TBT-fake dist')
#else
         call newDistribution(no_u,-1           ,fdist,name='TBT-fake dist')
#endif
          
         ! Create the data container
         call newzSpData1D(sp,fdist,delta%d,name='delta')
         
         call delete(sp)
         call delete(fdist)
         
      end if
      
      zM  => val(delta%d)
      nnz =  size(zM)

      ! If the dH file is not read by NF90_SHARE
      ! We must have the IO-node to read and distribute
      if ( is_real ) then
         allocate(rM(nnz))
      end if
      
      start(:)    =  1
      start(2)    =  delta%ispin
      select case ( lvl )
      case ( 2 )
         start(3) = ik
      case ( 3 )
         start(3) = iE
      case ( 4 )
         start(3) = iE
         start(4) = ik
      end select
         
#ifdef MPI
      if ( .not. cdf_r_parallel ) then
      if ( lvl == 1 .or. lvl == 2 ) then
         if ( is_real ) then
            call ncdf_get_var(grp,'delta',rM, start = start )
            call MPI_Bcast(rM,nnz,MPI_Double_Precision, 0, &
                MPI_Comm_World, MPIerror)
            zM(:) = rM(:)
            deallocate(rM)
         else
            call ncdf_get_var(grp,'delta',zM, start = start )
            call MPI_Bcast(zM,nnz,MPI_Double_Complex, 0, &
                 MPI_Comm_World, MPIerror)
         end if

         return

      else 
      if ( Node == 0 ) then
         do iN = 1 , Node - 1
            ! Retrieve the energy point (the k-point IS the same)
            call MPI_Recv(oE,1,MPI_Integer, iN, iN, &
                 MPI_Comm_World, status, MPIerror)
            start(3) = oE
            if ( is_real ) then
               call ncdf_get_var(grp,'delta',rM, start = start )
               call MPI_Send(rM,nnz,MPI_Double_Precision, iN, iN, &
                    MPI_Comm_World, MPIerror)
            else
               call ncdf_get_var(grp,'delta',zM, start = start )
               call MPI_Send(zM,nnz,MPI_Double_Complex, iN, iN, &
                    MPI_Comm_World, MPIerror)
            end if
         end do
         ! set back the original starting place
         start(3) = iE
      else
         call MPI_Send(iE,1,MPI_Integer, 0, Node, &
              MPI_Comm_World, MPIerror)
         if ( is_real ) then
            call MPI_Recv(rM,nnz,MPI_Double_Precision, 0, Node, &
                 MPI_Comm_World, status, MPIerror)
         else
            call MPI_Recv(zM,nnz,MPI_Double_Complex, 0, Node, &
                 MPI_Comm_World, status, MPIerror)
         end if
      end if
      end if
      end if
#endif
      
      if ( is_real ) then
         call ncdf_get_var(grp,'delta',rM, start = start )
         zM(:) = rM(:)
         deallocate(rM)
      else
         call ncdf_get_var(grp,'delta',zM, start = start )
      end if

      call ncdf_close(grp)

    end subroutine sub_read_delta
    
    function idx_k(bk,fbk) result(i)
      real(dp), intent(in) :: bk(3), fbk(:,:)
      integer :: i
      do i = 1 , size(fbk,dim=2)
         if ( abs( fbk(1,i) - bk(1) ) + &
              abs( fbk(2,i) - bk(2) ) + &
              abs( fbk(3,i) - bk(3) ) < 0.0001_dp ) then
            return
         end if
      end do
      i = 0
    end function idx_k

    function idx_E(E,fE) result(i)
      real(dp), intent(in) :: E, fE(:)
      integer :: i
      do i = 1 , size(fE)
         if ( abs(fE(i) - E) < 7.349806700083788e-06_dp ) then
            return
         end if
      end do
      i = 0
    end function idx_E

    subroutine correct_idx(ik,iE)
      integer, intent(inout), optional :: ik, iE
      if ( present(iE) .and. present(ik) ) then
         if ( iE == 0 ) then
            ik = 0
         else if ( ik == 0 ) then
            iE = 0
         end if
      end if
    end subroutine correct_idx

  end subroutine read_delta_next

  subroutine clean_delta( delta )
    type(tDelta), intent(inout) :: delta

    delta%lvl = -1
    delta%bkpt = 2.12425_dp
    call delete(delta%d)

  end subroutine clean_delta

  subroutine delete_delta( delta )
    type(tDelta), intent(inout) :: delta

    call clean_delta(delta)

    delta%ispin = 1
    delta%lvls(:) = 0
    if ( associated(delta%bkpt2) ) deallocate(delta%bkpt2)
    if ( associated(delta%E3) ) deallocate(delta%E3)
    if ( associated(delta%bkpt4) ) deallocate(delta%bkpt4)
    if ( associated(delta%E4) ) deallocate(delta%E4)
    
    nullify(delta%bkpt2, delta%E3, delta%bkpt4, delta%E4)

  end subroutine delete_delta

#endif

  ! Add the delta to the tri-diagonal matrix
  subroutine add_zdelta_TriMat( zd, GFinv_tri, r, pvt, sc_off, k)

    use class_zTriMat
    use class_Sparsity
    use m_region

    use intrinsic_missing, only : MODP

    type(zSpData1D), intent(inout) :: zd
    type(zTriMat), intent(inout) :: GFinv_tri
    type(tRgn), intent(in) :: r, pvt
    ! Super-cell offset and k-point
    real(dp), intent(in) :: sc_off(:,:), k(3)

    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), allocatable :: ph(:)
    complex(dp), pointer :: d(:), GFinv(:)

    integer :: idx, iu, ju, ind, jo, no

    sp => spar(zd)
    d => val(zd)

    call attach(sp, nrows_g=no, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    allocate( ph(0:size(sc_off,dim=2)-1) )
    do iu = 1 , size(sc_off, dim=2)
       ph(iu-1) = cdexp(dcmplx(0._dp, &
            k(1) * sc_off(1,iu) + &
            k(2) * sc_off(2,iu) + &
            k(3) * sc_off(3,iu)))
    end do

    Gfinv => val(Gfinv_tri)

!$OMP parallel do default(shared), private(iu,jo,ind,ju,idx)
    do ju = 1, r%n
       jo = r%r(ju) ! get the orbital in the big sparsity pattern
       if ( l_ncol(jo) /= 0 ) then
          
          ! Loop on entries here...
          do ind = l_ptr(jo) + 1 , l_ptr(jo) + l_ncol(jo)
             ! Look up in the pivoting array what the pivoted orbital is
             iu = pvt%r(MODP(l_col(ind), no))
             ! Check whether this element should be added
             if ( iu == 0 ) cycle
             
             idx = index(Gfinv_tri,ju,iu)
             
             GFinv(idx) = GFinv(idx) - d(ind) * ph( (l_col(ind)-1)/no )
          end do
          
       end if
    end do
!$OMP end parallel do

    deallocate(ph)
    
  end subroutine add_zdelta_TriMat
  
  ! Add the delta to the tri-diagonal matrix
  subroutine add_zdelta_Mat( zd , r, off1, n1, off2, n2, M, &
       sc_off, k)

    use class_Sparsity
    use m_region
    use intrinsic_missing, only : MODP

    type(zSpData1D), intent(inout) :: zd
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! The sizes and offsets of the matrix
    integer, intent(in) :: off1, n1, off2, n2
    complex(dp), intent(inout) :: M(n1,n2)
    ! Super-cell offset and k-point
    real(dp), intent(in) :: sc_off(:,:), k(3)

    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    complex(dp), allocatable :: ph(:)
    complex(dp), pointer :: d(:)

    integer :: iu, ju, ind, jo, no

    sp => spar(zd)
    d => val(zd)

    call attach(sp, nrows_g=no, &
         n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    ! Create the phases
    allocate( ph(0:size(sc_off,dim=2)-1) )
    do iu = 1 , size(sc_off, dim=2)
       ph(iu-1) = cdexp(dcmplx(0._dp, &
            k(1) * sc_off(1,iu) + &
            k(2) * sc_off(2,iu) + &
            k(3) * sc_off(3,iu)))
    end do
    
!$OMP parallel do default(shared), private(iu,jo,ind,ju)
    do ju = 1 , n1
       jo = r%r(off1+ju) ! get the orbital in the sparsity pattern
       
       if ( l_ncol(jo) /= 0 ) then
          
          do ind = l_ptr(jo) + 1 , l_ptr(jo) + l_ncol(jo)
             iu = rgn_pivot(r, MODP(l_col(ind), no)) - off2
             if ( iu < 1 .or. n2 < iu ) cycle
             
             M(ju,iu) = M(ju,iu) - d(ind) * ph( (l_col(ind)-1)/no )
          end do
          
       end if

    end do
!$OMP end parallel do

    deallocate(ph)
    
  end subroutine add_zdelta_Mat
  
end module m_tbt_delta
