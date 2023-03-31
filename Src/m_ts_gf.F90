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
module m_ts_GF
!
! Routines that are used for reading the Green function files.
! It is very closely related to m_ts_electrode (consider moving it to there)
! 
! This module constitutes the routines that are needed to ensure the
! correct format of the GF files. 
! In order to limit the size of the m_ts_electrode file this has been created.
! We also read in the next energy-points using routines in this module.
!
! A call to check_green is used to ensure the correct format of the GF file.
! It checks as much information as is available and dies if non-conforming
! entities exists. This routine should only be called by IONode!
!
! A call to read_Green will read in the header of the GF file and do "basic"
! checks against array sizes. Thus a call to check_green is advised before 
! read_Green!
! This routine will also distribute the arrays in an MPI run.

  use precision, only : dp

  implicit none

  public :: do_Green
  public :: read_Green, check_Green
  public :: reread_Gamma_Green
  public :: read_next_GS

  private

contains

  ! This method should only be called from transiesta (not tbtrans)
  subroutine do_Green(El, &
       ucell,nkpnt,kpoint,kweight, &
       xa_EPS, CalcDOS )
    
    use parallel  , only : IONode
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta
#endif
    use m_os, only : file_exist
    use m_ts_cctype
    use m_ts_electype
    use m_ts_electrode, only : create_Green

    use m_ts_contour_eq
    use m_ts_contour_neq

    implicit none
    
    ! ***********************
    ! * INPUT variables     *
    ! ***********************
    type(Elec), intent(inout) :: El
    integer, intent(in) :: nkpnt ! Number of k-points
    real(dp), intent(in) :: kpoint(3,nkpnt) ! k-points
    real(dp), intent(in) :: kweight(nkpnt) ! weights of kpoints
    real(dp), intent(in) :: xa_Eps ! coordinate precision check
    real(dp), dimension(3,3) :: ucell ! The unit cell of the CONTACT
    logical, intent(in) :: CalcDOS

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
    exist = file_exist(El%GFfile, Bcast = .true. )

    cReUseGF = El%ReUseGf
    ! If it does not find the file, calculate the GF
    if ( exist ) then
      if (IONode ) then
        write(*,*) 'Electrode Green function file: '//&
            trim(El%GFfile)//' already exist.'
        if ( .not. cReUseGF ) then
          write(*,*)'Green function file '//&
              trim(El%GFfile)//' is requested overwritten.'
        end if
      end if
    else
      cReUseGF = .false.
    end if

    errorGF = .false.

    ! we need to create all the contours
    NEn = N_Eq_E() + N_nEq_E()
    allocate(ce(NEn))
    iE = 0
    do i = 1 , N_Eq_E()
      c = Eq_E(i)
      ce(i) = c%e
    end do
    iE = N_Eq_E()
    if ( El%Eta > 0._dp ) then
      do i = 1 , N_nEq_E()
        c = nEq_E(i)
        ! We utilize the eta value for the electrode
        ce(iE+i) = cmplx(real(c%e,dp),El%Eta, dp)
      end do
    else
      ! Specified to use the device eta (ensure it is larger than 0)
      do i = 1 , N_nEq_E()
        c = nEq_E(i)
        ce(iE+i) = c%e
      end do
    end if

    ! We return if we should not calculate it
    if ( cReUseGF ) then
      
      ! Check that the Green functions are correct!
      ! This is needed as create_Green returns if user requests not to
      ! overwrite an already existing file.
      ! This check will read in the number of orbitals and atoms in the
      ! electrode surface Green function.
      ! Check the GF file
      if ( IONode ) then
        call io_assign(uGF)
        open(file=El%GFfile,unit=uGF,form='UNFORMATTED')
        
        call check_Green(uGF,El, &
            ucell, nkpnt, kpoint, kweight, NEn, ce, &
            xa_Eps, errorGF)
        
        write(*,'(/,4a,/)') "Using GF-file '",trim(El%GFfile),"'"
        
        call io_close(uGF)
      end if

#ifdef MPI
      call MPI_Bcast(errorGF,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif
    else
      
      call create_Green(El, ucell, nkpnt, kpoint, kweight, NEn, ce)

    end if

    ! Check the error in the GF file
    if ( errorGF ) &
        call die("Error in GFfile: "//trim(El%GFfile)//". Please move or delete")

    deallocate(ce)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS do_Green' )
#endif

  end subroutine do_Green


! ##################################################################
! ## Subroutine which read-in GS for the 1x1 surface cell         ##
! ##    The expansion of the arrays are performed else-where.     ##
! ##                                                              ##
! ##    Nick Papior Andersen, nickpapior@gmail.com                ##
! ## Fully recoded to conform with the memory reduced TranSIESTA  ##
! ##################################################################

  subroutine read_next_GS_Elec(uGF,NEReqs,ikpt,El,cE, &
       nwork,work, &
       forward)

    use parallel, only : IONode, Node, Nodes
    use units,    only : eV

    use m_ts_electype
    use m_ts_cctype

#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Send, MPI_Recv
    use mpi_siesta, only : MPI_Sum, MPI_Integer, MPI_Double_Complex
    use mpi_siesta, only : MPI_Status_Size, MPI_Comm_World
#endif

! *********************
! * INPUT variables   *
! *********************
    ! file-unit, and k-point index
    integer, intent(in) :: uGF, ikpt
    integer, intent(in) :: NEReqs
    ! The electrode also contains the arrays
    type(Elec), intent(in out) :: El
    type(ts_c_idx), intent(in)     :: cE
    ! The work array passed, this means we do not have
    ! to allocate anything down here.
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    ! If 'forward' is .true. we will read consecutively 
    ! and distribute from Node = 0, to Node = Nodes - 1
    ! (default)
    ! 
    ! if 'forward' is .false. we will read consecutively 
    ! and distribute from Node = Nodes - 1 to Node = 0
    logical, intent(in), optional :: forward

! *********************
! * LOCAL variables   *
! *********************
    real(dp), parameter :: EPS = 1.d-6
    integer :: read_Size, read_Size_HS

#ifdef MPI
    integer :: MPIerror, Status(MPI_Status_Size)
#endif

    complex(dp) :: ZE_cur
    integer :: iNode, ikGS, iEni
    integer :: iNodeS, iNodeE, iNodeStep
    logical :: lforward

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_next_GS' )
#endif

    lforward = .true.
    if ( present(forward) ) lforward = forward
    if ( lforward ) then
       iNodeS = 0
       iNodeE = NEReqs - 1
       iNodeStep = 1
    else
       iNodeS = Nodes - 1
       iNodeE = Nodes - NEReqs
       iNodeStep = -1
    end if

    ! Initial size of minimal size with expansion
    read_Size_HS = El%no_used ** 2 * El%Bloch%size() ! no_GS * no_GS * nq
    select case ( El%pre_expand )
    case ( 0 )
       ! Nothing is pre-expanded
       read_Size = read_Size_HS
    case ( 1 )
       ! Only the Green function is pre-expanded
       read_Size = read_Size_HS * El%Bloch%size()
    case ( 2 )
       ! Everything is pre-expanded
       read_Size = read_Size_HS * El%Bloch%size()
       read_Size_HS = read_Size
    end select

    ! Check if the number of energy points requested are 
    ! inconsistent
    if ( NEReqs <= 0 ) then
       if(IONode) &
            write(*,'(a,i0,a)') &
            'ERROR read_next_GS_Elec: Requested E-points=', &
            NEReqs,'< 0'
       call die('ERROR in reading GF file')
    else if ( Nodes < NEReqs ) then
       if(IONode) &
            write(*,'(2(a,i0))') &
            'ERROR read_next_GS_Elec: Requested E-points= ', &
            NEReqs,' > Nodes = ', Nodes
       call die('ERROR in reading GF file')
    end if

    if ( nwork < read_Size ) then
       write(*,*) 'Size of work array while reading GS was not large &
            &enough. Something went wrong.'
       call die('ERROR in reading GF file')
    end if
    
! Loop over nodes. If root send stuff out to relevant nodes, if
! not root receive information from root.

    ! We only loop over the requested energy points...
    do iNode = iNodeS, iNodeE, iNodeStep

       if ( IONode ) then
         ! read in header of GF-point
         read(uGF) ikGS, iEni, ZE_cur
         ! For gamma-only electrodes the read k-point is always correct
         ! I.e. there is only one!
         if ( El%is_gamma ) ikGS = ikpt
       end if

#ifdef MPI
       call MPI_Bcast(iEni,1,MPI_Integer, &
            0,MPI_Comm_World,MPIerror)

       ! distribute the current energy point
       if ( IONode .and. Node == iNode ) then
          ! do nothing
       else if ( Node == iNode ) then
          ! recieve from the host node
          call MPI_Recv(ZE_cur,1,MPI_Double_Complex,    0,iNode, &
               MPI_Comm_World,Status,MPIerror)
       else if ( IONode ) then
          call MPI_Send(ZE_cur,1,MPI_Double_Complex,iNode,iNode, &
               MPI_Comm_World,MPIerror)
       end if
#endif

       ! The test of the energy-point is performed on
       ! the calculating node...
       if ( Node == iNode ) then
          if ( cdabs(cE%e-ZE_cur) > 10._dp * EPS ) then
             write(*,*) 'GF-file: '//trim(El%GFfile)
             write(*,'(2(a,2(tr1,g20.13)))') 'Energies, TS / Gf:', &
                  cE%e / eV, ' /', ZE_cur / eV
             call die('Energy point in GF file does &
                  &not match the internal energy-point in transiesta. &
                  &Please correct your GF files.')
          end if
       end if
       
       ! If the k-point does not match what we expected...
       if ( IONode .and. ikpt /= ikGS ) then
          write(*,*) 'GF-file: '//trim(El%GFfile)
          write(*,'(2(a,i0))') 'k-point, TS / Gf: ', &
               ikpt, ' / ', ikGS
          call die('Read k-point in GF file does not match &
               &the requested k-point. Please correct your &
               &GF files.')
       end if

       if ( iEni == 1 ) then

          ! read in the electrode Hamiltonian and overlap...
          if ( IONode ) then
             if ( associated(El%HA) ) then
                read(uGF) El%HA
                read(uGF) El%SA
             else
                read(uGF) !El%HA
                read(uGF) !El%SA
             end if
          end if

#ifdef MPI
          if ( associated(El%HA) ) then
          call MPI_Bcast(El%HA(1,1,1),read_Size_HS,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(El%SA(1,1,1),read_Size_HS,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
          end if
#endif
       end if 


#ifdef MPI
       if ( IONode .and. iNode == Node ) then
#endif
          ! read in surface Green function
          read(uGF) El%GA
#ifdef MPI
       else if ( IONode ) then

          read(uGF) work(1:read_Size)
          
          call MPI_Send(work(1),read_Size,MPI_Double_Complex, &
               iNode,1,MPI_Comm_World,MPIerror) 
       else if ( Node ==  iNode ) then
          call MPI_Recv(El%GA(1),read_Size,MPI_Double_Complex, &
               0    ,1,MPI_Comm_World,Status,MPIerror)
       end if
#endif

    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS read_next_GS' )
#endif

  end subroutine read_next_GS_Elec


  ! Subroutine for reading in both the left and right next energy point
  subroutine read_next_GS(ispin,ikpt, bkpt, cE, &
       N_Elec, uGF, Elecs, &
       nzwork, zwork, &
       reread, forward, DOS, T )

    use parallel, only : IONode

#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum, MPI_Integer
    use mpi_siesta, only : MPI_Comm_World
#endif

    use m_ts_electype
    use m_ts_cctype
    use m_ts_electrode, only: calc_next_GS_Elec
      
    integer, intent(in) :: ispin, ikpt
    real(dp), intent(in) :: bkpt(3)
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: N_Elec, uGF(N_Elec)
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: zwork(nzwork)
    logical, intent(in), optional :: reread, forward
    ! When requesting density of states
    real(dp), intent(out), optional :: DOS(:,:), T(:)

    real(dp) :: kpt(3)
    integer :: NEReqs, i, j
#ifdef MPI
    integer :: MPIerror
#endif
    type(ts_c_idx) :: c

    c = cE
    ! save the current weight of the point
    ! This is where we include the factor-of-two for spin and
    ! and the (1/Pi) from DM = Im[G]/Pi
    ! Furthermore we include the weight of the k-point

    ! the number of points we wish to read in this segment
#ifdef MPI
    if ( any(Elecs(:)%out_of_core) ) then
      if ( cE%exist .and. .not. cE%fake) then
        i = 1
      else
        i = 0
      end if
      call MPI_AllReduce(i,NEReqs,1,MPI_Integer, MPI_Sum, &
          MPI_Comm_World, MPIerror)
    end if
#else
    NEReqs = 1
#endif

    ! TODO Move reading of the energy points
    ! directly into the subroutines which need them
    ! In this way we can save both GAA, Sigma AND Gamma arrays!!!!
    ! However, this will probably come at the expense 
    ! of doing the same "repetition" expansion twice, we can live with
    ! that!
    do i = 1 , N_Elec

       ! Transfer the k-point to the expanded supercell
       kpt(1) = bkpt(Elecs(i)%pvt(1)) / Elecs(i)%Bloch%B(1)
       kpt(2) = bkpt(Elecs(i)%pvt(2)) / Elecs(i)%Bloch%B(2)
       kpt(3) = bkpt(Elecs(i)%pvt(3)) / Elecs(i)%Bloch%B(3)

       ! Ensure zero in the semi-infinite direction
       select case ( Elecs(i)%t_dir )
       case ( 4 )
         kpt(2) = 0._dp
         kpt(3) = 0._dp
       case ( 5 )
         kpt(1) = 0._dp
         kpt(3) = 0._dp
       case ( 6 )
         kpt(1) = 0._dp
         kpt(2) = 0._dp
       case ( 7 )
         kpt(:) = 0._dp
       case default
         kpt(Elecs(i)%t_dir) = 0._dp
       end select
       
       ! If the index for the contour is negative
       ! It means that we are dealing with a Fermi
       ! charge correction
       ! If it is different from 1 we do not have an equilibrium contour
       if ( c%idx(1) /= 1 ) then
         ! In this case the energy is the eta value of the electrode
         if ( Elecs(i)%Eta > 0._dp ) then
#ifdef TBT_PHONON
           c%e = cmplx(real(cE%e,dp)**2,Elecs(i)%Eta, dp)
#else
           c%e = cmplx(real(cE%e,dp),Elecs(i)%Eta, dp)
#endif
         else
#ifdef TBT_PHONON
           c%e = cmplx(real(cE%e,dp)**2,aimag(cE%e)**2, dp)
#else
           c%e = cE%e
#endif
         end if
       end if
       if ( Elecs(i)%out_of_core ) then

          ! Backspace the file if needed
          if ( present(reread) ) then
            ! Currently the equilibrium energy points are just after
            ! the k-point, hence we will never need to backspace behind the
            ! HAA and SAA reads
            ! however, when we add kpoints for the bias contour, then it might be
            ! necessary!
            if ( IONode .and. reread ) then
              do j = 1 , NEReqs * 2
                backspace(unit=uGF(j))
              end do
              ! if ( new_kpt ) then
              !  do j = 1 , 2
              !   backspace(unit=uGFL)
              !   backspace(unit=uGFR)
              !  end do
              ! end if
            end if
          end if

          ! Set k-point for calculating expansion
          Elecs(i)%bkpt_cur = kpt
          
          call read_next_GS_Elec(uGF(i), NEReqs, &
               ikpt, Elecs(i), c, &
               nzwork, zwork, forward = forward)
       else
          ! This routine will automatically check
          ! (and SET) the k-point for the electrode.
          ! This is necessary for the expansion to work.
          if ( present(DOS) ) then
            j = Elecs(i)%no_u
            call calc_next_GS_Elec(Elecs(i),ispin,kpt,c%e, &
                nzwork, zwork, DOS(1:j,i) , T(i) )
#ifdef TBT_PHONON
             ! For phonons, we also require a factor of 2 * omega
             DOS(1:j,i) = 2._dp * real(cE%e,dp) * DOS(1:j,i)
#endif
          else
             call calc_next_GS_Elec(Elecs(i),ispin,kpt,c%e, &
                  nzwork, zwork)
          end if
       end if
    end do

  end subroutine read_next_GS


! This routine requires a call to check_Green before.
! It will return the k-points of the electrode Green function file

! ##################################################################
! ##            Read-in header of Green function file             ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ## Changed to F90 by Jose-Luis Mozos, jlm@icmab.es              ##
! ## Changed to only read header by Nick P. Andersen              ##
! ##    will only check against integer information and Ef shift. ##
! ##################################################################
  subroutine read_Green(funit,El,c_nkpar,c_NEn)
    
    use parallel,  only : IONode
    use sys ,      only : die
#ifdef MPI
    use mpi_siesta, only: MPI_Double_Precision
    use mpi_siesta, only: MPI_logical, MPI_Bcast
#endif
    use m_ts_electype
    real(dp) , parameter :: EPS = 1.e-6_dp
    
! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)  :: funit ! unit of gf-file
    type(Elec), intent(in) :: El
    integer, intent(in)  :: c_nkpar, c_NEn

! ***********************
! * LOCAL variables     *
! ***********************
    character(len=FILE_LEN) :: curGFfile  ! Name of the GF file

    integer :: nspin,nkpar,na,no,Bloch(3),NEn, pre_expand
    logical :: repeat
    integer :: na_u, no_u
    real(dp) :: mu ! The Fermi energy shift due to a voltage
    real(dp) :: ucell(3,3)
    logical :: errorGf

    ! we should only read if the GF-should exist
    if ( .not. El%out_of_core ) return

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_Green' )
#endif

    errorGF = .false.
    
    io_read: if ( IONode ) then

       ! Retrieve name of file currently reading
       inquire(unit=funit,name=curGFfile)

       ! read electrode information
       read(funit) nspin, ucell
       read(funit) na_u, no_u ! total atoms and total orbitals
       read(funit) na, no ! used atoms and used orbitals
       read(funit) ! xa, lasto
       read(funit) repeat, Bloch(:), pre_expand
       read(funit) mu
       ! read contour information
       read(funit) nkpar
       read(funit) ! kpoints, kweight
       read(funit) NEn
       read(funit)! ce

       ! Check Fermi shift
       if ( dabs(El%mu%mu-mu) > EPS ) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*)"The chemical shift in the electrode does not match the &
               &required shift!"
          write(*,'(2(a,f12.6))')"Found: ",mu,", expected: ",El%mu%mu
          errorGF = .true.
       end if

       ! Check # of energy points
       if (NEn .ne. c_NEn) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: NEn=',NEn,' expected:', c_NEn
          errorGF = .true.
       end if

       ! Check # of atoms
       if (na_u .ne. El%na_u) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: na_u=',na_u,' expected:', El%na_u
          errorGF = .true.
       end if

       ! Check # of atoms
       if (na .ne. El%na_used) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: na=',na,' expected:', El%na_used
          errorGF = .true.
       end if

       ! Check # of Bloch k-points
       if ( any(El%Bloch%B /= Bloch) ) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: unexpected no. Bloch expansion k-points'
          errorGF = .true.
       end if
       if ( El%repeat .neqv. repeat ) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: Ordering of Bloch repetitions is not the same (repeat or tile)'
          errorGF = .true.
       end if

       ! Check # of k-points
       if ( (El%is_gamma .and. nkpar /= 1) .or. &
           (.not. El%is_gamma .and. nkpar /= c_nkpar) ) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: Unexpected number of k-points'
          write(*,*) 'read_Green: ERROR: nkpt=',nkpar,' expected:', c_nkpar
          errorGF = .true.
       end if

       ! Check # of spin
       if (nspin .ne. El%nspin ) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: nspin=',nspin,' expected:', El%nspin
          errorGF = .true.
       end if

       ! Check # of orbitals
       if (no_u .ne. El%no_u) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: no_u=',no_u,' expected:', El%no_u
          errorGF = .true.
       end if

       ! Check # of orbitals
       if (no .ne. El%no_used) then
          write(*,*)"ERROR: Green function file: "//trim(curGFfile)
          write(*,*) 'read_Green: ERROR: no=',no,' expected:', El%no_used
          errorGF = .true.
       end if

       if ( El%Bloch%size() > 1 ) then
          if ( pre_expand /= El%pre_expand ) then
             write(*,*)"ERROR: Green function file: "//trim(curGFfile)
             write(*,*) 'read_Green: ERROR: Bloch pre-expansion not consistent'
             errorGF = .true.
          end if
       end if

    end if io_read

    if ( errorGF ) then
       call die("Error in reading GFfile: "//trim(curGFfile))
    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS read_Green' )
#endif

  end subroutine read_Green


  !> Re-opens a Gamma-only TSGF for a k-point sampled run
  !>
  !> When running transiesta calculations where one electrode is made of
  !> only periodicty along the semi-infinite direction we can skip creating
  !> the same self-energy and re-use the file content for a single k-point.
  !> This routine is a simple method to close/open/read-header of the
  !> GF files.
  !> This routine should *only* be called before incrementing the k-point
  !> in the calculation.
  !> This routine will do nothing if the electrode cannot re-use the Gamma-content.
  !> This routine does not check *any* content in the file as that *must* be done prior.
  subroutine reread_Gamma_Green(El, funit, NE, ispin)
    use parallel, only: IONode
    use m_ts_electype

    !< Electrode which is to be re-read.
    type(Elec), intent(in) :: El
    !< Unit that has the TSGF file currently open
    integer, intent(in) :: funit
    !< Total number of energy-points stored in the GF file
    integer, intent(in) :: NE
    !< Current segment of the spin components. I.e. if ispin == 2 the first spin-component will be read past
    integer, intent(in) :: ispin

    ! Local variables
    integer :: nk
    integer :: ie

    ! Only reread if io-node, the GF file is in use AND it is Gamma-only
    if ( .not. IONode ) return
    if ( .not. El%out_of_core ) return
    if ( .not. El%is_gamma ) return

    ! Close file
    close(funit)
    open(FILE=El%GFfile,UNIT=funit,FORM='UNFORMATTED')
    
    read(funit) ! nspin, ucell
    read(funit) ! na_u, no_u ! total atoms and total orbitals
    read(funit) ! na, no ! used atoms and used orbitals
    read(funit) ! xa, lasto
    read(funit) ! repeat, Bloch(:), pre_expand
    read(funit) ! mu
    ! read contour information
    read(funit) nk
    if ( nk /= 1 ) call die('reread_Green: nk must be 1, please delete the GF file and rerun')
    read(funit) ! kpoints, kweight
    read(funit) ! NEn
    read(funit) ! ce

    if ( ispin == 2 ) then

      ! We have to read past the first spin-component
      read(funit) ! H
      read(funit) ! S
      do iE = 1, NE
        read(funit) ! ik, ie, E
        ! We know H, S is actually just after iE == 1, however, since we are discarding the blocks
        read(funit) ! SE
      end do
 
    end if
    
  end subroutine reread_Gamma_Green


! ##################################################################
! ##    Check header of Green function file against settings      ##
! ##                            By                                ##
! ##      Nick Papior Andersen, nickpapior@gmail.com              ##
! ##                                                              ##
! ## Checks information an returns number of atoms and orbitals   ##
! ##################################################################
  subroutine check_Green(funit,El, &
      c_ucell,c_nkpar,c_kpar,c_wkpar, &
      c_NEn,c_ce, &
      xa_Eps, errorGF)

    use fdf, only: fdf_convfac
    use units, only: Ang
    use m_ts_cctype
    use m_ts_electype

    real(dp) , parameter :: EPS = 1d-6

! ***********************
! * INPUT variables     *
! ***********************
! file for reading, Green function file
    integer, intent(in)        :: funit
    type(Elec), intent(in)     :: El
    real(dp), intent(in)       :: c_ucell(3,3) ! Unit cell of the CONTACT
    ! k-point information
    integer, intent(in)        :: c_nkpar
    real(dp), intent(in)       :: c_kpar(3,c_nkpar) , c_wkpar(c_nkpar)
! Energy point on the contour used 
    integer, intent(in)        :: c_NEn
    complex(dp), intent(in)    :: c_ce(c_NEn)
    real(dp), intent(in)       :: xa_Eps
! ***********************
! * OUTPUT variables    *
! ***********************
! Return whether it is a correct Green function file
    logical, intent(out)       :: errorGF

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: mu ! The energy shift in the Fermi energy
    integer :: na_u, no_u ! total atoms and orbitals in unit-cell
    integer :: nspin, na, no, nkpar ! spin, # of atoms, # of orbs, # k-points
    integer :: Bloch(3) ! # repetitions in x, # repetitions in y
    real(dp), allocatable :: xa(:,:) ! electrode atomic coordinates
    integer, allocatable :: lasto(:) ! the electrode orbitals of the atoms
    real(dp), allocatable :: kpar(:,:) ! k-points
    real(dp), allocatable :: wkpar(:) ! k-point weights

    integer :: NEn ! # energy points on the contour
    complex(dp), allocatable :: ce(:)

! Helpers..
    character(FILE_LEN) :: curGFfile
    real(dp) :: ucell(3,3)
    integer :: iEn, pre_expand
    integer :: i, j, ia
    real(dp) :: kpt(3)
    logical :: localErrorGf, eXa, repeat
    logical :: showed_warn
    real(dp) :: Ry2eV

    ! we should only read if the GF-should exist
    errorGF = .false.
    if ( .not. El%out_of_core ) return

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE check_Green' )
#endif

    ! Conversion
    Ry2eV = fdf_convfac('Ry', 'eV')

    ! Initialize it to be an error unless the entire routine is runned through.
    ! Upon normal exit it will be changed to .FALSE.
    localErrorGf = errorGF
    
    ! Retrieve name of file currently reading
    inquire(unit=funit,name=curGFfile)

    ! Read in electrode information
    read(funit) nspin, ucell
    read(funit) na_u,no_u
    read(funit) na,no
    allocate(xa(3,na),lasto(na+1))
    read(funit) xa,lasto
    read(funit) repeat, Bloch(:),pre_expand
    read(funit) mu

    if ( El%nspin /= nspin ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of spin is wrong!"
      write(*,'(2(a,i0))') "Found: ",nspin,", expected: ",El%nspin
      localErrorGf = .true.
    end if
    if ( any(abs(El%cell-ucell) > EPS) ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Unit-cell is not consistent!"
      write(*,*) "Found (Ang):"
      write(*,'(3(3(tr1,f14.8),/))') ucell/Ang
      write(*,*) "Expected (Ang):"
      write(*,'(3(3(tr1,f14.8),/))') El%cell/Ang
      localErrorGf = .true.
    end if
    if ( El%na_u /= na_u ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of atoms is wrong!"
      write(*,'(2(a,i0))') "Found: ",na_u,", expected: ",El%na_u
      localErrorGf = .true.
    end if
    if ( El%na_used /= na ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of used atoms is wrong!"
      write(*,'(2(a,i0))') "Found: ",na,", expected: ",El%na_used
      localErrorGf = .true.
    end if
    if ( El%no_u /= no_u ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of orbitals is wrong!"
      write(*,'(2(a,i0))') "Found: ",no_u,", expected: ",El%no_u
      localErrorGf = .true.
    end if
    if ( El%no_used /= no ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of used orbitals is wrong!"
      write(*,'(2(a,i0))') "Found: ",no,", expected: ",El%no_used
      localErrorGf = .true.
    end if

    ! Initialize error parameter
    eXa = .false.
    ! We check elsewhere that the electrode is consistent with
    ! the FDF input
    kpt(:) = xa(:, 1) - El%xa_used(:, 1)
    do ia = 1 , min(na,El%na_used) ! in case it is completely wrong
      do i = 1 , 3
        eXa = eXa .or. &
            abs(xa(i,ia) - El%xa_used(i,ia) - kpt(i)) > xa_Eps
      end do
    end do
    if ( eXa ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Atomic coordinates are wrong:"
      write(*,'(1x,a,t35,a)') "Structure of GF electrode","| Electrode:"
      write(*,'(t3,3a15,''  |'',3a15)') &
          "X (Ang)","Y (Ang)","Z (Ang)", &
          "X (Ang)","Y (Ang)","Z (Ang)"
      do ia = 1, na
        write(*,'(t3,3(tr1,f14.8),''  |'',3(tr1,f14.8))') &
            xa(:,ia)/Ang, El%xa_used(:,ia)/Ang
      end do
      localErrorGf = .true.
    end if
    deallocate(xa)

    if ( any(lasto - El%lasto_used /= 0) ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of orbitals on used atoms is wrong!"
      write(*,'(a,1000i3)') "Found    lasto: ",lasto
      write(*,'(a,1000i3)') "Expected lasto: ",El%lasto_used
      localErrorGf = .true.
    end if
    deallocate(lasto)

    if ( any(El%Bloch%B(:)/=Bloch(:)) ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of Bloch expansions is wrong!"
      write(*,'(2(a,i3))') "Found Bloch-A1: ",Bloch(1),", expected: ",El%Bloch%B(1)
      write(*,'(2(a,i3))') "Found Bloch-A2: ",Bloch(2),", expected: ",El%Bloch%B(2)
      write(*,'(2(a,i3))') "Found Bloch-A3: ",Bloch(3),", expected: ",El%Bloch%B(3)
      localErrorGf = .true.
    end if
    if ( El%repeat .neqv. repeat ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*) 'Ordering of Bloch repetitions is not the same (repeat or tile)'
      localErrorGf = .true.
    end if
    if ( El%Bloch%size() > 1 ) then
      if ( pre_expand /= El%pre_expand ) then
        write(*,*)"ERROR: Green function file: "//trim(curGFfile)
        write(*,*)"Expecting a pre-expanded self-energy!"
        localErrorGf = .true.
      end if
    end if
    if ( abs(El%mu%mu-mu) > EPS ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"The chemical shift in the electrode does not match the &
          &required shift! [eV]"
      write(*,'(2(a,f12.6))')"Found: ",mu * Ry2eV,", expected: ",El%mu%mu * Ry2eV
      localErrorGf = .true.
    end if

    ! Read in general information about the context
    read(funit) nkpar
    allocate(kpar(3,nkpar),wkpar(nkpar))
    read(funit) kpar,wkpar

    if ( c_nkpar /= nkpar ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of k-points is wrong!"
      write(*,'(2(a,i4))') "Found: ",nkpar,", expected: ",c_nkpar
      localErrorGf = .true.
    end if

    ! Check k-points
    do i = 1 , min(c_nkpar,nkpar)
      ! As the k-points are in the Electrode unit cell
      ! we need to compare that with those of the CONTACT cell!
      ! The advantage of this is that the GF files can be re-used for
      ! the same system with different lengths between the electrode layers.
      call Elec_kpt(El,c_ucell,c_kpar(:,i),kpt, opt = 2)
      if ( abs(kpar(1,i)-kpt(1)) > EPS .or. &
          abs(kpar(2,i)-kpt(2)) > EPS .or. &
          abs(kpar(3,i)-kpt(3)) > EPS ) then
        write(*,*)"k-points are not the same:"
        do j = 1 , min(c_nkpar,nkpar)
          call Elec_kpt(El,c_ucell,c_kpar(:,j),kpt, opt = 2)
          write(*,'(3(tr1,f14.10),a,3(tr1,f14.10))') kpt(:),'  :  ',kpar(:,j)
        end do
        localErrorGf = .true.
        exit
      end if
      if ( dabs(c_wkpar(i)-wkpar(i)) > EPS ) then
        write(*,*)"k-point weights are not the same:"
        do j = 1 , c_nkpar
          write(*,'(f14.10,a,f14.10)') c_wkpar(j),'  :  ',wkpar(j)
        end do
        localErrorGf = .true.
        exit
      end if
    end do
    deallocate(kpar,wkpar)

    ! Read in information about the contour
    read(funit) NEn
    allocate(ce(NEn))
    read(funit) ce

    ! Check energy points
    if ( c_NEn /= NEn ) then
      write(*,*)"ERROR: Green function file: "//trim(curGFfile)
      write(*,*)"Number of energy points is not as expected!"
      write(*,'(2(a,i4))') "Found: ",NEn,", expected: ",c_NEn
      localErrorGf = .true.
    end if
    eXa = .false.
    showed_warn = .false.
    do iEn = 1 , min(c_NEn, NEn)
      if ( cdabs(ce(iEn)-c_ce(iEn)) > EPS ) then
        if ( showed_warn ) then
          ! do nothing
        else
          write(*,*) ' Warning: contours differ by >', EPS
          showed_warn = .true.
        end if

        if ( cdabs(ce(iEn)-c_ce(iEn)) > 10.d0*EPS ) then
          if ( eXa ) then
            ! do nothing
          else
            write(*,*) ' ERROR  : contours differ by >', 10.d0*EPS
            eXa = .true.
          end if
        end if

      end if
    end do
    if ( eXa ) then
      localErrorGf = .true.
      write(*,'(2(2(tr1,a25),tr1))') 'TS E.real [eV]', 'TS E.imag [eV]', &
          'File E.real [eV]', 'File E.imag [eV]'
      do iEn = 1 , min(c_NEn, NEn)
        write(*,'(2(2(tr1,e25.17),tr1))') c_ce(iEn) * Ry2eV, ce(iEn) * Ry2eV
      end do
    end if
    deallocate(ce)

    errorGF = localErrorGf
    
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS check_Green' )
#endif
    
  end subroutine check_Green

end module m_ts_GF
