module m_ts_io

  use precision, only : dp
  use parallel, only : Node

  implicit none

  public :: ts_read_TSHS_opt
  public :: ts_read_TSHS
  public :: ts_write_TSHS
  public :: fname_TSHS
  public :: TSHS_version
  public :: FC_index

  private

contains

  subroutine ts_read_TSHS_opt(TSHS,DUMMY,na_u,no_u,no_s,nspin,n_nzs, &
       xa, ucell, nsc, Qtot, Temp, Ef, &
       Gamma,Gamma_TS,kscell,kdispl,onlyS,lasto, &
       Bcast)

    use m_os, only : file_exist
#ifdef MPI
    use mpi_siesta
#endif

! ***********************
! * INPUT variables     *
! ***********************
    character(len=*), intent(in) :: TSHS

! ***********************
! * OUTPUT variables    *
! *********************** 
    integer, optional :: DUMMY ! MUST NEVER BE PASSED
    integer, intent(out), optional :: na_u, no_u, no_s, nspin, n_nzs, lasto(:), kscell(3,3), nsc(3)
    real(dp), intent(out), optional :: xa(:,:), ucell(3,3), Qtot, Temp, Ef, kdispl(3)
    logical, intent(out), optional :: Gamma, Gamma_TS, onlyS
    logical, intent(in), optional :: Bcast

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: lna_u, lno_u, lno_s, lnspin, ln_nzs
    integer :: uTSHS, version, tkscell(3,3), tnsc(3)
    real(dp), allocatable :: txa(:,:)
    real(dp) :: rtmp(3), tucell(3,3)
    logical :: fGamma, bool(3), lonlyS
#ifdef MPI
    integer :: buffer_size, ipos
    character(len=1), allocatable :: buffer(:)
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( present(DUMMY) ) call die('ts_read_TSHS_opt: Arguments has to be &
         &named. Please correct sources.')

    if ( .not. file_exist(TSHS, Bcast = Bcast) ) then
       call die('ERROR: Could not read '//trim(TSHS)//'.')
    end if

#ifdef NCDF_4
    ! If it is a netcdf file we read from that
    ! instead
    lna_u = len_trim(TSHS)
    if ( TSHS(lna_u-1:lna_u) == 'nc' ) then
       call ts_read_TSHS_opt_nc(TSHS,na_u=na_u,no_u=no_u,no_s=no_s, &
            nspin=nspin,n_nzs=n_nzs,xa=xa, ucell=ucell, &
            nsc=nsc, Qtot=Qtot, Temp=Temp, Ef=Ef, &
            Gamma=Gamma,Gamma_TS=Gamma_TS, &
            kscell=kscell,kdispl=kdispl,onlyS=onlyS,lasto=lasto, &
            Bcast = Bcast)
       return
    end if
#endif

    if ( Node == 0 ) then

       version = tshs_version(tshs)

       select case ( version )
       case ( 0 , 1 )
          ! do nothing
       case default
          call die('Unsupported TSHS version file [0,1]')
       end select

       call io_assign(uTSHS)
       open(file=trim(TSHS),unit=uTSHS,form='unformatted')
       if ( version /= 0 ) then
          read(uTSHS) version
       end if
       read(uTSHS) lna_u,lno_u,lno_s,lnspin,ln_nzs
       allocate(txa(3,lna_u))
       if ( present(na_u) ) na_u = lna_u
       if ( present(no_u) ) no_u = lno_u
       if ( present(no_s) ) no_s = lno_s
       if ( present(nspin) ) nspin = lnspin
       if ( present(n_nzs) ) n_nzs = ln_nzs
       if ( version == 0 ) then
          read(uTSHS) txa
          read(uTSHS) ! iza
          read(uTSHS) tucell
          tnsc = 0
       else if ( version == 1 ) then
          read(uTSHS) tnsc ! corrected nsc
          read(uTSHS) tucell, txa
       end if
       if ( present(nsc) ) nsc = tnsc
       if ( present(ucell) ) ucell = tucell
       if ( present(xa) )    xa    = txa
       if ( version == 0 ) then
          read(uTSHS) fGamma ! SIESTA_Gamma
          if ( present(Gamma) ) Gamma = fGamma
          read(uTSHS) lonlyS
          if ( present(OnlyS) ) onlys = lonlyS
          
          if ( present(Gamma_TS) ) then
             read(uTSHS) Gamma_TS
          else
             read(uTSHS) ! Gamma_TS
          end if
          if ( present(kscell) ) then
             read(uTSHS) kscell
          else
             read(uTSHS) ! ts_kscell_file
          end if
          if ( present(kdispl) ) then
             read(uTSHS) kdispl
          else
             read(uTSHS) ! ts_kdispl_file  
          end if
       else if ( version == 1 ) then
          read(uTSHS) fGamma, bool(1), lonlyS
          if ( present(Gamma) ) Gamma = fGamma
          if ( present(Gamma_TS) ) Gamma_TS = bool(1)
          if ( present(OnlyS) ) OnlyS = lonlyS
          read(uTSHS) tkscell, rtmp
          if ( present(kscell) ) kscell = tkscell
          if ( present(kdispl) ) kdispl = rtmp
          read(uTSHS) rtmp(1:3)
          if ( present(Ef) )   Ef   = rtmp(1)
          if ( present(Qtot) ) Qtot = rtmp(2)
          if ( present(Temp) ) Temp = rtmp(3)
       end if

       ! These quantities are not needed
       read(uTSHS) ! istep, ia1

       if ( present(lasto) ) then
          if ( size(lasto) /= lna_u+1 ) call die('ts_read_TSHS: Wrong size of lasto')
          read(uTSHS) lasto
       else
          read(uTSHS) ! lasto
       end if

       if ( version == 0 .and. .not. fGamma ) then
          read(uTSHS) ! indxuo
       end if

       read(uTSHS) ! numh

       if ( version == 0 ) then
          read(uTSHS) rtmp(1:2)
          if ( present(Qtot) ) Qtot = rtmp(1)
          if ( present(Temp) ) Temp = rtmp(2)

          if ( present(Ef) ) then
             read(uTSHS) Ef
          else
             read(uTSHS) ! Ef
          end if
       end if

       deallocate(txa)
       call io_close(uTSHS)

    end if

#ifdef MPI
    if ( present(Bcast) ) then
       ! if we do not request broadcasting, then return...
       if ( .not. Bcast ) return
    end if

    ! Broadcast na_u (for easy reference)
    call MPI_Bcast(lna_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)

#ifdef MPI_OLD
    if ( present(na_u) ) &
         call MPI_Bcast(na_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(no_u) ) &
         call MPI_Bcast(no_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(no_s) ) &
         call MPI_Bcast(no_s,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(nspin) ) &
         call MPI_Bcast(nspin,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(n_nzs) ) &
         call MPI_Bcast(n_nzs,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(nsc) ) &
         call MPI_Bcast(nsc,3,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(xa) ) &
         call MPI_Bcast(xa(1,1),3*lna_u,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(ucell) ) &
         call MPI_Bcast(ucell(1,1),9,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Gamma) ) &
         call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(OnlyS) ) &
         call MPI_Bcast(OnlyS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(Gamma_TS) ) &
         call MPI_Bcast(Gamma_TS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
    if ( present(kscell) ) &
         call MPI_Bcast(kscell,9,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(kdispl) ) &
         call MPI_Bcast(kdispl,3,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(lasto) ) &
         call MPI_Bcast(lasto(1),lna_u+1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    if ( present(Qtot) ) &
         call MPI_Bcast(Qtot,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Temp) ) &
         call MPI_Bcast(Temp,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    if ( present(Ef) ) &
         call MPI_Bcast(Ef,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#else

    ! this should be more than enough...
    buffer_size = 8 * (lna_u * 6 + 103)
    allocate(buffer(buffer_size))
    ! position of data in buffer...
    ipos = 0

    if ( Node == 0 ) then
       if ( present(na_u) ) & !  4
            call MPI_Pack(na_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_u) ) & !  4
            call MPI_Pack(no_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_s) ) & !  4
            call MPI_Pack(no_s,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nspin) ) & !  4
            call MPI_Pack(nspin,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(n_nzs) ) & !  4
            call MPI_Pack(n_nzs,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(xa) ) &    ! 8 * 3 * na_u
            call MPI_Pack(xa(1,1),3*lna_u,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(ucell) ) & ! 8 * 3 * 3
            call MPI_Pack(ucell(1,1),9,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &   ! 4
            call MPI_Pack(nsc(1),3,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) & ! 4
            call MPI_Pack(Gamma,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &! 4 * 3 * 3
            call MPI_Pack(kscell(1,1),9,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &! 8 * 3
            call MPI_Pack(kdispl(1),3,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(OnlyS) ) & ! 4
            call MPI_Pack(OnlyS,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &! 4
            call MPI_Pack(Gamma_TS,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(lasto) ) & ! 4 * (na_u+1)
            call MPI_Pack(lasto(1),lna_u+1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &  ! 8
            call MPI_Pack(Qtot,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &  ! 8
            call MPI_Pack(Temp,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &    ! 8
            call MPI_Pack(Ef,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)

       if ( ipos >= buffer_size .or. ipos < 0 .or. MPIerror /= MPI_Success ) then
          call die('Error in estimating the buffer-size for the &
               &TSHS reading. Please contact the developers')
       end if

    end if

    call MPI_Bcast(buffer,buffer_size,MPI_Packed, &
         0, MPI_Comm_World, MPIerror)

    if ( Node /= 0 ) then
       if ( present(na_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            na_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_s) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_s,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(nspin) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nspin,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(n_nzs) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            n_nzs,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(xa) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            xa(1,1),3*lna_u,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(ucell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            ucell(1,1),9,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nsc(1),3,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kscell(1,1),9,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kdispl(1),3,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(OnlyS) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            OnlyS,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma_TS,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(lasto) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            lasto(1),lna_u+1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Qtot,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Temp,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Ef,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
    end if       

    deallocate(buffer)

#endif

#endif

  end subroutine ts_read_TSHS_opt

#ifdef NCDF_4

  ! This is a routine for reading information from a 
  ! siesta NetCDF-4 format file
  subroutine ts_read_TSHS_opt_nc(TSHS,na_u,no_u,no_s,nspin,n_nzs, &
       xa, ucell, nsc, Qtot, Temp, Ef, &
       Gamma,Gamma_TS,kscell,kdispl,onlyS,lasto, &
       Bcast)

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta
#endif
! ***********************
! * INPUT variables     *
! ***********************
    character(len=*), intent(in) :: TSHS

! ***********************
! * OUTPUT variables    *
! *********************** 
    integer, intent(out), optional :: na_u, no_u, no_s, nspin, n_nzs, lasto(:), kscell(3,3), nsc(3)
    real(dp), intent(out), optional :: xa(:,:), ucell(3,3), Qtot, Temp, Ef, kdispl(3)
    logical, intent(out), optional :: Gamma, Gamma_TS, onlyS
    logical, intent(in), optional :: Bcast

! ***********************
! * LOCAL variables     *
! ***********************
    type(hNCDF) :: ncdf, grp
    integer :: lna_u, tnsc(3)
#ifdef MPI
    integer :: buffer_size, ipos
    character(len=1), allocatable :: buffer(:)
    integer :: MPIerror
#endif

    if ( present(onlyS) ) onlyS = .false. 

    call ncdf_open(ncdf,trim(TSHS),mode=NF90_NOWRITE)

    ! Now we read in the things
    call ncdf_inq_dim(ncdf,'na_u',len=lna_u)
    if ( present(na_u) ) na_u = lna_u
    if ( present(no_u) ) &
         call ncdf_inq_dim(ncdf,'no_u',len=no_u)
    if ( present(no_s) ) &
         call ncdf_inq_dim(ncdf,'no_s',len=no_s)
    if ( present(nspin) ) &
         call ncdf_inq_dim(ncdf,'spin',len=nspin)
    if ( present(xa) ) &
         call ncdf_get_var(ncdf,'xa',xa)
    if ( present(ucell) ) &
         call ncdf_get_var(ncdf,'cell',ucell)
    call ncdf_get_var(ncdf,'nsc',tnsc)
    if ( present(nsc) ) nsc = tnsc
    if ( present(Gamma) ) then
       Gamma = sum(tnsc) == 1
       if ( present(Gamma_TS) ) Gamma_TS = Gamma
    end if
    if ( present(Ef) ) &
         call ncdf_get_var(ncdf,'Ef',Ef)
    if ( present(Qtot) ) &
         call ncdf_get_var(ncdf,'Qtot',Qtot)

    if ( Node == 0 ) then
       if ( present(lasto) ) then
          if ( size(lasto) /= lna_u+1 ) call die('ts_read_TSHS: Wrong size of lasto')
          lasto(1) = 0
          call ncdf_get_var(ncdf,'lasto',lasto(2:lna_u+1))
       end if
    end if

    if ( present(n_nzs) ) then
       call ncdf_open_grp(ncdf,'SPARSE',grp)
       call ncdf_get_var(grp,'nnzs',n_nzs)
    end if

    call ncdf_open_grp(ncdf,'SETTINGS',grp)
    if ( present(Temp) ) &
         call ncdf_get_var(grp,'ElectronicTemperature',Temp)
    if ( present(kscell) ) &
         call ncdf_get_var(grp,'BZ',kscell)
    if ( present(kdispl) ) &
         call ncdf_get_var(grp,'BZ_displ',kdispl)

    call ncdf_close(ncdf)

#ifdef MPI
    if ( present(Bcast) ) then
       ! if we do not request broadcasting, then return...
       if ( .not. Bcast ) return
    end if

    ! Broadcast na_u (for easy reference)
    call MPI_Bcast(lna_u,1,MPI_Integer,0,MPI_Comm_World,MPIerror)

    ! this should be more than enough...
    buffer_size = 8 * (lna_u * 6 + 103)
    allocate(buffer(buffer_size))
    ! position of data in buffer...
    ipos = 0

    if ( Node == 0 ) then
       if ( present(na_u) ) & !  4
            call MPI_Pack(na_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_u) ) & !  4
            call MPI_Pack(no_u,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(no_s) ) & !  4
            call MPI_Pack(no_s,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nspin) ) & !  4
            call MPI_Pack(nspin,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(n_nzs) ) & !  4
            call MPI_Pack(n_nzs,1,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(xa) ) &    ! 8 * 3 * na_u
            call MPI_Pack(xa(1,1),3*lna_u,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(ucell) ) & ! 8 * 3 * 3
            call MPI_Pack(ucell(1,1),9,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &   ! 4
            call MPI_Pack(nsc(1),3,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) & ! 4
            call MPI_Pack(Gamma,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &! 4 * 3 * 3
            call MPI_Pack(kscell(1,1),9,MPI_Integer, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &! 8 * 3
            call MPI_Pack(kdispl(1),3,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &! 4
            call MPI_Pack(Gamma_TS,1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(lasto) ) & ! 4 * (na_u+1)
            call MPI_Pack(lasto(1),lna_u+1,MPI_Logical, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &  ! 8
            call MPI_Pack(Qtot,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &  ! 8
            call MPI_Pack(Temp,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &    ! 8
            call MPI_Pack(Ef,1,MPI_Double_Precision, &
            buffer,buffer_size, ipos, MPI_Comm_World, MPIerror)

       if ( ipos >= buffer_size .or. ipos < 0 .or. MPIerror /= MPI_Success ) then
          call die('Error in estimating the buffer-size for the &
               &TSHS reading. Please contact the developers')
       end if

    end if

    call MPI_Bcast(buffer,buffer_size,MPI_Packed, &
         0, MPI_Comm_World, MPIerror)

    if ( Node /= 0 ) then
       if ( present(na_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            na_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_u) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_u,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(no_s) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            no_s,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(nspin) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nspin,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(n_nzs) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            n_nzs,1,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(xa) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            xa(1,1),3*lna_u,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(ucell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            ucell(1,1),9,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(nsc) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            nsc(1),3,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(kscell) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kscell(1,1),9,MPI_Integer, &
            MPI_Comm_World, MPIerror)
       if ( present(kdispl) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            kdispl(1),3,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Gamma_TS) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Gamma_TS,1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(lasto) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            lasto(1),lna_u+1,MPI_Logical, &
            MPI_Comm_World, MPIerror)
       if ( present(Qtot) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Qtot,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Temp) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Temp,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
       if ( present(Ef) ) &
            call MPI_UnPack(buffer,buffer_size,ipos, &
            Ef,1,MPI_Double_Precision, &
            MPI_Comm_World, MPIerror)
    end if       

    deallocate(buffer)

#endif

  end subroutine ts_read_TSHS_opt_nc
#endif

  subroutine ts_read_TSHS(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, nsc, na_u, no_u, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       sp, H, S, isc_off, &
       Ef, Qtot, Temp, &
       istep, ia1, tag, &
       Bcast)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! logical       Gamma         : Is only gamma point used?
! logical       TSGamma       : Is only TS gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer n_nzs               : First dimension of listh, H, S and
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(n_nzs)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(n_nzs,Enspin)     : Hamiltonian in sparse form
! real*8  S(n_nzs)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,n_nzs),isc_off : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

    use sys,          only : die
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc
    use geom_helper,  only : iaorb, ucorb
    use m_io_s
    use m_os, only : file_exist
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData1D
    use class_dSpData2D
    use m_sparse, only : list_col_correct, xij_offset, calc_nsc

    implicit none

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(out) :: onlyS, Gamma, TSGamma
    real(dp), intent(out) :: ucell(3,3)
    integer, intent(out) :: nsc(3), na_u, no_u, nspin
    integer, intent(out) :: kscell(3,3)
    real(dp), intent(out) :: kdispl(3)
    real(dp), pointer :: xa(:,:)
    integer, pointer :: lasto(:) ! (0:na_u) 
    type(Sparsity), intent(inout) :: sp
    type(dSpData2D), intent(inout) :: H
    type(dSpData1D), intent(inout) :: S
    integer, pointer :: isc_off(:,:)
    real(dp), intent(out) :: Ef, Qtot,Temp
    ! These have to be set before entrance (makes it possible to read
    ! in FCrun TSHS files...)
    integer, intent(out) :: istep, ia1
    character(len=*), intent(in), optional :: tag
    ! If true it will broadcast every information within the code...
    logical, intent(in), optional :: Bcast
    
! ************************
! * LOCAL variables      *
! ************************
    type(OrbitalDistribution), pointer :: dit
    type(dSpData2D) :: xij
    integer :: iu, version, no_s, n_s, n_nzs
    integer :: i, j, ind
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer, allocatable :: indxuo(:)
    real(dp), pointer :: lxij(:,:)
    character(len=250) :: ltag
    logical :: lBcast, exist
#ifdef MPI
    integer :: all_I(0:9)
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( .not. file_exist(filename, Bcast = Bcast ) ) then
       call die('ERROR: Could not read '//trim(filename)//'.')
    end if

#ifdef NCDF_4
    ! If it is a NetCDF file, we call the netCDF
    ! routine
    i = len_trim(filename)
    if ( filename(i-1:i) == 'nc' ) then
       call ts_read_TSHS_nc(filename, &
       Gamma, TSGamma, &
       ucell, nsc, na_u, no_u, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       sp, H, S, isc_off, &
       Ef, Qtot, Temp, &
       tag = tag, Bcast = Bcast)

       OnlyS = .false.
       Temp = 0._dp
       istep = 0
       ia1 = 0

       return

    end if
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_read' )
#endif

    nullify(xa,lasto,ncol,l_ptr,l_col,isc_off)

    ltag = trim(filename)
    if ( present(tag) ) ltag = trim(tag)

    ! Determine whether to broadcast afterwards
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    ! Many variables are needed initialization due
    ! to compiletime unawareness
    nsc(:) = 0

    exist = file_exist(filename, Bcast = Bcast )

    if ( Node == 0 ) then

       if ( .not. exist ) then
          call die('ERROR: Could not read '//trim(filename)//'.')
       end if

       ! Get file version
       version = tshs_version(filename)

       select case ( version )
       case ( 0 , 1 )
          ! do nothing
       case default
          call die('Unsupported TSHS version file [0,1]')
       end select

       ! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='old' )
       
       ! Read Dimensions Information
       if ( version /= 0 ) then
          read(iu) ! version
       end if
       read(iu) na_u, no_u, no_s, nspin, n_nzs
       
       ! Read Geometry information
       allocate(xa(3,na_u)) 
       call memory('A','D',3*na_u,'iohs')
       call memory('A','I',na_u,'iohs')
       if ( version == 0 ) then
          read(iu) xa
          read(iu) ! iza
          read(iu) ucell
       else if ( version == 1 ) then
          read(iu) nsc ! original nsc, and the corrected nsc
          read(iu) ucell, xa
       end if

       ! Read k-point sampling information
       if ( version == 0 ) then
          read(iu) Gamma
          read(iu) onlyS
          read(iu) TSGamma
          read(iu) kscell
          read(iu) kdispl
       else if ( version == 1 ) then
          read(iu) Gamma, TSGamma, onlyS
          read(iu) kscell, kdispl
          read(iu) Ef, Qtot, Temp
       end if
       read(iu) istep, ia1

       allocate(lasto(0:na_u))
       call memory('A','I',1+na_u,'iohs')
       read(iu) lasto

       if ( version == 0 ) then
          ! We check that the format of the supercell indices
          ! are consistent with the general formalism
          allocate(indxuo(no_s))
          read(iu) indxuo
          do i = 1 , no_s
             if ( indxuo(i) /= ucorb(i,no_u) ) &
                  call die('Error in indxuo, not recognized')
          end do
          deallocate(indxuo)
       end if

    end if

#ifdef MPI
    if ( lBcast ) then
       ! Bcast initial sizes
       if ( Node == 0 ) then
          all_I(0:9) = (/version,na_u,no_u,nspin,n_nzs,nsc(1),nsc(2),nsc(3),istep,ia1/)
       end if
       call MPI_Bcast(all_I(0),10,MPI_Integer,0,MPI_Comm_World,MPIerror)
       version = all_I(0)
       na_u = all_I(1)
       no_u = all_I(2)
       nspin = all_I(3)
       n_nzs = all_I(4)
       nsc(1) = all_I(5)
       nsc(2) = all_I(6)
       nsc(3) = all_I(7)
       istep = all_I(8)
       ia1 = all_I(9)
       call MPI_Bcast(ucell(1,1),9,MPI_Double_Precision,0, &
            MPI_Comm_World,MPIerror)
       call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(TSGamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(onlyS,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kscell(1,1),9,MPI_Integer,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kdispl(1),3,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)

       if ( Node /= 0 ) then
          allocate(xa(3,na_u))
          call memory('A','D',3*na_u,'iohs')
          allocate(lasto(0:na_u))
          call memory('A','I',1+na_u,'iohs')
       end if
       call MPI_Bcast(xa(1,1),3*na_u,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(lasto(0),1+na_u,MPI_Integer,0, MPI_Comm_World,MPIerror)

    end if
#endif

    if ( version == 0 ) then

       nullify(ncol,l_col)
       allocate(ncol(no_u),l_col(n_nzs))

       if ( Node == 0 ) then
          
          read(iu) ncol

          ! Read Electronic Structure Information
          read(iu) Qtot,Temp
          read(iu) Ef
          
          j = 0
          do i = 1 , no_u
             read(iu) l_col(j+1:j+ncol(i))
             j = j + ncol(i)
          end do

       end if

#ifdef MPI
       if ( lBcast ) then
          call MPI_Bcast(ncol,no_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
          call MPI_Bcast(l_col,n_nzs,MPI_Integer,0, MPI_Comm_World,MPIerror)
       end if
#endif

       nullify(l_ptr)
       allocate(l_ptr(no_u))
       
       l_ptr(1) = 0
       do i = 2 , no_u
          l_ptr(i) = l_ptr(i-1) + ncol(i-1)
       end do

       ! Create the sparsity pattern
       call newSparsity(sp,no_u,no_u, &
            n_nzs, ncol, l_ptr, l_col, trim(ltag))

       deallocate(ncol,l_ptr,l_col)
       nullify(ncol,l_ptr,l_col)

    else if ( version == 1 ) then
       call io_read_Sp(iu,no_u,sp,tag=trim(ltag),Bcast=Bcast)
    end if

    ! Read in S
    call io_read_d1D(iu,sp,S,trim(ltag)//': S',Bcast=Bcast)

    ! Utilize the same distribution for the others
    dit => dist(S)

    if ( .not. onlyS ) then
       ! Read in H
       call io_read_d2D(iu,sp,H,nspin,trim(ltag)//': H',Bcast=Bcast,dit=dit)
    end if

    if ( .not. Gamma ) then

       if ( version == 0 ) then
          call io_read_d2D(iu,sp,xij,3,'xij',sparsity_dim=2, &
               dit=dit, Bcast=Bcast)

          call attach(sp,n_col=ncol)
          lxij => val(xij)

          ! A stupid "bug" in the very old TRANSIESTA code
          ! was transposing the xij array before writing... sigh...
          ind = 0
          do i = 1 , no_u
             lxij(:,ind+1:ind+ncol(i)) = transpose( &
                  reshape(lxij(:,ind+1:ind+ncol(i)), (/ncol(i),3/)))
             ind = ind + ncol(i)
          end do

          ! We do not need (MUST NOT) do bcast the routines:
          !   calc_nsc, list_col_correct, xij_offset
          ! They are duplicated on all nodes, hence b-casting
          ! will introduce wrong columns...
          call calc_nsc(ucell,na_u,xa,lasto,xij,nsc) 

          ! Ensure that list_col is correctly formatted (not always
          ! needed, but for consistency)
          call list_col_correct(ucell,nsc,na_u,xa,lasto,xij)
          
          call xij_offset(ucell,nsc,na_u,xa,lasto,xij,isc_off)

          ! We do not need the xij array anymore... :)
          call delete(xij)

       else if ( version == 1 ) then
             
          ! Number of supercells
          n_s = product(nsc)
          call re_alloc(isc_off,1,3,1,n_s)

          if ( Node == 0 ) then
             read(iu) isc_off
          end if
#ifdef MPI
          if ( lBcast ) then
             call MPI_Bcast(isc_off(1,1),3*n_s, MPI_Integer, 0, &
                  MPI_Comm_World,MPIerror)
          end if
#endif

       end if
    else

       ! Number of supercells
       n_s = product(nsc)
       call re_alloc(isc_off,1,3,1,n_s)
       isc_off(:,:) = 0

    end if

    if ( Node == 0 ) then
       ! Close file
       call io_close( iu )
    end if

#ifdef MPI
    if ( lBcast ) then
       call MPI_Bcast(Qtot,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Temp,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Ef,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
    end if
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_read' )
#endif

  end subroutine ts_read_TSHS

#ifdef NCDF_4
  subroutine ts_read_TSHS_nc(filename, &
       Gamma, TSGamma, &
       ucell, nsc, na_u, no_u, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       sp, H, S, isc_off, &
       Ef, Qtot, Temp, &
       tag, Bcast)

    use m_ncdf_io, only : cdf_r_Sp, cdf_r_d1D, cdf_r_d2D

    use netcdf_ncdf, ncdf_parallel => parallel
    use sys,          only : die
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData1D
    use class_dSpData2D

    implicit none

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(out) :: Gamma, TSGamma
    real(dp), intent(out) :: ucell(3,3)
    integer, intent(out) :: nsc(3), na_u, no_u, nspin
    integer, intent(out) :: kscell(3,3)
    real(dp), intent(out) :: kdispl(3)
    real(dp), pointer :: xa(:,:)
    integer, pointer :: lasto(:) ! (0:na_u) 
    type(Sparsity), intent(inout) :: sp
    type(dSpData2D), intent(inout) :: H
    type(dSpData1D), intent(inout) :: S
    integer, pointer :: isc_off(:,:)
    real(dp), intent(out) :: Ef, Qtot,Temp
    character(len=*), intent(in), optional :: tag
    ! If true it will broadcast every information within the code...
    logical, intent(in), optional :: Bcast
    
! ************************
! * LOCAL variables      *
! ************************
    type(hNCDF) :: ncdf, grp
    type(OrbitalDistribution), pointer :: dit
    integer :: n_s
    character(len=250) :: ltag
    logical :: lBcast
#ifdef MPI
    integer :: all_I(6)
    integer :: MPIerror
#endif

    nullify(xa,lasto,isc_off)

    ltag = trim(filename)
    if ( present(tag) ) ltag = trim(tag)

    ! Determine whether to broadcast afterwards
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    call ncdf_open(ncdf,filename,mode=NF90_NOWRITE)
    
    call ncdf_inq_dim(ncdf,'na_u',len=na_u)
    call ncdf_inq_dim(ncdf,'no_u',len=no_u)
    call ncdf_inq_dim(ncdf,'spin',len=nspin)
    
    call ncdf_get_var(ncdf,'nsc',nsc)
    call ncdf_get_var(ncdf,'cell',ucell)
    call ncdf_get_var(ncdf,'Ef',Ef)
    call ncdf_get_var(ncdf,'Qtot',Qtot)

    ! The Brillouin zone sampling for siesta
    call ncdf_open_grp(ncdf,'SETTINGS',grp)
    call ncdf_get_var(grp,'BZ',kscell)
    call ncdf_get_var(grp,'BZ_displ',kdispl)
    call ncdf_get_var(grp,'ElectronicTemperature',Temp)

#ifdef MPI
    if ( lBcast ) then
       ! Bcast initial sizes
       if ( Node == 0 ) then
          all_I(1:6) = (/na_u,no_u,nspin,nsc(1),nsc(2),nsc(3)/)
       end if
       call MPI_Bcast(all_I(1),6,MPI_Integer,0,MPI_Comm_World,MPIerror)
       na_u = all_I(1)
       no_u = all_I(2)
       nspin = all_I(3)
       nsc(1) = all_I(4)
       nsc(2) = all_I(5)
       nsc(3) = all_I(6)
       call MPI_Bcast(ucell(1,1),9,MPI_Double_Precision,0, &
            MPI_Comm_World,MPIerror)
       call MPI_Bcast(Gamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(TSGamma,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kscell(1,1),9,MPI_Integer,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(kdispl(1),3,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(Ef,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Qtot,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(Temp,1,MPI_Double_Precision,0, MPI_Comm_World,MPIerror)

    end if
#endif

    if ( Node == 0 ) then
       ! Read Geometry information
       allocate(xa(3,na_u))
       allocate(lasto(0:na_u))
       call ncdf_get_var(ncdf,'xa',xa)
       call ncdf_get_var(ncdf,'lasto',lasto(1:na_u))
       lasto(0) = 0
    end if

    call ncdf_open_grp(ncdf,'SPARSE',grp)

    ! Number of supercells
    n_s = product(nsc)
    call re_alloc(isc_off,1,3,1,n_s)
    call ncdf_get_var(grp,'isc_off',isc_off)

    ! Calculate the Gamma and TSGamma
    Gamma   = ( n_s == 1 )
    TSGamma = sum(kscell) == 1

#ifdef MPI
    if ( lBcast ) then
       if ( Node /= 0 ) then
          allocate(xa(3,na_u))
          allocate(lasto(0:na_u))
       end if
       call MPI_Bcast(xa(1,1),3*na_u,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(lasto(0),1+na_u,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(isc_off(1,1),3*n_s, MPI_Integer, 0, &
            MPI_Comm_World,MPIerror)
    end if
#endif

    call cdf_r_Sp(grp, no_u, sp, tag = trim(filename), Bcast = Bcast )
    call cdf_r_d1D(grp, 'S', sp, S, tag = trim(filename)//': S', &
         Bcast = Bcast )
    dit => dist(S)
    call cdf_r_d2D(grp, 'H', sp, H, nspin, tag = trim(filename)//': H', &
         Bcast = Bcast, dit = dit )
             
    call ncdf_close(ncdf)

  end subroutine ts_read_TSHS_nc
#endif

  subroutine ts_write_TSHS(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, nsc, isc_off, na_u, no_s, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       H2D, S1D, indxuo, &
       Ef, Qtot, Temp, &
       istep, ia1)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! logical       Gamma         : Is only gamma point used?
! logical       TSGamma       : Is only TS gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer n_nzs               : First dimension of l_col, H, S and
!                               second of xij
! integer ncol(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer l_ptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer l_col(n_nzs)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(n_nzs,Enspin)     : Hamiltonian in sparse form
! real*8  S(n_nzs)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,n_nzs)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *********************************************************************

    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData1D
    use class_dSpData2D
    use geom_helper, only : ucorb
    use m_sparse, only : xij_offset

    use m_io_s, only: io_write_Sp, io_write_d1D, io_write_d2D
#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(in) :: onlyS
    logical, intent(in) :: Gamma, TSGamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)
    integer, intent(in) :: nsc(3), na_u, no_s, nspin
    integer, intent(in) :: isc_off(3,product(nsc))
    real(dp), intent(in) :: xa(3,na_u)
    type(dSpData2D), intent(inout) :: H2D
    type(dSpData1D), intent(inout) :: S1D
    integer, intent(in) :: indxuo(no_s)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: Ef
    real(dp), intent(in) :: Qtot, Temp
    integer, intent(in) :: istep, ia1
    
! ************************
! * LOCAL variables      *
! ************************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer :: iu, no_l, no_u, n_nzs
    integer :: i, n_s
    integer :: n_nzsg
    integer, allocatable, target :: gncol(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE ts_io_write' )
#endif

    ! Gather sparse pattern
    dit => dist(H2D)
    sp => spar(H2D)
    call attach(sp,nrows=no_l,nrows_g=no_u,nnzs=n_nzs, &
         n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    ! This setting is (I am afraid) not constant
    ! with system. I had suspected that n_s ALWAYS would 
    ! equal no_s / no_u, but this seems not to be the case.
    n_s = no_s / no_u 
    if ( mod(no_s,no_u) /= 0 ) call die('Error in supercell orbitals, no_s')
    if ( n_s /= product(nsc) ) call die('Error in supercell orbitals, nsc')

    ! Check that indxuo is constructed in a sane format
    do i = 1 , no_s
       if ( indxuo(i) /= ucorb(i,no_u) ) &
            call die('Error in indxuo, could not understand the format, &
            &please consult the developers.')
    end do

#ifdef MPI
    ! get total number of non-zero elements
    call MPI_Reduce(n_nzs,n_nzsg,1,MPI_Integer, MPI_SUM, 0, &
         MPI_Comm_World,MPIerror)
#else
    n_nzsg = n_nzs
#endif

    if ( Node == 0 ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=filename, form='unformatted', status='unknown' )

       ! Write file version
       write(iu) 1 ! This is version ONE of the file format

       ! Write Dimensions Information
       write(iu) na_u, no_u, no_s, nspin, n_nzsg
       write(iu) nsc
       
       ! Write Geometry information
       write(iu) ucell, xa

       ! Write k-point sampling information
       write(iu) Gamma, TSGamma, onlyS
       write(iu) kscell, kdispl

       ! Write Electronic Structure Information
       write(iu) Ef, Qtot, Temp

       ! The phonon-calculation
       write(iu) istep, ia1

       write(iu) lasto

    end if

    allocate(gncol(no_u))
    gncol(1) = -1

    ! Write out sparsity pattern
    call io_write_Sp(iu,sp,dit,gncol=gncol)

    ! Write overlap matrix
    call io_write_d1D(iu,S1D,gncol=gncol)

    if ( .not. onlyS ) then

       ! Write Hamiltonian
       call io_write_d2D(iu,H2D,gncol=gncol)

    end if

    deallocate(gncol)

    if ( Node == 0 ) then

       if ( .not. Gamma ) then
          write(iu) isc_off
       end if

       ! Close file
       call io_close( iu )

    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS ts_io_write' )
#endif

  end subroutine ts_write_TSHS

  function fname_TSHS(slabel,istep,onlyS,ia1) result(fname)
    character(len=*), intent(in) :: slabel
    integer, intent(in), optional :: istep
    integer, intent(in), optional :: ia1
    logical, intent(in), optional :: onlyS
    character(len=255) :: fname
    integer :: fL
    logical :: lonlyS
    integer :: listep, lia1

    ! Initialize...
    fname = ' '

    lonlyS = .false.
    if ( present(onlyS) ) lonlyS = onlyS
    listep = -1
    if ( present(istep) ) listep = istep
    lia1 = 0
    if ( present(ia1) ) lia1 = ia1

    ! This is the criteria for not doing an FCrun
    ! There will never be an atom denoted by index 0
    if ( lia1 /= 0 .and. listep >= 0 ) then
       if ( listep == 0 ) then
          ! We need the extra 00000's to let python 
          ! sort in a stringent way... :(
          write(fname,'(a,i5.5)') '.',0
       else
          ! Be sure to have the python sorting working by adding
          ! zeroes
          write(fname,'(a,i5.5,a,i1)') '.',lia1,'-',listep
       end if
       fname = trim(slabel)//trim(fname)
    else if ( listep >= 0 ) then
       ! Simply a consecutive number increase.
       write(fname,'(a,i0)') '.',listep
       fname = trim(slabel)//trim(fname)
    else
       fname = slabel
    end if
    
    fL = len_trim(fname)
    if ( lonlyS ) then
       fname = trim(fname)//'.onlyS'
    else
       fname = trim(fname)//'.TSHS'
    end if

  end function fname_TSHS

  subroutine FC_index(istep,ia1,ostep,oa1) 
    integer, intent(in) :: istep
    integer, intent(in) :: ia1
    integer, intent(out) :: ostep
    integer, intent(out) :: oa1

    if ( istep == 0 ) then
       ostep = 0
       oa1 = ia1
    else
       ostep = mod(istep-1,6) + 1
       oa1 = (istep - mod(istep-1,6))/6 + ia1
    end if

  end subroutine FC_index


  function TSHS_version(fname) result(version)
    character(len=*), intent(in) :: fname
    integer :: version
    integer :: iu
    integer :: na_u, no_u, no_s, nspin, n_nzs, err
    
    external :: io_assign, io_close

    ! Initialize
    version = -1
    if ( Node /= 0 ) return

    ! Open file
    call io_assign( iu )
    open( iu, file=fname, form='unformatted', status='unknown' )

    read(iu,iostat=err) na_u, no_u, no_s, nspin, n_nzs
    if ( err == 0 ) then
       ! we can successfully read 5 integers
       version = 0
    else
       backspace(iu)
       read(iu,iostat=err) version
    end if

    call io_close(iu)

  end function TSHS_version

end module m_ts_io
