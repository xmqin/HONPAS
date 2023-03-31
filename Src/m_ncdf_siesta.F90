! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_siesta

  use precision, only : sp, dp, grid_p
  use parallel, only : Node, Nodes
#ifdef NCDF_4
  use siesta_options, only: cdf_comp_lvl, cdf_w_parallel
  use variable
  use dictionary
  use netcdf_ncdf, ncdf_parallel => parallel
  use m_ncdf_io
#endif
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  private

#ifdef NCDF_4

  public :: cdf_init_file
  public :: cdf_init_grid

  public :: cdf_save_settings
  public :: cdf_save_basis
  public :: cdf_save_fc
  public :: cdf_save_state
  public :: cdf_save_grid

#endif

contains

#ifdef NCDF_4

  subroutine cdf_init_file(fname, is_fc)
    use class_Sparsity
    use files, only : slabel
    use atomlist, only: no_u, no_s, lasto, Qtot
    use siesta_geom, only: na_u, nsc
    use sparse_matrices, only: sparse_pattern
    use m_spin, only : spin
    use siesta_options, only: sname, isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN
    use siesta_options, only: SOLVE_MINIM, SOLVE_TRANSI
    use siesta_options, only: SOLVE_PEXSI
    use siesta_options, only: savehs
    use siesta_options, only: fixspin, total_spin
    use siesta_options, only: ia1, ia2, dx ! FC information
    use m_timestamp, only: datestring
    use m_ts_electype, elec_name => name
    use m_ts_options, only : Volt, N_Elec, Elecs
    use m_ts_options, only : TS_HS_Save

    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: is_fc

    ! Local variables
    type(hNCDF) :: ncdf, grp, grp2
    type(dictionary_t) :: dic, d
    integer :: n_nzs, tmp, i, chks(3), iEl
    integer, allocatable :: ibuf(:)
    logical :: lis_fc
#ifdef MPI
    integer :: MPIerror
#endif

    lis_fc = .false.
    if ( present(is_fc) ) lis_fc = is_fc

    ! We always re-write the file...
    call ncdf_create(ncdf,fname, &
         mode=ior(NF90_WRITE,NF90_NETCDF4), overwrite=.true.)

#ifdef MPI
    tmp = nnzs(sparse_pattern)
    call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, &
         MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
    n_nzs = nnzs(sparse_pattern)
#endif

    ! First we create all the dimensions
    ! necessary
    d = ('DIMna_u'.kv.na_u)//('DIMno_u'.kv.no_u)
    d = d//('DIMno_s'.kv.no_s)//('DIMspin'.kv.spin%H)
    d = d//('DIMxyz'.kv. 3) // ('DIMn_s'.kv. product(nsc))
    d = d// ('DIMone'.kv. 1)
    call ncdf_crt(ncdf,d)
    call delete(d)

    dic = dic//('info'.kv.'Last orbital of equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
        compress_lvl=0, atts=dic)
    call ncdf_put_var(ncdf,'lasto',lasto(1:na_u))
    
    dic = dic//('info'.kv.'Total charge')
    call ncdf_def_var(ncdf,'Qtot',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)
    call ncdf_put_var(ncdf,'Qtot',Qtot)

    if ( fixspin ) then
      dic = dic//('info'.kv.'Total spin')
      call ncdf_def_var(ncdf,'Qspin',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
      call ncdf_put_var(ncdf,'Qspin',total_spin)
    end if
    
    ! Create all necessary containers...
    dic = dic//('info'.kv.'Number of supercells in each unit-cell direction') 
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Fermi level')//('unit'.kv.'Ry')
    if ( fixspin ) then
      call ncdf_def_var(ncdf,'Ef',NF90_DOUBLE,(/'spin'/), &
          compress_lvl=0, atts=dic)
    else
      call ncdf_def_var(ncdf,'Ef',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
    end if

    dic = dic//('info'.kv.'Atomic coordinates')//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Unit cell')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Atomic forces')//('unit'.kv.'Ry/Bohr')
    call ncdf_def_var(ncdf,'fa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Cell stress')//('unit'.kv.'Ry/Bohr**3')
    call ncdf_def_var(ncdf,'stress',NF90_DOUBLE,(/'xyz','xyz'/), &
        compress_lvl=0, atts=dic)

    call delete(dic)

    ! Create matrix group
    call ncdf_def_grp(ncdf,'SPARSE',grp)

    call ncdf_def_dim(grp,'nnzs',n_nzs)

    ! EDM is required to have its own spin (because of spin-orbit coupling
    ! where the spin == 4, and not 8)
    call ncdf_def_dim(grp,'spin_EDM',spin%EDM)


    dic = dic//('info'.kv.'Index of supercell coordinates')
    call ncdf_def_var(grp,'isc_off',NF90_INT,(/'xyz','n_s'/), &
        compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Number of non-zero elements per row')
    call ncdf_def_var(grp,'n_col',NF90_INT,(/'no_u'/), &
        compress_lvl=0, atts=dic)

    chks = (/n_nzs,1,1/)

    dic = dic//('info'.kv. 'Supercell column indices in the sparse format')
    call ncdf_def_var(grp,'list_col',NF90_INT,(/'nnzs'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    ! Create the overlap matrix (we know it will not change)
    dic = ('info'.kv.'Overlap matrix')
    call ncdf_def_var(grp,'S',NF90_DOUBLE,(/'nnzs'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    dic = dic//('info'.kv.'Density matrix')
    call ncdf_def_var(grp,'DM',NF90_DOUBLE,(/'nnzs','spin'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    dic = dic//('info'.kv.'Energy density matrix')
    dic = dic//('unit'.kv.'Ry')
    call ncdf_def_var(grp,'EDM',NF90_DOUBLE,(/'nnzs    ','spin_EDM'/), &
        compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)

    if ( savehs .or. TS_HS_save ) then
      ! Unit is already present in dictionary
      dic = dic//('info'.kv.'Hamiltonian')
      call ncdf_def_var(grp,'H',NF90_DOUBLE,(/'nnzs','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    ! Even though I think we could do without, I add the
    ! xij array to the file
    ! Note that xij can be re-created using
    !    nsc, isc_off and xa
    !    if ( .not. Gamma ) then
    !       dic = dic//('info'.kv. &
    !            'Distance between orbital i and j')
    !       dic = dic//('unit'.kv.'Bohr')
    !       call ncdf_def_var(grp,'xij',NF90_DOUBLE,(/'xyz ','nnzs'/), &
    !            compress_lvl=cdf_comp_lvl,atts=dic)
    !    end if

    ! Delete the dictionary
    call delete(dic)

    ! Create grid group
    call ncdf_def_grp(ncdf,'GRID',grp)

    ! Note that there are now 2 spin dimensions
    ! 1. in the top level NC file and one in the GRID group
    ! This is because for spin-orbit coupling the grid and
    ! matrix dimensions are not the same (4 vs 8, respectively)
    d = ('DIMspin' .kv. spin%Grid)
    call ncdf_crt(grp, d)

    call delete(d)
    call delete(dic)

    call ncdf_def_grp(ncdf,'SETTINGS',grp)

    dic = ('info'.kv.'Tolerance for converging the density matrix')
    call ncdf_def_var(grp,'DMTolerance',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Tolerance for converging the Hamiltonian')
    call ncdf_def_var(grp,'HTolerance',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Net charge of the system')
    call ncdf_def_var(grp,'NetCharge',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Mixing weight')
    call ncdf_def_var(grp,'MixingWeight',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Grid used for the Brillouin zone integration')
    call ncdf_def_var(grp,'BZ',NF90_INT,(/'xyz','xyz'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Grid displacement used in Brillouin zone') &
         //('unit'.kv.'b**-1')
    call ncdf_def_var(grp,'BZ_displ',NF90_DOUBLE,(/'xyz'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Temperature for electrons')// &
         ('unit'.kv.'Ry')
    call ncdf_def_var(grp,'ElectronicTemperature',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    dic = dic//('info'.kv.'Mesh cutoff for real space grid')
    call ncdf_def_var(grp,'MeshCutoff',NF90_DOUBLE,(/'one'/), &
         compress_lvl=0, atts=dic)

    ! Create FC group
    if ( lis_FC ) then
      
      call ncdf_def_grp(ncdf,'FC',grp)

      ! Create dimension arrays
      ! Create arrays for containers
      allocate(ibuf(ia2-ia1+1))
      do i = ia1, ia2
        ibuf(i-ia1+1) = i
      end do
      call ncdf_def_dim(grp,'na_fc',ia2-ia1+1)
      call ncdf_def_dim(grp,'m_p',2)

      dic = dic//('info'.kv.'Displacement length')//('unit'.kv.'Bohr')
      call ncdf_def_var(grp,'disp',NF90_DOUBLE,(/'one'/), &
          compress_lvl=0, atts=dic)
      call ncdf_put_var(grp,'disp',dx)

      call delete(dic)
      dic = dic//('info'.kv.'Displaced atoms')
      call ncdf_def_var(grp,'ia_fc',NF90_INT,(/'na_fc'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)

      call ncdf_put_var(grp,'ia_fc',ibuf)
      deallocate(ibuf)

      dic = dic//('info'.kv.'Undisplaced atomic forces')//('unit'.kv.'Ry/Bohr')
      call ncdf_def_var(grp,'fa0',NF90_DOUBLE,(/'xyz ','na_u'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)

      dic = dic//('info'.kv.'Displaced atomic forces')//('unit'.kv.'Ry/Bohr')
      call ncdf_def_var(grp,'fa',NF90_DOUBLE,(/'xyz  ','na_u ','m_p  ','xyz  ', 'na_fc'/), &
          compress_lvl=cdf_comp_lvl, atts=dic)
      
    end if

    
    call delete(dic)

    if ( isolve == SOLVE_TRANSI ) then

       ! Save all information about the transiesta method
       call ncdf_def_grp(ncdf,'TRANSIESTA',grp)

       dic = ('info'.kv.'Grid used for the Brillouin zone integration')
       call ncdf_def_var(grp,'BZ',NF90_INT,(/'xyz','xyz'/), &
            compress_lvl=0, atts=dic)

       dic = dic//('info'.kv.'Grid displacement used in Brillouin zone') &
            //('unit'.kv.'b**-1')
       call ncdf_def_var(grp,'BZ_displ',NF90_DOUBLE,(/'xyz'/), &
            compress_lvl=0, atts=dic)

       dic = dic//('info'.kv.'Applied voltage')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'Volt',NF90_DOUBLE,(/'one'/), &
            compress_lvl=0,atts=dic)
       call delete(dic)

       ! Add all the electrodes
       do iEl = 1 , N_Elec
          
          call ncdf_def_grp(grp,trim(Elecs(iEl)%name),grp2)

          tmp = TotUsedAtoms(Elecs(iEl))
          call ncdf_def_dim(grp2,'na',tmp)

          dic = ('info'.kv.'Atoms belonging to electrode')
          call ncdf_def_var(grp2,'a_idx',NF90_INT,(/'na'/), &
               compress_lvl=0,atts=dic)
          
          dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
          call ncdf_def_var(grp2,'mu',NF90_DOUBLE,(/'one'/), &
               compress_lvl=0,atts=dic)

          dic = dic//('info'.kv.'Electronic temperature')//('unit'.kv.'Ry')
          call ncdf_def_var(grp2,'kT',NF90_DOUBLE,(/'one'/), &
               compress_lvl=0,atts=dic)

          call delete(dic)

       end do

    end if

    ! Save all things necessary here
    dic = ('time'.kv.datestring())
    dic = dic//('name'.kv.trim(sname))
    dic = dic//('label'.kv.trim(slabel))

    if ( isolve == SOLVE_DIAGON ) then
       dic = dic//('method'.kv.'diagon')
    else if ( isolve == SOLVE_ORDERN ) then
       dic = dic//('method'.kv.'order-n')
    else if ( isolve == SOLVE_MINIM ) then
       dic = dic//('method'.kv.'omm')
    else if ( isolve == SOLVE_TRANSI ) then
       dic = dic//('method'.kv.'transiesta')
    else if ( isolve == SOLVE_PEXSI ) then
       dic = dic//('method'.kv.'pexsi')
    end if

    ! Attributes are collective
    call ncdf_put_gatt(ncdf,atts=dic)

    call delete(dic)

    call ncdf_put_var(ncdf,'lasto',lasto(1:na_u))

    if ( isolve == SOLVE_TRANSI ) then

       ! Save all information about the transiesta method
       call ncdf_open_grp(ncdf,'TRANSIESTA',grp)

       call ncdf_put_var(grp,'Volt',Volt)
       
       ! Add all the electrodes
       do iEl = 1 , N_Elec
          
          call ncdf_open_grp(grp,trim(Elecs(iEl)%name),grp2)

          tmp = TotUsedAtoms(Elecs(iEl))

          allocate(ibuf(tmp))
          do i = 1 , tmp
             ibuf(i) = Elecs(iEl)%idx_a + i - 1
          end do
          call ncdf_put_var(grp2,'a_idx',ibuf)
          deallocate(ibuf)
          
          call ncdf_put_var(grp2,'mu',Elecs(iEl)%mu%mu)
          call ncdf_put_var(grp2,'kT',Elecs(iEl)%mu%kT)

       end do

    end if

    ! Close the file
    call ncdf_close(ncdf)

  end subroutine cdf_init_file

  subroutine cdf_init_grid(fname, ntm)
    use fdf, only : fdf_get, leqi
    use siesta_options, only: saverho, savedrho, savevh, savevna
    use siesta_options, only: savevt, savepsch, savetoch, saverhoxc
    use siesta_options, only: savebader
    use siesta_options, only: save_initial_charge_density

    character(len=*), intent(in) :: fname
    integer, intent(in) :: ntm(3)

    ! Local variables
    type(hNCDF) :: ncdf, grp
    integer :: prec, chks(3)
    character(len=64) :: key
    type(dictionary_t) :: dic

    ! We always re-write the file...
    call ncdf_open(ncdf,fname, &
        mode=ior(NF90_WRITE,NF90_NETCDF4))

    ! Create grid group
    call ncdf_open_grp(ncdf,'GRID', grp)

    call ncdf_def_dim(grp,'nx',ntm(1))
    call ncdf_def_dim(grp,'ny',ntm(2))
    call ncdf_def_dim(grp,'nz',ntm(3))

    ! Create the grid functions...

    ! all grids are using the grid_p precision
    if ( grid_p == dp ) then
      prec = NF90_DOUBLE
      ! In case the user thinks the double precision
      ! is needed, but only single precision
      ! is needed saving, we allow that.
      key = fdf_get('CDF.Grid.Precision','double')
      if ( leqi(key,'single') ) prec = NF90_FLOAT
      if ( leqi(key,'float') )  prec = NF90_FLOAT
    else
      ! The grid is in single precision, so
      ! we save it in that precision.
      prec = NF90_FLOAT
    end if
    
    chks = (/ntm(1),ntm(2),1/)

    ! Define the units for all charge variables
    dic = ('unit'.kv.'e/Bohr**3')

    if ( save_initial_charge_density ) then
      dic = dic//('info'.kv.'Initial charge density')
      call ncdf_def_var(grp,'RhoInit',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( saverho ) then
      dic = dic//('info'.kv.'Charge density')
      call ncdf_def_var(grp,'Rho',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savepsch ) then
      dic = dic//('info'.kv.'Diffuse ionic charge')
      call ncdf_def_var(grp,'Chlocal',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savetoch ) then
      dic = dic//('info'.kv.'Total charge')
      call ncdf_def_var(grp,'RhoTot',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savedrho ) then
      dic = dic//('info'.kv.'Density difference from atomic densities')
      call ncdf_def_var(grp,'RhoDelta',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( saverhoxc ) then
      dic = dic//('info'.kv.'Density used to calculate XC functional')
      call ncdf_def_var(grp,'RhoXC',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savebader ) then
      dic = dic//('info'.kv.'Bader charge')
      call ncdf_def_var(grp,'RhoBader',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    ! Define the units for the potential grids
    dic = dic//('unit'.kv.'Ry')

    if ( savevna ) then
      dic = dic//('info'.kv.'Neutral atom potential')
      call ncdf_def_var(grp,'Vna',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savevh ) then
      dic = dic//('info'.kv.'Electrostatic potential')
      call ncdf_def_var(grp,'Vh',prec,(/'nx','ny','nz'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    if ( savevt ) then
      dic = dic//('info'.kv.'Total potential')
      call ncdf_def_var(grp,'Vt',prec,(/'nx  ','ny  ','nz  ','spin'/), &
          compress_lvl=cdf_comp_lvl,atts=dic,chunks=chks)
    end if

    call delete(dic)

    ! Close the file
    call ncdf_close(ncdf)

  end subroutine cdf_init_grid

  subroutine cdf_save_settings(fname)

    use kpoint_grid,    only: kscell, kdispl
    use siesta_options, only: dDtol, dHtol, charnet, wmix, temp, g2cut
    use siesta_options, only: isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN, SOLVE_TRANSI
    use m_ts_kpoints,   only: ts_kscell, ts_kdispl

    character(len=*), intent(in) :: fname
    
    type(hNCDF) :: ncdf, grp

    ! We just open it (prepending)
    call ncdf_open(ncdf,fname, &
        mode=ior(NF90_WRITE,NF90_NETCDF4))

    call ncdf_open_grp(ncdf,'SETTINGS',grp)

    ! Save settings
    call ncdf_put_var(grp,'BZ',kscell)
    call ncdf_put_var(grp,'BZ_displ',kdispl)
    call ncdf_put_var(grp,'DMTolerance',dDtol)
    call ncdf_put_var(grp,'HTolerance',dHtol)
    call ncdf_put_var(grp,'NetCharge',charnet)
    call ncdf_put_var(grp,'MixingWeight',wmix)
    call ncdf_put_var(grp,'ElectronicTemperature',Temp)
    call ncdf_put_var(grp,'MeshCutoff',g2cut)

    ! I suggest that all MD settings are 
    ! put in the MD-group

    if ( isolve == SOLVE_TRANSI ) then
       call ncdf_open_grp(ncdf,'TRANSIESTA',grp)

       call ncdf_put_var(grp,'BZ',ts_kscell)
       call ncdf_put_var(grp,'BZ_displ',ts_kdispl)

    end if

    call ncdf_close(ncdf)

  end subroutine cdf_save_settings


  subroutine cdf_save_state(fname,dic_save)
    !    use m_gamma, only : Gamma
    use m_energies, only: Ef, Efs
    use atomlist, only : Qtot
    use siesta_options, only: fixspin, total_spin
    use siesta_geom, only: na_u, ucell, xa, va
    use siesta_geom, only: nsc, isc_off
    use sparse_matrices, only: sparse_pattern, block_dist
    use sparse_matrices, only: S_1D
    use sparse_matrices, only: DM_2D, EDM_2D, H_2D
!    use sparse_matrices, only: xij_2D
    use m_stress, only : stress
    use m_forces, only: fa

    character(len=*), intent(in) :: fname
    ! Dictionary containing keys that we will save
    type(dictionary_t), intent(in) :: dic_save
    type(hNCDF) :: ncdf, grp

    call timer('CDF',1)

    ! We just open it (prepending)
#ifdef MPI
    if ( cdf_w_parallel ) then
       call ncdf_open(ncdf,fname, &
            mode=ior(NF90_WRITE,NF90_MPIIO), comm=MPI_Comm_World)
    else
#endif
       call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4))
#ifdef MPI
    end if
#endif

    ! Add attribute of current Fermi-level
    if ( fixspin ) then
      if ( ('Ef' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Ef',Efs)
      if ( ('Qspin' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Qspin',total_spin)
    else
      if ( ('Ef' .in. dic_save) .and. Node == 0 ) &
          call ncdf_put_var(ncdf,'Ef',Ef)
    end if
    if ( ('Qtot' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'Qtot',Qtot)
    ! Save nsc, xa, fa, lasto, ucell
    if ( ('nsc' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'nsc',nsc)
    if ( ('cell' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'cell',ucell)
    if ( ('xa' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'xa',xa(:,1:na_u))
    if ( ('va' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'va',va(:,1:na_u))
    if ( ('fa' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'fa',fa(:,1:na_u))
    if ( ('stress' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'stress',stress)

    ! Sparsity format
    call ncdf_open_grp(ncdf,'SPARSE',grp)

    if ( 'isc_off'.in. dic_save) &
         call ncdf_put_var(grp,'isc_off',isc_off)
    if ( 'sp' .in. dic_save ) &
         call cdf_w_Sp(grp,block_dist,sparse_pattern)
    if ( 'S' .in. dic_save ) &
         call cdf_w_d1D(grp,'S',S_1D)
!    if ( .not. Gamma .and. ('xij' .in. dic_save) ) then
       ! Write the xij array, it will not change during SCF
!       call cdf_w_d2D(grp,'xij',xij_2D)
    !    end if
    if ( 'H' .in. dic_save ) &
         call cdf_w_d2D(grp,'H',H_2D)
    if ( 'DM' .in. dic_save ) &
         call cdf_w_d2D(grp,'DM',DM_2D)
    if ( 'EDM' .in. dic_save ) &
         call cdf_w_d2D(grp,'EDM',EDM_2D)

    call ncdf_close(ncdf)
    
    call timer('CDF',2)

  end subroutine cdf_save_state

  subroutine cdf_save_grid(fname,vname,nspin,nmeshl,grid)

    character(len=*), intent(in) :: fname, vname
    integer, intent(in) :: nspin, nmeshl(3)
    real(grid_p), intent(in) :: grid(product(nmeshl),nspin)

    type(hNCDF) :: ncdf
    integer :: is

    call timer('CDF-grid',1)

    ! We just open it (prepending)
#ifdef MPI
    if ( cdf_w_parallel ) then
       call ncdf_open(ncdf,fname, group='GRID', &
            mode=ior(NF90_WRITE,NF90_MPIIO), &
            comm=MPI_Comm_World)
    else
#endif
       call ncdf_open(ncdf,fname, group='GRID', &
            mode=ior(NF90_WRITE,NF90_NETCDF4))
#ifdef MPI
    end if
#endif

    ! Save the grid
    if ( nspin > 1 ) then
       do is = 1 , nspin 
          call cdf_w_grid(ncdf,vname,nmeshl,grid(:,is),idx=is)
       end do
    else
       call cdf_w_grid(ncdf,vname,nmeshl,grid(:,1))
    end if

    call ncdf_close(ncdf)

    call timer('CDF-grid',2)

  end subroutine cdf_save_grid

  subroutine cdf_save_basis(fname)

    use siesta_geom, only: na_u, isa

    use atmparams, only : nt => NTBMAX
    use atm_types, only : species_info, species, nspecies
    use radial, only : rad_func

    character(len=*), intent(in) :: fname

    ! Local variables
    type(species_info), pointer :: spp
    type(rad_func), pointer :: p

    type(hNCDF) :: nf, ncdf, grp
    type(dictionary_t) :: dic, d
    character(len=DICTIONARY_KEY_LENGTH) :: key
    type(variable_t) :: v
    integer :: is, i

    ! Used for saving variables
    integer :: no, nk
    integer, allocatable :: iaux(:)
    real(dp), allocatable :: raux(:)


    ! Unluckily the new basis saves only
    ! saved on the IO node. 
    ! No other node must therefore access this routine
    if ( Node /= 0 ) return

    call timer('CDF-basis',1)

    call ncdf_open(nf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4))

    ! create the BASIS group
    call ncdf_def_grp(nf,'BASIS',ncdf)

    ! Create a list of the species associated with each atom
    dic = ('info'.kv.'Basis of each atom by ID')
    call ncdf_def_var(ncdf,'basis',NF90_INT,(/'na_u'/), atts= dic)
    call delete(dic)
    call ncdf_put_var(ncdf,'basis',isa(1:na_u))

    do is = 1 , nspecies

       ! Get current specie
       spp => species(is)

       ! Get array sizes
       no = spp%n_orbnl
       nk = spp%n_pjnl
       if ( no == 0 ) cycle

       ! Create the group
       call ncdf_def_grp(ncdf,trim(spp%label),grp)

       call ncdf_def_dim(grp,'norbs',no)
       if ( nk > 0 ) then
          call ncdf_def_dim(grp,'nkbs',nk)
       end if
       call ncdf_def_dim(grp,'ntb',nt)

       ! Save the orbital global attributes
       dic = ('Element'.kv.trim(spp%symbol))
       dic = dic//('Label'.kv.trim(spp%label))
       dic = dic//('Atomic_number'.kv.spp%z)
       dic = dic//('Valence_charge'.kv.spp%zval)
       dic = dic//('Mass'.kv.spp%mass)
       dic = dic//('Self_energy'.kv.spp%self_energy)
       dic = dic//('Number_of_orbitals'.kv.spp%norbs)
       dic = dic//('L_max_basis'.kv.spp%lmax_basis)
       dic = dic//('Number_of_projectors'.kv.spp%nprojs)
       dic = dic//('L_max_projs'.kv.spp%lmax_projs)
       dic = dic//('ID'.kv.is)
       call ncdf_put_gatt(grp,atts=dic)
       call delete(dic)

       ! Create all orbital variables...
       dic = ('orbnl_l'.kv.NF90_INT)//('orbnl_n'.kv.NF90_INT)
       dic = dic//('orbnl_z'.kv.NF90_INT)//('orbnl_ispol'.kv.NF90_INT)
       dic = dic//('orbnl_pop'.kv.NF90_DOUBLE)//('cutoff'.kv.NF90_DOUBLE)
       dic = dic//('delta'.kv.NF90_DOUBLE)
       d = .first. dic
       do while ( .not. (.empty. d) )
          key = .key. d
          v = .val. d
          call assign(i,v)
          call ncdf_def_var(grp,trim(key),i,(/'norbs'/), &
               compress_lvl=0,chunks=(/no/))
          d = .next. d
       end do
       call delete(dic)

       ! Create all projector variables...
       if ( nk > 0 ) then
          call ncdf_def_var(grp,'proj',NF90_DOUBLE, (/'ntb ','nkbs'/), &
             compress_lvl=cdf_comp_lvl,chunks=(/nt,1/))

          dic = ('pjnl_l'.kv.NF90_INT)//('pjnl_j'.kv.NF90_DOUBLE)
          dic = dic//('pjnl_n'.kv.NF90_INT)//('pjnl_ekb'.kv.NF90_DOUBLE)
          dic = dic//('kbcutoff'.kv.NF90_DOUBLE)//('kbdelta'.kv.NF90_DOUBLE)
          d = .first. dic
          do while ( .not. (.empty. d) )
             key = .key. d
             v = .val. d
             call assign(i,v)
             call ncdf_def_var(grp,trim(key),i,(/'nkbs'/), &
                  compress_lvl=0,chunks=(/nk/))
             d = .next. d
          end do
          call delete(dic)
       end if
       call delete(v)

       ! Create orbital projector
       call ncdf_def_var(grp,'orb',NF90_DOUBLE,(/'ntb  ','norbs'/), &
           compress_lvl=cdf_comp_lvl,chunks=(/nt,1/))

       if ( spp%z > 0 ) then ! negative are floating orbitals
       
         ! Local potential
         call save_rad_func('vna', spp%vna)
         
         ! Local potential charge density
         call save_rad_func('chlocal', spp%chlocal)
         
         ! Reduced local potential (rV+2*Zval)
         call save_rad_func('reduced_vlocal', spp%reduced_vlocal)

         if ( spp%there_is_core ) then
           ! Core charge, if it exists, the variable will be created
           ! Hence, the old way of designating whether core is present
           ! or not is removed.
           ! I.e. no Core_flag [1|0] will be saved, Core_flag == .true. is 
           ! apt if 'core' variable exists.
           call save_rad_func('core', spp%core)
         end if

       end if

       call delete(dic)

       ! Save all variables to the group

       ! Save orbital
       call ncdf_put_var(grp,'orbnl_l',spp%orbnl_l(1:no))
       call ncdf_put_var(grp,'orbnl_n',spp%orbnl_n(1:no))
       call ncdf_put_var(grp,'orbnl_z',spp%orbnl_z(1:no))

       allocate(iaux(no))
       do i = 1, no
          if ( spp%orbnl_ispol(i) ) then
             iaux(i) = 1
          else
             iaux(i) = 0 
          end if
       end do
       call ncdf_put_var(grp,'orbnl_ispol',iaux(1:no))
       deallocate(iaux)
       call ncdf_put_var(grp,'orbnl_pop',spp%orbnl_pop(1:no))

       ! Save projector
       if ( nk > 0 ) then
          call ncdf_put_var(grp,'pjnl_l',spp%pjnl_l(1:nk))
          allocate(raux(nk))
          do i = 1, nk
            raux(i) = 0._dp
          end do
          call ncdf_put_var(grp,'pjnl_j',raux)
          deallocate(raux)
          call ncdf_put_var(grp,'pjnl_n',spp%pjnl_n(1:nk))
          call ncdf_put_var(grp,'pjnl_ekb',spp%pjnl_ekb(1:nk))
       end if

       do i = 1, nk
          p => spp%pjnl(i)
          call ncdf_put_var(grp,'proj',p%f(1:nt),start=(/1,i/))
          call ncdf_put_var(grp,'kbcutoff',p%cutoff,start=(/i/))
          call ncdf_put_var(grp,'kbdelta',p%delta,start=(/i/))
       end do

       do i = 1, no
          p => spp%orbnl(i)
          call ncdf_put_var(grp,'orb',p%f(1:nt),start=(/1,i/))
          call ncdf_put_var(grp,'cutoff',p%cutoff,start=(/i/))
          call ncdf_put_var(grp,'delta',p%delta,start=(/i/))
       end do

    end do

    call ncdf_close(nf)

    call timer('CDF-basis',2)

  contains

    subroutine save_rad_func(name, rfunc)
      use radial, only: rad_func

      character(len=*), intent(in) :: name
      type(rad_func), intent(in) :: rfunc
      type(dictionary_t) :: dic

      ! Only create it if it exists in the pseudo
      if ( rfunc%n <= 0 ) return

      dic = ('cutoff'.kv.rfunc%cutoff) // ('delta'.kv.rfunc%delta)
      call ncdf_def_var(grp, name, NF90_DOUBLE, (/'ntb'/), atts=dic, &
          compress_lvl=cdf_comp_lvl,chunks=(/nt/))
      ! Save immediately
      call ncdf_put_var(grp, name, rfunc%f(1:nt))
      call delete(dic)

    end subroutine save_rad_func

  end subroutine cdf_save_basis

  subroutine cdf_save_fc(fname,istep)
    use siesta_geom, only: na_u
    use m_forces, only: fa

    character(len=*), intent(in) :: fname
    integer, intent(in) :: istep

    type(hNCDF) :: ncdf
    integer :: ia, dir, pm

    if ( Node /= 0 ) return

    ! open the file...
    call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4),group='FC')
    
    if ( istep == 0 ) then
      call ncdf_put_var(ncdf,'fa0',fa)
      
    else
      ! Atomic index (with respect to ia1)
      ia = mod(istep - 1, 6) + 1
      dir = (ia - 1) / 2 + 1
      pm = mod(ia-1, 2) + 1
      ia = (istep - 1) / 6 + 1
      
      call ncdf_put_var(ncdf,'fa',fa,count=(/3,na_u/),start=(/1,1,pm,dir,ia/))
    end if
    
    call ncdf_close(ncdf)

  end subroutine cdf_save_fc
  
#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_siesta
