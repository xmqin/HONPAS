! This module has been fully created by:
! Nick Papior Andersen, nickpapior@gmail.com
! It is to be considered copyrighted and ONLY to be used for scientific purposes.
! Every dimension ID, variable ID and NetCDF file ID will now be tracked by
! this module. This will let you track all variables by name.
! Furthermore, the module will be developed to handle parallel writes using
! the MPI-communication layer without the user "knowing". One simply opens/creates
! the NetCDF file in parallel mode.
!
! Furthermore the inquiries, and returning variables has been abstracted to their
! type and size. Also all references should be done with names of variables, dimensions
! etc.
! Example. You have created a NetCDF file with the following variable:
! Energy(MD): double
! If you need to read the variables content in 10:13 you will read by:
! ncdf_get_var(ncdf,"Energy",E(1:4),start=(/10/))
!
! A specialized ncdf_[put|get]_var has been developed for complex variables.
! In ordinary NetCDF contents one has to define two variables, do a conversion
! for real and imaginary parts, etc.
! In this case we have created:
! ncdf_get_var(ncdf,"ComplexEnergy",CE(1:4),start=(/10/))
! to read in the variables:
! "ReComplexEnergy" and "ImComplexEnergy".
! This heavily relieves the burden of tracking content and imaginary double IDs for the
! NetCDF interface.
!
! NOTE: this limits the fetch of a COMPLEX variable, automatically converted to an
! integer. The programmer needs to do this by other means.
!
! The routines for variables/attributes/dimensions are the following.
! Here an asterisk (*) denotes an optional argument. ALWAYS refer to
! optional arguments by keywords. ALWAYS.
! - ncdf_inq_var
! (ncdf,name,*exist,*id,*size)
! - ncdf_get_var
! (ncdf,name,*exist,*start,*count,*stride)
! - ncdf_put_var
! (ncdf,name,var,*start,*count,*stride)
! - ncdf_inq_dim
! (ncdf,name,*exist,*id,*len)
! - ncdf_inq_att
! (ncdf,var,name,*exist)
! - ncdf_get_att
! (ncdf,var,name,att<int,real,char>)
! - ncdf_inq_gatt
! (ncdf,name,*exist)
! - ncdf_get_gatt
! (ncdf,var,name,att<int,real,char>)
! - ncdf_def_var
! (ncdf,var,type,dims,
! *atts=(dictionary),*compress_lvl,*fill) (the fill variable is not implemented in NetCDF for now)
! A wrapper module for doing netcdf operations
! The idea is that this module should be able to do parallel IO when needed
! Currently it does not have this implemented, but it provides a wrapper basis
! which makes it easy.
module netcdf_ncdf
  ! Globalize the variables needed for generating the requested features
  use netcdf
  implicit none
  public
  ! Privatize used variables
  integer, private, parameter :: ih = selected_int_kind(4)
  integer, private, parameter :: is = selected_int_kind(9)
  integer, private, parameter :: il = selected_int_kind(18)
  integer, private, parameter :: sp = selected_real_kind(p=6)
  integer, private, parameter :: dp = selected_real_kind(p=15)
  ! The IONode setting
  logical, save :: IONode = .true.
  private :: IONode
  ! Local routines
  private :: ncdf_def_var_generic
  private :: cat_char_ncdf
  private :: cat_ncdf_char
  private :: ncdf_die
  private :: h_reset
  ! We add a specific NetCDF handle for dealing with files
  type :: hNCDF
     ! The top file id (in case Groups are used)
     ! The file-handle for the netCDF-file
     integer :: f_id = -1
     integer :: id = -1
     ! Whether the file is handled parallely
     logical :: parallel = .false.
     ! The mode of the file
     integer :: mode
     ! If define < 0, then it is a netCDF-4 file, enddef will only
     ! be called if ncdf_enddef is called
     ! If define == 0 then it is in define mode (needed for netCDF-3)
     ! If define == 1 then it is in data mode (needed for netCDF-3)
     integer :: define
     ! The name of the netCDF-file
     character(len=256) :: name = " "
     ! the group of the netCDF-file (i.e. a file within a file)
     character(len=NF90_MAX_NAME) :: grp = " "
     ! The communicator describing the parallel activity
     integer :: comm = - 1
     ! Default compression level
     integer :: comp_lvl = 0
  end type hNCDF
  ! Interface the concatenation
  interface operator (//)
     module procedure cat_char_ncdf
     module procedure cat_ncdf_char
  end interface operator (//)
  ! Interface
  interface parallel
     module procedure parallel_
  end interface parallel
  private :: parallel_
  interface parallel_io
     module procedure parallel_io_
  end interface parallel_io
  private :: parallel_io_
  ! Interface for creating variables
  interface ncdf_def_var
     module procedure ncdf_def_var_integer
     module procedure ncdf_def_var_logical
  end interface ncdf_def_var
  private :: ncdf_def_var_integer
  private :: ncdf_def_var_logical
  ! Interface for acquiring information about a file...
  interface ncdf_inq
    module procedure ncdf_inq_ncdf
    module procedure ncdf_inq_name
  end interface ncdf_inq
  private :: ncdf_inq_name
  private :: ncdf_inq_ncdf
  ! Add new data types
  ! We need them to be logical due to the interface of the def_var.
  ! Otherwise they would have the same interface (due to the optional argument var_id)
  logical, parameter :: NF90_DOUBLE_COMPLEX = .true. ! for true it is double
  logical, parameter :: NF90_FLOAT_COMPLEX = .false. ! for false it is float
  ! This variable is used to enable the NOFILL on certain variables (for the interface to ncdf_def_var)
  logical, parameter :: NF90_VAR_NOFILL = .true. ! for false it is float
  ! This variable is used to re-enable the FILL on certain variables (for the interface to ncdf_def_var)
  logical, parameter :: NF90_VAR_FILL = .false. ! for false it is float
  ! Added interface
interface ncdf_put_var
module procedure put_var_h0_name
module procedure put_var_h1_name
module procedure put_var_h2_name
module procedure put_var_h3_name
module procedure put_var_s0_name
module procedure put_var_s1_name
module procedure put_var_s2_name
module procedure put_var_s3_name
module procedure put_var_d0_name
module procedure put_var_d1_name
module procedure put_var_d2_name
module procedure put_var_d3_name
module procedure put_var_c0_name
module procedure put_var_c1_name
module procedure put_var_c2_name
module procedure put_var_c3_name
module procedure put_var_z0_name
module procedure put_var_z1_name
module procedure put_var_z2_name
module procedure put_var_z3_name
module procedure put_var_i0_name
module procedure put_var_i1_name
module procedure put_var_i2_name
module procedure put_var_i3_name
end interface ncdf_put_var
interface ncdf_get_var
module procedure get_var_h0_name
module procedure get_var_h1_name
module procedure get_var_h2_name
module procedure get_var_h3_name
module procedure get_var_s0_name
module procedure get_var_s1_name
module procedure get_var_s2_name
module procedure get_var_s3_name
module procedure get_var_d0_name
module procedure get_var_d1_name
module procedure get_var_d2_name
module procedure get_var_d3_name
module procedure get_var_c0_name
module procedure get_var_c1_name
module procedure get_var_c2_name
module procedure get_var_c3_name
module procedure get_var_z0_name
module procedure get_var_z1_name
module procedure get_var_z2_name
module procedure get_var_z3_name
module procedure get_var_i0_name
module procedure get_var_i1_name
module procedure get_var_i2_name
module procedure get_var_i3_name
end interface ncdf_get_var
interface ncdf_put_gatt
module procedure put_gatt
module procedure put_gatt_s0
module procedure put_gatt_s1
module procedure put_gatt_d0
module procedure put_gatt_d1
module procedure put_gatt_i0
module procedure put_gatt_i1
end interface ncdf_put_gatt
interface ncdf_put_att
module procedure put_att
module procedure put_att_s0
module procedure put_att_s1
module procedure put_att_d0
module procedure put_att_d1
module procedure put_att_i0
module procedure put_att_i1
end interface ncdf_put_att
interface ncdf_get_gatt
module procedure get_gatt
module procedure get_gatt_s0
module procedure get_gatt_s1
module procedure get_gatt_d0
module procedure get_gatt_d1
module procedure get_gatt_i0
module procedure get_gatt_i1
end interface ncdf_get_gatt
interface ncdf_get_att
module procedure get_att
module procedure get_att_s0
module procedure get_att_s1
module procedure get_att_d0
module procedure get_att_d1
module procedure get_att_i0
module procedure get_att_i1
end interface ncdf_get_att
interface ncdf_inq_var
module procedure ncdf_inq_var_def
module procedure inq_var_h0
module procedure inq_var_s0
module procedure inq_var_d0
module procedure inq_var_c0
module procedure inq_var_z0
module procedure inq_var_i0
end interface ncdf_inq_var
interface ncdf_def_fill
module procedure def_fill_h0
module procedure def_fill_s0
module procedure def_fill_d0
module procedure def_fill_c0
module procedure def_fill_z0
module procedure def_fill_i0
end interface ncdf_def_fill
contains
  ! Every routine in this module needs NetCDF
  ! So it is sourrounded by this...
  function parallel_(this) result(par)
    type(hNCDF), intent(in) :: this
    logical :: par
    par = this%parallel
  end function parallel_
  function parallel_io_(this) result(par)
    type(hNCDF), intent(in) :: this
    logical :: par
    par = iand(NF90_MPIIO,this%mode) == NF90_MPIIO
    if ( par ) return
    par = iand(NF90_MPIPOSIX,this%mode) == NF90_MPIPOSIX
    if ( par ) return
    par = iand(NF90_PNETCDF,this%mode) == NF90_PNETCDF
  end function parallel_io_
  subroutine ncdf_copy(this,copy)
    type(hNCDF), intent(in) :: this
    type(hNCDF), intent(out) :: copy
    copy = this
  end subroutine ncdf_copy
  ! Local subroutine to reset the ncdf handle
  subroutine h_reset(this)
    type(hNCDF), intent(inout) :: this
    this%f_id = -1
    this%id = -1
    this%parallel = .false.
    this%define = 0
    this%mode = 0
    this%name = " "
    this%grp = " "
    this%comm = -1
    this%comp_lvl = 0
  end subroutine h_reset
  subroutine ncdf_init(this,name,mode,parallel,comm,overwrite,compress_lvl)
    type(hNCDF), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    integer, optional, intent(in) :: mode
    logical, optional, intent(in) :: parallel
    integer, optional, intent(in) :: comm
    logical, optional, intent(in) :: overwrite
    integer, optional, intent(in) :: compress_lvl
    integer :: format
    logical :: exist
    call h_reset(this)
    if ( present(name) ) this%name = name
    if ( present(compress_lvl) ) this%comp_lvl = compress_lvl
    exist = .false.
    if ( present(mode) ) exist = .true.
    if ( present(comm) ) then
      if ( comm >= 0 ) exist = .true.
    end if
    if ( present(parallel) ) then
      if ( parallel ) exist = .false.
    end if
    if ( .not. exist ) then
      ! Our default is the 64 bit offset files...
      ! The best backwards compatibility format
      this%mode = IOR(this%mode,NF90_64BIT_OFFSET)
    end if
    ! set the mode
    if ( present(mode) ) then
      this%mode = mode
    end if
    ! This will create the correct order
    exist = .false.
    if ( present(parallel) ) then
      if ( parallel ) exist = .true.
    end if
    if ( present(comm) ) then
      if ( comm >= 0 ) exist = .false.
    end if
    if ( exist ) then
      this%parallel = parallel
      ! The parallel flag is for the sequential parallel access !
      ! This will be reset to a zero mode if a communicator is supplied
      ! In this way we can have "parallel" access for reading purposes...
      ! Check that the mode is not existing in the passed mode
      this%mode = ior(this%mode, NF90_SHARE)
    end if
    if ( present(comm) ) then
      if ( comm < 0 ) then
        ! If the communicator is negative we
        ! must assume that it is not parallel
        this%comm = comm
      end if
    end if
    ! Define whether it should not be clobbered
    inquire(file=""//this,exist=exist)
    if ( present(overwrite) ) then
      if ( exist .and. .not. overwrite ) then
        this%mode = ior(this%mode, NF90_NOCLOBBER)
      end if
    end if
  end subroutine ncdf_init
  subroutine ncdf_create(this,filename,mode,overwrite,parallel,comm, &
      compress_lvl)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, optional, intent(in) :: mode
    logical, optional, intent(in) :: overwrite
    logical, optional, intent(in) :: parallel
    integer, optional, intent(in) :: comm
    integer, optional, intent(in) :: compress_lvl
    integer :: file_format
    logical :: exist
    call ncdf_init(this,name=filename, &
        mode=mode,parallel=parallel,comm=comm, &
        overwrite=overwrite, &
        compress_lvl=compress_lvl)
    ! We need to correct the definition for netCDF-3 files
    inquire(file=this%name, exist=exist)
    if ( present(overwrite) ) then
      if ( overwrite ) then
        exist = .false.
      end if
    end if
    if ( iand(NF90_NETCDF4,this%mode) == NF90_NETCDF4 ) then
      this%define = -1
    end if
    if ( .not. ncdf_participate(this) ) return
    if ( exist ) then
      call ncdf_die("File: "//this//" already exists! "//&
          "Please delete the file (or request overwriting).")
    end if
    if ( this%parallel .and. this%comm >= 0 ) then
      call ncdf_err(-100,"Not compiled with communicater parallel")
    else if ( this%parallel ) then
      call ncdf_err(nf90_create(filename, this%mode , this%f_id), &
          "Creating file: "//this//" in parallel")
    else
      call ncdf_err(nf90_create(filename, this%mode , this%f_id), &
          "Creating file: "//this)
    end if
    ! We could check for mode == NF90_SHARE in case of parallel...
    ! However, it does not make sense as the other is still correct, just slow
    this%id = this%f_id
  end subroutine ncdf_create
  subroutine ncdf_open(this,filename,group,mode,parallel,comm,compress_lvl)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: filename
    character(len=*), optional, intent(in) :: group
    integer, optional, intent(in) :: mode
    logical, optional, intent(in) :: parallel
    integer, optional, intent(in) :: comm
    integer, optional, intent(in) :: compress_lvl
    integer :: file_format, i
    logical :: exist
    ! Save the general information which should be accesible to all processors
    call ncdf_init(this,name=filename,mode=mode, &
        parallel=parallel,comm=comm, &
        compress_lvl=compress_lvl)
    ! When we open a file, it will always be in data mode...
    if ( iand(NF90_NETCDF4,this%mode) == NF90_NETCDF4 ) then
      this%define = -1
    else
      this%define = 1
    end if
    if ( .not. ncdf_participate(this) ) return
    inquire(file=filename, exist=exist)
    if ( .not. exist ) then
      call ncdf_die("File: "//trim(filename)//" does not exist! "//&
          "Please check your inqueries.")
    end if
    ! If we have not added a mode it must be for non-writing purposes
    if ( .not. present(mode) ) then
      this%mode = IOR(this%mode,NF90_NOWRITE)
    end if
    if ( this%parallel .and. this%comm >= 0 ) then
      call ncdf_err(-100,"Code not compiled with NCDF_PARALLEL")
    else if ( this%parallel ) then
      call ncdf_err(nf90_open(filename, this%mode , this%f_id), &
          "Opening file: "//this//" in parallel")
    else
      call ncdf_err(nf90_open(filename, this%mode , this%f_id), &
          "Opening file: "//this)
    end if
    ! Copy so that we can create inquiry
    this%id = this%f_id
    if ( present(group) ) then
      this%grp = '/'//trim(group)
      call ncdf_err(nf90_inq_grp_full_ncid(this%f_id, trim(this%grp), this%id))
    end if
  end subroutine ncdf_open
  ! A specific routine for generating a NetCDF file
  ! group based on a dictionary input.
  recursive subroutine ncdf_crt(this,dic)
    use variable
    use dictionary
    type(hNCDF), intent(inout) :: this
    type(dictionary_t), intent(inout) :: dic
    ! Local dictionary keys and variables
    type(hNCDF) :: grp
    type(dictionary_t) :: d, d_var, dv, atts
    type(variable_t) :: v
    character(len=DICTIONARY_KEY_LENGTH) :: key
    character(len=NF90_MAX_NAME) :: name, char
    character(len=64), allocatable :: dims(:)
    integer, pointer :: chunks(:) => null()
    integer :: d_n, type, i, n_d, j
    ! Type declarations
    integer :: itype
    logical :: ltype
    if ( .not. ncdf_participate(this) ) return
    ! Create all dimensions
    d = .first. dic
    do while ( .not. (.empty. d) )
      key = .key. d
      if ( key(1:3) == 'DIM' ) then
        key = key(4:)
        ! We have a dimension
        ! { DIMname : <size> }
        ! 1. Get value
        call associate(v,d)
        ! 2. Assign value
        call assign(d_n,v)
        call ncdf_def_dim(this,key,d_n)
        ! 3. just nullify the variable, let the user delete
        ! the dictionary, we cannot know whether the
        ! variable is a pointer.
        call nullify(v)
      end if
      d = .next. d
    end do
    ! Create all variables
    d = .first. dic
    do while ( .not. (.empty. d) )
      key = .key. d
      if ( key(1:3) == 'VAR' ) then
        name = key(4:)
        !print *,'Creating variable: ',trim(key)
        ! We have a variable
        ! { VARname : <dict> } dict ->
        ! {
        ! dims : '1st,2nd,3rd', ! Comma-seperated list
        ! type : NF90_INT|NF90_FLOAT|..., ! Integer, or logical for complex
        ! atts : <dict>, ! Dictionary of attributes.
        ! chunks : <array of int>, ! For specifying chunk-sizes
        ! }
        call associate(d_var,d)
        !call print(d_var)
        ! We have the dictionary of the variable
        ! 0. in case 'name' exists, it must be the name
        if ( 'name'.in.d_var ) then
          call associate(v,d_var,'name')
          if ( which(v) /= 'a1' ) then
            call ncdf_err(-200, &
                'Name of variable is not a character variable.')
          end if
          name = ' '
          call assign(name,v)
        end if
        ! 1. Get the dimensions
        if ( 'dims'.nin. d_var ) then
          call ncdf_err(-200, &
              'Unable to retrieve the dimension &
              &key from a variable dictionary. A variable &
              &MUST be defined with a comma separated dimension.')
        end if
        !print *,'Retrieve dims (1): ',trim(key)
        call associate(v,d_var,'dims')
        ! The dimensions has to be given in a comma separated list
        if ( which(v) /= 'a1' ) then
          call ncdf_err(-200, &
              'Dimension variable is not a character variable.')
        end if
        ! Ensure we have a completely empty character.
        !print *,'Retrieve dims (2): ',trim(key)
        key = ' '
        call assign(key,v)
        ! Count number of dimensions
        n_d = 1
        do i = 1 , len_trim(key)
          if ( key(i:i) == ',' ) then
            n_d = n_d + 1
          end if
        end do
        allocate(dims(n_d))
        ! Copy over the dimensions
        j = 1
        n_d = 1
        do i = 2 , len_trim(key)
          if ( key(i:i) == ',' ) then
            dims(n_d) = ' '
            dims(n_d) = trim(adjustl(key(j:i-1)))
            n_d = n_d + 1
            j = i + 1
          end if
        end do
        ! Grab the last dimension
        dims(n_d) = trim(adjustl(key(j:)))
        ! Figure out the type
        if ( 'type' .nin. d_var ) then
          call ncdf_err(-200, &
              'Unable to retrieve the type &
              &key from a variable dictionary. A variable &
              &MUST have a clear type.')
        end if
        !print *,'Retrieve type: ',trim(key)
        call associate(v,d_var,'type')
        char = which(v)
        if ( trim(char) /= 'b0' .and. &
            trim(char) /= 'i0' ) then
          call ncdf_err(-200, &
              'Type of variable is not defined with a &
              &proper variable designator.')
        end if
        if ( trim(char) == 'i0' ) then
          call assign(itype,v)
        else
          call assign(ltype,v)
        end if
        ! Check if there are any attributes in the dictionary
        if ( 'atts' .in. d_var ) then
          call associate(atts,d_var,'atts')
        end if
        ! We have now gathered all information.
        ! Lets create that variable, alrighty! :)
        if ( trim(char) == 'i0' ) then
          call ncdf_def_var(this,name,itype, &
              dims=dims,atts=atts)
        else
          call ncdf_def_var(this,name,ltype, &
              dims=dims,atts=atts)
        end if
        deallocate(dims)
        call nullify(atts)
        call nullify(v)
        call nullify(d_var)
      else if ( key(1:5) == 'GROUP' ) then
        name = key(6:) ! default name
        if ( 'name' .in. d_var ) then
          call associate(v,d_var,'name')
          name = ' '
          call assign(name,v)
        end if
        call ncdf_def_grp(this,name,grp)
        call associate(d_var,d)
        call ncdf_crt(grp,d_var)
        call nullify(d_var)
      end if
      d = .next. d
    end do
  end subroutine ncdf_crt
  ! Handy routine for completely deleting a
  ! dictionary containing a "construction" dictionary
  ! for NetCDF files.
  recursive subroutine ncdf_crt_delete(dic)
    use dictionary
    type(dictionary_t), intent(inout) :: dic
    type(dictionary_t) :: ld, v_dic, att_dic
    character(len=DICTIONARY_KEY_LENGTH) :: key
    ! Delete all entries
    ld = .first. dic
    do while ( .not. (.empty.ld) )
      key = .key. ld
      if ( key(1:3) == 'VAR' ) then
        call associate(v_dic,ld)
        if ( 'atts'.in.v_dic ) then
          call associate(att_dic,v_dic)
          call delete(att_dic)
        end if
        call delete(v_dic)
      else if ( key(1:5) == 'GROUP' ) then
        call associate(v_dic,ld)
        ld = .next. ld
        call ncdf_crt_delete(v_dic)
        cycle
      end if
      ld = .next. ld
    end do
    call delete(dic)
  end subroutine ncdf_crt_delete
  subroutine ncdf_par_access(this,name,access)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: access
  end subroutine ncdf_par_access
  subroutine ncdf_close(this)
    type(hNCDF), intent(inout) :: this
    if ( .not. ncdf_participate(this) ) return
    if ( this%f_id < 0 ) return
    call ncdf_err(nf90_close(this%f_id),"Closing NetCDF file: "//this)
    ! Reset
    call h_reset(this)
  end subroutine ncdf_close
  subroutine ncdf_inq_ncdf(this,dims,vars,atts,format,grps,exist, &
      dict_dim, dict_att)
    use dictionary
    type(hNCDF), intent(inout) :: this
    integer, optional, intent(out) :: dims, vars, atts, format, grps
    logical, optional, intent(out) :: exist
    ! possibly obtain all attributes, dimensions
    type(dictionary_t), optional, intent(inout) :: dict_dim, dict_att
    integer :: ldims, lvars, latts, lformat, lgrps, val, i
    integer, allocatable :: grp_id(:)
    character(len=NF90_MAX_NAME) :: key
    if ( .not. ncdf_participate(this) ) return
    ! A file-check has been requested...
    if ( present(exist) ) then
      inquire(file=this%name,exist=exist)
      ! if it does not exist we simply return
      ! this ensures that the user can request all the information
      ! at once
      if ( .not. exist ) then
        return
      end if
    end if
    call ncdf_err(nf90_inquire(this%id,ldims,lvars,latts,formatNum=lformat), &
        "Inquiring file information "//this)
    ! Copy over requested information...
    if ( present(dims) ) dims = ldims
    if ( present(vars) ) vars = lvars
    if ( present(atts) ) atts = latts
    if ( present(format) ) format = lformat
    if ( present(dict_dim) ) then
      call delete(dict_dim)
      do i = 1 , ldims
        call ncdf_err(nf90_inquire_dimension(this%id,i,name=key))
        call ncdf_inq_dim(this,key,len=val)
        dict_dim = dict_dim // (trim(key).kv.val)
      end do
    end if
    if ( present(dict_att) ) then
      call delete(dict_att)
      call get_atts_id(this,NF90_GLOBAL,dict_att)
    end if
    if ( present(grps) ) then
      if ( IAND(this%mode,NF90_NETCDF4) == NF90_NETCDF4 ) then
        allocate(grp_id(50))
        call ncdf_err(nf90_inq_grps(this%id,grps,grp_id), &
            "Inquiring file information "//this)
        if ( grps > size(grp_id) ) then
          deallocate(grp_id)
          allocate(grp_id(grps))
          call ncdf_err(nf90_inq_grps(this%id,grps,grp_id), &
              "Inquiring file information "//this)
          deallocate(grp_id)
        end if
      else
        grps = -1
      end if
    end if
  end subroutine ncdf_inq_ncdf
  subroutine ncdf_inq_name(name,dims,vars,atts,format,grps,exist, &
      dict_dim, dict_att)
    use dictionary
    character(len=*), intent(in) :: name
    integer, optional, intent(out) :: dims, vars, atts, format, grps
    logical, optional, intent(out) :: exist
    ! possibly obtain all attributes, dimensions
    type(dictionary_t), optional, intent(inout) :: dict_dim, dict_att
    type(hNCDF) :: this
    ! A file-check has been requested...
    if ( present(exist) ) then
      inquire(file=name,exist=exist)
      if ( .not. exist ) then
        return
      end if
    end if
    ! Open the file...
    call ncdf_open(this,name,parallel=.false.)
    ! Do the inquiry...
    call ncdf_inq_ncdf(this,dims=dims,vars=vars,atts=atts,grps=grps,&
        format=format,dict_dim=dict_dim,dict_att=dict_att)
    ! Close the file
    call ncdf_close(this)
  end subroutine ncdf_inq_name
  subroutine ncdf_inq_grp(this,group,exist,dims,vars,atts,format,grps, &
      dict_dim, dict_att)
    use dictionary
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: group
    logical, optional, intent(out) :: exist
    integer, optional, intent(out) :: dims, vars, atts, format, grps
    ! possibly obtain all attributes, dimensions
    type(dictionary_t), optional, intent(inout) :: dict_dim, dict_att
    integer :: ldims, lvars, latts, lformat, lgrps, val, i
    integer, allocatable :: grp_id(:)
    character(len=NF90_MAX_NAME) :: key
    type(hNCDF) :: grp
    if ( .not. ncdf_participate(this) ) return
    if ( present(exist) ) then
      i = nf90_inq_grp_ncid(this%id,group,val)
      if ( i == NF90_NOERR ) then
        exist = .true.
      else if ( i == NF90_ENOGRP ) then
        exist = .false.
      else
        call ncdf_err(i, &
            'Inquiring group information'//this)
      end if
      if ( exist ) then
        ! We can fetch the other information
        call ncdf_open_grp(this,group,grp)
        call ncdf_inq(grp,dims=dims,vars=vars,atts=atts,format=format, &
            grps=grps,dict_dim=dict_dim,dict_att=dict_att)
      end if
      return
    end if
    call ncdf_open_grp(this,group,grp)
    call ncdf_inq(grp,dims=dims,vars=vars,atts=atts,format=format, &
        grps=grps,dict_dim=dict_dim,dict_att=dict_att)
  end subroutine ncdf_inq_grp
  ! Routine to assert that a NetCDF file has
  ! certain dimensions, attributes, and variables.
  ! There exists two interfaces for figuring out the
  ! dims/vars
  ! and
  ! has_dims/has_vars
  ! The former checks dimensions and their values.
  ! The latter only checks if they exist in the
  ! file.
  subroutine ncdf_assert(this,assert,dims,vars, &
      has_dims,has_vars,s_EPS,d_EPS)
    use variable
    use dictionary
    type(hNCDF), intent(inout) :: this
    logical, intent(out) :: assert
    type(dictionary_t), intent(in), optional :: dims, vars
    type(dictionary_t), intent(in), optional :: has_dims, has_vars
    real(sp), intent(in), optional :: s_EPS
    real(dp), intent(in), optional :: d_EPS
    ! We currently do not check attributes.
    ! This is a little more tricky as strings, chars, etc... :(
    ! It just needs to be done...
    ! We can currently only check integers :(
    character(len=DICTIONARY_KEY_LENGTH) :: key
    character(len=VARIABLE_TYPE_LENGTH) :: t
    type(dictionary_t) :: dic ! local loop dictionary...
    type(variable_t) :: ivar
    logical :: success
    integer, pointer :: i1(:), i2(:,:)
    integer, allocatable :: i1a(:), i2a(:,:)
    ! We will also allow comparison of single/doubles
    real(sp) :: ls_EPS
    real(sp), pointer :: s1(:), s2(:,:)
    real(sp), allocatable :: s1a(:), s2a(:,:)
    real(dp) :: ld_EPS
    real(dp), pointer :: d1(:), d2(:,:)
    real(dp), allocatable :: d1a(:), d2a(:,:)
    integer :: i, i0
    assert = .true.
    if ( .not. ncdf_participate(this) ) return
    if ( present(dims) ) then
      ! We check the dimensions of the file
      dic = .first. dims
      do while ( .not. (.empty. dic) )
        key = .key. dic
        ! Check that the dimension exists
        call ncdf_inq_dim(this,key,exist=assert)
        if ( .not. assert ) exit
        ! Get the value in the dictionary
        ivar = .valp. dic
        call assign(i0,ivar,success=success)
        if ( .not. success ) then
          ! Error in type of dictionary...
          call ncdf_err(-100, &
              'Request of dimension in ncdf_assert &
              &went wrong, the dimension is not an integer.')
        end if
        call nullify(ivar)
        ! Now find the dimension size
        call ncdf_inq_dim(this,key,len=i)
        ! Now we can actually check it...
        assert = ( i == i0 )
        if ( .not. assert ) exit
        dic = .next. dic
      end do
      ! Clean up... (the variable allocates the "enc")
      call nullify(ivar)
      if ( .not. assert ) return
    end if
    if ( present(vars) ) then
      ! We retrieve the epsilon for check
      ls_EPS = 1.e-6_sp
      if ( present(s_EPS) ) ls_EPS = s_EPS
      ld_EPS = 1.e-6_dp
      if ( present(d_EPS) ) ld_EPS = d_EPS
      ! We check the dimensions of the file
      dic = .first. vars
      do while ( .not. (.empty. dic) )
        key = .key. dic
        ! Check that the variable exists
        call ncdf_inq_var(this,key,exist=assert)
        if ( .not. assert ) exit
        ! Get the value in the dictionary
        ivar = .valp. dic
        t = which(ivar)
        select case ( t )
        case ( 'i0' )
          call assign(i0,ivar,success=success)
        case ( 'i1' )
          call associate(i1,ivar,success=success)
        case ( 'i2' )
          call associate(i2,ivar,success=success)
        case ( 's1' )
          call associate(s1,ivar,success=success)
        case ( 's2' )
          call associate(s2,ivar,success=success)
        case ( 'd1' )
          call associate(d1,ivar,success=success)
        case ( 'd2' )
          call associate(d2,ivar,success=success)
        case default
          success = .false.
        end select
        if ( .not. success ) then
          ! Error in type of dictionary...
          call ncdf_err(-100, &
              'Request of variable in input vars in ncdf_assert &
              &went wrong, the variable is not [i0,i1,i2,s1,s2,d1,d2].')
        end if
        ! Do not deallocate the dictionary stuff
        call nullify(ivar)
        ! Now grab the first elements of the variable
        ! First we need to allocate the read in data
        ! array
        select case ( t )
        case ( 'i0' )
          call ncdf_get_var(this,key,i)
          assert = ( i == i0 )
        case ( 'i1' )
          allocate(i1a(size(i1)))
          call ncdf_get_var(this,key,i1a)
          assert = all( i1 == i1a )
          deallocate(i1a)
        case ( 'i2' )
          allocate(i2a(size(i2,dim=1),size(i2,dim=2)))
          call ncdf_get_var(this,key,i2a)
          assert = all( i2 == i2a )
          deallocate(i2a)
        case ( 's1' )
          allocate(s1a(size(s1)))
          call ncdf_get_var(this,key,s1a)
          assert = all( abs(s1 - s1a) <= ls_EPS )
          deallocate(s1a)
        case ( 's2' )
          allocate(s2a(size(s2,dim=1),size(s2,dim=2)))
          call ncdf_get_var(this,key,s2a)
          assert = all( abs(s2 - s2a) <= ls_EPS )
          deallocate(s2a)
        case ( 'd1' )
          allocate(d1a(size(d1)))
          call ncdf_get_var(this,key,d1a)
          assert = all( abs(d1 - d1a) <= ld_EPS )
          deallocate(d1a)
        case ( 'd2' )
          allocate(d2a(size(d2,dim=1),size(d2,dim=2)))
          call ncdf_get_var(this,key,d2a)
          assert = all( abs(d2 - d2a) <= ld_EPS )
          deallocate(d2a)
        end select
        ! no success...
        if ( .not. assert ) exit
        dic = .next. dic
      end do
      ! Clean up...
      call nullify(ivar)
      if ( .not. assert ) return
    end if
    if ( present(has_dims) ) then
      ! We check the dimensions of the file
      dic = .first. has_dims
      do while ( .not. (.empty. dic) )
        key = .key. dic
        call ncdf_inq_dim(this,key,exist=assert)
        if ( .not. assert ) return
        dic = .next. dic
      end do
    end if
    if ( present(has_vars) ) then
      dic = .first. has_vars
      do while ( .not. (.empty. dic) )
        key = .key. dic
        call ncdf_inq_var(this,key,exist=assert)
        if ( .not. assert ) return
        dic = .next. dic
      end do
    end if
  end subroutine ncdf_assert
  ! Simplify the addition of any dimension
  subroutine ncdf_def_dim(this,name,size)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: size
    integer :: id
    if ( .not. ncdf_participate(this) ) return
    ! ensure definition step
    ! in case of netCDF-3 this will not change anything
    call ncdf_redef(this)
    call ncdf_err(nf90_def_dim(this%id, name, size, id),&
        "Defining dimension: "//trim(name)//" in file: "//this)
  end subroutine ncdf_def_dim
  ! Simplify the renaming of any dimension
  subroutine ncdf_rename_var(this,old_name,new_name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: old_name, new_name
    integer :: id
    if ( .not. ncdf_participate(this) ) return
    call ncdf_redef(this)
    call ncdf_inq_var(this,old_name,id=id)
    call ncdf_err(nf90_rename_var(this%id, id, new_name),&
        "Renaming variable: "//trim(old_name)//" to "//&
        trim(new_name)//" in file: "//this)
  end subroutine ncdf_rename_var
  ! Simplify the renaming of any dimension
  subroutine ncdf_rename_dim(this,old_name,new_name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: old_name, new_name
    integer :: id
    if ( .not. ncdf_participate(this) ) return
    call ncdf_redef(this)
    call ncdf_inq_dim(this,old_name,id=id)
    call ncdf_err(nf90_rename_dim(this%id, id, new_name),&
        "Renaming dimension: "//trim(old_name)//" to "//&
        trim(new_name)//" in file: "//this)
  end subroutine ncdf_rename_dim
  subroutine ncdf_rename_att(this,var,old_name,new_name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: var, old_name, new_name
    integer :: id
    if ( .not. ncdf_participate(this) ) return
    call ncdf_redef(this)
    call ncdf_inq_var(this,var,id=id)
    call ncdf_err(nf90_rename_att(this%id, id, old_name, new_name),&
        "Renaming variable ("//trim(var)//") attribute: "//trim(old_name)&
        //" to "//trim(new_name)//" in file: "//this)
  end subroutine ncdf_rename_att
  subroutine ncdf_rename_gatt(this,old_name,new_name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: old_name, new_name
    if ( .not. ncdf_participate(this) ) return
    call ncdf_redef(this)
    call ncdf_err(nf90_rename_att(this%id, NF90_GLOBAL, old_name, new_name),&
        "Renaming global attribute: "//trim(old_name)&
        //" to "//trim(new_name)//" in file: "//this)
  end subroutine ncdf_rename_gatt
  ! Simplify the addition of a variable...
  ! This routine *MUST* be called after ncdf_participate
  ! (however, as it is a local routine the burden is ours, not the programmers)
  subroutine ncdf_def_var_generic(this,name,type,dims,id,atts, &
      compress_lvl,shuffle, access, chunks)
    use variable
    use dictionary
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    character(len=*), intent(in) :: dims(:)
    integer, intent(out) :: id
    type(dictionary_t), optional :: atts
    integer, intent(in), optional :: compress_lvl
    logical, intent(in), optional :: shuffle
    integer, intent(in), optional :: access, chunks(:)
    integer :: iret, i, ldims(size(dims))
    call ncdf_redef(this)
    do i = 1 , size(dims)
      call ncdf_inq_dim(this,trim(dims(i)),id=ldims(i))
    end do
    ! Determine whether we have NetCDF 4 enabled, in that case do compression if asked for
    ! In case of NetCDF 3
    iret = nf90_def_var(this%id, name, type, ldims, id)
    call ncdf_err(iret,"Defining variable: "//trim(name)//" in file: "//this)
    if ( present(atts) ) then
      call put_atts_id(this,id,atts)
    end if
    if ( present(access) ) then
      call ncdf_par_access(this,name=name,access=access)
    end if
  end subroutine ncdf_def_var_generic
  subroutine ncdf_def_var_integer(this, name, type, dims, &
      atts, compress_lvl, shuffle, fill, &
      access, chunks)
    use dictionary
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    character(len=*), intent(in) :: dims(:)
    type(dictionary_t), optional :: atts
    integer, intent(in), optional :: compress_lvl
    logical, intent(in), optional :: shuffle, fill
    integer, intent(in), optional :: access, chunks(:)
    integer :: id
    if ( .not. ncdf_participate(this) ) then
      ! in case the attributes are present, we
      ! still need to clean-up if asked
      if ( present(atts) ) then
        if ( 'ATT_DELETE' .in. atts ) then
          call delete(atts)
        end if
      end if
      return
    end if
    call ncdf_def_var_generic(this, name, type, dims, id, &
        atts=atts, compress_lvl=compress_lvl, shuffle=shuffle, &
        access=access, chunks=chunks)
    if ( present(fill) ) then
      if ( fill ) then
        call ncdf_err(nf90_def_var_fill(this%id,id, 1, 0), &
            "Setting the variable "//trim(name)//" to NOFILL in file "//this)
      else
        call ncdf_err(nf90_def_var_fill(this%id,id, 0, 0), &
            "Setting the variable "//trim(name)//" to FILL in file "//this)
      end if
    end if
  end subroutine ncdf_def_var_integer
  subroutine ncdf_def_var_logical(this, name, type, dims, &
      atts, compress_lvl, shuffle, fill, &
      access, chunks)
    use dictionary
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    logical, intent(in) :: type
    character(len=*), intent(in) :: dims(:)
    type(dictionary_t), optional :: atts
    integer, intent(in), optional :: compress_lvl
    logical, intent(in), optional :: shuffle, fill
    integer, intent(in), optional :: access, chunks(:)
    integer :: id
    integer :: ltype
    if ( .not. ncdf_participate(this) ) then
      ! in case the attributes are present, we
      ! still need to clean-up if asked
      if ( present(atts) ) then
        if ( 'ATT_DELETE' .in. atts ) then
          call delete(atts)
        end if
      end if
      return
    end if
    if ( type .eqv. NF90_DOUBLE_COMPLEX ) then
      ltype = NF90_DOUBLE
    else
      ltype = NF90_FLOAT
    end if
    if ( type .eqv. NF90_DOUBLE_COMPLEX ) then
      ltype = NF90_DOUBLE
    else
      ltype = NF90_FLOAT
    end if
    call ncdf_def_var_generic(this, "Re"//trim(name), ltype, dims, id, &
        atts=atts, compress_lvl=compress_lvl, shuffle=shuffle, &
        access=access, chunks=chunks)
    if ( present(fill) ) then
      if ( fill ) then
        call ncdf_err(nf90_def_var_fill(this%id,id, 1, 0), &
            "Setting the variable "//trim(name)//" to NOFILL in file "//this)
      else
        call ncdf_err(nf90_def_var_fill(this%id,id, 0, 0), &
            "Setting the variable "//trim(name)//" to FILL in file "//this)
      end if
    end if
    call ncdf_def_var_generic(this, "Im"//trim(name), ltype, dims, id, &
        atts=atts, compress_lvl=compress_lvl, shuffle=shuffle, &
        access=access, chunks=chunks)
    if ( present(fill) ) then
      if ( fill ) then
        call ncdf_err(nf90_def_var_fill(this%id,id, 1, 0), &
            "Setting the variable "//trim(name)//" to NOFILL in file "//this)
      else
        call ncdf_err(nf90_def_var_fill(this%id,id, 0, 0), &
            "Setting the variable "//trim(name)//" to FILL in file "//this)
      end if
    end if
  end subroutine ncdf_def_var_logical
  subroutine ncdf_default(this,access,compress_lvl)
    type(hNCDF), intent(inout) :: this
    integer, intent(in), optional :: access, compress_lvl
    if ( .not. ncdf_participate(this) ) return
    if ( present(access) ) call ncdf_par_access(this,access=access)
  end subroutine ncdf_default
  subroutine ncdf_inq_var_def(this,name,exist,id,size,atts)
    use dictionary
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    logical, optional, intent(out) :: exist
    integer, optional, intent(out) :: id
    integer, optional, intent(out) :: size(:)
    type(dictionary_t), optional, intent(inout) :: atts
    integer :: iret ! We need to retain any error message...
    integer :: lid, nids, i
    integer :: ldids(10) ! In case the user only wishes to read a sub-part of the size
    character(len=NF90_MAX_NAME) :: dim
    logical :: lexist
    if ( .not. ncdf_participate(this) ) return
    ! Figure out if the dimension exists
    iret = nf90_inq_varid(this%id, trim(name), lid)
    ! The variable must exist
    lexist = iret == NF90_NOERR
    if ( present(exist) ) then
      exist = lexist
    else if ( .not. lexist ) then
      call ncdf_err(iret,"Retrieving information about: "//trim(name)//" in file: "//this)
    end if
    ! If there is nothing to inquire: return
    if ( .not. lexist ) return
    if ( present(id) ) id = lid
    ! If the user has requested information about the size of the variable...
    if ( present(size) ) then
      call ncdf_err(nf90_inquire_variable(this%id, lid, ndims=nids, dimids=ldids))
      do i = 1 , min(nids,ubound(size,1))
        call ncdf_err(nf90_inquire_dimension(this%id,ldids(i),name=dim), &
            "Retrieving dimension name in inq_var for file: "//this)
        ! Save the dimension size in array "size"
        call ncdf_inq_dim(this,trim(dim),len=size(i))
      end do
    end if
    ! The user has requested information about the attributes associated...
    if ( present(atts) ) then
      call get_atts_id(this,lid,atts=atts)
    end if
  end subroutine ncdf_inq_var_def
  subroutine ncdf_inq_dim(this,name,exist,id,len)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    logical, optional, intent(out) :: exist
    integer, optional, intent(out) :: id
    integer, optional, intent(out) :: len
    integer :: iret ! We need to retain any error message...
    integer :: lid
    logical :: lexist
    if ( .not. ncdf_participate(this) ) return
    ! Figure out if the dimension exists
    iret = nf90_inq_dimid(this%id, trim(name), lid)
    lexist = iret == NF90_NOERR
    if ( present(exist) ) then
      exist = lexist
    else if ( .not. lexist ) then
      call ncdf_err(iret,"Retrieving information about: "//trim(name)//" in file: "//this)
    end if
    ! If there is nothing to inquire: return
    if ( .not. lexist ) return
    if ( present(id) ) id = lid
    if ( present(len) ) then
      call ncdf_err(nf90_inquire_dimension(this%id, lid, len=len), &
          "Retrieving length of dimension: "//trim(name)//" in file: "//this)
    end if
  end subroutine ncdf_inq_dim
  subroutine ncdf_inq_gatt(this,name,exist,len,xtype)
    use variable
    type(hNCDF), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    logical, optional, intent(out) :: exist
    integer, optional, intent(out) :: len, xtype
    integer :: iret ! We need to retain any error message...
    logical :: lexist
    if ( .not. ncdf_participate(this) ) return
    ! Figure out if the dimension exists
    iret = nf90_inquire_attribute(this%id, NF90_GLOBAL, trim(name), &
        len = len, xtype = xtype)
    lexist = iret == NF90_NOERR
    if ( present(exist) ) then
      exist = lexist
    else if ( .not. lexist ) then
      call ncdf_err(iret,"Retrieving information about: "//trim(name)//" in file: "//this)
    end if
  end subroutine ncdf_inq_gatt
  subroutine ncdf_inq_att(this,var,name,exist,len,xtype)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: name
    logical, optional, intent(out) :: exist
    integer, optional, intent(out) :: len, xtype
    integer :: iret ! We need to retain any error message...
    integer :: id
    logical :: lexist
    if ( .not. ncdf_participate(this) ) return
    call ncdf_inq_var(this,var,id=id)
    ! Figure out if the dimension exists
    iret = nf90_inquire_attribute(this%id, id, trim(name), &
        len = len, xtype = xtype)
    lexist = iret == NF90_NOERR
    if ( present(exist) ) then
      exist = lexist
    else if ( .not. lexist ) then
      call ncdf_err(iret,"Retrieving information about: "//trim(name)//" in file: "//this)
    end if
  end subroutine ncdf_inq_att
  subroutine put_gatt(this,name,att,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    type(variable_t), optional, intent(inout) :: att
    type(dictionary_t), optional, intent(inout) :: atts
    if ( .not. ncdf_participate(this) ) return
    if ( present(name) .and. present(att) ) then
      call put_att_id(this,NF90_GLOBAL,trim(name),att)
    else if ( present(atts) ) then
      call put_atts_id(this,NF90_GLOBAL,atts)
    else
      call ncdf_err(-100, &
          'Programming error: put_gatt interface not properly populated')
    end if
  end subroutine put_gatt
  subroutine get_gatt(this,name,att,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    type(variable_t), optional, intent(inout) :: att
    type(dictionary_t), optional, intent(inout) :: atts
    if ( .not. ncdf_participate(this) ) return
    if ( present(name) .and. present(att) ) then
      call get_att_id(this,NF90_GLOBAL,trim(name),att)
    else if ( present(atts) ) then
      call get_atts_id(this,NF90_GLOBAL,atts)
    else
      call ncdf_err(-100, &
          'Programming error: get_gatt interface not properly populated')
    end if
  end subroutine get_gatt
  subroutine put_att(this,var,name,att,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: var
    character(len=*), optional, intent(in) :: name
    type(variable_t), optional, intent(inout) :: att
    type(dictionary_t), optional, intent(inout) :: atts
    integer :: ID
    if ( .not. ncdf_participate(this) ) return
    call ncdf_inq_var(this,var,id=ID)
    if ( present(name) .and. present(att) ) then
      call put_att_id(this,ID,trim(name),att)
    else if ( present(atts) ) then
      call put_atts_id(this,ID,atts)
    else
      call ncdf_err(-100, &
          'Programming error: put_att interface not properly populated')
    end if
  end subroutine put_att
  subroutine get_att(this,var,name,att,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: var
    character(len=*), optional, intent(in) :: name
    type(variable_t), optional, intent(inout) :: att
    type(dictionary_t), optional, intent(inout) :: atts
    integer :: ID
    if ( .not. ncdf_participate(this) ) return
    call ncdf_inq_var(this,var,id=ID)
    if ( present(name) .and. present(att) ) then
      call get_att_id(this,ID,trim(name),att)
    else if ( present(atts) ) then
      call get_atts_id(this,ID,atts)
    else
      call ncdf_err(-100, &
          'Programming error: get_att interface not properly populated')
    end if
  end subroutine get_att
  subroutine put_atts_id(this,id,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    integer, intent(in) :: ID
    type(dictionary_t), intent(inout) :: atts
    type(dictionary_t) :: att
    type(variable_t) :: at_var
    character(len=NF90_MAX_NAME) :: key
    if ( len(atts) == 0 ) return
    att = .first. atts
    att_loop: do
      if ( .empty. att ) exit att_loop
      key = .key. att
      if ( key == "ATT_DELETE" ) then
        att = .next. att
        cycle
      end if
      ! we do not copy any data
      call associate(at_var,att)
      call put_att_id(this,id,trim(key),at_var)
      att = .next. att
    end do att_loop
    ! clean-up encoding
    call nullify(at_var)
    ! If the user adds this key, the dictionary will be deleted
    ! after usage...
    if ( "ATT_DELETE" .in. atts ) then
      call delete(atts)
    end if
  end subroutine put_atts_id
  subroutine put_att_id(this,id,name,att)
    use variable
    type(hNCDF), intent(inout) :: this
    integer, intent(in) :: ID
    character(len=*), intent(in) :: name
    type(variable_t), intent(inout) :: att
    integer :: iret
    character(len=NF90_MAX_NAME) :: tmp
    integer(ih), pointer :: h0, h1(:)
    integer(is), pointer :: i0, i1(:)
    real(sp), pointer :: s0, s1(:)
    real(dp), pointer :: d0, d1(:)
    call ncdf_redef(this)
    select case ( which(att) )
    case ( 'a1' ) ! character array
      call assign(tmp,att)
      iret = nf90_put_att(this%id, id, trim(name), tmp)
    case ( 'h0' )
      call associate(h0,att)
      iret = nf90_put_att(this%id, id, trim(name), h0)
    case ( 'h1' )
      call associate(h1,att)
      iret = nf90_put_att(this%id, id, trim(name), h1)
    case ( 'i0' )
      call associate(i0,att)
      iret = nf90_put_att(this%id, id, trim(name), i0)
    case ( 'i1' )
      call associate(i1,att)
      iret = nf90_put_att(this%id, id, trim(name), i1)
    case ( 's0' )
      call associate(s0,att)
      iret = nf90_put_att(this%id, id, trim(name), s0)
    case ( 's1' )
      call associate(s1,att)
      iret = nf90_put_att(this%id, id, trim(name), s1)
    case ( 'd0' )
      call associate(d0,att)
      iret = nf90_put_att(this%id, id, trim(name), d0)
    case ( 'd1' )
      call associate(d1,att)
      iret = nf90_put_att(this%id, id, trim(name), d1)
    case default
      iret = -100
    end select
    call ncdf_err(iret, &
        "Saving attribute: "//trim(name)// &
        " in file: "//this)
  end subroutine put_att_id
  subroutine get_atts_id(this,id,atts)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    integer, intent(in) :: ID
    type(dictionary_t), intent(inout) :: atts
    integer :: i, nAtts
    character(len=NF90_MAX_NAME) :: name
    type(variable_t) :: att
    if ( id == NF90_GLOBAL ) then
      call ncdf_err(nf90_inquire(this%id, nAttributes=nAtts), &
          "Retrieving number of associated attributes in inquire for file: "//this)
    else
      call ncdf_err(nf90_inquire_variable(this%id, id, nAtts=nAtts), &
          "Retrieving number of associated attributes in inq_var for file: "//this)
    end if
    do i = 1 , nAtts
      name = ' '
      call ncdf_err(nf90_inq_attname(this%id, id, i, name), &
          "Retrieving the attribute name for file: "//this)
      call get_att_id(this,id,name,att)
      call extend(atts,(trim(name).kv.att))
    end do
    call delete(att)
  end subroutine get_atts_id
  subroutine get_att_id(this,ID,name,att)
    use dictionary
    use variable
    type(hNCDF), intent(inout) :: this
    integer, intent(in) :: ID
    character(len=*), intent(in) :: name
    type(variable_t), intent(inout) :: att
    integer :: xtype, att_len
    character(len=512) :: att_char
    real(sp), allocatable :: a_sp(:)
    real(dp), allocatable :: a_dp(:)
    integer(ih), allocatable :: a_ih(:)
    integer(is), allocatable :: a_is(:)
    ! retrieve the attribute length and data-type
    call ncdf_err(nf90_inquire_attribute(this%id,id,trim(name), &
        xtype=xtype,len=att_len),'Retriving inquire_attribute: '//this)
    select case ( xtype )
    case ( NF90_CHAR )
      att_char = ' '
      call ncdf_err(nf90_get_att(this%id, id, trim(name), att_char), &
          "Retrieving the attribute value for file: "//this)
      call assign(att,trim(att_char))
    case ( NF90_SHORT )
      allocate(a_ih(att_len))
      call ncdf_err(nf90_get_att(this%id, id, trim(name), a_ih), &
          "Retrieving the attribute value for file: "//this)
      if ( att_len == 1 ) then
        call assign(att,a_ih(1))
      else
        call assign(att,a_ih)
      end if
      deallocate(a_ih)
    case ( NF90_INT )
      allocate(a_is(att_len))
      call ncdf_err(nf90_get_att(this%id, id, trim(name), a_is), &
          "Retrieving the attribute value for file: "//this)
      if ( att_len == 1 ) then
        call assign(att,a_is(1))
      else
        call assign(att,a_is)
      end if
      deallocate(a_is)
    case ( NF90_FLOAT )
      allocate(a_sp(att_len))
      call ncdf_err(nf90_get_att(this%id, id, trim(name), a_sp), &
          "Retrieving the attribute value for file: "//this)
      if ( att_len == 1 ) then
        call assign(att,a_sp(1))
      else
        call assign(att,a_sp)
      end if
      deallocate(a_sp)
    case ( NF90_DOUBLE )
      allocate(a_dp(att_len))
      call ncdf_err(nf90_get_att(this%id, id, trim(name), a_dp), &
          "Retrieving the attribute value for file: "//this)
      if ( att_len == 1 ) then
        call assign(att,a_dp(1))
      else
        call assign(att,a_dp)
      end if
      deallocate(a_dp)
    end select
  end subroutine get_att_id
  ! Delete attributes
  subroutine ncdf_del_att(this,var,name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: name
    integer :: iret ! We need to retain any error message...
    integer :: id
    call ncdf_redef(this)
    if ( .not. ncdf_participate(this) ) return
    call ncdf_inq_var(this,var,id=id)
    ! Figure out if the dimension exists
    iret = nf90_inquire_attribute(this%id, id, trim(name))
    if ( iret == NF90_NOERR ) then
      call ncdf_err(nf90_del_att(this%id, id, trim(name)), &
          "Deleting attribute: "//trim(name)//" for variable "//&
          trim(var)//" in file: "//this)
    end if
  end subroutine ncdf_del_att
  subroutine ncdf_del_gatt(this,name)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer :: iret ! We need to retain any error message...
    integer :: id
    call ncdf_redef(this)
    if ( .not. ncdf_participate(this) ) return
    ! Figure out if the dimension exists
    iret = nf90_inquire_attribute(this%id, NF90_GLOBAL, trim(name))
    if ( iret == NF90_NOERR ) then
      call ncdf_err(nf90_del_att(this%id, NF90_GLOBAL, trim(name)), &
          "Deleting global attribute: "//trim(name)//" in file: "//this)
    end if
  end subroutine ncdf_del_gatt
  subroutine ncdf_fill(this,fill,old_fill)
    type(hNCDF), intent(inout) :: this
    integer, optional, intent(in) :: fill
    integer, optional, intent(out) :: old_fill
    integer :: lf, lof
    ! option collect
    lf = NF90_NOFILL
    if ( present(fill) ) lf = fill
    call ncdf_err(nf90_set_fill(this%id,lf, lof), &
        "Setting fill mode in file: "//this)
    if ( present(old_fill) ) old_fill = lof
    if ( .not. present(fill) ) &
        call ncdf_err(nf90_set_fill(this%id,lof, lf), &
        "Re-setting fill mode in file: "//this)
  end subroutine ncdf_fill
  ! Use the ncdf.sh script to generate the needed code...
subroutine put_var_h0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_h0_name
subroutine get_var_h0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer(ih) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_h0_name
subroutine def_fill_h0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  integer(ih) :: lfill_val
  integer :: tmp_lr
  integer(ih) :: lr
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  call ncdf_err(nf90_def_var_fill(this%id, id, lfill, tmp_lr), &
       'Setting fill (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine def_fill_h0
subroutine inq_var_h0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  integer(ih) :: lfill_val
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) '//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  !if ( present(fill_val) ) fill_val = lfill_val
  fill_val = 0.
end subroutine inq_var_h0
subroutine put_var_h1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_h1_name
subroutine get_var_h1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_h1_name
subroutine put_var_h2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_h2_name
subroutine get_var_h2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_h2_name
subroutine put_var_h3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_h3_name
subroutine get_var_h3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(ih), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_h3_name
subroutine put_gatt_s0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: att
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_s0
subroutine put_att_s0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(sp), intent(in) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_s0
subroutine get_att_s0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(sp), intent(out) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_s0
subroutine get_gatt_s0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: att
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_s0
subroutine put_var_s0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_s0_name
subroutine get_var_s0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  real(sp) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_s0_name
subroutine def_fill_s0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  real(sp) :: lfill_val
  integer :: tmp_lr
  real(sp) :: lr
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  call ncdf_err(nf90_def_var_fill(this%id, id, lfill, tmp_lr), &
       'Setting fill (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine def_fill_s0
subroutine inq_var_s0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  real(sp) :: lfill_val
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) '//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  !if ( present(fill_val) ) fill_val = lfill_val
  fill_val = 0.
end subroutine inq_var_s0
subroutine put_gatt_s1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_s1
subroutine put_att_s1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(sp), intent(in) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_s1
subroutine get_att_s1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(sp), intent(out) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_s1
subroutine get_gatt_s1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_s1
subroutine put_var_s1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_s1_name
subroutine get_var_s1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_s1_name
subroutine put_var_s2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_s2_name
subroutine get_var_s2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_s2_name
subroutine put_var_s3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_s3_name
subroutine get_var_s3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(sp), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_s3_name
subroutine put_gatt_d0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: att
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_d0
subroutine put_att_d0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(dp), intent(in) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_d0
subroutine get_att_d0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(dp), intent(out) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_d0
subroutine get_gatt_d0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: att
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_d0
subroutine put_var_d0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_d0_name
subroutine get_var_d0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  real(dp) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_d0_name
subroutine def_fill_d0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  real(dp) :: lfill_val
  integer :: tmp_lr
  real(dp) :: lr
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  call ncdf_err(nf90_def_var_fill(this%id, id, lfill, tmp_lr), &
       'Setting fill (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine def_fill_d0
subroutine inq_var_d0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  real(dp) :: lfill_val
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) '//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  !if ( present(fill_val) ) fill_val = lfill_val
  fill_val = 0.
end subroutine inq_var_d0
subroutine put_gatt_d1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_d1
subroutine put_att_d1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(dp), intent(in) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_d1
subroutine get_att_d1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  real(dp), intent(out) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_d1
subroutine get_gatt_d1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_d1
subroutine put_var_d1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_d1_name
subroutine get_var_d1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_d1_name
subroutine put_var_d2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_d2_name
subroutine get_var_d2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_d2_name
subroutine put_var_d3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_d3_name
subroutine get_var_d3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_d3_name
subroutine put_var_c0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(sp), allocatable :: r
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_c0_name
subroutine get_var_c0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  complex(sp) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_c0_name
subroutine def_fill_c0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  complex(sp) :: lfill_val
  integer :: tmp_lr
  real(sp) :: lr
  integer :: id(2)
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  lr = real(lfill_val, sp)
  call ncdf_err(nf90_def_var_fill(this%id, id(1), lfill, tmp_lr), &
       'Setting fill (VAR) Re'//trim(name)//' in file: '//this)
  lr = aimag(lfill_val)
  call ncdf_err(nf90_def_var_fill(this%id, id(2), lfill, tmp_lr), &
       'Setting fill (VAR) Im'//trim(name)//' in file: '//this)
end subroutine def_fill_c0
subroutine inq_var_c0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id(2)
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  real(sp) :: lfill_valr,lfill_valc
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,'Re'//name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id(1) = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) Re-'//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  call ncdf_inq_var_def(this,'Im'//name,id=lid,size=size,atts=atts)
  if ( present(id) ) id(2) = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) Im-'//trim(name)//' in file: '//this)
  if ( present(fill) ) then
     if ( fill /= lfill ) then
      call ncdf_err(-100,'Fill-value for real and imaginary part'//&
   ' are not the same. This is not allowed.')
     end if
  end if
! if ( present(fill_val) ) then
     !fill_val = cmplx(lfill_valr,lfill_valc, sp)
     fill_val = cmplx(0,0, sp)
! end if
end subroutine inq_var_c0
subroutine put_var_c1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(sp), allocatable :: r (:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var)))
  r (:) = real(var, sp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_c1_name
subroutine get_var_c1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(sp), allocatable :: r (:) , i (:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var)))
  allocate(i(size(var)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:) = cmplx(r,i, sp)
  deallocate(r,i)
end subroutine get_var_c1_name
subroutine put_var_c2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(sp), allocatable :: r (:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2)))
  r (:,:) = real(var, sp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:,:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_c2_name
subroutine get_var_c2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(sp), allocatable :: r (:,:) , i (:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2)))
  allocate(i(size(var,1),size(var,2)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:,:) = cmplx(r,i, sp)
  deallocate(r,i)
end subroutine get_var_c2_name
subroutine put_var_c3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(sp), allocatable :: r (:,:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2),size(var,3)))
  r (:,:,:) = real(var, sp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:,:,:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_c3_name
subroutine get_var_c3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(sp), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(sp), allocatable :: r (:,:,:) , i (:,:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2),size(var,3)))
  allocate(i(size(var,1),size(var,2),size(var,3)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:,:,:) = cmplx(r,i, sp)
  deallocate(r,i)
end subroutine get_var_c3_name
subroutine put_var_z0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(dp), allocatable :: r
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_z0_name
subroutine get_var_z0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  complex(dp) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_z0_name
subroutine def_fill_z0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  complex(dp) :: lfill_val
  integer :: tmp_lr
  real(dp) :: lr
  integer :: id(2)
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  lr = real(lfill_val, dp)
  call ncdf_err(nf90_def_var_fill(this%id, id(1), lfill, tmp_lr), &
       'Setting fill (VAR) Re'//trim(name)//' in file: '//this)
  lr = aimag(lfill_val)
  call ncdf_err(nf90_def_var_fill(this%id, id(2), lfill, tmp_lr), &
       'Setting fill (VAR) Im'//trim(name)//' in file: '//this)
end subroutine def_fill_z0
subroutine inq_var_z0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id(2)
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  real(dp) :: lfill_valr,lfill_valc
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,'Re'//name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id(1) = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) Re-'//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  call ncdf_inq_var_def(this,'Im'//name,id=lid,size=size,atts=atts)
  if ( present(id) ) id(2) = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) Im-'//trim(name)//' in file: '//this)
  if ( present(fill) ) then
     if ( fill /= lfill ) then
      call ncdf_err(-100,'Fill-value for real and imaginary part'//&
   ' are not the same. This is not allowed.')
     end if
  end if
! if ( present(fill_val) ) then
     !fill_val = cmplx(lfill_valr,lfill_valc, dp)
     fill_val = cmplx(0,0, dp)
! end if
end subroutine inq_var_z0
subroutine put_var_z1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(dp), allocatable :: r (:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var)))
  r (:) = real(var, dp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_z1_name
subroutine get_var_z1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(dp), allocatable :: r (:) , i (:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var)))
  allocate(i(size(var)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:) = cmplx(r,i, dp)
  deallocate(r,i)
end subroutine get_var_z1_name
subroutine put_var_z2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(dp), allocatable :: r (:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2)))
  r (:,:) = real(var, dp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:,:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_z2_name
subroutine get_var_z2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(dp), allocatable :: r (:,:) , i (:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2)))
  allocate(i(size(var,1),size(var,2)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:,:) = cmplx(r,i, dp)
  deallocate(r,i)
end subroutine get_var_z2_name
subroutine put_var_z3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  real(dp), allocatable :: r (:,:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2),size(var,3)))
  r (:,:,:) = real(var, dp)
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
     "Saving variable (VAR) Re"//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  r (:,:,:) = aimag(var)
  call ncdf_err(nf90_put_var(this%id, id, r,start=start,count=count), &
       'Saving variable (VAR) Im'//trim(name)//' in file: '//this)
  deallocate(r)
end subroutine put_var_z3_name
subroutine get_var_z3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  complex(dp), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  real(dp), allocatable :: r (:,:,:) , i (:,:,:)
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  allocate(r(size(var,1),size(var,2),size(var,3)))
  allocate(i(size(var,1),size(var,2),size(var,3)))
  call ncdf_inq_var(this,'Re'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, r,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Re'//trim(name)//' in file: '//this)
  call ncdf_inq_var(this,'Im'//name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, i,start=start,count=count,stride=stride), &
       'Retrieving variable (VAR) Im'//trim(name)//' in file: '//this)
  var (:,:,:) = cmplx(r,i, dp)
  deallocate(r,i)
end subroutine get_var_z3_name
subroutine put_gatt_i0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: att
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_i0
subroutine put_att_i0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  integer(is), intent(in) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_i0
subroutine get_att_i0(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  integer(is), intent(out) :: att
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_i0
subroutine get_gatt_i0(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: att
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_i0
subroutine put_var_i0_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  call ncdf_put_var(this,name,(/var/),start=start,count=count)
end subroutine put_var_i0_name
subroutine get_var_i0_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: var
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer(is) :: v(1)
  call ncdf_get_var(this,name,v,start=start,count=count,stride=stride)
  var = v(1)
end subroutine get_var_i0_name
subroutine def_fill_i0(this, name, fill_val, fill)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: fill_val ! non-optional to allow interfacing
  integer, intent(in), optional :: fill
  integer :: lfill
  integer(is) :: lfill_val
  integer :: tmp_lr
  integer(is) :: lr
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_redef(this)
  call ncdf_inq_var(this,name,id=id,fill=lfill,fill_val=lfill_val)
  if ( present(fill) ) lfill = fill
  lfill_val = fill_val
  call ncdf_err(nf90_def_var_fill(this%id, id, lfill, tmp_lr), &
       'Setting fill (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine def_fill_i0
subroutine inq_var_i0(this, name,fill_val,exist,id,size,atts,fill)
  use dictionary
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: fill_val ! non-optional to allow interfacing
  logical, intent(out), optional :: exist
  integer, intent(out), optional :: id
  integer, intent(out), optional :: size(:)
  type(dictionary_t), intent(inout), optional :: atts
  integer, intent(out), optional :: fill
  integer :: lid, lfill
  integer :: tmp_lr
  integer(is) :: lfill_val
  if ( .not. ncdf_participate(this) ) return
  call ncdf_inq_var_def(this,name,id=lid,exist=exist,size=size,atts=atts)
  if ( present(exist) ) then
     if ( .not. exist ) return
  end if
  if ( present(id) ) id = lid
  call ncdf_err(nf90_inq_var_fill(this%id, lid, lfill, tmp_lr), &
       'Retrieving variable-fill (VAR) '//trim(name)//' in file: '//this)
  if ( present(fill) ) fill = lfill
  !if ( present(fill_val) ) fill_val = lfill_val
  fill_val = 0.
end subroutine inq_var_i0
subroutine put_gatt_i1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine put_gatt_i1
subroutine put_att_i1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  integer(is), intent(in) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine put_att_i1
subroutine get_att_i1(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  integer(is), intent(out) :: att (:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '//trim(var)//' : '//trim(name)//' in file: '//this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '//trim(var)//' attribute: '//trim(name)//' in file: '//this)
end subroutine get_att_i1
subroutine get_gatt_i1(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: att (:)
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '//trim(name)//' in file: '//this)
end subroutine get_gatt_i1
subroutine put_var_i1_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_i1_name
subroutine get_var_i1_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: var (:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_i1_name
subroutine put_var_i2_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_i2_name
subroutine get_var_i2_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: var (:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_i2_name
subroutine put_var_i3_name(this,name,var,start,count)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_put_var(this%id, id, var,start=start,count=count), &
       'Saving variable (VAR) '//trim(name)//' in file: '//this)
end subroutine put_var_i3_name
subroutine get_var_i3_name(this, name, var, start, count, stride)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(is), intent(out) :: var (:,:,:)
  integer, intent(in), optional :: start(:)
  integer, intent(in), optional :: count(:)
  integer, intent(in), optional :: stride(:)
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if ( this%define > -1 ) call ncdf_enddef(this)
  call ncdf_inq_var(this,name,id=id)
  call ncdf_err(nf90_get_var(this%id, id, var,start=start,count=count,stride=stride), &
       'Retrieving (VAR) variable, '//trim(name)//' in file: '//this)
end subroutine get_var_i3_name
  subroutine ncdf_enddef(this)
    type(hNCDF), intent(inout) :: this
    integer :: i
    ! A NetCDF4 file still needs to define/redefine
    ! in collective manner, yet it is taken care of
    ! internally.
    if ( this%define == 1 ) return
    if ( this%define == 0 ) this%define = 1
    if ( .not. ncdf_participate(this) ) return
    i = nf90_enddef(this%id)
    if ( i == nf90_noerr ) return
    if ( i == nf90_enotindefine ) then
      ! we pass, the not-in-define mode
      ! just tells us that we are already in data-mode
      return
    end if
    call ncdf_err(i,"End definition segment of file: "//this)
  end subroutine ncdf_enddef
  subroutine ncdf_sync(this)
    type(hNCDF), intent(in) :: this
    ! We need (must) not sync when in definition mode...
    if ( this%define == 0 ) return
    if ( .not. ncdf_participate(this) ) return
    call ncdf_err(nf90_sync(this%f_id), &
        "File syncronization for file"//this)
  end subroutine ncdf_sync
  subroutine ncdf_redef(this)
    type(hNCDF), intent(inout) :: this
    integer :: i
    ! Already in define mode:
    if ( this%define == 0 ) return
    if ( this%define == 1 ) this%define = 0
    if ( .not. ncdf_participate(this) ) return
    i = nf90_redef(this%id)
    if ( i == nf90_noerr ) return
    if ( i == nf90_eindefine ) then
      ! we pass, the in-define mode
      ! just tells us that we are already in define-mode
      return
    end if
    call ncdf_err(i,"Redef definition segment in file: "//this)
  end subroutine ncdf_redef
  !"#############################################################"
  !"##############" Specialized routines for handling"#########"
  !"#########" the different aspects of this module"###########"
  !"#############################################################"
  ! A simple error checker for NetCDF
  subroutine ncdf_err(status,msg)
    integer, intent(in) :: status
    character(len=*), optional, intent(in) :: msg
    ! We should never encounter this region of the code with out having enabled some
    ! reading or writing of the NetCDF handle...
    ! if ( .not. ncdf_participate() ) return
    if (status .ne. nf90_noerr) then
      if (present(msg)) write(*,"(a)") trim(msg)
      write(*,*)
      write(*,"(a)") "Error occured in NCDF:"
      write(0,"(a)") "Error occured in NCDF:"
      select case ( status )
      case ( -130 ) ! include/netcdf.h:"define" NC_ECANTEXTEND (-130)
        ! < Attempt to extend dataset during ind. I/O operation. >!
        write(*,"(a)") 'Attempt to extend (UNLIMITED dimension) a parallel file in INDEPENDENT ACCESS mode'
        write(0,"(a)") 'Attempt to extend (UNLIMITED dimension) a parallel file in INDEPENDENT ACCESS mode'
      case default
        write(*,"(a)") trim(nf90_strerror(status))
        write(0,"(a)") trim(nf90_strerror(status))
      end select
      write(*,"(a,tr1,i0)") "Status number:",status
      write(0,"(a,tr1,i0)") "Status number:",status
      call ncdf_die("Stopped due to error in NetCDF file")
    endif
  end subroutine ncdf_err
  !"#############################################################"
  !"#############" Routines for handling the groups"############"
  !"#########" of the NetCDF files. We allow this to"##########"
  !"#########" always be present due to reduction of"##########"
  !"###############" preprocessor flags"#######################"
  !"#############################################################"
  ! Create groups in a NetCDF4 file
  subroutine ncdf_def_grp(this,name,grp)
    type(hNCDF), intent(in out) :: this
    character(len=*), intent(in) :: name
    type(hNCDF), intent(out) :: grp
    ! Copy the information regarding the parent ncdf
    ! We need to do this out-side (as the information
    ! required to denote the owner of the file)
    call ncdf_copy(this,grp)
    if ( .not. ncdf_participate(grp) ) return
    ! Save the group name... (we save it with hiercharal notice /"grp1"/"grp2")
    grp%grp = trim(grp%grp)//"/"//trim(name)
    ! Create the group and return
    call ncdf_err(nf90_def_grp(this%id,name,grp%id), &
        "Creating group "//trim(name)//" in file "//this)
  end subroutine ncdf_def_grp
  ! Open a group from an existing file
  subroutine ncdf_open_grp(this,name,grp)
    type(hNCDF), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(hNCDF), intent(out) :: grp
    ! Copy the information regarding the parent ncdf
    ! We need to do this out-side (as the information
    ! required to denote the owner of the file)
    call ncdf_copy(this,grp)
    if ( .not. ncdf_participate(grp) ) return
    ! Save the group name... (we save it with hiercharal notice /"grp1"/"grp2")
    grp%grp = trim(this%grp)//"/"//trim(name)
    ! Find the group and return
    call ncdf_err(nf90_inq_grp_full_ncid(grp%f_id, trim(grp%grp), grp%id))
  end subroutine ncdf_open_grp
  !"#############################################################"
  !"###############" End of group routines"####################"
  !"#############################################################"
  !" Returns" a logical determining the participation of the node
  function ncdf_participate(this) result(participate)
    type(hNCDF), intent(in) :: this
    logical :: participate
    ! In all cases this should be the correct way to do it
    ! If wire is attached, the correct parallel setting is created
    ! If a communicator is attached the parallel flag is also set
    participate = this%parallel .or. IONode
  end function ncdf_participate
  ! These routines or functions are global available even if the NetCDF is not used...
  ! functions for concatenating strings and ncdf handles.
  function cat_char_ncdf(char,this) result(cat)
    character(len=*), intent(in) :: char
    type(hNCDF), intent(in) :: this
    character(len=len(char)+len_trim(this%name)) :: cat
    cat = char//trim(this%name)
  end function cat_char_ncdf
  function cat_ncdf_char(this,char) result(cat)
    type(hNCDF), intent(in) :: this
    character(len=*), intent(in) :: char
    character(len=len(char)+len_trim(this%name)) :: cat
    cat = trim(this%name)//char
  end function cat_ncdf_char
  subroutine ncdf_IONode(IO_Node)
    logical, intent(in) :: IO_Node
    IONode = IO_Node
  end subroutine ncdf_IONode
  subroutine ncdf_print(this)
    type(hNCDF), intent(inout) :: this
    integer :: ndims,nvars,ngatts,file_format,ngrps
    integer :: Node,Nodes
    if ( .not. ncdf_participate(this) ) return
    Node = 0
    Nodes = 1
    ! This will fail if it is not the 0th Node in the communicator
    ! For instance a subgroup in the Comm_World...
    if ( Node == 0 ) then
      write(*,"(a20,a)") "NetCDF filename:    ",trim(this%name)
      if ( len_trim(this%grp) /= 0 ) then
        write(*,"(a20,a)") "NetCDF group name:  ",trim(this%grp)
      end if
      write(*,"(a20,i7)") "NetCDF ID:          ",this%id
      if ( this%parallel ) then
        write(*,"(a20,a)") "Parallel access:    ","True"
        write(*,"(a20,tr1,i0)") "Parallel processors:",Nodes
      else
        write(*,"(a20,a)") "Parallel access:    ","False"
      end if
      if ( this%define == 0 ) then
        write(*,"(a20,a)") "In define-mode:     ","True"
      else if ( this%define == 1 ) then
        write(*,"(a20,a)") "In define-mode:     ","False"
      end if
      if ( this%f_id >= 0 ) then
        call ncdf_inq(this, dims=ndims, vars=nvars, atts=ngatts, &
            grps=ngrps, format=file_format)
        select case ( file_format )
        case ( NF90_FORMAT_CLASSIC )
          write(*,"(a20,a)") "File format:        ","Classic"
        case ( NF90_FORMAT_64BIT )
          write(*,"(a20,a)") "File format:        ","Classic 64Bit"
        case ( NF90_FORMAT_NETCDF4 )
          write(*,"(a20,a)") "File format:        ","NetCDF4"
          write(*,"(a20,i7)")"Default compression:",this%comp_lvl
        case ( NF90_FORMAT_NETCDF4_CLASSIC )
          write(*,"(a20,a)") "File format:        ","NetCDF4 Classic format"
          write(*,"(a22,i7)")"Default compression:  ",this%comp_lvl
        case default
          write(*,"(a20,a)") "File format:        ","Could not be determined"
        end select
        write(*,"(a20,i7)") "Number of dimensions:  ",ndims
        write(*,"(a20,i7)") "Number of variables:   ",nvars
        write(*,"(a20,i7)") "Number of attributes:  ",ngatts
        if ( ngrps >= 0 ) then
          write(*,"(a20,i7)") "Number of groups:      ",ngrps
        end if
      end if
        if ( iand(NF90_WRITE,this%mode) == NF90_WRITE ) &
            write(*,"(a20,a)") "NetCDF mode:        ","NF90_WRITE"
        if ( iand(NF90_NOCLOBBER,this%mode) == NF90_NOCLOBBER ) then
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_NOCLOBBER"
        else
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_CLOBBER"
      end if
      if ( iand(NF90_NOFILL,this%mode) == NF90_NOFILL ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_NOFILL"
      if ( iand(NF90_64BIT_OFFSET,this%mode) == NF90_64BIT_OFFSET ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_64BIT_OFFSET"
      if ( iand(NF90_LOCK,this%mode) == NF90_LOCK ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_LOCK"
      if ( iand(NF90_SHARE,this%mode) == NF90_SHARE ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_SHARE"
      if ( iand(NF90_NETCDF4,this%mode) == NF90_NETCDF4 ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_NETCDF4"
      if ( iand(NF90_CLASSIC_MODEL,this%mode) == NF90_CLASSIC_MODEL ) &
          write(*,"(a20,a)") "NetCDF mode:        ","NF90_CLASSIC_MODEL"
    end if
  end subroutine ncdf_print
  ! A standard die routine... It is not pretty... But it works...
  ! Recommended to be adapted!
  subroutine ncdf_die(str)
    character(len=*), intent(in) :: str
    write(0,"(2a)") 'ncdf: ',trim(str)
    write(6,"(2a)") 'ncdf: ',trim(str)
    stop
  end subroutine ncdf_die
end module netcdf_ncdf
