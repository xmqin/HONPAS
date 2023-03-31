module flook_siesta

#ifdef SIESTA__FLOOK
  use flook

  use siesta_options
  use siesta_geom
  use parallel, only : IONode, Node, Nodes

  use sys, only: bye, die

#endif

  implicit none

  private

  ! Signals to LUA
  ! Right after reading initial options 
  integer, parameter, public :: LUA_INITIALIZE = 1
  ! Right before SCF step starts, but at each MD step
  integer, parameter, public :: LUA_INIT_MD = 2
  ! at the start of each SCF step
  integer, parameter, public :: LUA_SCF_LOOP = 3
  ! after each SCF has finished
  integer, parameter, public :: LUA_FORCES = 4
  ! when moving the atoms, right after the FORCES step
  integer, parameter, public :: LUA_MOVE = 5
  ! when Siesta is about to do analysis
  integer, parameter, public :: LUA_ANALYSIS = 6
  ! when SIESTA is complete, just before it exists
  integer, parameter, public :: LUA_FINALIZE = 7

#ifdef SIESTA__FLOOK

  public :: slua_init, slua_call, slua_close

  ! Internal parameters
  logical, save :: slua_run = .false.
  character(len=512), save, public :: slua_file = ' '
  ! Debugging flag for both parallel and serial debugging
  logical, save, public :: slua_debug = .false.
  ! Interactive Lua?
  logical, save :: slua_interactive = .false.

contains

  subroutine slua_init(LUA, md_lua)

    use fdf, only : fdf_get
    use m_os, only : file_exist
    
    type(luaState), intent(inout) :: LUA
    logical, intent(in) :: md_lua

    character(len=30) :: fortran_msg

    character(*), parameter :: fortran_static_lua = '&
siesta = { &
    Node = 0, &
    Nodes = 1, &
    INITIALIZE = 1, &
    INIT_MD = 2, &
    SCF_LOOP = 3, &
    FORCES = 4, &
    MOVE = 5, &
    ANALYSIS = 6, &
    FINALIZE = 7, &
    state = 0, &
    IOprint = function(self, ...) &
       if self.IONode then &
          print(...) &
       end &
    end, &
    print = function(self, ...) &
       print(...) &
    end, &
} &
IOprint = function(...) &
   siesta:IOprint(...) &
end &
siesta_comm = function(...) end'

    character(*), parameter :: unit_static_lua = '&
siesta.Units = { &
    Ang    = 1. / 0.529177, &
    eV     = 1. / 13.60580, &
    kBar   = 1. / 1.47108e5, &
    Debye  = 0.393430, &
    amu    = 2.133107, &
} &
siesta.Units.GPa = siesta.Units.kBar * 10 &
siesta.Units.Kelvin = siesta.Units.eV / 11604.45'

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg

    ! Whether the interactive lua stuff should be runned
    slua_interactive = fdf_get('LUA.Interactive', .false.)
    slua_interactive = slua_interactive .AND. IONode
    
    ! First retrieve lua file
    slua_file = fdf_get('LUA.Script',' ')
    ! Immediately return if the file is not specified...
    if ( len_trim(slua_file) == 0 ) return

    ! Default debugging only on the io-node.
    slua_debug = fdf_get('LUA.Debug',.false.)
    slua_debug = slua_debug .and. IONode 
    if ( fdf_get('LUA.Debug.MPI',.false.) ) then
       ! Only if requesting parallel debug should all processors
       ! use the debugging.
       slua_debug = .true.
    end if

    ! Check that all sees the files
    ! Currently this is a limitation of this simple
    ! flook hook.
    ! This is because the calling of the fortran
    ! function from Lua is not following the
    ! same path. One could easily provide
    ! signals to the other routines by send-catch
    ! messages, but for now this is not the scope
    slua_run = file_exist(slua_file, all = .true. )

    ! If not all processors sees the file, 
    ! we return immediately...
    if ( .not. slua_run ) then
       ! Let the user know if their system
       ! does not allow running flook
       if ( IONode .and. file_exist(slua_file) ) then
          write(*,'(a)') 'siesta-lua: WARNING'
          write(*,'(a)') 'siesta-lua: The file-system does not &
               &support file-views on all running processors.'
          write(*,'(a)') 'siesta-lua: siesta-lua CANNOT be runned on &
               &this file-system. Try running with only one node.'
          write(*,'(a)') 'siesta-lua: WARNING'
       else if ( IONode ) then
          write(*,'(a)') 'siesta-lua: WARNING'
          write(*,'(3a)') 'siesta-lua: File ', trim(slua_file), &
               ' could not be found!'
          write(*,'(a)') 'siesta-lua: WARNING'
       end if

       if ( md_lua ) then
          call die('Requested Lua-MD, however no Lua ' // &
               'script could be found !')
       end if
       
       return
    end if

    call timer('LUA', 1)
    call timer('LUA-init', 1)

    ! Initialize the Lua state
    call lua_init(LUA)

    ! Create LUA table for data container
    call lua_run(LUA, code = fortran_static_lua )
    ! Append the unit table for SIESTA unit conversion
    call lua_run(LUA, code = unit_static_lua )

    ! Register siesta calls to communicate to the lua layer
    ! Old names for backwards compatibility
    call lua_register(LUA,'siesta_get', slua_receive_siesta)
    call lua_register(LUA,'siesta_return', slua_send_siesta)

    call lua_register(LUA,'siesta_receive', slua_receive_siesta)
    call lua_register(LUA,'siesta_send', slua_send_siesta)
    ! Make local siesta.receive and siesta.send
    call lua_run(LUA, code = 'siesta.receive = siesta_receive' )
    call lua_run(LUA, code = 'siesta.send = siesta_send' )

    ! Only used for printing information about
    ! what can be retrieved
    ! This function will return different things
    ! dependen on where in the routine it is called
    call lua_register(LUA,'_internal_print_allowed', slua_siesta_print_objects)
    call lua_run(LUA, code = 'siesta.print_allowed = _internal_print_allowed' )

    write(fortran_msg,'(a,i0)') 'siesta.Node = ',Node
    call lua_run(LUA, code = fortran_msg )
    write(fortran_msg,'(a,i0)') 'siesta.Nodes = ',Nodes
    call lua_run(LUA, code = fortran_msg )
    if ( IONode ) then
       call lua_run(LUA, code = 'siesta.IONode = true' )
    else
       call lua_run(LUA, code = 'siesta.IONode = false' )
    end if

    ! Run the requested lua-script
    err_msg = " "
    call lua_run(LUA, slua_file, error = err, message=err_msg)
    if ( err /= 0 ) then
      write(*,'(a)') trim(err_msg)
      call die('LUA initialization failed, please check your Lua script!!!')
    else if ( IONode ) then
      write(*,'(/2a/)') 'LUA successfully initialized: ',trim(slua_file)
    end if

    call timer('LUA-init', 2)
    call timer('LUA', 2)

  end subroutine slua_init

  subroutine slua_call(LUA, state)
    type(luaState), intent(inout) :: LUA
    integer, intent(in) :: state
    character(len=30) :: tmp

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg

    ! To finally build the user-defined command
    integer :: n_chars
    type ll_line
      character, allocatable :: line(:)
      type(ll_line), pointer :: next => null()
    end type ll_line
    type(ll_line), target :: lines
    integer :: iostat
    
    ! Return immediately if we should not run
    if ( .not. (slua_run .or. slua_interactive) ) return

    call timer('LUA', 1)
    call timer('LUA-call', 1)

    ! Transfer the state to the lua interpreter such
    ! that decisions can be made as to which steps
    ! to take
    write(tmp,'(a,i0)') 'siesta.state = ',state
    call lua_run(LUA, code = tmp )

    if ( slua_debug ) then
      write(*,'(a,i0)') 'siesta-lua: calling siesta_comm() @ ',state
    end if

   if ( slua_interactive ) then

     ! Start an interactive session.
     ! This allows the user to inspect things but *NOT* send things
     ! Although we don't prohibit the user from doing so, it may end up
     ! breaking things as the interactive session is currently only available
     ! from the IO node.
     select case ( state )
     case ( LUA_INITIALIZE )
       tmp = 'INITIALIZE'
     case ( LUA_INIT_MD )
       tmp = 'INIT_MD'
     case ( LUA_SCF_LOOP )
       tmp = 'SCF_LOOP'
     case ( LUA_FORCES )
       tmp = 'FORCES'
     case ( LUA_MOVE )
       tmp = 'MOVE'
     case ( LUA_ANALYSIS )
       tmp = 'ANALYSIS'
     case ( LUA_FINALIZE )
       tmp = 'FINALIZE'
     end select
       
     write(*,'(/2a)') 'Entering Lua-interactive @ siesta.state = ', trim(tmp)
     write(*,'(a)') 'The following commands are available in the interactive session:'
     write(*,*) ! newline
     write(*,'(a)') '   /debug    turn on/off debugging information'
     write(*,'(a)') '   /show     show the currently collected lines of code'
     write(*,'(a)') '   /clear    clears the currently collected lines of code'
     write(*,'(a)') '   ;         ending a line in ; will also run all collected lines of code'
     write(*,'(a)') '   /run      execute the currently collected lines of code'
     write(*,'(a)') '   /cont     execute and continue Siesta'
     write(*,'(a)') '   /stop     execute and stop using the interactive session'
     write(*,'(a)') '   ^D (EOF)  closes interactive console'
     write(*,*) ! newline

     ! Run interactive Lua-shell
     call interactive_run()

   end if
   
   ! Call communicator
   if ( slua_run ) then
     call lua_run(LUA, code = 'siesta_comm()', error = err, message=err_msg )
     if ( err /= 0 ) then
       write(*,'(a)') trim(err_msg)
       call die('LUA could not run siesta_comm() without an error, please &
           &check your Lua script')
     end if
   end if

   call timer('LUA-call', 2)
   call timer('LUA', 2)
   
 contains

   subroutine interactive_run()
     character(len=512) :: line
     integer :: i_chars, i
     type(ll_line), pointer :: next
     logical :: first

     ! Total number of characters assembled
     n_chars = 0

     ! Start assembling the lines
     next => lines
     first = .true.
     interactive_loop: do 

       ! Read line
       if ( first ) then
         write(*,"(a)", advance="no") "LUA> "
       else
         write(*,"(a)", advance="no") "   > "
       end if
       read(*,"(a)",iostat=iostat) line
       if ( iostat /= 0 ) exit
       
       select case ( trim(line) )
       case ( "/debug" )
         
         if ( slua_debug ) then
           write(*,'(a)') 'Debugging OFF!'
           slua_debug = .false.
         else
           write(*,'(a)') 'Debugging ON!'
           slua_debug = .true.
         end if

         cycle interactive_loop
         
       case ( "/show" )
         
         call interactive_show()
         cycle interactive_loop

       case ( "/clear" )
         
         call interactive_clean()
         next => lines
         first = .true.
         cycle interactive_loop

       case ( "/run", ";" )
         
         ! Run and continue
         call interactive_execute()
         next => lines
         first = .true.
         cycle interactive_loop

       case ( "/cont", "/continue" )
         
         ! Run and exit interactive session
         call interactive_execute()
         exit interactive_loop

       case ( "/stop" )
         
         ! Run and exit interactive session
         call interactive_execute()
         write(*,'(a)') 'Stopping future interactive sessions!'
         slua_interactive = .false.
         
         exit interactive_loop

       end select

       ! Not first line
       first = .false.
       
       ! Add line to linked-list
       i_chars = len_trim(line)
       allocate(next%line(i_chars + 1))
       do i = 1, i_chars
         next%line(i) = line(i:i)
       end do
       next%line(i_chars+1) = char(10)
       n_chars = n_chars + i_chars + 1
       allocate(next%next)
       next => next%next

       ! In case the user asks for immediate run
       if ( line(i_chars:i_chars) == ';' ) then
         ! Run and continue
         call interactive_execute()
         next => lines
         first = .true.
       end if

     end do interactive_loop

     ! Ensure the linked-list is clean
     call interactive_clean()
     
   end subroutine interactive_run
   
   subroutine interactive_execute()
     use variable, only: cunpack
     character, allocatable :: interactive(:)
     type(ll_line), pointer :: next, nnext
     integer :: i_chars

     call interactive_collect(interactive)
     ! Clean the linked-list
     call interactive_clean()

     if ( size(interactive) > 0 ) then
       call lua_run(lua, code = cunpack(interactive), error=err, message=err_msg )
       if ( err /= 0 ) then
         write(*,'(a)') trim(err_msg)
       end if
     end if
     
     deallocate(interactive)

   end subroutine interactive_execute

   subroutine interactive_collect(interactive)
     character, intent(inout), allocatable :: interactive(:)
     type(ll_line), pointer :: next
     integer :: i_chars

     ! Now we can build the single input
     allocate(interactive(n_chars))

     ! Reset n_chars to count stuff to execute
     n_chars = 1
     
     next => lines
     do while ( allocated(next%line) )
       i_chars = size(next%line)
       ! Append line
       interactive(n_chars:n_chars+i_chars-1) = next%line(:)
       n_chars = n_chars + i_chars
       next => next%next
     end do

   end subroutine interactive_collect

   subroutine interactive_show()
     use variable, only: cunpack
     character, allocatable :: interactive(:)
     integer :: old_n_chars

     ! We have to store the current size counter
     old_n_chars = n_chars
     call interactive_collect(interactive)
     
     if ( size(interactive) > 0 ) then
       write(*,'(a)') cunpack(interactive)
     end if
     deallocate(interactive)

     n_chars = old_n_chars

   end subroutine interactive_show

   subroutine interactive_clean()
     type(ll_line), pointer :: next, nnext

     ! Clean up the linked-list and ensure we reset everything
     next => lines%next
     do while ( associated(next) )
       nnext => next%next
       nullify(next%next)
       if ( allocated(next%line) ) deallocate(next%line)
       deallocate(next)
       next => nnext
     end do
     nullify(lines%next)
     if ( allocated(lines%line) ) deallocate(lines%line)

     ! Reset!
     n_chars = 0

   end subroutine interactive_clean
   
  end subroutine slua_call

  subroutine slua_close(LUA)
    type(luaState), intent(inout) :: LUA
    ! Return immediately if we should not run
    if ( .not. slua_run ) return
    call lua_close(LUA)
  end subroutine slua_close


  ! ! ! ! ! ! ! 
  ! The remaining functions/routines are private
  ! methods.
  ! ! ! ! ! ! !


  !> Lua-exposed function to retrieve data from Siesta
  !!
  !! This function is called from Lua for asking for data transfer
  !! from Siesta.
  function slua_receive_siesta(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    type(dictionary_t) :: keys

    if ( slua_debug ) then
       write(*,'(a,i0)') '  lua: siesta_receive, Node = ',Node + 1
    end if

    call lua_init(LUA, state)

    ! Retrieve information
    call slua_get_tbl_to_dict(lua, keys)
    call slua_expand_tbl_dict(keys, options)
    call slua_expand_tbl_dict(keys, variables)

    ! open global siesta table
    tbl = lua_table(LUA, 'siesta')
    
    ! Expose the dictionary
    call slua_put_dict(tbl, options, keys)
    call slua_put_dict(tbl, variables, keys)

    call lua_close_tree(tbl)

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function slua_receive_siesta

  !> Lua-exposed function to send data to Siesta
  !!
  !! This function is called from Lua to transfer data from Lua
  !! to Siesta.
  function slua_send_siesta(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl
    type(dictionary_t) :: keys

    if ( slua_debug ) then
       write(*,'(a,i0)') '  lua: siesta_send, Node = ',Node + 1
    end if

    call lua_init(LUA, state)

    ! Retrieve information
    call slua_get_tbl_to_dict(lua, keys)
    call slua_expand_tbl_dict(keys, options)
    call slua_expand_tbl_dict(keys, variables)

    ! open global siesta table
    tbl = lua_table(LUA, 'siesta')
    
    ! Expose the dictionary
    call slua_get_dict(tbl, options, keys)
    call slua_get_dict(tbl, variables, keys)

    call lua_close_tree(tbl)

    ! Check whether the user has requested to abort siesta
    if ( 'Stop' .in. keys ) then
      call bye('LUA code has asked to stop execution!') ! Exit siesta
    end if

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function slua_send_siesta

  subroutine slua_get_tbl_to_dict(lua,keys)

    use variable
    use dictionary

    type(luaState), intent(inout) :: lua
    type(dictionary_t), intent(inout) :: keys

    type(luaTbl) :: tbl
    character(len=255) :: name
    integer :: i, N

    ! Clean the dictionary
    call delete(keys)

    ! Retrieve the table @ the top
    tbl = lua_table(LUA)

    ! Traverse all elements in the table
    N = len(tbl)
    do i = 1 , N
       call lua_get(tbl, i, name)
       keys = keys // (trim(name).kv.1)
    end do

    call lua_close(tbl)

  end subroutine slua_get_tbl_to_dict

  subroutine slua_expand_tbl_dict(keys, d)
    use dictionary

    type(dictionary_t), intent(inout) :: keys
    type(dictionary_t), intent(inout) :: d

    character(len=DICTIONARY_KEY_LENGTH) :: d_key, key
    type(dictionary_t) :: pd, pk ! pointer to the dictionary
    type(dictionary_t) :: added

    pd = .first. d
    do while ( .not. (.empty. pd) )
      ! This is the key in the dictionary
      d_key = trim(.key. pd)

      ! Loop on the keys
      ! This loop is somewhat superflous.
      ! Instead we should partition the d_key into
      ! <>.<>.<>.!
      ! where <> are build up as acceptable table calls.
      ! However, this requires extra work and it shouldn't be
      ! a bottle-neck for these small number of entries.
      pk = .first. keys
      do while ( .not. (.empty. pk) )
        
        ! This is the key in the keys
        key = trim(.key. pk)
        if ( index(trim(d_key), trim(key) // '.') == 1 ) then
          added = added // (trim(d_key).kv.1)
        end if
        pk = .next. pk
      end do
      
      pd = .next. pd
    end do

    ! Append new keys
    ! Since appending "merges" the two dictionaries
    ! we should not de-allocate 'added' since that will
    ! be done when keys are de-allocated.
    keys = keys // added

  end subroutine slua_expand_tbl_dict

  subroutine slua_put_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dictionary_t), intent(inout) :: dic
    type(dictionary_t), intent(inout), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICTIONARY_KEY_LENGTH) :: key
    type(dictionary_t) :: pd ! pointer to the dictionary
    type(variable_t) :: v

    if ( present(keys) ) then

       ! Loop over all entries in the keys dictionary
       pd = .first. keys
       do while ( .not. (.empty. pd) )
          key = trim(.key. pd)
          if ( key .in. dic ) then
             call associate(v,dic,key)
             call slua_put_var(key)
             call nullify(v) ! do not delete, simply nullify
          end if
          pd = .next. pd
       end do
       
    else
       
       ! Loop over all entries
       pd = .first. dic
       do while ( .not. (.empty. pd) )
          key = .key. pd
          call associate(v,dic,trim(key))
          call slua_put_var(key)
          call nullify(v)
          pd = .next. pd
       end do

    end if

  contains

    subroutine slua_put_var(key)
      character(len=*), intent(in) :: key
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=VARIABLE_TYPE_LENGTH) :: t
      character(len=255) :: lkey, rkey
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( slua_debug ) then
         write(*,'(4a)') '    siesta2lua; dtype = ',t,', var = ',trim(key)
      end if

      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
         ! We need to handle the variables secondly
         call key_split(key,lkey,rkey)
         if ( len_trim(rkey) == 0 ) then
            lvls = 0
            rkey = lkey
         else
            call lua_open(tbl,lkey,lvls = lvls)
         end if
      end select
      select case ( t )
      case ( 'a1' )
         call associate(a1,v)
         lkey = cunpack(a1)
         call lua_set(tbl,rkey,lkey(1:len_trim(lkey)))
      case ( 'b0' )
         call associate(b0,v)
         call lua_set(tbl,rkey,b0)
      case ( 'b1' )
         call associate(b1,v)
         call lua_set(tbl,key,b1)
      case ( 'b2' ) 
         call associate(b2,v)
         call lua_set(tbl,key,b2)
      case ( 'i0' ) 
         call associate(i0,v)
         call lua_set(tbl,rkey,i0)
      case ( 'i1' ) 
         call associate(i1,v)
         call lua_set(tbl,key,i1)
      case ( 'i2' ) 
         call associate(i2,v)
         call lua_set(tbl,key,i2)
      case ( 's0' ) 
         call associate(s0,v)
         call lua_set(tbl,rkey,s0)
      case ( 's1' ) 
         call associate(s1,v)
         call lua_set(tbl,key,s1)
      case ( 's2' ) 
         call associate(s2,v)
         call lua_set(tbl,key,s2)
      case ( 'd0' ) 
         call associate(d0,v)
         call lua_set(tbl,rkey,d0)
      case ( 'd1' ) 
         call associate(d1,v)
         call lua_set(tbl,key,d1)
      case ( 'd2' ) 
         call associate(d2,v)
         call lua_set(tbl,key,d2)
      end select
       
      call lua_close(tbl,lvls = lvls)
    end subroutine slua_put_var

  end subroutine slua_put_dict

  subroutine slua_get_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dictionary_t), intent(inout) :: dic
    type(dictionary_t), intent(in), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICTIONARY_KEY_LENGTH) :: key
    type(dictionary_t) :: pd ! pointer to the dictionary
    type(variable_t) :: v

    if ( present(keys) ) then

       ! Loop over all entries in the keys dictionary
       pd = .first. keys
       do while ( .not. (.empty. pd) )
          key = .key. pd
          if ( key .in. dic ) then
             call associate(v,dic,trim(key))
             call slua_get_var(key)
             call nullify(v)
          end if
          pd = .next. pd
       end do
       
    else
       
       ! Loop over all entries
       pd = .first. dic
       do while ( .not. (.empty. pd) )
          key = .key. pd
          call associate(v,dic,trim(key))
          call slua_get_var(key)
          call nullify(v)
          pd = .next. pd
       end do

    end if

  contains

    subroutine slua_get_var(key)
      character(len=*), intent(in) :: key
      character(len=256) :: V0
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=VARIABLE_TYPE_LENGTH) :: t
      character(len=255) :: lkey, rkey
      integer :: na1
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( slua_debug ) then
         write(*,'(4a)') '    lua2siesta; dtype = ',t,', var = ',trim(key)
      end if
      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
         ! We need to handle the variables secondly
         call key_split(key,lkey,rkey)
         if ( len_trim(rkey) == 0 ) then
            lvls = 0
            rkey = lkey
         else
            call lua_open(tbl,lkey,lvls = lvls)
         end if
      end select

      select case ( t )
      case ( 'a1' )
         call associate(a1,v)
         call lua_get(tbl,rkey,V0)
         na1 = len_trim(V0)
         a1 = ' '
         a1(1:na1) = cpack(V0(1:na1))
      case ( 'b0' ) 
         call associate(b0,v)
         call lua_get(tbl,rkey,b0)
      case ( 'b1' ) 
         call associate(b1,v)
         call lua_get(tbl,key,b1)
      case ( 'b2' ) 
         call associate(b2,v)
         call lua_get(tbl,key,b2)
      case ( 'i0' ) 
         call associate(i0,v)
         call lua_get(tbl,rkey,i0)
      case ( 'i1' ) 
         call associate(i1,v)
         call lua_get(tbl,key,i1)
      case ( 'i2' ) 
         call associate(i2,v)
         call lua_get(tbl,key,i2)
      case ( 's0' ) 
         call associate(s0,v)
         call lua_get(tbl,rkey,s0)
      case ( 's1' ) 
         call associate(s1,v)
         call lua_get(tbl,key,s1)
      case ( 's2' ) 
         call associate(s2,v)
         call lua_get(tbl,key,s2)
      case ( 'd0' ) 
         call associate(d0,v)
         call lua_get(tbl,rkey,d0)
      case ( 'd1' ) 
         call associate(d1,v)
         call lua_get(tbl,key,d1)
      case ( 'd2' ) 
         call associate(d2,v)
         call lua_get(tbl,key,d2)
      end select

      call lua_close(tbl,lvls = lvls)
    end subroutine slua_get_var

  end subroutine slua_get_dict
  
  subroutine key_split(key,lkey,rkey)
    character(len=*), intent(in) :: key
    character(len=*), intent(inout) :: lkey, rkey
    integer :: i
    i = index(key,'.',back=.true.)
    if ( i > 0 ) then
       lkey = trim(adjustl(key(1:i-1)))
       rkey = key(i+1:)
    else
       lkey = trim(adjustl(key))
       rkey = ' '
    end if
  end subroutine key_split


  function slua_siesta_print_objects(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    use parallel, only : Node

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(dictionary_t) :: et
    character(len=DICTIONARY_KEY_LENGTH) :: key
    character(len=*), parameter :: fmt = '(tr2,a,'','')'

    ! Currently we only let the current io-node
    ! print out information.
    nret = 0
    if ( .not. slua_debug .and. .not. IONode ) return

    ! Print out information
    write(*,'(a)') '-- siesta table structure available in LUA'
    write(*,'(a)') 'siesta = {'
    write(*,'(tr2,a,i0,'','')') 'Node = ',Node + 1

    ! Loop across all keys in the dictionaries
    et = .first. options
    do while ( .not. (.empty. et) )
       key = .key. et
       write(*,fmt) trim(key)
       et = .next. et
    end do

    ! Loop across all keys in the dictionaries
    et = .first. variables
    do while ( .not. (.empty. et) )
       key = .key. et
       write(*,fmt) trim(key)
       et = .next. et
    end do
    
    write(*,'(a)') '}'

  end function slua_siesta_print_objects

#endif

end module flook_siesta
