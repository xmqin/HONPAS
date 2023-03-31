module m_iterator

  implicit none

  private

  public :: itt1, itt2, itt3, itt4, itt5

  type :: itt1
     private
     ! beginning
     integer, pointer :: start => null()
     ! end
     integer, pointer :: end => null()
     ! integer current
     integer, pointer :: cur => null()
     ! stepping sequence
     integer, pointer :: step => null()
     ! whether or not it has stepped
     logical, pointer :: stepped => null()
  end type itt1

  type :: itt2
     private
     ! outer loop iterator
     type(itt1), pointer :: it1 => null()
     ! inner loop iterator
     type(itt1), pointer :: it2 => null()
  end type itt2

  type :: itt3
     private
     ! outer loop iterator
     type(itt1), pointer :: it1 => null()
     ! middle loop iterator
     type(itt2), pointer :: it2 => null()
  end type itt3

  type :: itt4
     private
     ! outer loop iterator
     type(itt2), pointer :: it1 => null()
     ! middle loop iterator
     type(itt2), pointer :: it3 => null()
  end type itt4

  type :: itt5
     private
     ! outer loop iterator
     type(itt2), pointer :: it1 => null()
     ! middle loop iterator
     type(itt3), pointer :: it3 => null()
  end type itt5

  public :: itt_init, itt_attach, itt_reset
  public :: itt_last, itt_end
  public :: itt_step, fitt_step
  public :: sitt_step
  public :: itt_done
  public :: itt_stepped
  public :: itt_destroy

  public :: itt_cur_step
  public :: itt_steps

  interface itt_init
     module procedure init_itt1
     module procedure init_itt2, init_itt2_of_itt1
     module procedure init_itt3, init_itt3_of_itt1
     module procedure init_itt4, init_itt4_of_itt1
     module procedure init_itt5, init_itt5_of_itt1
  end interface

  interface itt_stepped
     module procedure stepped_itt1, stepped_itt2
     module procedure stepped_itt3, stepped_itt4
     module procedure stepped_itt5
  end interface

  interface itt_last
     module procedure last_itt1, last_itt2
     module procedure last_itt3, last_itt4
     module procedure last_itt5
  end interface

  interface itt_end
     module procedure end_itt1, end_itt2
     module procedure end_itt3, end_itt4
     module procedure end_itt5
  end interface

  interface itt_attach
     module procedure attach_itt1, attach_itt2
     module procedure attach_itt3, attach_itt4
     module procedure attach_itt5
  end interface

  interface itt_reset
     module procedure reset_itt1, reset_itt2
     module procedure reset_itt3, reset_itt4
     module procedure reset_itt5
  end interface

  interface itt_done
     module procedure done_itt1, done_itt2
     module procedure done_itt3, done_itt4
     module procedure done_itt5
  end interface

  interface itt_step
     module procedure f_step_itt1, f_step_itt2
     module procedure f_step_itt3, f_step_itt4
     module procedure f_step_itt5
  end interface

  interface itt_destroy
     module procedure destroy_itt1, destroy_itt2
     module procedure destroy_itt3, destroy_itt4
     module procedure destroy_itt5
  end interface

  interface itt_cur_step
     module procedure cur_step_itt1, cur_step_itt2
     module procedure cur_step_itt3, cur_step_itt4
     module procedure cur_step_itt5
  end interface

  interface itt_steps
     module procedure steps_itt1, steps_itt2
     module procedure steps_itt3, steps_itt4
     module procedure steps_itt5
  end interface

  interface fitt_step
     module procedure f_step_itt1, f_step_itt2
     module procedure f_step_itt3, f_step_itt4
     module procedure f_step_itt5
  end interface

  interface sitt_step
     module procedure s_step_itt1, s_step_itt2
     module procedure s_step_itt3, s_step_itt4
     module procedure s_step_itt5
  end interface

contains

  pure subroutine init_itt1(this,start,end,step)
    type(itt1), intent(inout) :: this
    integer, intent(in), optional :: start , end, step
    if ( .not. associated(this%start) ) then
       allocate(this%start)
       allocate(this%end)
       allocate(this%cur)
       allocate(this%step)
       allocate(this%stepped)
    end if
    this%start = 1
    if ( present(start) ) this%start = start
    this%end = 0
    if ( present(end) ) this%end = end
    this%step = 1
    if ( present(step) ) this%step = step
    if ( this%step == 0 ) then
       !call die('iterator: stepping must be non-zero')
    end if
    call itt_reset(this)
  end subroutine init_itt1

  pure subroutine destroy_itt1(this)
    type(itt1), intent(inout) :: this
    if ( associated(this%start) ) then
       deallocate(this%start)
       nullify(this%start)
       deallocate(this%end)
       nullify(this%end)
       deallocate(this%cur)
       nullify(this%cur)
       deallocate(this%step)
       nullify(this%step)
       deallocate(this%stepped)
       nullify(this%stepped)
    end if
  end subroutine destroy_itt1

  pure subroutine reset_itt1(this)
    type(itt1), intent(inout) :: this
    this%stepped = .false.
    ! ensures that the first step is always valid
    this%cur = this%start - this%step
  end subroutine reset_itt1

  pure function last_itt1(this) result(last)
    type(itt1), intent(in) :: this
    logical :: last
    if ( this%step > 0 ) then
       last = this%cur + this%step > this%end
    else
       last = this%cur + this%step < this%end
    end if
  end function last_itt1

  pure subroutine end_itt1(this)
    type(itt1), intent(inout) :: this
    this%cur = this%end + this%step
    this%stepped = .false.
  end subroutine end_itt1

  pure subroutine attach_itt1(this,i,cur)
    type(itt1), intent(inout) :: this
    integer, pointer, optional :: i,cur
    if ( present(i) ) i => this%cur
    if ( present(cur) ) cur => this%cur
  end subroutine attach_itt1

  pure subroutine s_step_itt1(this)
    type(itt1), intent(inout) :: this

    if ( .not. itt_done(this) ) then
       ! if it ain't done, we surely will step
       this%stepped = .true.
    else
       ! we don't need to step the current iterator,
       this%stepped = .false.
       return
    end if

    this%cur = this%cur + this%step
  end subroutine s_step_itt1

  pure function done_itt1(this) result(done)
    type(itt1), intent(in) :: this
    logical :: done
    if ( this%step > 0 ) then
       done = this%cur > this%end
    else
       done = this%cur < this%end
    end if
  end function done_itt1

  pure function stepped_itt1(this) result(stepped)
    type(itt1), intent(in) :: this
    logical :: stepped
    stepped = this%stepped
  end function stepped_itt1

  function f_step_itt1(this) result(done)
    type(itt1), intent(inout) :: this
    logical :: done
    call sitt_step(this)
    done = itt_done(this)
  end function f_step_itt1

  pure function cur_step_itt1(this) result(cur_step)
    type(itt1), intent(in) :: this
    integer :: cur_step
    if ( itt_done(this) ) then
       cur_step = this%end - this%start + 1
       if ( cur_step > this%start ) &
            cur_step = cur_step / this%step
    else if ( this%stepped ) then
       cur_step = this%cur - this%start + 1
       if ( cur_step > this%start ) &
            cur_step = cur_step / this%step
    else
       cur_step = 0
    end if
  end function cur_step_itt1

  pure function steps_itt1(this) result(steps)
    type(itt1), intent(in) :: this
    integer :: steps
    steps = this%end - this%start + 1
    if ( steps > 1 ) then ! correct number of steps by the stepping
       steps = steps / this%step
    end if
  end function steps_itt1

  


  pure subroutine init_itt2(this,start1,end1,step1,start2,end2,step2)
    type(itt2), intent(inout) :: this
    integer, intent(in), optional :: start1, end1, step1, start2, end2, step2
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
       allocate(this%it2)
    end if
    call itt_init(this%it1,start=start1,end=end1,step=step1)
    call itt_init(this%it2,start=start2,end=end2,step=step2)
  end subroutine init_itt2

  subroutine init_itt2_of_itt1(this,it1,it2)
    type(itt2), intent(inout) :: this
    type(itt1), intent(inout), target :: it1, it2
    this%it1 => it1
    this%it2 => it2
    call itt_reset(this)
  end subroutine init_itt2_of_itt1

  function last_itt2(this,i) result(last)
    type(itt2), intent(in) :: this
    integer, intent(in), optional :: i
    logical :: last
    if ( present(i) ) then
       if ( i < 1 .or. 2 < i ) then
          !call die('iterator: index of iterator not existing')
       end if
       if ( i == 1 ) then
          last = itt_last(this%it1)
       else if ( i == 2 ) then
          last = itt_last(this%it2)
       end if
    else
       last = itt_last(this%it1) .and. itt_last(this%it2)
    end if
  end function last_itt2

  pure subroutine end_itt2(this)
    type(itt2), intent(inout) :: this
    call itt_end(this%it1)
    call itt_end(this%it2)
  end subroutine end_itt2

  pure subroutine reset_itt2(this)
    type(itt2), intent(inout) :: this
    call itt_reset(this%it1)
    call itt_reset(this%it2)
  end subroutine reset_itt2

  pure subroutine destroy_itt2(this)
    type(itt2), intent(inout) :: this
    if ( associated(this%it1) ) then
       call itt_destroy(this%it1)
       call itt_destroy(this%it2)
       deallocate(this%it1)
       nullify(this%it1)
       deallocate(this%it2)
       nullify(this%it2)
    end if
  end subroutine destroy_itt2

  pure subroutine attach_itt2(this,it1,i1,cur1,it2,i2,cur2)
    type(itt2), intent(inout) :: this
    type(itt1), pointer, optional :: it1, it2
    integer, pointer, optional :: i1, cur1, i2, cur2
    if ( present(it1) ) it1 => this%it1
    if ( present(it2) ) it2 => this%it2
    call itt_attach(this%it1,i=i1,cur=cur1)
    call itt_attach(this%it2,i=i2,cur=cur2)
  end subroutine attach_itt2

  pure subroutine s_step_itt2(this)
    type(itt2), intent(inout) :: this
    logical :: fs

    fs = .false.
    ! we need to ensure that the first step
    ! will also step the first counter
    if ( this%it1%step > 0 ) then
       fs = ( this%it1%cur < this%it1%start ) 
    else
       fs = ( this%it1%cur > this%it1%start ) 
    end if
    if ( fs ) call sitt_step(this%it1)
    ! try to step the outermost loop
    call sitt_step(this%it2)
    if ( itt_done(this%it2) ) then
       if ( .not. fs ) &
            call sitt_step(this%it1)
       if ( itt_done(this%it1) ) then
          ! now both it1 and it2 are done
       else
          call itt_reset(this%it2)
          call sitt_step(this%it2)
       end if
    else
       if ( .not. fs ) &
            this%it1%stepped = .false.
    end if

    if ( itt_done(this) ) then
       ! if we are done after a step
       ! it must meen that every loop is done
       ! we immediately set everything to the end-value
       ! to ensure that a step will not be caught
       call itt_end(this)
    end if
  end subroutine s_step_itt2

  pure function done_itt2(this) result(done)
    type(itt2), intent(in) :: this
    logical :: done
    ! if the outer loop is done
    done = itt_done(this%it1)
  end function done_itt2

  function f_step_itt2(this) result(done)
    type(itt2), intent(inout) :: this
    logical :: done
    call sitt_step(this)
    done = itt_done(this)
  end function f_step_itt2

  pure function stepped_itt2(this,it) result(stepped)
    type(itt2), intent(in) :: this
    integer, intent(in) :: it
    logical :: stepped
    if ( it < 1 .or. 2 < it ) then
       !call die('iterator: requested level does not exist in iterator')
    end if
    if ( it == 1 ) then
       stepped = this%it1%stepped
    else if ( it == 2 ) then
       stepped = this%it2%stepped
    end if
  end function stepped_itt2

  pure function cur_step_itt2(this) result(cur_step)
    type(itt2), intent(in) :: this
    integer :: cur_step
    cur_step = itt_cur_step(this%it1)
    if ( cur_step > 1 ) then
       cur_step = itt_steps(this%it2) * (cur_step-1) + itt_cur_step(this%it2)
    else
       cur_step = itt_cur_step(this%it2)
    end if
  end function cur_step_itt2

  pure function steps_itt2(this) result(steps)
    type(itt2), intent(in) :: this
    integer :: steps
    steps = itt_steps(this%it2) * itt_steps(this%it1)
  end function steps_itt2




  pure subroutine init_itt3(this,&
       start1,end1,step1, &
       start2,end2,step2, &
       start3,end3,step3)
    type(itt3), intent(inout) :: this
    integer, intent(in), optional :: start1, end1, step1
    integer, intent(in), optional :: start2, end2, step2
    integer, intent(in), optional :: start3, end3, step3
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
       allocate(this%it2)
    end if
    call itt_init(this%it1,start=start1,end=end1,step=step1)
    call itt_init(this%it2,start1=start2,end1=end2,step1=step2, &
         start2=start3,end2=end3,step2=step3)
  end subroutine init_itt3

  subroutine init_itt3_of_itt1(this,it1,it2,it3)
    type(itt3), intent(inout) :: this
    type(itt1), intent(inout), target :: it1, it2, it3
    this%it1 => it1
    if ( .not. associated(this%it2) ) then
       allocate(this%it2)
    end if
    this%it2%it1 => it2
    this%it2%it2 => it3
    call itt_reset(this)
  end subroutine init_itt3_of_itt1

  function last_itt3(this,i) result(last)
    type(itt3), intent(in) :: this
    integer, intent(in), optional :: i
    logical :: last
    if ( present(i) ) then
       if ( i < 1 .or. 3 < i ) then
          !call die('iterator: index of iterator not existing')
       end if
       if ( i == 1 ) then
          last = itt_last(this%it1)
       else if ( i == 2 ) then
          last = itt_last(this%it2,1)
       else if ( i == 3 ) then
          last = itt_last(this%it2,2)
       end if
    else
       last = itt_last(this%it1) .and. itt_last(this%it2)
    end if
  end function last_itt3

  pure subroutine end_itt3(this)
    type(itt3), intent(inout) :: this
    call itt_end(this%it1)
    call itt_end(this%it2)
  end subroutine end_itt3

  pure subroutine destroy_itt3(this)
    type(itt3), intent(inout) :: this
    if ( associated(this%it1) ) then
       call itt_destroy(this%it1)
       call itt_destroy(this%it2)
       deallocate(this%it1)
       nullify(this%it1)
       deallocate(this%it2)
       nullify(this%it2)
    end if
  end subroutine destroy_itt3

  pure subroutine reset_itt3(this)
    type(itt3), intent(inout) :: this
    call itt_reset(this%it1)
    call itt_reset(this%it2)
  end subroutine reset_itt3

  pure subroutine attach_itt3(this,it1,i1,cur1,it2,i2,cur2,it3,i3,cur3)
    type(itt3), intent(inout) :: this
    type(itt1), pointer, optional :: it1, it2, it3
    integer, pointer, optional :: i1, cur1, i2, cur2, i3, cur3
    if ( present(it1) ) it1 => this%it1
    if ( present(it2) ) it2 => this%it2%it1
    if ( present(it3) ) it3 => this%it2%it2
    call itt_attach(this%it1,i=i1,cur=cur1)
    call itt_attach(this%it2,i1=i2,cur1=cur2,i2=i3,cur2=cur3)
  end subroutine attach_itt3

  pure subroutine s_step_itt3(this)
    type(itt3), intent(inout) :: this
    logical :: fs

    fs = .false.
    ! we need to ensure that the first step
    ! will also step the first counter
    if ( this%it1%step > 0 ) then
       fs = ( this%it1%cur < this%it1%start ) 
    else
       fs = ( this%it1%cur > this%it1%start ) 
    end if
    if ( fs ) call sitt_step(this%it1)
    ! try to step the outermost loop
    call sitt_step(this%it2)
    if ( itt_done(this%it2) ) then
       if ( .not. fs ) &
            call sitt_step(this%it1)
       if ( itt_done(this%it1) ) then
          ! now both it1 and it2 are done
       else
          call itt_reset(this%it2)
          call sitt_step(this%it2)
       end if
    else
       if ( .not. fs ) &
            this%it1%stepped = .false.
    end if
  end subroutine s_step_itt3

  pure function done_itt3(this) result(done)
    type(itt3), intent(in) :: this
    logical :: done
    ! if the outer loop is done
    done = itt_done(this%it1)
  end function done_itt3

  function f_step_itt3(this) result(done)
    type(itt3), intent(inout) :: this
    logical :: done
    call sitt_step(this)
    done = itt_done(this)
  end function f_step_itt3

  pure function stepped_itt3(this,it) result(stepped)
    type(itt3), intent(in) :: this
    integer, intent(in) :: it
    logical :: stepped
    if ( it < 1 .or. 3 < it ) then
       !call die('iterator: requested level does not exist in iterator')
    end if
    if ( it == 1 ) then
       stepped = this%it1%stepped
    else if ( it == 2 ) then
       stepped = this%it2%it1%stepped
    else if ( it == 3 ) then
       stepped = this%it2%it2%stepped
    end if
  end function stepped_itt3

  pure function cur_step_itt3(this) result(cur_step)
    type(itt3), intent(in) :: this
    integer :: cur_step
    cur_step = itt_cur_step(this%it1)
    if ( cur_step > 1 ) then
       cur_step = itt_steps(this%it2) * (cur_step-1) + itt_cur_step(this%it2)
    else
       cur_step = itt_cur_step(this%it2)
    end if
  end function cur_step_itt3

  pure function steps_itt3(this) result(steps)
    type(itt3), intent(in) :: this
    integer :: steps
    steps = itt_steps(this%it2) * itt_steps(this%it1)
  end function steps_itt3




  pure subroutine init_itt4(this,start1,end1,step1, &
       start2,end2,step2,start3,end3,step3, &
       start4,end4,step4)
    type(itt4), intent(inout) :: this
    integer, intent(in), optional :: start1, end1, step1
    integer, intent(in), optional :: start2, end2, step2
    integer, intent(in), optional :: start3, end3, step3
    integer, intent(in), optional :: start4, end4, step4
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
       allocate(this%it3)
    end if
    call itt_init(this%it1,start1=start1,end1=end1,step1=step1, &
         start2=start2,end2=end2,step2=step2)
    call itt_init(this%it3,start1=start3,end1=end3,step1=step3, &
         start2=start4,end2=end4,step2=step4)
  end subroutine init_itt4

  subroutine init_itt4_of_itt1(this,it1,it2,it3,it4)
    type(itt4), intent(inout) :: this
    type(itt1), intent(inout), target :: it1, it2, it3, it4
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
    end if
    this%it1%it1 => it1
    this%it1%it2 => it2
    if ( .not. associated(this%it3) ) then
       allocate(this%it3)
    end if
    this%it3%it1 => it3
    this%it3%it2 => it4
    call itt_reset(this)
  end subroutine init_itt4_of_itt1

  function last_itt4(this,i) result(last)
    type(itt4), intent(in) :: this
    integer, intent(in), optional :: i
    logical :: last
    if ( present(i) ) then
       if ( i < 1 .or. 4 < i ) then
          !call die('iterator: index of iterator not existing')
       end if
       if ( i == 1 ) then
          last = itt_last(this%it1,1)
       else if ( i == 2 ) then
          last = itt_last(this%it1,2)
       else if ( i == 3 ) then
          last = itt_last(this%it3,1)
       else if ( i == 4 ) then
          last = itt_last(this%it3,2)
       end if
    else
       last = itt_last(this%it1) .and. itt_last(this%it3)
    end if
  end function last_itt4

  pure subroutine end_itt4(this)
    type(itt4), intent(inout) :: this
    call itt_end(this%it1)
    call itt_end(this%it3)
  end subroutine end_itt4

  pure subroutine destroy_itt4(this)
    type(itt4), intent(inout) :: this
    if ( associated(this%it1) ) then
       call itt_destroy(this%it1)
       call itt_destroy(this%it3)
       deallocate(this%it1)
       nullify(this%it1)
       deallocate(this%it3)
       nullify(this%it3)
    end if
  end subroutine destroy_itt4

  pure subroutine reset_itt4(this)
    type(itt4), intent(inout) :: this
    call itt_reset(this%it1)
    call itt_reset(this%it3)
  end subroutine reset_itt4

  pure subroutine attach_itt4(this,it1,i1,cur1,it2,i2,cur2, &
       it3,i3,cur3,it4,i4,cur4)
    type(itt4), intent(inout) :: this
    type(itt1), pointer, optional :: it1, it2, it3, it4
    integer, pointer, optional :: i1, cur1, i2, cur2 
    integer, pointer, optional :: i3, cur3, i4, cur4
    if ( present(it1) ) it1 => this%it1%it1
    if ( present(it2) ) it2 => this%it1%it2
    if ( present(it3) ) it3 => this%it3%it1
    if ( present(it4) ) it4 => this%it3%it2
    call itt_attach(this%it1,i1=i1,cur1=cur1,i2=i2,cur2=cur2)
    call itt_attach(this%it3,i1=i3,cur1=cur3,i2=i4,cur2=cur4)
  end subroutine attach_itt4

  pure subroutine s_step_itt4(this)
    type(itt4), intent(inout) :: this
    logical :: fs

    fs = .false.
    ! we need to ensure that the first step
    ! will also step the first counter
    if ( this%it1%it1%step > 0 ) then
       fs = ( this%it1%it1%cur < this%it1%it1%start ) 
    else
       fs = ( this%it1%it1%cur > this%it1%it1%start ) 
    end if
    if ( fs ) call sitt_step(this%it1)
    ! try to step the outermost loop
    call sitt_step(this%it3)
    if ( itt_done(this%it3) ) then
       if ( .not. fs ) &
            call sitt_step(this%it1)
       if ( itt_done(this%it1) ) then
          ! now both it1 and it3 are done
       else
          call itt_reset(this%it3)
          call sitt_step(this%it3)
       end if
    else
       if ( .not. fs ) then
          this%it1%it1%stepped = .false.
          this%it1%it2%stepped = .false.
       end if
    end if
  end subroutine s_step_itt4

  pure function done_itt4(this) result(done)
    type(itt4), intent(in) :: this
    logical :: done
    ! if the outer loop is done
    done = itt_done(this%it1)
  end function done_itt4

  function f_step_itt4(this) result(done)
    type(itt4), intent(inout) :: this
    logical :: done
    call sitt_step(this)
    done = itt_done(this)
  end function f_step_itt4

  pure function stepped_itt4(this,it) result(stepped)
    type(itt4), intent(in) :: this
    integer, intent(in) :: it
    logical :: stepped
    if ( it < 1 .or. 4 < it ) then
       !call die('iterator: requested level does not exist in iterator')
    end if
    if ( it == 1 ) then
       stepped = this%it1%it1%stepped
    else if ( it == 2 ) then
       stepped = this%it1%it2%stepped
    else if ( it == 3 ) then
       stepped = this%it3%it1%stepped
    else if ( it == 4 ) then
       stepped = this%it3%it2%stepped
    end if
  end function stepped_itt4

  pure function cur_step_itt4(this) result(cur_step)
    type(itt4), intent(in) :: this
    integer :: cur_step
    cur_step = itt_cur_step(this%it1)
    if ( cur_step > 1 ) then
       cur_step = itt_steps(this%it3) * (cur_step-1) + itt_cur_step(this%it3)
    else
       cur_step = itt_cur_step(this%it3)
    end if
  end function cur_step_itt4

  pure function steps_itt4(this) result(steps)
    type(itt4), intent(in) :: this
    integer :: steps
    steps = itt_steps(this%it3) * itt_steps(this%it1)
  end function steps_itt4




  pure subroutine init_itt5(this,start1,end1,step1, &
       start2,end2,step2,start3,end3,step3, &
       start4,end4,step4,start5,end5,step5)
    type(itt5), intent(inout) :: this
    integer, intent(in), optional :: start1, end1, step1
    integer, intent(in), optional :: start2, end2, step2
    integer, intent(in), optional :: start3, end3, step3
    integer, intent(in), optional :: start4, end4, step4
    integer, intent(in), optional :: start5, end5, step5
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
       allocate(this%it3)
    end if
    call itt_init(this%it1,start1=start1,end1=end1,step1=step1, &
         start2=start2,end2=end2,step2=step2)
    call itt_init(this%it3,start1=start3,end1=end3,step1=step3, &
         start2=start4,end2=end4,step2=step4, &
         start3=start5,end3=end5,step3=step5)
  end subroutine init_itt5

  subroutine init_itt5_of_itt1(this,it1,it2,it3,it4,it5)
    type(itt5), intent(inout) :: this
    type(itt1), intent(inout), target :: it1, it2, it3, it4, it5
    if ( .not. associated(this%it1) ) then
       allocate(this%it1)
    end if
    this%it1%it1 => it1
    this%it1%it2 => it2
    if ( .not. associated(this%it3) ) then
       allocate(this%it3)
    end if
    this%it3%it1 => it3
    this%it3%it2%it1 => it4
    this%it3%it2%it2 => it5
    call itt_reset(this)
  end subroutine init_itt5_of_itt1

  function last_itt5(this,i) result(last)
    type(itt5), intent(in) :: this
    integer, intent(in), optional :: i
    logical :: last
    if ( present(i) ) then
       if ( i < 1 .or. 5 < i ) then
          !call die('iterator: index of iterator not existing')
       end if
       if ( i == 1 ) then
          last = itt_last(this%it1,1)
       else if ( i == 2 ) then
          last = itt_last(this%it1,2)
       else if ( i == 3 ) then
          last = itt_last(this%it3,1)
       else if ( i == 4 ) then
          last = itt_last(this%it3,2)
       else if ( i == 5 ) then
          last = itt_last(this%it3,3)
       end if
    else
       last = itt_last(this%it1) .and. itt_last(this%it3)
    end if
  end function last_itt5

  pure subroutine end_itt5(this)
    type(itt5), intent(inout) :: this
    call itt_end(this%it1)
    call itt_end(this%it3)
  end subroutine end_itt5

  pure subroutine destroy_itt5(this)
    type(itt5), intent(inout) :: this
    if ( associated(this%it1) ) then
       call itt_destroy(this%it1)
       call itt_destroy(this%it3)
       deallocate(this%it1)
       nullify(this%it1)
       deallocate(this%it3)
       nullify(this%it3)
    end if
  end subroutine destroy_itt5

  pure subroutine reset_itt5(this)
    type(itt5), intent(inout) :: this
    call itt_reset(this%it1)
    call itt_reset(this%it3)
  end subroutine reset_itt5

  pure subroutine attach_itt5(this,it1,i1,cur1,it2,i2,cur2, &
       it3,i3,cur3,it4,i4,cur4,it5,i5,cur5)
    type(itt5), intent(inout) :: this
    type(itt1), pointer, optional :: it1, it2, it3, it4, it5
    integer, pointer, optional :: i1, cur1, i2, cur2 
    integer, pointer, optional :: i3, cur3, i4, cur4, i5, cur5
    if ( present(it1) ) it1 => this%it1%it1
    if ( present(it2) ) it2 => this%it1%it2
    if ( present(it3) ) it3 => this%it3%it1
    if ( present(it4) ) it4 => this%it3%it2%it1
    if ( present(it5) ) it5 => this%it3%it2%it2
    call itt_attach(this%it1,i1=i1,cur1=cur1,i2=i2,cur2=cur2)
    call itt_attach(this%it3,i1=i3,cur1=cur3,i2=i4,cur2=cur4, &
         i3=i5,cur3=cur5)
  end subroutine attach_itt5

  pure subroutine s_step_itt5(this)
    type(itt5), intent(inout) :: this
    logical :: fs

    fs = .false.
    ! we need to ensure that the first step
    ! will also step the first counter
    if ( this%it1%it1%step > 0 ) then
       fs = ( this%it1%it1%cur < this%it1%it1%start ) 
    else
       fs = ( this%it1%it1%cur > this%it1%it1%start ) 
    end if
    if ( fs ) call sitt_step(this%it1)
    ! try to step the outermost loop
    call sitt_step(this%it3)
    if ( itt_done(this%it3) ) then
       if ( .not. fs ) &
            call sitt_step(this%it1)
       if ( itt_done(this%it1) ) then
          ! now both it1 and it3 are done
       else
          call itt_reset(this%it3)
          call sitt_step(this%it3)
       end if
    else
       if ( .not. fs ) then
          this%it1%it1%stepped = .false.
          this%it1%it2%stepped = .false.
       end if
    end if
  end subroutine s_step_itt5

  pure function done_itt5(this) result(done)
    type(itt5), intent(in) :: this
    logical :: done
    ! if the outer loop is done
    done = itt_done(this%it1)
  end function done_itt5

  function f_step_itt5(this) result(done)
    type(itt5), intent(inout) :: this
    logical :: done
    call sitt_step(this)
    done = itt_done(this)
  end function f_step_itt5

  pure function stepped_itt5(this,it) result(stepped)
    type(itt5), intent(in) :: this
    integer, intent(in) :: it
    logical :: stepped
    if ( it < 1 .or. 5 < it ) then
       !call die('iterator: requested level does not exist in iterator')
    end if
    if ( it == 1 ) then
       stepped = this%it1%it1%stepped
    else if ( it == 2 ) then
       stepped = this%it1%it2%stepped
    else if ( it == 3 ) then
       stepped = this%it3%it1%stepped
    else if ( it == 4 ) then
       stepped = this%it3%it2%it1%stepped
    else if ( it == 5 ) then
       stepped = this%it3%it2%it2%stepped
    end if
  end function stepped_itt5

  pure function cur_step_itt5(this) result(cur_step)
    type(itt5), intent(in) :: this
    integer :: cur_step
    cur_step = itt_cur_step(this%it1)
    if ( cur_step > 1 ) then
       cur_step = itt_steps(this%it3) * (cur_step-1) + itt_cur_step(this%it3)
    else
       cur_step = itt_cur_step(this%it3)
    end if
  end function cur_step_itt5

  pure function steps_itt5(this) result(steps)
    type(itt5), intent(in) :: this
    integer :: steps
    steps = itt_steps(this%it3) * itt_steps(this%it1)
  end function steps_itt5

end module m_iterator
