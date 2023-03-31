!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This module enables the pivoting of orbitals such that a
! smallest bandwidth solution can be created.

module m_ts_pivot

  use precision, only : dp, i8b
  use m_region

  use m_ts_electype

  implicit none

  public :: ts_pivot
  public :: crt_El_priority
  public :: consecutive_Elec_orb
  
  private
  
contains

  subroutine ts_pivot( dit, sp, &
       N_Elec, Elecs, &
       cell, na_u, xa, lasto, &
       r_pvt, pvt_str, extend)
    
    use alloc, only : re_alloc, de_alloc
    use parallel, only : IONode
    use fdf, only : fdf_get, leqi
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use intrinsic_missing, only: VNORM, VEC_PROJ, VEC_PROJ_SCA

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use m_sparsity_handling

    use geom_helper, only: iaorb
    
    use m_pivot

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    type(OrbitalDistribution), intent(in) :: dit
    ! The sparse pattern we wish to pivot
    ! This sparsity pattern *MUST* be in a unit-cell format
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    ! Initially the region that we want to pivot
    type(tRgn), intent(inout) :: r_pvt
    ! The pivoting string
    character(len=*), intent(inout) :: pvt_str
    ! Whether the region should extend beyond that
    ! found by connectivity
    logical, intent(in), optional :: extend
    
    ! The temporary sparsity
    logical :: lextend
    type(Sparsity) :: tmp_Sp
    character(len=len(pvt_str)) :: str_tmp, from_elec

    integer :: i, iEl, n, n_pvt, j
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp, r_tmp2, c_pvt
    type(tRgn) :: r_Els, priority

    ! For certain pivoting schemes it can be good to make
    ! dependencies for the tight-binding model
    ! This pvt_option decides whether a
    !   fan, fan2d/front
    ! pivoting scheme is used.
    integer, parameter :: PVT_NONE = 0
    integer, parameter :: PVT_FAN = 1
    integer, parameter :: PVT_FRONT = 2
    integer, parameter :: PVT_FAN_NONE = -2
    integer, parameter :: PVT_FAN_MEAN = 0
    integer, parameter :: PVT_FAN_MIN = -1
    integer, parameter :: PVT_FAN_MAX = 1
    integer :: pvt_option
    logical, allocatable :: r_logical(:)
    logical :: pvt_orb, orb_1, is_rev, is_priority
    integer :: fan_option
    integer :: fan1, fan2
    real(dp) :: p_n(3,2), p_cr(3), p_center(3,2), tmp3(3)
    real(dp) :: p_cc(3,2)
    ! For linear-least-squares problem
    real(dp) :: llsA(3,2), llsB(3), work(20)

    call timer('pivot', 1)

    ! We keep this copy and edit this to
    ! create the correct pivoting string
    lextend = .true.
    if ( present(extend) ) lextend = extend

    ! If there is a fan option, grab it
    pvt_option = PVT_NONE
    fan_option = -2
    if ( str_contain(pvt_str,'fan-mean') ) then
       pvt_option = PVT_FAN
       fan_option = PVT_FAN_MEAN
    else if ( str_contain(pvt_str,'fan-min') ) then
       pvt_option = PVT_FAN
       fan_option = PVT_FAN_MIN
    else if ( str_contain(pvt_str,'fan-max') ) then
       pvt_option = PVT_FAN
       fan_option = PVT_FAN_MAX
    else if ( str_contain(pvt_str,'front') ) then
       ! This forces the pvt_option
       ! to only consider the normal vector to the
       ! first electrode plane
       pvt_option = PVT_FRONT
       fan_option = PVT_FAN_MAX
    else if ( str_contain(pvt_str, 'fan') ) then
       pvt_option = PVT_FAN
    end if

    ! If the electrode is specified, then add it
    from_elec = ' '
    do i = 1, N_Elec
       if ( str_contain(pvt_str, '+' // trim(Elecs(i)%name), &
            delete=.false.) ) then
          from_elec = trim(from_elec) // '+' // trim(Elecs(i)%name)
          ! Electrode orbitals in the device
          call rgn_append(r_Els, Elecs(i)%o_inD, r_Els)
       end if
    end do
    if ( rgn_size(r_Els) == 0 ) then
       from_elec = '+' // trim(Elecs(1)%name)
       call rgn_copy(Elecs(1)%o_inD, r_Els)
    end if

    pvt_orb = str_contain(pvt_str, 'orb')
    ! variable to check for single orbital cases
    orb_1 = na_u == lasto(na_u)
    
    if ( .not. pvt_orb ) then
       ! Simply pop the atom keyword
       pvt_orb = str_contain(pvt_str,'atom') ! just to pop 'atom'
       
       ! In case we have a tight-binding scheme with 1 orbital
       ! per atom, no need to act on orbital sparsity
       ! pattern (they are equivalent).
       pvt_orb = orb_1

    end if

    if ( pvt_orb ) then
       str_tmp = 'orb'

       ! Copy pivoting region
       call rgn_copy(r_pvt,c_pvt)

       tmp_Sp = sp

    else
       str_tmp = 'atom'

       ! Convert pivoting region to atomic 
       call rgn_Orb2Atom(r_pvt,na_u,lasto,c_pvt)
       
       pvt_orb = str_contain(pvt_str,'atom') ! just to pop 'atom'
       pvt_orb = .false. ! it isn't orbital sorted

       call rgn_Orb2Atom(r_Els,na_u,lasto,r_tmp)
       call rgn_copy(r_tmp,r_Els)
       call rgn_delete(r_tmp)
       call SpOrb_to_SpAtom(dit,sp,na_u,lasto,tmp_Sp)
       ! *** the distribution will always
       !     be bigger than for the atoms, hence we need
       !     not re-construct it
       ! ***

    end if

    ! The size of the sparsity pattern
    ! Note this is WITH buffer atoms...
    n = nrows_g(tmp_Sp)
    ! This is without buffer atoms
    n_pvt = c_pvt%n

    ! Create priority list for electrodes
    call rgn_init(priority,n)
    call crt_El_priority(N_Elec, Elecs, priority, &
         na_u, lasto, is_orb = pvt_orb )

    ! Sort the pivoting region
    call rgn_sort(c_pvt)

    ! Sort the electrode region to make it faster
    call rgn_sort(r_Els)

    if ( str_contain(pvt_str,'none') ) then
       str_tmp = 'orb+none'

       ! do nothing, the pivoting array
       ! is already arranged correctly
       ! We simply copy the r_pvt to c_pvt
       ! a 'none' pivoting is equivalent
       ! to an 'orb+none'
       call rgn_copy(r_pvt, c_pvt)
       pvt_orb = .true.

       ! We have to have the order of largest string to be able to distinguish them.
       
    else if ( str_contain(pvt_str,'rev-CM+priority') ) then
       str_tmp = trim(str_tmp)//'+rev-CM+priority' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_CUTHILL_MCKEE, c_pvt, start = r_Els , &
            priority = priority%r, only_sub = .not. lextend)

    else if ( str_contain(pvt_str,'rev-CM') ) then
       str_tmp = trim(str_tmp)//'+rev-CM' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_CUTHILL_MCKEE, c_pvt, start = r_Els, &
            only_sub = .not. lextend)

    else if ( str_contain(pvt_str,'CM+priority') ) then
       str_tmp = trim(str_tmp)//'+CM+priority' // trim(from_elec)
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_CUTHILL_MCKEE, c_pvt, start = r_Els, &
            priority = priority%r, only_sub = .not. lextend)

    else if ( str_contain(pvt_str,'CM') ) then
       str_tmp = trim(str_tmp)//'+CM' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_CUTHILL_MCKEE, c_pvt, start = r_Els, &
            only_sub = .not. lextend)

    else if ( str_contain(pvt_str,'rev-GGPS+priority') ) then
       str_tmp = trim(str_tmp)//'+rev-GGPS+priority'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_GGPS, c_pvt, &
            priority = priority%r )

    else if ( str_contain(pvt_str,'rev-GGPS') ) then
       str_tmp = trim(str_tmp)//'+rev-GGPS'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_GGPS, c_pvt)

    else if ( str_contain(pvt_str,'GGPS+priority') ) then
       str_tmp = trim(str_tmp)//'+GGPS+priority'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_GGPS, c_pvt, &
            priority = priority%r )
       
    else if ( str_contain(pvt_str,'GGPS') ) then
       str_tmp = trim(str_tmp)//'+GGPS'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_GGPS, c_pvt)

    else if ( str_contain(pvt_str,'rev-GPS+priority') ) then
       str_tmp = trim(str_tmp)//'+rev-GPS+priority'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_GPS, c_pvt, &
            priority = priority%r )
       
    else if ( str_contain(pvt_str,'rev-GPS') ) then
       str_tmp = trim(str_tmp)//'+rev-GPS'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_GPS, c_pvt)
       
    else if ( str_contain(pvt_str,'GPS+priority') ) then
       str_tmp = trim(str_tmp)//'+GPS+priority'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_GPS, c_pvt, &
            priority = priority%r )
       
    else if ( str_contain(pvt_str,'GPS') ) then
       str_tmp = trim(str_tmp)//'+GPS'
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_GPS, c_pvt)

    else if ( str_contain(pvt_str,'rev-PCG+priority') ) then
       str_tmp = trim(str_tmp)//'+rev-PCG+priority' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_PCG, c_pvt, &
            start = r_Els, priority = priority%r, only_sub = .not. lextend)

    else if ( str_contain(pvt_str,'rev-PCG') ) then
       str_tmp = trim(str_tmp)//'+rev-PCG' // trim(from_elec)
       
       call sp_pvt(n,tmp_Sp,r_pvt, PVT_REV_PCG, c_pvt, start = r_Els, &
            only_sub = .not. lextend)
       
    else if ( str_contain(pvt_str,'PCG+priority') ) then
       str_tmp = trim(str_tmp)//'+PCG+priority' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_PCG, c_pvt, &
            start = r_Els, priority = priority%r, only_sub = .not. lextend )

    else if ( str_contain(pvt_str,'PCG') ) then
       str_tmp = trim(str_tmp)//'+PCG' // trim(from_elec)

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_PCG, c_pvt, start = r_Els, only_sub = .not. lextend )

#ifdef SIESTA__METIS
    else if ( str_contain(pvt_str,'metis+priority') .or. &
        str_contain(pvt_str,'nodend+priority') ) then
       str_tmp = trim(str_tmp)//'+nodend+priority'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_NODEND, c_pvt, &
            priority = priority%r)
    else if ( str_contain(pvt_str,'metis') .or. &
        str_contain(pvt_str,'nodend') ) then
       str_tmp = trim(str_tmp)//'+nodend'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_NODEND, c_pvt)

     else if ( str_contain(pvt_str,'partgraphkway+priority') ) then
       str_tmp = trim(str_tmp)//'+partgraphkway+priority'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_PARTGRAPHKWAY, c_pvt, &
           priority = priority%r)
       
    else if ( str_contain(pvt_str,'partgraphkway') ) then
       str_tmp = trim(str_tmp)//'+partgraphkway'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_PARTGRAPHKWAY, c_pvt)

    else if ( str_contain(pvt_str,'partgraphrecursive+priority') ) then
       str_tmp = trim(str_tmp)//'+partgraphrecursive+priority'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_PARTGRAPHRECURSIVE, c_pvt, &
           priority = priority%r)
       
    else if ( str_contain(pvt_str,'partgraphrecursive') ) then
       str_tmp = trim(str_tmp)//'+partgraphrecursive'

       call sp_pvt(n,tmp_Sp,r_pvt, PVT_METIS_PARTGRAPHRECURSIVE, c_pvt)

#endif
    else if ( str_contain(pvt_str,'CG') .or. pvt_option == PVT_NONE ) then

       is_rev = str_contain(pvt_str, 'rev-')
       is_priority = str_contain(pvt_str, '+priority')
       if ( is_rev ) &
            str_tmp = trim(str_tmp)//'+rev-CG'
       if ( is_priority ) &
            str_tmp = trim(str_tmp)//'+priority'

       call rgn_delete(r_Els, r_tmp)

       ! Allocate logical for all the elements
       allocate(r_logical(n))
       r_logical = .false.
       
       ! Collect the electrode orbitals
       iEl = 0
       if ( N_Elec >= 1 ) then
         
         do i = 1 , N_Elec
           if ( .not. str_contain(pvt_str,Elecs(i)%name) ) cycle
           str_tmp = trim(str_tmp)//'+'//trim(Elecs(i)%name)

           if ( pvt_orb ) then
             call rgn_assoc(r_tmp, Elecs(i)%o_inD)
           else
             call rgn_Orb2Atom(Elecs(i)%o_inD,na_u,lasto,r_tmp)
           end if

           do j = 1, r_tmp%n
             r_logical(r_tmp%r(j)) = .true.
           end do

           call rgn_append(r_Els, r_tmp, r_Els)

           ! Sort this region
           call rgn_sp_sort(r_Els, tmp_Sp, r_tmp, R_SORT_MAX_FRONT, r_logical=r_logical)

           iEl = iEl + 1
         end do

         if ( pvt_orb ) then
           call rgn_nullify(r_tmp)
         else
           call rgn_delete(r_tmp)
         end if

       end if
       
       if ( iEl == 0 ) then
         ! Default to the first electrode
         str_tmp = trim(str_tmp)//'+'//trim(Elecs(1)%name)

         if ( pvt_orb ) then
           call rgn_copy(Elecs(1)%o_inD, r_Els)
         else
           call rgn_Orb2Atom(Elecs(1)%o_inD,na_u,lasto,r_Els)
         end if

         do j = 1, r_Els%n
           r_logical(r_Els%r(j)) = .true.
         end do

         ! Sort this region
         call rgn_sp_sort(r_Els, tmp_Sp, r_Els, R_SORT_MAX_FRONT, r_logical=r_logical)

       end if

       deallocate(r_logical)

       ! do pivoting
       i = PVT_CONNECT
       if ( is_rev ) i = PVT_REV_CONNECT

       if ( is_priority ) then
          call sp_pvt(n,tmp_Sp,r_pvt, i, c_pvt, start=r_Els, &
               priority = priority%r, only_sub = .not. lextend )
       else
          call sp_pvt(n,tmp_Sp,r_pvt, i, c_pvt, start=r_Els, only_sub = .not. lextend )
       end if


    else ! the user *must* have supplied an electrode with a specific
       ! pivoting option (fan/front/...)

       ! prepare the initial pivoting region
       call rgn_init(r_pvt, c_pvt%n)
       r_pvt%n = 0
       
       ! Create exception region
       call rgn_range(priority,1,n)
       call rgn_complement(c_pvt,priority,priority)

       ! Allocate logical for all the elements
       allocate(r_logical(n))
       r_logical(:) = .false.

       call rgn_delete(r_tmp)
       
       ! Figure out which electrode(s) has been given
       ! as a starting point
       ! Note that for several electrodes one could
       ! possibly get a better sparsity pattern by 
       ! using several as a starting point
       iEl = 0
       fan1 = 0
       fan2 = 0
       do i = 1 , N_Elec
          if ( str_contain(pvt_str,Elecs(i)%name) .or. N_Elec == 1 ) then
             str_tmp = trim(str_tmp)//'+'//trim(Elecs(i)%name)
             if ( fan1 == 0 ) fan1 = i

             if ( pvt_orb ) then
                call rgn_assoc(r_tmp, Elecs(i)%o_inD)
             else
                call rgn_Orb2Atom(Elecs(i)%o_inD,na_u,lasto,r_tmp)
             end if

             if ( .not. rgn_push(r_pvt, r_tmp) ) then
                call die('ts_pivot: programming error -- 1')
             end if
             
             do j = 1, r_tmp%n
               r_logical(r_tmp%r(j)) = .true.
             end do

             ! Sort this region
             call rgn_sp_sort(r_pvt, tmp_Sp, r_tmp, R_SORT_MAX_FRONT, r_logical=r_logical)
             
             iEl = iEl + 1
          else
             if ( fan2 == 0 .and. fan1 /= 0 ) fan2 = i
          end if
       end do

       if ( pvt_orb ) then
         call rgn_nullify(r_tmp)
       else
         call rgn_delete(r_tmp)
       end if

       if ( iEl == 0 ) then
          print *,'Currently assembled options: ',trim(str_tmp)
          print *,'Remaining option string: ',trim(pvt_str)
          call die('Could not find any electrode in &
               &BTD.Pivot in the list of electrodes, &
               &please correct sorting method.')
       end if

       ! The user has explicitly requested a fan-pivoting scheme
       
       ! Calculate the normal vector along the semi-infinite direction
       ! of the first electrode.
       if ( Elecs(fan1)%t_dir > 3 ) then
         call die('Pivoting fan/front/* with real-space self-energies will not work')
       end if
       p_n(:,1) = Elecs(fan1)%cell(:,Elecs(fan1)%t_dir)
       p_n(:,1) = p_n(:,1) / VNORM(p_n(:,1))
       
       ! In case the semi-infinite direction is positive,
       ! we flip the vector
       if ( Elecs(fan1)%inf_dir == INF_POSITIVE ) then
          p_n(:,1) = -p_n(:,1)
       end if
       
       ! Find the middle of the electrode atoms (in the device region)!
       call rgn_Orb2Atom(Elecs(fan1)%o_inD,na_u,lasto,r_tmp)
       llsB = 0._dp
       do i = 1 , r_tmp%n
          llsB = llsB + xa(:,r_tmp%r(i))
       end do
       llsB = llsB / r_tmp%n ! take average
       ! llsB is the middle of the electrode fan1
       
       ! This locates the atom which has the largest projection vector
       ! on the normal vector to the atomic plane
       ! Note that we subtract the middle of the electrode to
       ! have the origo at the middle of the atoms.
       iEl = 0
       work(1) = -huge(0._dp)
       do i = 1 , r_tmp%n
          ! we take the largest one by projecting onto the vectors
          tmp3 = xa(:,r_tmp%r(i)) - llsB
          work(2) = VEC_PROJ_SCA(p_n(:,1),tmp3)
          ! this means that they point in the same direction
          if ( work(2) > work(1) ) then
             iEl = r_tmp%r(i)
             work(1) = work(2)
          end if
       end do
       if ( iEl == 0 ) then
          print *,'Electrode: ',fan1,trim(elecs(fan1)%name)
          call die('Using the fan method for BTD matrices seems like &
               &a bad choice for this system.')
       end if
       ! Get the actual top-atom
       work(1:3) = xa(:,iEl)

       ! Calculate the electrode plane middle coordinate
       p_center(:,1) = 0._dp
       do i = 1 , r_tmp%n
          tmp3 = xa(:,r_tmp%r(i))
          tmp3 = tmp3 - VEC_PROJ(p_n(:,1),tmp3-work(1:3))
          p_center(:,1) = p_center(:,1) + tmp3
       end do
       p_center(:,1) = p_center(:,1) / r_tmp%n

       ! Temporary clean-up
       call rgn_delete(r_tmp)


       ! Set the second fan value in case it hasn't been set
       if ( fan2 == 0 ) then
          ! Correct fan2 index
          if ( fan1 < N_Elec ) then
             fan2 = fan1 + 1
          else
             fan2 = 1
          end if
       end if

       ! Calculate the plane-crossing between the first two electrodes
       if ( Elecs(fan2)%t_dir > 3 ) then
         call die('Pivoting fan/front/* with real-space self-energies will not work')
       end if
       p_n(:,2) = Elecs(fan2)%cell(:,Elecs(fan2)%t_dir)
       p_n(:,2) = p_n(:,2) / VNORM(p_n(:,2))
       ! 2nd electrode points outwards
       if ( Elecs(fan2)%inf_dir == INF_POSITIVE ) then
          p_n(:,2) = -p_n(:,2)
       end if

       ! Create vector pointing along the line intersecting the
       ! two planes.
       call cross(p_n(:,1),p_n(:,2),p_cr)

       ! Now we check for 
       ! If the length of the cross-product vector is too small,
       ! they *must* be parallel planes.
       ! Then we remove the pivoting option.
       ! We do the pivoting check later...

       if ( VNORM(p_cr) > 1.e-6_dp .and. pvt_option == PVT_FAN ) then

          ! Default to mean
          if ( fan_option == PVT_FAN_NONE ) fan_option = PVT_FAN_MEAN

       else

          ! The electrodes are parallel
          ! so we might as well only use the frontal search
          pvt_option = PVT_FRONT

          ! Default to max
          if ( fan_option == PVT_FAN_NONE ) fan_option = PVT_FAN_MAX

       end if

       select case ( pvt_option )
       case ( PVT_FAN )
          select case ( fan_option )
          case ( PVT_FAN_MIN )
             str_tmp = trim(str_tmp)//'+fan-min'
          case (  PVT_FAN_MEAN ) 
             str_tmp = trim(str_tmp)//'+fan-mean'
          case ( PVT_FAN_MAX )
             str_tmp = trim(str_tmp)//'+fan-max'
          end select
       case ( PVT_FRONT )
          str_tmp = trim(str_tmp)//'+front'
       end select
       
       if ( pvt_option == PVT_FAN ) then

          ! Normalize vector for line
          p_cr = p_cr / VNORM(p_cr)

          ! do end-fan electrode
          call rgn_Orb2Atom(Elecs(fan2)%o_inD,na_u,lasto,r_tmp)
          llsB = 0._dp
          do i = 1 , r_tmp%n
             llsB = llsB + xa(:,r_tmp%r(i))
          end do
          llsB = llsB / r_tmp%n

          ! Calculate the bottom atom
          iEl = 0
          work(1) = huge(0._dp)
          do i = 1 , r_tmp%n
             ! we take the largest one by projecting onto 
             tmp3 = xa(:,r_tmp%r(i)) - llsB
             work(2) = VEC_PROJ_SCA(p_n(:,2),tmp3)
             ! this means that they point in the same direction
             if ( work(2) < work(1) ) then
                iEl = r_tmp%r(i)
                work(1) = work(2)
             end if
          end do
          if ( iEl == 0 ) then
             print *,'Electrode: ',fan2,trim(elecs(fan2)%name)
             call die('Using the fan method for BTD matrices seems like &
                  &a bad choice for this system.')
          end if
          ! Get the actual top-atom
          work(1:3) = xa(:,iEl)

          ! Calculate the electrode plane middle coordinate
          p_center(:,2) = 0._dp
          do i = 1 , r_tmp%n
             tmp3 = xa(:,r_tmp%r(i))
             tmp3 = tmp3 - VEC_PROJ(p_n(:,2),tmp3 - work(1:3))
             p_center(:,2) = p_center(:,2) + tmp3
          end do
          p_center(:,2) = p_center(:,2) / r_tmp%n

          ! Temporary clean-up
          call rgn_delete(r_tmp)

          ! Now we need to find the closests point that intersects
          ! the line that touches both planes

          ! Solve the linear least square problem
          ! to find the x_0 point on the line that intersects
          ! the planes
          llsA = p_n
          llsB(1) = sum( p_n(:,1) * p_center(:,1) )
          llsB(2) = sum( p_n(:,2) * p_center(:,2) )
          llsB(3) = 0._dp
          call dgels('T',3,2,1,llsA,3,llsB,3,work,20,i)
          if ( i /= 0 ) then
             write(*,'(a,i0)')'dgels: ERROR: ',i
             call die('Could not solve dgels problem, &
                  &add nofan to BTD.Pivot')
          end if

          ! We now have the point x_0 on the line
          tmp3 = llsB

          ! Calculate t from minimizing:
          !  \sum_i l_i^2 = \sum_i|x0 + t * p_cr - c_i|^2
          llsA(:,1) = p_center(:,1) - tmp3
          llsA(:,2) = p_center(:,2) - tmp3
          work(2) = sum( llsA(:,1) * p_cr )
          work(3) = sum( llsA(:,2) * p_cr )
          work(1) = .25_dp * (work(2)+work(3)) / sum( p_cr**2 )

          ! llsB is now the point that minimizes the distance between
          ! c_1, c_2 and the line that crosses the planes
          llsB = tmp3 + work(1) * p_cr

          ! Create the two vectors pointing from
          ! the closest point to the center of the surface
          ! plane
          p_cc(:,1) = p_center(:,1) - llsB
          p_cc(:,1) = p_cc(:,1) / VNORM(p_cc(:,1))
          p_cc(:,2) = p_center(:,2) - llsB
          p_cc(:,2) = p_cc(:,2) / VNORM(p_cc(:,2))

!!$          print '(a,3(tr1,e10.4))','x_0:',tmp3
!!$          print '(a,3(tr1,e10.4))','C:',llsB
!!$          print '(a,3(tr1,e10.4))','t:',work(1)
!!$          print '(a,3(tr1,e10.4))','p_c1:',p_center(:,1)
!!$          print '(a,3(tr1,e10.4))','p_c2:',p_center(:,2)
!!$          print '(a,3(tr1,e10.4))','p_c:',p_cr
!!$          print '(a,e10.4)','Angle between plane:', &
!!$               acos( sum( p_cc(:,1)*p_cc(:,2) ) ) * 180._dp / 3.1415926353_dp
          
       end if

       
       ! Loop the connectivity graph
       do 

          ! Create attached region starting from electrodes specified
          call rgn_sp_connect(r_pvt, dit, tmp_Sp, r_tmp, except=priority )

          if ( r_tmp%n == 0 .and. r_pvt%n /= n_pvt .and. lextend ) then
             ! This check ensures that in case there is vacuum
             ! (i.e. the previous segment does not connect to anything)
             ! we will still form the full sparsity pattern.
             ! However, we select a "random" orbital/atom which may be
             ! very sub-optimal!

             ! we need to ensure that all (even non-connected)
             ! orbitals are taken into account
             do i = 1 , n_pvt
                ! if it already exists, skip it
                if ( in_rgn(r_pvt,c_pvt%r(i)) ) cycle

                call rgn_init(r_tmp,1,val=c_pvt%r(i))
                exit

             end do

             if ( r_tmp%n /= 0 ) then

                if ( pvt_orb .and. .not. orb_1 ) then
                  ! get the atomic index to create atomic range of all orbitals
                  ! on the given atom
                  i = iaorb(r_tmp%r(1),lasto)
                  call rgn_range(r_tmp,lasto(i-1)+1,lasto(i))
                else
                  i = r_tmp%r(1)
                end if
                if ( IONode ) then
                   write(*,'(/,a,i0,a)')'WARNING: &
                        &Random atom ',i,' has been added &
                        &due to non-completeness of the connectivity &
                        &graph.'
                   write(*,'(a)')'WARNING: Expect sub-optimal &
                        &BTD format.'
                end if

                ! Assert that none of the orbitals exist in the, region
                do i = 1 , r_tmp%n
                   if ( in_rgn(r_pvt,r_tmp%r(i)) ) then
                      call die('This is extremely difficult. &
                           &Please do not sort the BTD format as &
                           &it cannot figure out what to do.')
                   end if
                end do
                
                do j = 1, r_tmp%n
                  r_logical(r_tmp%r(j)) = .true.
                end do

                if ( .not. rgn_push(r_pvt, r_tmp) ) then
                   call die('ts_pivot: programming error -- 2')
                end if
                cycle

             end if
             
          end if
          
          ! If no additional orbitals are found, exit
          if ( r_tmp%n == 0 ) exit

          ! We want to add the 'fan' atoms
          
          if ( pvt_orb .and. .not. orb_1 ) then
            call rgn_Orb2Atom(r_tmp,na_u,lasto,r_Els)
          else
            call rgn_assoc(r_Els, r_tmp)
          end if

          if ( pvt_option == PVT_FAN ) then

             ! Calculate max angle between plane and vector 'C'-xa
             select case ( fan_option ) 
             case ( PVT_FAN_MIN )
                ! minimum
                work(1) = huge(1._dp)
                work(2) = -huge(1._dp)
             case ( PVT_FAN_MEAN )
                ! mean
                work(1:2) = 0._dp
             case ( PVT_FAN_MAX )
                ! maximum
                work(1) = -huge(1._dp)
                work(2) = huge(1._dp)
             end select
             do i = 1 , r_Els%n

                ! Get C->xa vector, and normalize |a.b| / (|a|*|b|)
                tmp3 = xa(:,r_Els%r(i)) - llsB
                tmp3 = tmp3 / VNORM(tmp3)

                ! There exists two planes, and two angles

                ! 1. the planar angle only between the first plane and the vector
                work(3) = sum( tmp3 * p_cc(:,1) )
                !work(3) = sum( tmp3 * p_n(:,1) )


                ! 2. the planar angle only between the second plane and the vector
                work(4) = sum( tmp3 * p_cc(:,2) )
                !work(4) = sum( tmp3 * p_n(:,2) )


                ! Convert to radians
                ! We take the absolute as sin differs sign on
                ! positive and negative angles, we always count one way
                work(3:4) = acos(work(3:4))
                !work(3:4) = abs( asin(work(3:4)) )
                !print *,work(3:4)/3.1415926353_dp*180._dp

                select case( fan_option )
                case ( PVT_FAN_MIN )
                   work(1) = min(work(1),work(3))
                   work(2) = max(work(2),work(4))
                case ( PVT_FAN_MEAN )
                   work(1) = work(1) + work(3)
                   work(2) = work(2) + work(4)
                case ( PVT_FAN_MAX )
                   work(1) = max(work(1),work(3))
                   work(2) = min(work(2),work(4))
                end select

             end do
             if ( fan_option == PVT_FAN_MEAN ) then
               work(1:2) = work(1:2) / r_Els%n
             end if

             ! Now we need to find all atoms below those angles
             if ( .not. (pvt_orb .and. .not. orb_1) ) call rgn_nullify(r_Els)

             ! Create the list of atoms that we should search as 
             ! possible candidates for adding based on the 'fan'
             call rgn_append(r_pvt, r_tmp, r_tmp2)
             call rgn_complement(r_tmp2, c_pvt, r_Els)
             if ( pvt_orb .and. .not. orb_1 ) &
                  call rgn_Orb2Atom(r_Els,na_u,lasto,r_Els)

             ! Pre-allocate maximum size of the added atoms
             call rgn_init(r_tmp2, r_Els%n)
             r_tmp2%n = 0 ! allows doing rgn_push

             ! 'r_Els' is now the list of atoms that hasn't been processed yet
             do i = 1 , r_Els%n

                ! Get C->xa vector
                tmp3 = xa(:,r_Els%r(i)) - llsB
                tmp3 = tmp3 / VNORM(tmp3)

                ! 1. the planar angle only between the first plane and the vector
                work(3) = sum( tmp3 * p_cc(:,1) )
                work(3) = acos(work(3))
                !work(3) = sum( tmp3 * p_n(:,1) )
                !work(3) = abs( asin(work(3)) )
                if ( work(3) <= work(1) ) then
                   if ( .not. rgn_push(r_tmp2,r_Els%r(i)) ) then
                      call die('Error in programming')
                   end if
                   ! Immediately search the next one (to not add it twice)
                   cycle
                end if

                ! 2. the planar angle only between the second plane and the vector
                work(3) = sum( tmp3 * p_cc(:,2) )
                work(3) = acos(work(3))
                !work(3) = sum( tmp3 * p_n(:,2) )
                !work(3) = abs( asin(work(3)) )
                if ( work(3) >= work(2) ) then
                   if ( .not. rgn_push(r_tmp2,r_Els%r(i)) ) then
                      call die('Error in programming')
                   end if
                end if
                
             end do

!!$             ! Debugging
!!$             print *,'Currently applying: ',work(1:2) / 3.1415926353_dp * 180._dp,r_tmp2%n

          else if ( pvt_option == PVT_FRONT ) then

             ! Calculate the plane that we wish to fan about
             select case ( fan_option ) 
             case ( PVT_FAN_MIN )
               ! minimum
               work(1) = huge(1._dp)
             case ( PVT_FAN_MEAN )
               ! mean
               work(1) = 0._dp
             case ( PVT_FAN_MAX )
               ! maximum
               work(1) = -huge(1._dp)
             end select
              
             do i = 1 , r_Els%n

                ! Get elec-center -> xa vector,
                tmp3 = xa(:,r_Els%r(i)) - p_center(:,1)
                ! Scalar projection onto vector
                work(2) = VEC_PROJ_SCA(p_n(:,1), tmp3)
                
                select case ( fan_option ) 
                case ( PVT_FAN_MIN )
                   ! minimum
                   work(1) = min(work(1),work(2))
                case ( PVT_FAN_MEAN )
                   ! mean
                   work(1) = work(1) + work(2)
                case ( PVT_FAN_MAX )
                   ! maximum
                   work(1) = max(work(1),work(2))
                end select
                
             end do
             if ( fan_option == PVT_FAN_MEAN ) work(1) = work(1) / r_Els%n

             ! Now we need to find all atoms below the coordinate
             if ( .not. (pvt_orb .and. .not. orb_1) ) call rgn_nullify(r_Els)

             ! Create the list of atoms that we should search as 
             ! possible candidates for adding based on the 'fan'
             call rgn_append(r_pvt, r_tmp, r_tmp2)
             call rgn_complement(r_tmp2,c_pvt,r_Els)
             if ( pvt_orb .and. .not. orb_1 ) &
                  call rgn_Orb2Atom(r_Els,na_u,lasto,r_Els)

             ! Pre-allocate maximum size of the added atoms
             call rgn_init(r_tmp2, r_Els%n)
             r_tmp2%n = 0

             ! 'r_Els' is now the list of atoms that hasn't been processed yet
             do i = 1 , r_Els%n

                ! Get elec-center -> xa vector,
                tmp3 = xa(:,r_Els%r(i)) - p_center(:,1)
                ! Project onto vector
                work(2) = VEC_PROJ_SCA(p_n(:,1), tmp3)

                ! Check whether the atom should be added
                if ( work(1) >= work(2) ) then
                   if ( .not. rgn_push(r_tmp2,r_Els%r(i)) ) then
                      call die('Error in programming')
                   end if
                end if
                
             end do

          end if

          ! Add the new atoms
          if ( r_tmp2%n > 0 ) then
             ! Extend the added region to have the new fanned region
             if ( pvt_orb .and. .not. orb_1 ) then
                call rgn_atom2orb(r_tmp2,na_u,lasto,r_Els)
                call rgn_append(r_tmp,r_Els,r_tmp)
             else
                call rgn_append(r_tmp,r_tmp2,r_tmp)
             end if
          end if
          
          ! Clean-up
          call rgn_delete(r_Els,r_tmp2)

          do j = 1, r_tmp%n
            r_logical(r_tmp%r(j)) = .true.
          end do
          
          ! Append the newly found region that is connecting out to the
          ! full region
          if ( .not. rgn_push(r_pvt, r_tmp) ) then
             call die('ts_pivot: programming error -- 3')
          end if
          
          ! we sort the newly attached region
          call rgn_sp_sort(r_pvt, tmp_Sp, r_tmp, R_SORT_MAX_BACK, r_logical=r_logical)

       end do

       deallocate(r_logical)

    end if

    ! Clean-up
    call delete(tmp_Sp)
    
    if ( r_pvt%n /= c_pvt%n .and. lextend ) then
       if ( IONode ) then
          call rgn_print(r_pvt)
          call rgn_print(c_pvt)
          call rgn_print(priority)
       end if
       call die('ts_pivot: Error in size estimation, the sparse pattern &
            &removal is erroneous')
    end if

    if ( .not. pvt_orb ) then

       ! re-create the orbital pivoting
       call rgn_copy(r_pvt,r_tmp)
       call rgn_atom2orb(r_tmp,na_u,lasto,r_pvt)

    end if

    ! Reduce the pivoting array to its size
    if ( r_pvt%n /= size(r_pvt%r) ) then
       call rgn_copy(r_pvt, r_tmp)
       call rgn_delete(r_pvt)
       call rgn_assoc(r_pvt, r_tmp)
       call rgn_nullify(r_tmp)
    end if

    ! Clean-up
    call rgn_delete(r_tmp,r_tmp2,c_pvt)
    call rgn_delete(r_Els,priority)

    pvt_str = str_tmp

    call timer('pivot', 2)

  contains

    function str_contain(str,name, delete) result(contain)
      use m_char, only : lcase
      character(len=*), intent(inout) :: str
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: delete
      logical :: contain

      character(len=len(str)) :: lstr
      character(len=len(name)) :: lname
      
      integer :: i
      
      contain = .false.
      
      lstr = lcase(str)
      lname = lcase(name)
      
      ! Get the index of the str
      i = index(lstr,trim(lname))
      contain = i > 0
      if ( present(delete) ) then
         if ( .not. delete ) return
      end if
      if ( contain ) then
         str(i:i+len_trim(name)-1) = ' '
      end if
      
      ! Remove any ' +', '+ '
      i = index(str,' +')
      if ( i > 0 ) str(i:i+1) = ' '
      i = index(str,'+ ')
      if ( i > 0 ) str(i:i+1) = ' '
      
    end function str_contain
    
  end subroutine ts_pivot

  subroutine crt_El_priority(N_Elec,Elecs,pr,na_u,lasto,is_orb)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(tRgn), intent(inout) :: pr ! Needs to be pre-allocated
    integer, intent(in) :: na_u, lasto(0:na_u)
    logical, intent(in) :: is_orb
    integer :: i, iEl, j
    type(tRgn) :: rgn

    ! Initialize
    pr%r = 0
    do iEl = 1 , N_Elec
       if ( is_orb ) then
          call rgn_copy(Elecs(iEl)%o_inD,rgn)
       else
          call rgn_Orb2Atom(Elecs(iEl)%o_inD,na_u,lasto,rgn)
       end if

       i = rgn%n
       do j = 1 , rgn%n
          pr%r(rgn%r(j)) = i
          i = i - 1
       end do
    end do
    call rgn_delete(rgn)

  end subroutine crt_El_priority

  function consecutive_Elec_orb(El,r) result(con)
    type(Elec), intent(in) :: El
    type(tRgn), intent(in) :: r
    type(tRgn) :: o_inD, r_tmp
    integer :: con
    integer :: i_Elec, io, idx_Elec, idx, no

    idx = El%idx_o - 1
    no = TotUsedOrbs(El)

    ! Create the pivoting table for the electrodes
    call rgn_init(r_tmp,no)
    do io = 1 , no
       ! equals the io'th orbital index in the 
       !           TS region.  io == ts_i
       r_tmp%r(io) = rgn_pivot(r,idx+io)
    end do
       
    ! Sort it to be able to gather the indices
    ! in the correct order
    call rgn_sort(r_tmp)
    ! pivot the o_inD back
    call rgn_init(o_inD,no)
    do io = 1 , no 
       o_inD%r(io) = r%r(r_tmp%r(io))
    end do

    con    = 0
    i_Elec = 1
    do while ( i_Elec <= o_inD%n ) 

       idx_Elec = rgn_pivot(r,o_inD%r(i_Elec))
       io = 1
       do while ( i_Elec + io <= o_inD%n ) 
          idx = rgn_pivot(r,o_inD%r(i_Elec+io))
          ! In case it is not consecutive
          if ( idx - idx_Elec /= io ) exit
          io = io + 1
       end do
       con = con + 1
       
       i_Elec = i_Elec + io

    end do

    call rgn_delete(r_tmp,o_inD)

  end function consecutive_Elec_orb
  
end module m_ts_pivot
