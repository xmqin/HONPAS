! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior, 2014.
module m_tbt_kpoint
  
  use precision, only : dp

  implicit none

  private
  save

  logical, public :: Gamma
  integer, public :: nkpnt ! Total number of k-points

  real(dp), pointer, public :: kweight(:) 
  real(dp), pointer, public :: kpoint(:,:)

  integer,  public :: kscell(3,3) = 0
  real(dp), public :: kdispl(3) = 0.0_dp

  public :: setup_kpoint_grid
  public :: write_k_points
  public :: read_kgrid
  public :: tbt_iokp

  ! The local k-point method
  integer, public :: K_METHOD = 1
  integer, parameter, public :: K_METHOD_MONKHORST_PACK = 1
  integer, parameter, public :: K_METHOD_SIMP_MIX = 2
  integer, parameter, public :: K_METHOD_BOOLE_MIX = 3
  integer, parameter, public :: K_METHOD_GAUSS_LEGENDRE = 4
  integer, parameter, public :: K_METHOD_TANH_SINH = 5
  integer, parameter, public :: K_METHOD_PATH = 6
  integer, parameter, public :: K_METHOD_LIST = 7

contains

  subroutine read_kgrid(bName, TRS, cell, kpt, wkpt, &
       is_b, &
       kcell, kdispl)

    use parallel, only : IONode
    use fdf
    use fdf_extra, only : fdf_bnext
    use intrinsic_missing, only : EYE, VNORM, SPC_PROJ, VEC_PROJ
    use units, only : Pi
    use m_find_kgrid, only : find_kgrid, trim_kpoint_list

    use m_ts_tdir, only: ts_tidx

    use m_integrate
    use m_gauss_quad

    ! INPUT
    character(len=*), intent(in) :: bName
    ! Whether time-reversal symmetry applies
    logical, intent(in) :: TRS
    real(dp), intent(in) :: cell(3,3)
    ! OUTPUT
    real(dp), pointer :: kpt(:,:), wkpt(:)
    logical, intent(in), optional :: is_b
    ! Optional OUTPUT
    integer, intent(out), optional :: kcell(3,3)
    real(dp), intent(out), optional :: kdispl(3)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer :: i, ik, j, k, nkpt
    real(dp) :: rcell(3,3), displ(3), ksize(3), rtmp, p(3), q(3)
    real(dp) :: prev_k(3), next_k(3), k_path_length
    real(dp) :: contrib
    integer :: kscell(3,3), inkpt(3)
    real(dp), allocatable :: k3_1(:,:), k3_2(:,:), k3_3(:,:)
    real(dp), allocatable :: tmp3(:,:)

    logical :: is_block
    character(len=50) :: ctmp
    logical :: even_path

    ! Initialize values
    ksize(:) = 1._dp
    displ(:) = 0._dp
    kscell(:,:) = 0
    kscell(1,1) = 1
    kscell(2,2) = 1
    kscell(3,3) = 1
    K_METHOD = -1
    is_block = .true.

    ! If the block does not exist, simply 
    ! create the Gamma-point
    nullify(kpt,wkpt)

    if ( fdf_islist(bName) ) then
       
       ! This is an easy read...
       K_METHOD = K_METHOD_MONKHORST_PACK
       ! Try and read the list information
       call fdf_list(bName, 3, kscell(:,1))
       kscell(2,2) = kscell(2,1)
       kscell(2,1) = 0
       kscell(3,3) = kscell(3,1)
       kscell(3,1) = 0

       is_block = .false.
       
    else if ( .not. fdf_block(bName,bfdf) ) then
       K_METHOD = K_METHOD_MONKHORST_PACK

       ! the block does not exist, hence the user
       ! requests a Gamma-point.
       allocate(kpt(3,1),wkpt(1))
       kpt(:,:) = 0._dp
       wkpt(:) = 1._dp
       if ( present(kcell) ) kcell = kscell
       if ( present(kdispl) ) kdispl = displ
       
       return
       
    end if

    even_path = .false.
    k_path_length = 0._dp
    nkpt = 0
    i = 0
    ! Pre-read paths if requested (we need to pre-allocate
    ! number of k-points along path)

    ! We need this to calculate the correct length in the 
    ! Brillouin zone
    call reclat(cell,rcell,1)

    if ( is_block ) then
    do while ( fdf_bnext(bfdf,pline) )
       if ( fdf_bnnames(pline) > 0 ) then
          ctmp = fdf_bnames(pline,1)

          if ( leqi(ctmp,'path') ) then
             if ( K_METHOD /= -1 ) call die('Do not mix different &
                  &k-grid methods. Only specify one method.')
             K_METHOD = K_METHOD_PATH
             i = i + 1
             call read_path(bfdf,pline,rcell,.false.,0._dp,prev_k,next_k,j)
             call kpoint_convert(rcell,prev_k,p,-2)
             call kpoint_convert(rcell,next_k,q,-2)
             p = q - p
             k_path_length = k_path_length + VNORM(p)
             nkpt = nkpt + j
             
          else if ( leqi(ctmp,'path-even') .or. &
               leqi(ctmp,'even-path') ) then
             ! Check if we should make as even a spacing as possible
             even_path = .true.
             
          else if ( leqi(ctmp,'list') ) then
             if ( K_METHOD /= -1 ) call die('Do not mix different &
                  &k-grid methods. Only specify one method.')

             K_METHOD = K_METHOD_LIST

          else if ( leqi(ctmp,'method') ) then
             
             if ( K_METHOD /= -1 ) call die('Do not mix different &
                  &k-grid methods. Only specify one method.')
             
             ! read the method
             ctmp = fdf_bnames(pline,2)
             if ( leqi(ctmp,'monkhorst-pack') .or. &
                  leqi(ctmp,'MP') ) then
                K_METHOD = K_METHOD_MONKHORST_PACK
             else if ( leqi(ctmp,'gauss-legendre') .or. &
                  leqi(ctmp,'g-legendre') ) then
                K_METHOD = K_METHOD_GAUSS_LEGENDRE
             else if ( leqi(ctmp,'tanh-sinh') ) then
                K_METHOD = K_METHOD_TANH_SINH
             else if ( leqi(ctmp,'simpson-mix') .or. &
                  leqi(ctmp,'simp-mix') ) then
                K_METHOD = K_METHOD_SIMP_MIX
             else if ( leqi(ctmp,'boole-mix') ) then
                K_METHOD = K_METHOD_BOOLE_MIX
             else
                call die('Could not recognize the k-grid method, &
                     &please assert you have written something &
                     &available. Check the manual.')
             end if

          end if
          
       end if
    end do
    ! Reset the method to a pristine Monkhorst-Pack grid
    ! if none specified
    if ( K_METHOD == -1 ) K_METHOD = K_METHOD_MONKHORST_PACK
    if ( K_METHOD == K_METHOD_PATH ) then
       j = nkpt
       if ( even_path ) then
          ! For even paths we cannot assert full range, add a few 
          ! numbers
          j = nkpt + i * 2
       end if
       allocate(kpt(3,j))
    end if

    ! Rewind the block to read again
    call fdf_brewind(bfdf)

    ! Read in the blocks
    ik = 0
    do while ( fdf_bnext(bfdf,pline) )
       
       if ( fdf_bnnames(pline) > 0 ) then
          
          ctmp = fdf_bnames(pline,1)
          
          ! We have some kind of designation
          if ( K_METHOD == K_METHOD_PATH ) then
             
             if ( leqi(ctmp,'path') ) then

                ! Read in options (for even path, j is overwritten)
                j = nkpt
                call read_path(bfdf,pline,rcell,even_path, &
                     k_path_length,prev_k,next_k,j)

                ! Get dk (the next path curve has the end-point, if prev is used)
                p = (next_k - prev_k) / real(j,dp)
                do i = 1 , j
                   ik = ik + 1
                   kpt(:,ik) = prev_k + p * (i-1)
                   call check_zero(kpt(:,ik),ts_tidx)
                end do
                
             end if

          else if ( K_METHOD == K_METHOD_LIST ) then

             if ( leqi(ctmp,'list') ) then
                
                ! Get number of k-points
                nkpt = fdf_bintegers(pline,1)
                ! allocate for the k-points
                allocate(kpt(3,nkpt),wkpt(nkpt))
                ! reset all weights to be equal
                wkpt(1) = -1._dp

                do ik = 1 , nkpt
                   if ( .not. fdf_bnext(bfdf,pline) ) &
                        call die('Could not read correct number of k-points in list')
                   
                   kpt(1,ik) = fdf_bvalues(pline,1)
                   kpt(2,ik) = fdf_bvalues(pline,2)
                   kpt(3,ik) = fdf_bvalues(pline,3)
                   call check_zero(kpt(:,ik),ts_tidx)

                   if ( fdf_bnvalues(pline) > 3 ) then
                      ! fine, read in the weight
                      wkpt(ik) = fdf_bvalues(pline,4)
                   else if ( ik > 1 .and. wkpt(1) > 0._dp ) then
                      call die('Could not read weight for k-point, either &
                           &supply all, or none (which means equal weight)')
                   end if

                end do

                if ( wkpt(1) < 0._dp ) then
                   wkpt(:) = 1._dp / real(nkpt,dp)
                end if

             end if
             
          else if ( leqi(ctmp(1:5),'diag-') .or. &
               leqi(ctmp(1:9),'diagonal-') ) then

             if ( fdf_bnintegers(pline) /= 1 ) then
                call die('Please correct your input, you have not supplied a number in diag-.')
             end if
             
             if ( ctmp(5:5) == '-' ) ctmp = ctmp(6:)
             if ( ctmp(9:9) == '-' ) ctmp = ctmp(10:)

             if ( leqi(ctmp,'A1') .or. leqi(ctmp,'a') ) then
                ik = 1
             else if ( leqi(ctmp,'A2') .or. leqi(ctmp,'b') ) then
                ik = 2
             else if ( leqi(ctmp,'A3') .or. leqi(ctmp,'c') ) then
                ik = 3
             else
                call die('Tbt: could not figure out the diagonal direction specified')
             end if

             ! Set the diagonal
             kscell(ik,ik) = max(1,fdf_bnintegers(pline,1))
             
          else if ( leqi(ctmp,'diagonal') .or. &
               leqi(ctmp,'diag') ) then
             
             ik = fdf_bnintegers(pline)
             
             if ( ik < 3 .and. IONode ) then
                write(*,'(/,a)') 'tbt: POSSIBLE WARNING'
                write(*,'(a,i0,a)') 'tbt: You have only supplied ', &
                     ik,' of the 3 diagonal k-cell elements.'
                write(*,'(a)') 'tbt: Will assume this order A1-A2-A3'
             end if
             
             ! Set the diagonal
             do i = 1 , ik
                kscell(i,i) = max(1,fdf_bintegers(pline,i))
             end do

          else if ( leqi(ctmp,'displacement') .or. &
               leqi(ctmp,'displ') ) then

             displ(1) = mod(fdf_bvalues(pline,1), 1._dp)
             displ(2) = mod(fdf_bvalues(pline,2), 1._dp)
             displ(3) = mod(fdf_bvalues(pline,3), 1._dp)
             
          else if ( leqi(ctmp,'size') ) then

             ksize(1) = fdf_bvalues(pline,1)
             ksize(2) = fdf_bvalues(pline,2)
             ksize(3) = fdf_bvalues(pline,3)

             if ( any(ksize > 1._dp) .or. any(ksize <= 0._dp) ) then
                call die('The size of the Brillouin zone MUST be &
                     &less than or equal to 1.')
             end if

          end if

       else if ( fdf_bnvalues(pline) > 3 .and. &
            K_METHOD /= K_METHOD_PATH .and. &
            K_METHOD /= K_METHOD_LIST ) then

          do ik = 1 , 3
             if ( fdf_bnintegers(pline) < 3 ) then
                call die('Could not read three integers from the &
                     &tbtrans Monkhorst-Pack grid.')
             end if
             kscell(1,ik) = fdf_bintegers(pline,1)
             kscell(2,ik) = fdf_bintegers(pline,2)
             kscell(3,ik) = fdf_bintegers(pline,3)
             if ( fdf_bnvalues(pline) > 3 ) then
                displ(ik) = mod(fdf_bvalues(pline,4), 1._dp)
             end if

             ! To not error out of only 3 lines grids
             if ( ik == 3 ) cycle
             if ( .not. fdf_bnext(bfdf,pline) ) &
                  call die('Could not read kgrid from block: '//trim(bName))
          end do

       end if

    end do
    end if ! whether the input is a block

    ! We do not allow different methods if we do not have a diagonal
    ! size matrix
    select case ( K_METHOD )
    case ( K_METHOD_MONKHORST_PACK )
       ! do nothing
    case default
       if ( kscell(2,1) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
       if ( kscell(3,1) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
       if ( kscell(1,2) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
       if ( kscell(3,2) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
       if ( kscell(1,3) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
       if ( kscell(2,3) /= 0 ) K_METHOD = K_METHOD_MONKHORST_PACK
    end select

    if ( K_METHOD == K_METHOD_PATH ) then

       if ( IONode ) then
          write(*,'(a)')'tbt: k-points are following paths in the Brillouin zone.'
          write(*,'(a)')'WARNING: The averaged transmission will not necessarily &
               &reflect the total transmission!'
       end if

       ! Correct number of k-points
       ! ik is our loop counter in the above loop
       nkpt = ik

       if ( nkpt /= size(kpt,dim=2) ) then

          allocate(k3_1(3,nkpt))
          k3_1(:,:) = kpt(:,1:nkpt)
          deallocate(kpt)
          nullify(kpt)
          allocate(kpt(3,nkpt))
          kpt = k3_1
          deallocate(k3_1)

       end if
       
       ! Allocate weights
       allocate(wkpt(nkpt))
       wkpt = 1._dp / nkpt

    else if ( K_METHOD /= K_METHOD_LIST ) then

       if ( ts_tidx > 0 ) then
          i = ts_tidx
          ! the contribution along this vector is too much
          ! to disregard the elongation along this
          ! direction.
          ! We *MUST* kill all k-points in this direction
          kscell(:,i) = 0
          kscell(i,:) = 0
          kscell(i,i) = 1
          displ(i) = 0._dp
          ksize(i) = 1._dp
       end if

       if ( present(kcell) ) kcell = kscell
       if ( present(kdispl) ) kdispl = displ

       select case ( K_METHOD )
       case ( K_METHOD_MONKHORST_PACK )

          if ( IONode ) write(*,*) ! new-line
          
          call EYE(3, rcell, 2._dp * Pi)
          call find_kgrid(rcell, kscell, displ, .true., &
               TRS , &
               nkpt, kpt, wkpt, rtmp)

       case default

          allocate( k3_1(kscell(1,1),2) )
          allocate( k3_2(kscell(2,2),2) )
          allocate( k3_3(kscell(3,3),2) )

          call method_MP(K_METHOD,kscell(1,1),k3_1(:,1),k3_1(:,2))
          call method_MP(K_METHOD,kscell(2,2),k3_2(:,1),k3_2(:,2))
          call method_MP(K_METHOD,kscell(3,3),k3_3(:,1),k3_3(:,2))

          if ( TRS ) then
             
             ! Cut in half the largest one
             i = 3
             if ( kscell(2,2) > kscell(i,i) ) i = 2
             if ( kscell(1,1) > kscell(i,i) ) i = 1
             j = kscell(i,i)
             allocate( tmp3(j/2 + mod(j,2), 2) )
             j = j/2 + 1

             select case ( i )
             case ( 1 )
                tmp3(:,:) = k3_1(j:,:)
                deallocate(k3_1)
             case ( 2 )
                tmp3(:,:) = k3_2(j:,:)
                deallocate(k3_2)
             case ( 3 )
                tmp3(:,:) = k3_3(j:,:)
                deallocate(k3_3)
             end select
             
             if ( mod(j,2) == 1 ) then
                tmp3(2:,2) = tmp3(2:,2) * 2._dp
             else
                tmp3(:,2)  = tmp3(:,2) * 2._dp
             end if
             
             select case ( i )
             case ( 1 )
                allocate( k3_1(size(tmp3,1),2) )
                k3_1 = tmp3
             case ( 2 )
                allocate( k3_2(size(tmp3,1),2) )
                k3_2 = tmp3
             case ( 3 )
                allocate( k3_3(size(tmp3,1),2) )
                k3_3 = tmp3
             end select
             deallocate(tmp3)

          end if
          
          ! Create the k-points
          inkpt(1) = size(k3_1,1)
          inkpt(2) = size(k3_2,1)
          inkpt(3) = size(k3_3,1)
          nkpt = product(inkpt)
          allocate( kpt(3,nkpt) )
          allocate( wkpt(nkpt) )

          ik = 0
          do k = 1 , inkpt(3)
          do j = 1 , inkpt(2)
          do i = 1 , inkpt(1)
             ik = ik + 1
             kpt(1,ik) = k3_1(i,1)
             kpt(2,ik) = k3_2(j,1)
             kpt(3,ik) = k3_3(k,1)
             call check_zero(kpt(:,ik),ts_tidx)
             wkpt(ik)  = k3_1(i,2) * k3_2(j,2) * k3_3(k,2)
          end do
          end do
          end do

          ! Reduce number of k-points via explicit
          ! comparison
          call trim_kpoint_list(nkpt,kpt,wkpt)

          deallocate(k3_1)
          deallocate(k3_2)
          deallocate(k3_3)

       end select

    else

       ! Print out warning if the sum of the weights
       ! does not equal 1
       contrib = 0._dp
       do ik = 1 , nkpt
          contrib = contrib + wkpt(ik)
       end do

       if ( IONode .and. abs(contrib - 1._dp) > 1.e-7_dp ) then
          write(*,'(a)')'WARNING: Weights for k-points in &
               & %block '//trim(bName)//' does not sum to 1.'
       end if

    end if

    if ( K_METHOD /= K_METHOD_LIST .and. &
         K_METHOD /= K_METHOD_PATH ) then

       ! Re-scale the k-points to the correct size
       do ik = 1 , nkpt
          kpt(:,ik) = kpt(:,ik) * ksize(:)
          call check_zero(kpt(:,ik),ts_tidx)
       end do

       ! Rescale the weights if the size of the
       ! kgrid zone has been narrowed
       ! This is simple, the basic size is 1.
       ! hence if the size of the k-cell is 1/2
       ! then the weights must also be 1/2
       ! I.e. the size is the weight scale.
       do ik = 1 , 3
          ! we do not need to rescale if the number
          ! of k-points in the i'th direction is
          !  !! regardless of what the user says ;) !!
          if ( kscell(ik,ik) == 1 ) cycle
          wkpt(:) = wkpt(:) * ksize(ik)
       end do

    end if

    if ( present(is_b) ) then
       if ( is_b ) return
    end if

    ! Transform the reciprocal units into length
    call reclat(cell,rcell,1)
    do ik = 1 , nkpt
       p(:) = kpt(:,ik)
       call kpoint_convert(rcell,p(:),kpt(:,ik),-2)
    end do

  contains

    subroutine check_zero(k,tidx)
      real(dp), intent(in) :: k(3)
      integer, intent(in) :: tidx
      if ( tidx > 0 ) then
         if ( abs(k(tidx)) > 1.e-5_dp ) then
            write(*,'(a)') 'An univocal transport direction is existing'
            write(*,'(a)') 'You should not use k-points in that direction'
            write(*,'(a)') 'as it is redundant work.'
            call die('Check the output, you are using redundant k-points.')
         end if
      end if
      
    end subroutine check_zero

    subroutine read_path(bfdf,pline,rcell,even,kl,prev_k,next_k,nkpt)
      type(block_fdf), intent(inout) :: bfdf
      type(parsed_line), pointer :: pline
      real(dp), intent(in) :: rcell(3,3)
      logical, intent(in) :: even
      real(dp), intent(in) :: kl
      real(dp), intent(inout) :: prev_k(3), next_k(3)
      integer, intent(inout) :: nkpt

      character(len=50) :: ctmp
      integer :: i
      real(dp) :: p(3), q(3)

      ! Read number of k-points in this path
      i = fdf_bintegers(pline,1)

      ! Read from
      if ( .not. fdf_bnext(bfdf,pline) ) &
           call die('Could not read from in k-path')
      if ( fdf_bnnames(pline) > 1 ) then
         ctmp = fdf_bnames(pline,2)
         if ( leqi(ctmp,'prev') .or. leqi(ctmp,'previous') ) then
            ! perfect
            prev_k = next_k
         else
            call die('Could not find previous/prev in the from line')
         end if
      else
         prev_k(1) = fdf_bvalues(pline,1)
         prev_k(2) = fdf_bvalues(pline,2)
         prev_k(3) = fdf_bvalues(pline,3)
      end if

      if ( .not. fdf_bnext(bfdf,pline) ) &
           call die('Could not read to in k-path')
      next_k(1) = fdf_bvalues(pline,1)
      next_k(2) = fdf_bvalues(pline,2)
      next_k(3) = fdf_bvalues(pline,3)

      if ( even ) then
         call kpoint_convert(rcell,prev_k,p,-2)
         call kpoint_convert(rcell,next_k,q,-2)
         p = q - p
         ! Calculate the fraction of k-points in this
         ! segment (nkpt is then the full number of k-points)
         nkpt = nint(vnorm(p) / kl * nkpt)
      else
         nkpt = i
      end if
      
    end subroutine read_path

    subroutine method_MP(method,n,kpt,wkpt)
      integer, intent(in) :: method, n
      real(dp), intent(out) :: kpt(n), wkpt(n)
      
      if ( n == 1 ) then
         kpt(1) = 0._dp
         wkpt(1) = 1._dp
         return
      end if
      
      select case ( K_METHOD )
      case ( K_METHOD_GAUSS_LEGENDRE )
         
         call Gauss_Legendre_Rec(n,0,-0.5_dp,0.5_dp, kpt, wkpt)

       case ( K_METHOD_TANH_SINH )

          call TanhSinh_Exact(n, kpt, wkpt, -0.5_dp, 0.5_dp )

       case ( K_METHOD_SIMP_MIX )

          call Simpson_38_3_rule(n, kpt, wkpt, -0.5_dp, 0.5_dp )

       case ( K_METHOD_BOOLE_MIX )

          call Booles_Simpson_38_3_rule(n, kpt, wkpt, -0.5_dp, 0.5_dp )

       end select
       
    end subroutine method_MP
    
  end subroutine read_kgrid

  subroutine read_kdiag(fName, TRS, cell, kpt, wkpt, &
       is_b, &
       kcell, kdispl)

    use parallel, only : IONode
    use fdf
    use intrinsic_missing, only : EYE
    use units, only : Pi
    use m_find_kgrid, only : find_kgrid

    use m_ts_tdir, only: ts_tidx

    ! INPUT
    character(len=*), intent(in) :: fName
    ! Whether time-reversal symmetry applies
    logical, intent(in) :: TRS
    real(dp), intent(in) :: cell(3,3)
    ! OUTPUT
    real(dp), pointer :: kpt(:,:), wkpt(:)
    logical, intent(in), optional :: is_b
    ! Optional OUTPUT
    integer, intent(out), optional :: kcell(3,3)
    real(dp), intent(out), optional :: kdispl(3)

    ! Local variables
    ! Read diagonal elements
    integer :: kdiag(3)
    integer :: i, nkpt
    integer :: kscell(3,3)
    real(dp) :: rcell(3,3), p(3), rtmp
    real(dp) :: displ(3)

    ! Initialize constants
    ! This is the used method
    K_METHOD = K_METHOD_MONKHORST_PACK
    kscell = 0
    displ = 0._dp


    ! Initialize to 1 along all diagonal elements
    kdiag = 1
    ! Read in the diagonal kpoints
    call fdf_list(fName, 3, kdiag)

    ! Ensure a positive (>0) number
    do i = 1 , 3
       kscell(i,i) = max( kdiag(i) , 1 )
    end do
    ! Forcefully set the cell to 1 along the transport direction
    if ( ts_tidx > 0 ) then
       i = ts_tidx
       kscell(:,i) = 0
       kscell(i,:) = 0
       kscell(i,i) = 1
    end if

    ! Copy over read options (if requested)
    if ( present(kcell) ) kcell = kscell
    if ( present(kdispl) ) kdispl = displ

    ! Create the k-points
    if ( IONode ) write(*,*) ! new-line
    
    call EYE(3, rcell, 2._dp * Pi)
    call find_kgrid(rcell, kscell, displ, .true., &
         TRS , &
         nkpt, kpt, wkpt, rtmp)

    if ( present(is_b) ) then
       if ( is_b ) return
    end if

    ! Transform the reciprocal units into length
    call reclat(cell, rcell, 1)
    do i = 1 , nkpt
       p(:) = kpt(:,i)
       call kpoint_convert(rcell, p(:), kpt(:,i), -2)
    end do

  end subroutine read_kdiag
  
  subroutine setup_kpoint_grid( cell )
    
    use fdf, only : fdf_get, fdf_block, fdf_islist
    use fdf, only : leqi, block_fdf
    use parallel, only  : IONode

    ! Local Variables
    real(dp), intent(in)   :: cell(3,3)

    real(dp) :: bkpt(3)
    integer :: i
    type(block_fdf) :: bfdf
    ! Whether we should apply time-reversal symmetry
    logical :: TRS
    character(len=250) :: user_kfile

    nullify(kweight,kpoint)

    ! In case the user requests to utilize a provided
    ! k-point file
    user_kfile = fdf_get('TBT.k.File','NONE')
    if ( .not. leqi(user_kfile,'NONE') ) then

       if ( IONode ) then
          write(*,'(a)') 'tbt: Reading user specified k-points.'
          write(*,'(2a)')'tbt: k-points found in file: ',trim(user_kfile)
          write(*,'(2a)')'tbt: *** Responsibility is on your side! ***'
       end if

       K_METHOD = K_METHOD_LIST

       call tbt_iokp_read(user_kfile,nkpnt,kpoint,kweight)

       ! Convert the k-points to local units
       do i = 1 , nkpnt
          bkpt(:) = kpoint(:,i)
          call kpoint_convert(cell,bkpt(:),kpoint(:,i),-1)
       end do

    else    

       TRS = .not. fdf_get('SpinSpiral',.false.)
       TRS = fdf_get('TBT.Symmetry.TimeReversal',TRS)

       if ( fdf_block('TBT.k',bfdf) .or. fdf_islist('TBT.k') ) then
          call read_kgrid('TBT.k', &
               TRS, cell, kpoint, kweight, &
               kcell=kscell, kdispl=kdispl)
       else if ( fdf_block('TBT.kgrid.MonkhorstPack',bfdf) ) then
          call read_kgrid('TBT.kgrid.MonkhorstPack', &
               TRS, cell, kpoint, kweight, &
               kcell=kscell, kdispl=kdispl)
       else
          call read_kgrid('kgrid.MonkhorstPack', &
               TRS, cell, kpoint, kweight, &
               kcell=kscell, kdispl=kdispl)
       end if
       nkpnt = size(kweight)
       
    end if

    ! the user could do multiple Gamma-points (for no apparent reason)
    Gamma = .true.
    do i = 1, nkpnt
       Gamma = Gamma .and. &
            dot_product(kpoint(:,1),kpoint(:,1)) < 1.0e-20_dp
    end do
    
    call write_k_points()

  end subroutine setup_kpoint_grid
  
  subroutine write_k_points()
    use parallel, only : IONode
    use m_verbosity, only : verbosity
    
    integer :: ik, ix, i

    if ( .not. IONode ) return

    write(*,'(a,i0)')  'tbt: Number of transport k-points = ', nkpnt

    select case ( K_METHOD )
    case ( K_METHOD_MONKHORST_PACK )
       write(*,'(a)') 'tbt: Method = Monkhorst-Pack grid.'
    case ( K_METHOD_SIMP_MIX )
       write(*,'(a)') 'tbt: Method = Simpson grid.'
    case ( K_METHOD_BOOLE_MIX )
       write(*,'(a)') 'tbt: Method = Booles grid.'
    case ( K_METHOD_GAUSS_LEGENDRE )
       write(*,'(a)') 'tbt: Method = Gauss-Legendre grid.'
    case ( K_METHOD_TANH_SINH )
       write(*,'(a)') 'tbt: Method = Tanh-Sinh grid.'
    case ( K_METHOD_PATH )
       write(*,'(a)') 'tbt: Method = User path in Brillouin zone.'
    case ( K_METHOD_LIST )
       write(*,'(a)') 'tbt: Method = User specified k-points as list.'
    case default
       write(*,'(a)') 'tbt: Unknown k-point method.'
       call die('Unsupported, programming error.')
    end select

    if ( verbosity > 5 ) then
       write(*,'(/,a)') 'tbt: k-point coordinates (Bohr**-1) and weights:'
       write(*,'(a,i4,3f12.6,3x,f12.6)') &
            ('tbt: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik), &
            ik=1,nkpnt)
       write(*,*) ! new line
    end if

    ! Always write the TranSIESTA k-points
    call tbt_iokp( nkpnt, kpoint, kweight )

    ! Only print out the k-points if the grid is used
    if ( K_METHOD == K_METHOD_PATH ) return
    if ( K_METHOD == K_METHOD_LIST ) return

    write(*,'(a)') 'tbt: k-grid: Supercell and displacements'
    do ix = 1 , 3
       write(*,'(a,3i4,3x,f8.3)') 'tbt:         ', &
            (kscell(i,ix),i=1,3), kdispl(ix)
    end do

  end subroutine write_k_points
  
  subroutine tbt_iokp( nk, points, weight , fend)
! *******************************************************************
! Saves tbtrans k-points (only writing) Bohr^-1
! Emilio Artacho, Feb. 1999
! Modified by Nick Papior Andersen to not overwrite the SIESTA/TranSIESTA k-points
! ********** INPUT **************************************************
! integer nk           : Number of k-points
! real*8  points(3,nk) : k-point coordinates
! real*8  weight(3,nk) : k-point weight
! *******************************************************************
    use files, only : slabel
    use m_tbt_save, only : save_dir

    integer, intent(in) :: nk
    real(dp), intent(in) :: points(3,nk), weight(nk)
    character(len=*), intent(in), optional :: fend
    external :: io_assign, io_close

    ! Internal 
    integer :: iu, ik, ix

    call io_assign( iu )
    if ( present(fend) ) then
       open( iu, file=trim(save_dir)//trim(slabel)//'.'//trim(fend), &
            form='formatted', status='unknown' ) 
    else
#ifdef TBT_PHONON
       open( iu, file=trim(save_dir)//trim(slabel)//'.PHT.KP', &
            form='formatted', status='unknown' ) 
#else
       open( iu, file=trim(save_dir)//trim(slabel)//'.TBT.KP', &
            form='formatted', status='unknown' ) 
#endif
    end if

    write(iu,'(i6)') nk
    write(iu,'(i6,3f12.6,3x,f12.6)') &
         (ik, (points(ix,ik),ix=1,3), weight(ik), ik=1,nk)

    call io_close( iu )
    
  end subroutine tbt_iokp

  ! The user can specify their own k-points
  subroutine tbt_iokp_read(fname,nkpt,kpt,wkpt)
    use parallel, only : Node
    use m_os, only : file_exist
#ifdef MPI
    use mpi_siesta
#endif
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: wsum

#ifdef MPI
    integer :: MPIerror
#endif

    ! The user has requested to read in 
    ! k-points from a specific file...
    ! We will do that for him/her.

    integer :: iu, ik, ix, stat

    if ( Node == 0 ) then
       
       if ( .not. file_exist(trim(fname)) ) then
          call die('Could not locate file '//trim(fname)// &
               ' please ensure that the file exists.')
       end if

       call io_assign( iu )
       open( iu, file=trim(fname), form='formatted', status='old' ) 

       ! Read number of k-points
       read(iu,*,iostat=stat) nkpt
       call kill_iokp(stat,0)

    end if

#ifdef MPI
    call MPI_Bcast(nkpt,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

    nullify(kpt,wkpt)
    allocate(kpt(3,nkpt),wkpt(nkpt))

    if ( Node == 0 ) then
       ! Read in the k-points
       wsum = 0._dp
       do ik = 1 , nkpt
          ! read current k-point
          read(iu,*,iostat=stat) ix, kpt(:,ik), wkpt(ik) ! (i6,3f12.6,3x,f12.6)
          call kill_iokp(stat,ik)
          wsum = wsum + wkpt(ik)
       end do
      
       if ( abs(wsum - 1._dp) > 1.e-7_dp ) then
          write(*,'(a)')'WARNING: Weights for user specified k-points does &
               &not sum to 1.'
       end if

       call io_close( iu )

    end if

#ifdef MPI
    call MPI_Bcast(kpt(1,1),3*nkpt,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(wkpt(1),nkpt,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#endif

  contains
    
    subroutine kill_iokp(stat,line)
      integer, intent(in) :: stat, line
      if ( stat == 0 ) return
      write(*,*) 'TBtrans iokp could not read your input file'
      write(*,*) 'The k-points MUST be in units of reciprocal vectors!'
      write(*,*) 'TBtrans will convert the unit to correct units.'
      write(*,*) 'Also the sum of weights MUST equal 1.'
      write(*,*) !
      if ( line == 0 ) then
         write(*,*) 'Error occured on reading number of k-points (first line)'
      else
         write(*,'(a,i0,a)') 'Error occured on reading the ',line,' kpoint.'
      end if
      write(*,*) 'Please format your file like this:'
      write(*,*) ' $> cat '//trim(fname)
      write(*,*) ' <nkpt>'
      write(*,*) '     1  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) '     2  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) ' ....'
      write(*,*) ' <nkpt> <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      
      call die('TBT reading user specified k-points')
      
    end subroutine kill_iokp
    
  end subroutine tbt_iokp_read
  
end module m_tbt_kpoint
