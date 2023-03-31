!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for correcting the density matrix for retaining a constant charge density
! The idea is to introduce several different schemes of charge corrections.

module ts_charge_m

  use precision, only: dp
  
  implicit none

  public 

  ! Info parameters for obtaining charge calculations (mulliken charges in regions)
  integer, parameter :: TS_Q_INFO_FULL = 0
  integer, parameter :: TS_Q_INFO_SCF = 1

  private :: dp

contains

  ! Retrive the mulliken charges in each region of the transiesta setup
  subroutine ts_charge_get(N_Elec,dit, sp, nspin, n_nzs, DM, S, Q, Qtot)

    use m_ts_method
    use parallel, only : Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB
    use m_ts_electype

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: N_Elec
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The charge in the regions
    real(dp), intent(out), optional :: Q(0:1+1+N_Elec*2, nspin), Qtot
    
! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ir, jr, r
    real(dp) :: Qtmp(0:1+1+N_Elec*2, nspin)
#ifdef MPI
    real(dp) :: tmp
    integer :: MPIerror
#endif

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    ! Initialize charges
    Qtmp(:,:) = 0._dp

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ir,ind,jo,jr,r), &
!$OMP&reduction(+:Qtmp)
    do lio = 1 , no_lo

       ! obtain the global index of the orbital.
       io = index_local_to_global(dit,lio,Node)
       ir = orb_type(io)

       ! Loop number of entries in the row... (index frame)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          ! as the local sparsity pattern is a super-cell pattern,
          ! we need to check the unit-cell orbital
          ! The unit-cell column index
          jo = UCORB(l_col(ind),no_u)
          jr = orb_type(jo)

          if ( ir == TYP_BUFFER .and. jr == TYP_BUFFER ) then
            r = 1 ! buffer
          else if ( ir == TYP_BUFFER .or. jr == TYP_BUFFER ) then
            r = 0 ! other
          else if ( ir == TYP_DEVICE .and. jr == TYP_DEVICE ) then
            r = 2 ! device
          else if ( ir == TYP_DEVICE .or. jr == TYP_DEVICE ) then
            r = 4+(ir+jr-1)*2 ! device/electrode
          else if ( ir == jr ) then
            r = 3+(ir-1)*2 ! electrode/electrode
          else
            r = 0 ! other
          end if
          Qtmp(r,:) = Qtmp(r,:) + DM(ind,:) * S(ind)
       end do
    end do
!$OMP end parallel do

#ifdef MPI
    if ( present(Q) .and. present(Qtot) ) then
      call MPI_AllReduce(Qtmp(0,1),Q(0,1),size(Qtmp), &
          MPI_Double_Precision, MPI_SUM, MPI_Comm_World,MPIerror)
      Qtot = sum(Q)
    else if ( present(Q) ) then
      call MPI_AllReduce(Qtmp(0,1),Q(0,1),size(Qtmp), &
          MPI_Double_Precision,MPI_SUM, MPI_Comm_World,MPIerror)
    else if ( present(Qtot) ) then
      tmp = sum(Qtmp)
      call MPI_AllReduce(tmp, Qtot, 1, &
          MPI_Double_Precision,MPI_SUM, MPI_Comm_World,MPIerror)
    end if
#else
    if ( present(Q)    ) Q    = Qtmp
    if ( present(Qtot) ) Qtot = sum(Qtmp)
#endif

  end subroutine ts_charge_get

  ! A subroutine for printing out the charge distribution in the cell
  ! it will currently only handle the full charge distribution, and
  ! not per k-point.
  subroutine ts_charge_print(N_Elec,Elecs,Qtot,dit, sp, &
      nspin, n_nzs, DM, S, &
      method)
    use parallel, only : IONode
    use m_ts_electype
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

    ! **********************
    ! * INPUT variables    *
    ! **********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The requested number of electrons in the simulation
    real(dp), intent(in) :: Qtot
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    ! The method by which it should be printed out...
    integer, intent(in), optional :: method

    ! **********************
    ! * LOCAL variables    *
    ! **********************
    integer :: i
    real(dp), allocatable :: Q(:,:)
    integer :: ispin, lmethod
    logical :: has_buffer

    lmethod = TS_Q_INFO_FULL
    if ( present(method) ) lmethod = method

    allocate(Q(0:2+N_Elec*2,nspin))

    call ts_charge_get(N_Elec,dit, sp, nspin, n_nzs, DM, S, Q = Q)
    has_buffer = sum(Q(1,:)) > 0._dp

    ! it will only be the IONode which will write out...
    if ( .not. IONode ) then
      deallocate(Q)
      return
    end if

    if ( lmethod == TS_Q_INFO_FULL ) then
      write(*,'(/,a,f12.5)') 'transiesta: Charge distribution, target = ',Qtot
      if ( nspin > 1 ) then
        write(*,'(a,3(tr1,f12.5))') &
            'Total charge                  [Q]  :', &
            sum(Q(:,1)),sum(Q(:,2)),sum(Q)
        write(*,'(a,2(tr1,f12.5))') &
            'Device                        [D]  :',Q(2,1), Q(2,2)
        do i = 1 , N_Elec
          write(*,'(a,t31,a,i0,a,2(tr1,f12.5))') &
              trim(name(Elecs(i))),'[E',i,'] :', &
              Q(3+(i-1)*2,1), Q(3+(i-1)*2,2)
          write(*,'(a,t22,a,i0,a,2(tr1,f12.5))') &
              trim(name(Elecs(i))),'/ device [C',i,'] :', &
              Q(4+(i-1)*2,1), Q(4+(i-1)*2,2)
        end do
        write(*,'(a,2(tr1,f12.5))') &
            'Other                         [O]  :',Q(0,1), Q(0,2)
        if ( has_buffer ) then
          write(*,'(a,2(tr1,f12.5))') &
            'Buffer                        [B]  :',Q(1,1), Q(1,2)
        end if
      else
        write(*,'(a,tr1,f12.5)') &
            'Total charge                  [Q]  :', sum(Q(:,1))
        write(*,'(a,tr1,f12.5)') &
            'Device                        [D]  :',Q(2,1)
        do i = 1 , N_Elec
          write(*,'(a,t31,a,i0,a,tr1,f12.5)') &
              trim(name(Elecs(i)))         ,'[E',i,'] :',Q(3+(i-1)*2,1)
          write(*,'(a,t22,a,i0,a,tr1,f12.5)') &
              trim(name(Elecs(i))),'/ device [C',i,'] :',Q(4+(i-1)*2,1)
        end do
        if ( has_buffer ) then
          write(*,'(a,tr1,f12.5)') &
            'Buffer                        [B]  :',Q(1,1)
        end if
        write(*,'(a,tr1,f12.5)') &
            'Other                         [O]  :',Q(0,1)
      end if
      write(*,'(a,tr1,es12.5,/)') &
          'Excess charge                [dQ]  :',sum(Q) - Qtot


    else if ( lmethod == TS_Q_INFO_SCF ) then

      ! We write out the information from the SCF cycle...
      write(*,'(a,1x,a9)',advance='no') 'ts-q:','D'
      do i = 1 , N_Elec
        if ( i > 9 ) then
          write(*,'(1x,a7,i2,1x,a7,i2)',advance='no') 'E',i,'C',i
        else
          write(*,'(1x,a8,i1,1x,a8,i1)',advance='no') 'E',i,'C',i
        end if
      end do
      if ( has_buffer ) then
        write(*,'(1x,a9)',advance='no') 'B'
      end if
      if ( nspin > 1 ) then
        write(*,'(2(1x,a9))') 'dQ','Qup-Qdn'
      else
        write(*,'(1x,a9)') 'dQ'
      end if
      write(*,'(a,1x,f9.3)',advance='no') 'ts-q:', sum(Q(2,:))
      do i = 1 , N_Elec
        write(*,'(2(1x,f9.3))',advance='no') sum(Q(3+(i-1)*2,:)),sum(Q(4+(i-1)*2,:))
      end do
      if ( has_buffer ) then
        write(*,'(1x,f9.3)',advance='no') sum(Q(1,:))
      end if
      if ( nspin == 2 ) then
        write(*,'(2(1x,e9.3))') sum(Q) - Qtot, sum(Q(:,1)) - sum(Q(:,2))
      else
        write(*,'(1x,e9.3)') sum(Q) - Qtot
      end if

    end if

    deallocate(Q)

  end subroutine ts_charge_print

end module ts_charge_m
