subroutine print_spin(qspin)
  !
  ! Prints spin in output and CML files
  !
  use m_spin,          only: spin
  use atomlist,        only: no_l
  use sparse_matrices, only: listhptr, numh
  use sparse_matrices, only: S, Dscf   ! Dscf could have 1, 2, 4, or 8 components
  use siesta_cml
  use parallel,        only: IOnode
  use precision,       only: dp
#ifdef MPI
  use m_mpi_utils,     only: globalize_sum
#endif

  implicit none
  
  real(dp), intent(out)  :: qspin(spin%Grid)

  integer  :: ispin, io, j, ind
  real(dp) :: qaux
  real(dp) :: Stot        ! Total spin magnitude
  real(dp) :: Svec(3)     ! Total spin vector
#ifdef MPI
  real(dp) :: qtmp(spin%Grid)
#endif

  qspin(:) = 0.0_dp

  if ( spin%Grid < 2 ) return

  if ( spin%SO ) then
     
    do io = 1,no_l
      do j = 1,numh(io)
        ind = listhptr(io)+j
        ! In the SOC case, hermitify Dscf
        qspin(1:2) = qspin(1:2) + dscf(ind,1:2) * S(ind)
        qspin(3) = qspin(3) + 0.5_dp*(dscf(ind,3)+dscf(ind,7)) * S(ind)
        qspin(4) = qspin(4) + 0.5_dp*(dscf(ind,4)+dscf(ind,8)) * S(ind)
      end do
    end do

  else
     
    do io = 1,no_l
      do j = 1,numh(io)
        ind = listhptr(io)+j
        qspin(:) = qspin(:) + Dscf(ind,:) * S(ind)
      end do
    end do
     
  endif
#ifdef MPI
  ! Global reduction of spin components
  call globalize_sum(qspin(1:spin%Grid),qtmp(1:spin%Grid))
  qspin(1:spin%Grid) = qtmp(1:spin%Grid)
#endif

  ! We are left with printing out to stdout
  if ( .not. IONode ) return
  
  if ( spin%Grid == 2) then

    Svec(1) = 0.0_dp
    Svec(2) = 0.0_dp
    Svec(3) = qspin(1) - qspin(2)
    Stot = Svec(3)
    write(6,'(5x,a,f10.5,2f10.1,f10.5)') 'spin moment: S , {S} = ', Stot, Svec
    if (cml_p) call cmlAddProperty(xf=mainXML,            &
        value=qspin(1)-qspin(2), dictref='siesta:stot', &
        units='siestaUnits:spin')
    
  else if ( spin%Grid == 4 ) then

    call spnvec( spin%Grid, qspin, qaux, Stot, Svec )
    write(6,'(5x,a,4f10.5)') 'spin moment: S , {S} = ', Stot, Svec
    if (cml_p) then
      call cmlAddProperty(xf=mainXML, value=Stot,  &
          dictref='siesta:stot', units='siestaUnits:spin')
      call cmlAddProperty(xf=mainXML, value=Svec,  &
          dictref='siesta:svec', units='siestaUnits:spin')
    end if !cml_p

  end if

end subroutine print_spin
