module m_spin
  use precision, only: dp
  private

  ! Number of spin components
  integer, save, public           :: nspin

  ! If diagonalization mixes spins, we need to double the block
  ! size for the parallel distribution of matrices. 
  ! In that case, MColl=2,  else MColl=1
  integer, save, public           :: MColl     

  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  public :: init_spin

CONTAINS

  subroutine init_spin()
    use m_fdf_global
    use alloc, only: re_alloc
    use parallel, only: IOnode

    implicit none

    logical  :: SPpol, NonCol

    call fdf_global_get(SPpol,'SpinPolarized',.false.)
    call fdf_global_get(NonCol,'NonCollinearSpin',.false.)

!
!   NonCol takes precedence. The printout clarifies
!   whether Up/Down or non-collinear spin is used.
!
    if (NonCol) then
       nspin     = 4
       MColl     = 2
       SPpol     = .false.
    elseif (SPpol) then
       nspin     = 2
       MColl     = 1
    else 
       nspin     = 1
       MColl     = 1
    endif

    if (ionode) then
       write(6,'(a,4x,l1)') 'redata: Non-Collinear-spin run           = ',NonCol
       write(6,'(a,4x,l1)') 'redata: SpinPolarized (Up/Down) run      = ',SPpol
       write(6,'(a,4x,i1)') 'redata: Number of spin components        = ',nspin
    end if

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

  end subroutine init_spin

end module m_spin
