! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_normalize_dm

  implicit none

  private
  
  public :: normalize_dm

contains

  subroutine normalize_dm( first )
    
    use precision, only: dp
    use sparse_matrices, only: Dscf, Escf, maxnh, S
    use siesta_options,  only: dm_normalization_tol
    use siesta_options,  only: normalize_dm_during_scf
    use siesta_options,  only: mix_scf_first
    use atomlist, only: qtot
    use parallel, only: IOnode
    use sys,      only: die
    use m_spin,   only: spin
#ifdef MPI
    use m_mpi_utils, only: globalize_sum
#endif
    
    ! Whether the SCF step is the first
    ! In this case a different tolerance is used...
    ! Say if reading an old DM
    logical, intent(in) :: first

    integer :: io, is
    ! Total unnormalized electron charge
    real(dp):: qsol
    ! Scaling for normalization
    real(dp):: scale
    character(len=132) :: msg
#ifdef MPI
    real(dp):: buffer1
#endif

    ! Normalize density matrix to exact charge
    
    qsol = 0.0_dp
    do is = 1 , spin%spinor
       do io = 1 , maxnh
          qsol = qsol + Dscf(io,is) * S(io)
       end do
    end do
    
#ifdef MPI
    call globalize_sum(qsol, buffer1)
    qsol = buffer1
#endif

    ! Calculate the relative difference from the total charge
    scale = abs(qsol / qtot - 1._dp)
    
    ! We might consider using abs(qsol-qtot) instead of the fractional accuracy.
    ! For large numbers of electrons one might have a significant net charge
    ! before the fractional tolerance is reached.

    ! Two degrees of reporting:
    if ( first ) then
       ! If 1st scf step, it might be that a (reused) DM was not
       ! normalized properly
       if ( scale > 1.0e-3_dp ) then
          if (IOnode) then
             write(6,'(a,2f20.8)') &
                  'Note: For starting DM, Qtot, Tr[D*S] =', qtot, qsol
          end if
       end if
       ! We change whether one may mix the first SCF step
       ! in case there is missing some initial charge
       !   10e-4 electrons
       if ( abs(qsol - qtot) > 0.0001_dp ) then
          ! See discussion in m_new_dm.F90 (end of new_DM)
          mix_scf_first = .false.
       end if
    else
       ! In later steps, the lack of normalization is more serious
       ! The tolerance is tighter by default
       if ( scale > dm_normalization_tol ) then
          write(msg,'(a,2f20.8)') &
               'Bad DM normalization: Qtot, Tr[D*S] =', qtot, qsol
          call die(msg)
       end if
    end if
    
    ! Actually normalize only if not disallowed by the user
    if ( normalize_dm_during_scf ) then
       
       ! Figure out the scaling factor for normalization
       scale = qtot / qsol
       
!$OMP parallel default(shared), private(is,io)
       
       do is = 1 , spin%DM
!$OMP do
         do io = 1 , maxnh
             Dscf(io,is) = Dscf(io,is) * scale
          end do
!$OMP end do nowait
       end do

       do is = 1 , spin%EDM
!$OMP do
          do io = 1 , maxnh
             Escf(io,is) = Escf(io,is) * scale
          end do
!$OMP end do nowait
       end do

!$OMP end parallel

    else if ( IONode ) then
       write(6,'(a,g15.6)') &
            'Note: Not normalizing DM: log(Tr[D*S]/Qtot - 1) =', &
            log10(scale)
    end if
      
  end subroutine normalize_dm

end module m_normalize_dm
