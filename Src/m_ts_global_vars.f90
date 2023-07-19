!==========================================================================*
!                                                                          *
!  TRANSIESTA MODULE m_ts_global_vars : Declaration of the TS variables    *
!  that are accessed in different parts of the code by using a:            *
!      use m_ts_global_vars                                                *
!  declaration, instead of passing as dummy arguments                      *  
!                                                                          *
!  Written by F.D.Novaes, Apr'10                                           *
!==========================================================================*                                         


module m_ts_global_vars

USE precision, only : dp

integer :: TSiscf=0

! Used to hold the "input" H when mixing the Hamiltonian
real(dp), dimension(:,:), pointer :: Hold => null()   

 
logical :: TSinit=.false.,TSrun=.false.

integer :: ts_istep      ! FC step in phonon calculation

end module m_ts_global_vars
