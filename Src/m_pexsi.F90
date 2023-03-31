
module m_pexsi
#ifdef SIESTA__PEXSI
  use precision, only: dp
  use iso_c_binding

  integer(c_intptr_t), public :: plan

  public :: pexsi_initialize_scfloop
  public :: pexsi_finalize_scfloop

  private

  CONTAINS

    subroutine pexsi_initialize_scfloop(World_Comm,npPerPole,mpirank,info)
      use f_ppexsi_interface
      integer, intent(in) :: npPerPole, mpirank
      integer, intent(in) :: World_Comm
      integer, intent(out):: info
    
      integer :: numProcRow, numProcCol
      integer :: outputFileIndex
    
      call get_row_col(npPerPole,numProcRow,numProcCol)
    
      ! Set the outputFileIndex to be the pole index.
      ! Starting from PEXSI v0.8.0, the first processor for each pole outputs
      ! information
    
      if( mod( mpirank, npPerPole ) .eq. 0 ) then
        outputFileIndex = mpirank / npPerPole;
      else
        outputFileIndex = -1;
      endif
    
      plan = f_ppexsi_plan_initialize(&
        World_Comm,&
        numProcRow,&
        numProcCol,&
        outputFileIndex,&
        info) 
    
    end subroutine pexsi_initialize_scfloop
    
    
    subroutine pexsi_finalize_scfloop()
      use f_ppexsi_interface
    
      integer :: info
    
      call f_ppexsi_plan_finalize( plan, info )
    
    end subroutine pexsi_finalize_scfloop
    
    subroutine get_row_col(np,nrow,ncol)
      integer, intent(in)  :: np
      integer, intent(out) :: nrow, ncol
      !
      ! Finds the factors nrow and ncol such that nrow*ncol=np,
      ! are as similar as possible, and nrow>=ncol.
      ! For prime np, ncol=1, nrow=np.
    
      ncol  = floor(sqrt(dble(np)))
      do
        nrow = np/ncol
        if (nrow*ncol == np) exit
        ncol = ncol - 1
      enddo
    end subroutine get_row_col
#endif
end module m_pexsi
