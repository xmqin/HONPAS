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
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

module m_ts_mumps_init

  use precision, only : dp

  implicit none

#ifdef SIESTA__MUMPS 

  integer, public, save :: MUMPS_mem = 20
  integer, public, save :: MUMPS_ordering = 7
  integer, public, save :: MUMPS_block = -8 ! blocking factor

  public :: read_ts_mumps
  public :: mum_err
  public :: init_MUMPS
  public :: analyze_MUMPS
  public :: prep_LHS
  public :: prep_RHS_Eq
  public :: prep_RHS_nEq
  public :: insert_Self_Energies

  private
  
contains

  subroutine read_ts_mumps( )

    use fdf, only : fdf_get, leqi
    character(len=200) :: chars

    MUMPS_mem   = fdf_get('TS.MUMPS.Mem',20)
    MUMPS_block = fdf_get('TS.MUMPS.BlockingFactor',112)
    chars = fdf_get('TS.MUMPS.Ordering','auto')
    if ( leqi(chars,'auto') ) then
       MUMPS_ordering = 7
    else if ( leqi(chars,'amd') ) then
       MUMPS_ordering = 0
    else if ( leqi(chars,'amf') ) then
       MUMPS_ordering = 2
    else if ( leqi(chars,'scotch') ) then
       MUMPS_ordering = 3
    else if ( leqi(chars,'pord') ) then
       MUMPS_ordering = 4
    else if ( leqi(chars,'metis') ) then
       MUMPS_ordering = 5
    else if ( leqi(chars,'qamd') ) then
       MUMPS_ordering = 6
    else
       call die('Unknown MUMPS ordering.')
    end if


  end subroutine read_ts_mumps

  subroutine init_MUMPS(mum,ID)
#ifdef MPI
    use mpi_siesta
#endif
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: ID
    character(len=25) :: file

    ! define processors and MATRIX-type
    mum%SYM =  0 ! non-symmetric
    mum%PAR =  1 ! sequential
    mum%COMM = MPI_COMM_SELF ! only communicate with it-self

    mum%JOB = -1 ! initialise
    call zMUMPS(mum)
    call mum_err(mum, 'MUMPS initialization had an error.')

    ! print out control (verbosity)
    mum%ICNTL(4) = 2

    ! Setup the control parameters
    mum%ICNTL(5) = 0 ! assembled format
    ! The ordering of the matrix
    mum%ICNTL(7) = MUMPS_ordering

    ! We request a sparse right-hand side: 
    !   20 == 1, we need not initialize RHS_SPARSE!
    mum%ICNTL(20) = 1
    mum%ICNTL(21) = 0
    
    ! Allow memory increase handled by the user
    mum%ICNTL(14) = MUMPS_Mem

    ! Sets the blocking factor (opt for change in future versions)
    ! Current MUMPS version: 4.10.0
    mum%ICNTL(27) = MUMPS_block

    ! Request specific elements of the inverse matrix: 
    !   30 == 1 we MUST allocate it, but need not initialize it
    mum%ICNTL(30) = 1
    ! this however, posses some other problems.

    ! For each processor, add their own output files.
    ! As they become quite big we make them rewritten 
    ! in every SCF
    write(file,'(a,i0,a)') 'TS_MUMPS_',ID,'.dat'
    call io_assign(mum%ICNTL(1))
    mum%ICNTL(2) = mum%ICNTL(1)
    mum%ICNTL(3) = mum%ICNTL(1)
    open(mum%ICNTL(1),file=trim(file),status='replace', &
         action='write',form='formatted')

  end subroutine init_MUMPS

  subroutine analyze_MUMPS(mum)
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer :: iu

    ! We force the analysis to be done without
    ! knowing the values of A
    ! otherwise A contains "spurious" values
    mum%ICNTL(6) = 1

    ! analyse the MUMPS solver, this will determine
    ! factorization strategy
    mum%JOB = 1
    call zMUMPS(mum)
    call mum_err(mum, 'MUMPS analysis step had an error.')

    ! Write out estimated memory requirements
    iu = mum%ICNTL(1)
    write(iu,'(/,a)')'### Memory estimation from ANALYSIS step...'
    write(iu,'(a,i0,a)')'### Minimum memory required for calculation: ', &
         mum%INFO(15),' MB'
    write(iu,'(a,i0,a)')'### MUMPS is allocating: ', &
         nint(mum%INFO(15)*real(100+mum%ICNTL(14),dp)/100._dp),' MB'
    write(iu,'(a,/)')"### MUMPS memory can be altered using TS.MUMPS.Mem."

  end subroutine analyze_MUMPS

  subroutine prep_LHS(IsVolt, mum,N_Elec,Elecs)
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_elec_se, only: UC_minimum_worksize
    use m_ts_sparse, only : ts_sp_uc
    use m_ts_electype
    use m_ts_method, only : orb_offset
    use create_Sparsity_Union, only: crtSparsity_Union
    include 'zmumps_struc.h'

    logical, intent(in) :: IsVolt
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)

    type(OrbitalDistribution) :: dit
    type(Sparsity) :: tmpSp1, tmpSp2
    integer :: iEl, idx, no, io, ind, nr, ioff
    integer, pointer :: l_ptr(:),l_ncol(:), l_col(:)

#ifdef MPI
    call newDistribution(nrows_g(ts_sp_uc),MPI_Comm_Self,dit, &
         name='MUMPS UC distribution')
#else    
    call newDistribution(nrows_g(ts_sp_uc),-1,dit, &
         name='MUMPS UC distribution')
#endif

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    tmpSp2 = ts_sp_uc
    do iEl = 1 , N_Elec

       idx = Elecs(iEl)%idx_o
       no = TotUsedOrbs(Elecs(iEl))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(dit,tmpSp1, &
            idx,idx,no,no, tmpSp2)
    end do
    call delete(tmpSp1)
    call delete(dit)

    ! Create the index 
    call attach(tmpSp2,list_ptr=l_ptr, &
         n_col=l_ncol,list_col=l_col,nrows_g=nr, &
         nnzs=mum%NZ)
    call UC_minimum_worksize(IsVolt, N_Elec, Elecs, io)
    io = max(mum%NZ, io)

    ! Allocate LHS, first we need to
    ! calculate the actual size of the matrix
    ! We know that it must be the Hamiltonian sparsity
    ! pattern UNION DENSE-electrodes!
    call memory('A', 'I', mum%NZ * 2, 'prep_LHS')
    allocate( mum%IRN( mum%NZ ) )
    allocate( mum%JCN( mum%NZ ) )
    call memory('A', 'Z', io, 'prep_LHS')
    allocate( mum%A( io ) )

!$OMP parallel do default(shared), private(io,ioff,ind)
    do io = 1 , nr

       if ( l_ncol(io) /= 0 ) then
       
       ioff = io - orb_offset(io)
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 
             
          mum%IRN(ind) = l_col(ind) - orb_offset(l_col(ind))
          mum%JCN(ind) = ioff

       end do

       end if
    end do
!$OMP end parallel do

    call delete(tmpSp2)

  end subroutine prep_LHS

  subroutine allocate_mum(mum,nsize,N_Elec,Elecs,GF)
    use m_ts_electype
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: nsize, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)

    integer :: no, iEl, io

    ! We allocate for the equilibrium GF
    mum%NZ_RHS = nsize ! number of non-zero RHS elements

    if ( associated(mum%IRHS_PTR) ) then
       call memory('D','I',size(mum%IRHS_PTR), 'allocate_mum')
       deallocate(mum%IRHS_PTR)
       call memory('D','I',size(mum%IRHS_SPARSE), 'allocate_mum')
       deallocate(mum%IRHS_SPARSE)
       nullify(mum%IRHS_PTR,mum%IRHS_SPARSE,mum%RHS_SPARSE)
    end if
    call memory('A','I',mum%NRHS+1, 'allocate_mum')
    allocate( mum%IRHS_PTR(mum%NRHS+1) )
    call memory('A','I',mum%NZ_RHS, 'allocate_mum')
    allocate( mum%IRHS_SPARSE(mum%NZ_RHS) )

    no = 0
    do iEl = 1 , N_Elec
       io = TotUsedOrbs(Elecs(iEl))
       no = no + io ** 2
    end do
    ! Allocate maximum space available
    if ( associated(Gf) ) then
       call memory('D','Z',size(Gf), 'allocate_mum')
       deallocate(Gf)
       nullify(Gf)
    end if
    ! for large electrodes and not so large device
    ! we need to allocate more space
    allocate( Gf(max(no,mum%NZ_RHS)) )
    call memory('A','Z',size(Gf), 'allocate_mum')
    mum%RHS_SPARSE => Gf(1:mum%NZ_RHS)

    no = 0
    do iEl = 1 , N_Elec

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       io = TotUsedOrbs(Elecs(iEl)) ** 2
       Elecs(iEl)%Sigma => Gf(no+1:no+io)
       no = no + io

    end do

  end subroutine allocate_mum

  subroutine prep_RHS_Eq(mum,no_u_TS,nzs,N_Elec,Elecs,GF)
    use m_ts_electype
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_method, only : orb_offset, orb_type, TYP_BUFFER
    use class_Sparsity
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: no_u_TS, nzs, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)
    integer :: io, j, nr, ind
    integer, pointer :: l_ptr(:), l_ncol(:), l_col(:)

    ! Create the index 
    call attach(tsup_sp_uc,list_ptr=l_ptr, &
         n_col=l_ncol,list_col=l_col,nrows_g=nr)

    call allocate_mum(mum,nzs,N_Elec,Elecs,GF)

    ! TODO, this requires that the sparsity pattern is symmetric
    ! Which it always is!
    mum%IRHS_PTR(:) = 1 

!$OMP parallel do default(shared), private(io,j,ind)
    do io = 1 , nr
       if ( orb_type(io) /= TYP_BUFFER ) then
       mum%IRHS_PTR(io-orb_offset(io)) = l_ptr(io) + 1
       if ( l_ncol(io) /= 0 ) then ! no entries
       ! Create the row-index
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          mum%IRHS_SPARSE(ind) = l_col(ind) - orb_offset(l_col(ind))
       end do

       end if
       end if

    end do
!$OMP end parallel do

    mum%IRHS_PTR(no_u_TS+1) = nzs + 1

  end subroutine prep_RHS_Eq

  subroutine prep_RHS_nEq(mum,no_u_TS,N_Elec, Elecs,Gf)
    use m_ts_electype
    use m_ts_method, only : ts2s_orb
    use class_Sparsity
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum

    integer, intent(in) :: no_u_TS, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)
    integer :: iEl, no, i, j, io, ind

    ! We only need a partial size of the Green function
    no = sum(TotUsedOrbs(Elecs))

    call allocate_mum(mum,no*no_u_TS,N_Elec,Elecs,GF)

    ind = 0
    do i = 1 , no_u_TS
       mum%IRHS_PTR(i) = ind + 1
       ! get correct siesta-orbital
       io = ts2s_orb(i)
       iElec: do iEl = 1 , N_Elec
          if ( .not. OrbInElec(Elecs(iEl),io) ) cycle
          ! Create the row-index
          do j = 1 , no_u_TS
             ind = ind + 1
             mum%IRHS_SPARSE(ind) = j
          end do
          exit iElec
       end do iElec
    end do
    mum%IRHS_PTR(no_u_TS+1) = no*no_u_TS + 1

    if ( ind /= mum%NZ_RHS ) then
       call die('Error in sparsity pattern of non-equilibrium')
    end if

  end subroutine prep_RHS_nEq

  subroutine insert_Self_Energies(mum, El)
    use m_ts_electype
    use m_ts_method, only : orb_offset, ts2s_orb
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    type(Elec), intent(in) :: El

    integer :: off, ii, no, ind, iso, jso

    no = TotUsedOrbs(El)
    off = El%idx_o - 1

    if ( El%Bulk ) then
!$OMP do private(ind,jso,iso,ii)
       do ind = 1 , mum%NZ
          jso = ts2s_orb(mum%JCN(ind))
          if ( OrbInElec(El,jso) ) then 
          iso = ts2s_orb(mum%IRN(ind))
          if ( OrbInElec(El,iso) ) then
             ii = (jso - El%idx_o) * no + iso - off
             mum%A(ind) = El%Sigma(ii)
          end if
          end if
       end do
!$OMP end do nowait
    else
!$OMP do private(ind,jso,iso,ii)
       do ind = 1 , mum%NZ
          jso = ts2s_orb(mum%JCN(ind))
          if ( OrbInElec(El,jso) ) then
          iso = ts2s_orb(mum%IRN(ind))
          if ( OrbInElec(El,iso) ) then
             ii = (jso - El%idx_o) * no + iso - off
             mum%A(ind) = mum%A(ind) - El%Sigma(ii)
          end if
          end if
       end do
!$OMP end do nowait
    end if

  end subroutine insert_Self_Energies


  subroutine mum_err(mum,m)
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    character(len=*), intent(in) :: m

    ! We here check the error message that MUMPS
    ! has given
    if ( mum%INFOG(1) == 0 .and. mum%INFO(1) == 0 ) return ! no error

    ! write out the message
    call msg(m)

    select case ( mum%INFOG(1) ) 
    case ( -5 ) 
       call die('Analysis step could not allocate space. &
            &Do you have enough memory on your machine?')
    case ( -8 , -14 )
       call msg('Integer work-array too small, increase &
            &TS.MUMPS.Mem.')
       write(*,'(a,i0)') 'Additional elements are needed: ',mum%INFOG(2)
       write(0,'(a,i0)') 'Additional elements are needed: ',mum%INFOG(2)
    case ( -9 ) 
       call msg('Work-array S too small, increase &
            &TS.MUMPS.Mem.')
       write(*,'(a,i0)') 'Additional elements are needed: ',-mum%INFOG(2)*1e6
       write(0,'(a,i0)') 'Additional elements are needed: ',-mum%INFOG(2)*1e6
    case ( -17, -20 ) 
       call msg('Work-array S too small, increase &
            &TS.MUMPS.Mem.')
    end select

    call die('Check the MUMPS documentation for &
         &the error message as well as the output.')
    
  contains
    
    subroutine msg(m)
      character(len=*), intent(in) :: m
      write(*,'(a)') m
      write(0,'(a)') m
    end subroutine msg
    
  end subroutine mum_err
    

#endif
end module m_ts_mumps_init
