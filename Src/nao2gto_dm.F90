! *** Module: nao2gto_dm ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Density-matrix utilities for NAO2GTO
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_dm

  use precision,     only: dp      ! Double precision
  use parallel,      only: Node    ! Local node   
  use parallel,      only: Nodes   ! Total number of nodes
  use atmfuncs,      only: lofio   ! Returns angular momentum number
  use atmfuncs,      only: mofio   ! Returns magnetic quantum number
  use atomlist,      only: indxuo  ! Index of equivalent orbital in
                                   !   the unit cell
  use listsc_module, only: listsc  ! Orbital, equivalent to JUO,  
                                   !   which is related to IO
                                   !   exactly as JUO is related to IUO
  use nao2gto_data,  only: nso     ! Number of Slater orbitals per angular
                                   !   momentum shell

  implicit none

  private

  public :: sparse_dense, build_D2Sindx,  get_pmax_shell

contains

! *****************************************************************************
!> \brief Builds the dense density matrix from the sparse one.
!! 
!!   Here we apply the condition written right below
!!   Eq. (45) of the SIESTA paper
!!   J. M. Soler et al.,
!!   J. Phys.: Condens. Matter 14, 2745 (2002).
!!   \f$\rho_{\mu \nu} \equiv \rho_{\mu^\prime \nu^\prime}\f$
!!   if \f$ (\nu \mu) \f$ are equivalent by translation to
!!   \f$(\nu^\prime \mu^\prime)\f$
!!
!!   The output is written in dense format in the variable M_dense, whose
!!
!!   First  index: number of orbitals in the supercell
!!
!!   Second index: number of orbitals in the supercell
!!
!!   Third  index: spin
!! 
!!   As implemented here, all the variables are Node independent
!!   (i.e. all the nodes know everything, since all the relevant
!!   variable have been previously globalized)
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nuo: number of orbitals within the unit cell 
!! \param[in] norb: local number of orbitals within the supercell
!!              (norb=nuotot*ncells)
!! \param[in] maxnd: maximum size of the sparse matrix
!! \param[in] numd(nuo): number of nonzero elements of each row
!! \param[in] listd(maxnd): nonzero hamiltonian-matrix element
!! \param[in] listdptr(nuo): control vector of density matrix
!! \param[in] M_sparse(maxnd): sparse matrix
!! \param[out] M_dense(maxnd): dense matrix
! *****************************************************************************
  subroutine sparse_dense(nspin, nuo, nuotot, norb, maxnd, numd, listdptr, &
&                listd, M_sparse, M_dense)

    use precision, only: dp

    implicit none

    integer , intent(in)  :: nspin   ! Number of spin components 
    integer , intent(in)  :: nuo     ! Same as nuotot in this subroutine
    integer , intent(in)  :: nuotot  ! Total number of orbitals in the unit cell
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
    integer , intent(in)  :: norb    ! Total number of orbitals in the supercell
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
    integer , intent(in)  :: maxnd   ! Maximum number of elements of the
                                     !   density matrix
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: numd(nuo)
                                     ! Number of non-zero element
                                     !   per row in the matrices
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: listdptr(nuo)
                                     ! Pointer to the start of
                                     !   each row (-1) of the
                                     !   hamiltonian matrix 
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: listd(*)
                                     ! Nonzero hamiltonian-matrix
                                     !   element column indices
                                     !   for each matrix row
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    real(dp), intent(in)  :: M_sparse(maxnd, nspin)
                                     ! Density matrix in sparse form
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    real(dp), intent(out) :: M_dense(norb,norb,nspin)
                                     ! Density matrix in dense form
                                     ! NOTE: When running in parallel,
                                     !   this is core independent

    ! Local variables
    integer :: ispin, io, jo, j, iu, ju, ind
    real(dp), dimension(:,:), pointer :: M_aux => null()

    ! -------------------------------------------------------------------------

    allocate(M_aux(nuotot,norb))

    M_aux(:,:)     = 0.0_dp
    M_dense(:,:,:) = 0.0_dp

    do ispin=1,nspin
      do io=1,nuotot
        do j=1,numd(io)
          ind = listdptr(io) + j
          jo = listd(ind)
          M_aux(io,jo) = M_sparse(ind,ispin)
          M_dense(io,jo,ispin) = M_aux(io,jo)
        enddo
      enddo

      do io=nuotot+1,norb
        iu = indxuo(io)
        do j=1,numd(iu)
          ind=listdptr(iu) +j
          ju=listd(ind)
          jo=listsc(io,iu,ju)
          M_dense(io,jo,ispin)=M_aux(iu,ju)
        enddo
      enddo
    enddo

    deallocate(M_aux)

  end subroutine sparse_dense


  subroutine build_D2Sindx(nuotot, norb, maxnd, numd, listdptr,listd, D2Sindx)

       implicit none

       integer , intent(in)  :: nuotot,norb, maxnd
       integer , intent(in)  :: numd(nuotot), listdptr(nuotot),listd(maxnd)
       integer, intent(inout) :: D2Sindx(norb,norb)
      ! Local variables
       integer :: io, jo, j, iu, ju, ind

       do io = 1, norb
          iu = indxuo(io)
        do j=1,numd(iu)
           ind=listdptr(iu) + j
           ju=listd(ind)
           jo=listsc(io,iu,ju)
           D2Sindx(io,jo) = ind
         enddo
      enddo

   end subroutine build_D2Sindx

! *****************************************************************************
!> \brief Computes the largest value of the density matrix between two
!! atomic orbitals of given angular momenta in the supercell
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nspin:    Number of spin states
!! \param[in] norb:     Number of orbitals in the supercell
!! \param[in] iaorb:    Atomic index of each orbital
!! \param[in] iphorb:   Orbital index of each orbital in its atom
!! \param[in] nuo:      Number of orbitals in the unit cell 
!! \param[in] na:       Number of atoms in the supercell
!! \param[in] isa:      Species index of each atom
!! \param[in] maxnd:    Maximum number of elements of the density matrix
!! \param[in] numd:     Number of nonzero elements of each row
!! \param[in] listdptr: Control vector of density matrix
!! \param[in] listd:    Nonzero hamiltonian or density-matrix element
!! \param[in] M_dense:  Density matrix in dense format
!! \param[out] P_max:   Largest value of the density matrix between two 
!!                      atomic orbitals of given angular momentums 
!!                      in the supercell
! *****************************************************************************
  subroutine get_pmax_shell(nspin, norb, iaorb, iphorb, nuo, na, isa, &
&                maxnd, numd, listdptr, listd, D2Sindx, Dscf,  Dscf_max)

    ! Arguments
    integer , intent(in)  :: nspin   ! Number of spin components
    integer , intent(in)  :: norb    ! Total number of orbitals in the supercell
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
    integer , intent(in)  :: iaorb(norb)
                                     ! Atomic index of each orbital
    integer , intent(in)  :: iphorb(norb)
                                     ! Orbital index of each orbital
                                     !   in its atom
    integer , intent(in)  :: nuo     ! Total number of orbitals in the unit cell
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
    integer , intent(in)  :: na      ! Number of atoms in the supercell 
    integer , intent(in)  :: isa(na) ! Species index of each atom
    integer , intent(in)  :: maxnd   ! Maximum number of elements of the
                                     !   density matrix
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: numd(nuo)
                                     ! Number of non-zero element
                                     !   per row in the matrices
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: listdptr(nuo)
                                     ! Pointer to the start of
                                     !   each row (-1) of the
                                     !   hamiltonian matrix 
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: listd(*)
                                     ! Nonzero hamiltonian-matrix
                                     !   element column indices
                                     !   for each matrix row
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
                                     !   because it has been globalized
                                     !   before calling the subroutine
    integer , intent(in)  :: D2Sindx(norb, norb)
    real(dp), intent(in)  :: Dscf(maxnd, nspin)
                                     ! Density matrix in dense form
                                     ! NOTE: When running in parallel,
                                     !   this is core independent
    real(dp), intent(out) :: Dscf_max(maxnd)
                                     ! Largest value of the density matrix 
                                     !   between two atomic orbitals in the  
                                     !   supercell.
                                     !   Although it is defined as a square 
                                     !   matrix of size equal to the number
                                     !   of orbitals in the supercell,
                                     !   this is computed only the first time 
                                     !   that a given orbital with a defined 
                                     !   angular momentum appears.
                                     !   In that element of P_max, it stores
                                     !   the maximum value of the density 
                                     !   matrix between all the orbitals 
                                     !   for that atom with the same angular 
                                     !   momentum 
                                     !   Example: might be the maximum value
                                     !   of the DM would be between a 
                                     !   p_x (im=1) orbital within an atom
                                     !   and a d_yz (jm = -1) orbital of the
                                     !   neighbour.
                                     !   This value would be stored with the
                                     !   index corresponding to the p_y (im=-1)
                                     !   and d_xy (jm=-2) orbitals.
           
    ! Local variables
    integer :: ispin, ia, ja, io, ioa, is, jo, joa, js, j, iu, ju, ind, &
&     il, im, jl, jm

    integer :: io_in_shell, jo_in_shell, ind_in_shell

! ------------------------------------------------------------------------------

!   Loop on all the atomic orbitals in the supercell
!   all the nodes run the loop over all the atomic orbitals
    do ispin = 1, nspin
       do io  = 1, norb
!     Identify the atomic to which the orbital belongs to
          ia = iaorb(io)
!     Identify the index of the orbital within that atom
          ioa = iphorb(io)
!     Identify the species index of the atom to which the orbital belongs to
          is = isa(ia)
!     Identify the angular and magnetic quantum numbers
          il = lofio(is,ioa)
          im = mofio(is,ioa)
!     we are interested only in the first time a given il appears for that 
!     atom
          if ( im .ne. -il ) cycle

!     Identify the equivalent orbital within the unit cell
          iu = indxuo(io)

!     Run over all the neighbours
          do j=1, numd(iu)
!       Indentify the number of the neigbout in the list and the orbital index
!       in the supercell, jo
             ind = listdptr(iu) + j
             ju = listd(ind)
             jo = listsc(io,iu,ju)

!       Identify the atomic to which the neighbour orbital belongs to
             ja = iaorb(jo)

!       Identify the index of the neighbour orbital within that atom
             joa = iphorb(jo)

!       Identify the species index of the atom to which the neighbour 
!       orbital belongs to
             js = isa(ja)

!       Identify the angular and magnetic quantum numbers of the
!       neighbour orbital
             jl = lofio(js,joa)
             jm = mofio(js,joa)

!!       For debugging
!        write(6,'(a,7i5)')                                     &
! &        'get_pmax_shell: j, jo, ja, js, jl, jm, nso(jl) = ', &
! &         j, jo, ja, js, jl, jm, nso(jl)
!!       End debugging

!       jm runs between -jl <= jm <= jl
!       we are interested only in the first time a given il appears for that 
!       atom
            if ( jm .ne. -jl ) cycle

!       We select the largest value of the density matrix between two atomic
!       orbitals in the supercell

            do io_in_shell = io, io-1+nso(il)
               do jo_in_shell = jo, jo-1+nso(jl)
                  ind_in_shell = D2Sindx(io_in_shell,jo_in_shell)
                  if((ind_in_shell.ne.0) .and. (abs(Dscf(ind_in_shell,ispin)) .gt. Dscf_max(ind)) ) then
                     Dscf_max(ind)= abs(Dscf(ind_in_shell,ispin))
                  endif
               enddo
            enddo
!!       For debugging
!        write(6,'(a,5i5,f12.5)')                                &
! &        'get_pmax_shell: Node, io, jo, io-1+nso(il), jo:jo-1+nso(jl), P_max(io,jo) = ', &
! &         Node, io, jo, io-1+nso(il), jo-1+nso(jl), P_max(io,jo)
!!       End debugging
      enddo
    enddo
  enddo

  end subroutine get_pmax_shell

end module nao2gto_dm
