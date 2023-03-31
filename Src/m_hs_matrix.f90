!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_hs_matrix
!
! Routines that are used for the distribution of the Hamiltonian 
! and scattering matrices into a full size matrix instead of sparse matrices.
! It has several options for creating different kinds of matrices.
! This module is a serial version. It requires that the Node has the full
! sparse matrix available.
!
! Specifically it can be used to remove the z-direction connections.
! Also the inner-cell distances are in the xij arrays. Through routine calls
! these inner cell distances can be removed.
!
! The reason for having this is future use of routines which can utilize the
! Hamiltonian created in different ways.
! Say you want to test calculate a transmission from a Hamiltonian which is
! created by SIESTA. For that purpose you need to remove the cell connection 
! in the z-direction.
!
! The usage of this module is highly encouraged in future utilities where
! the need for the Hamiltonian and/or overlap matrix is needed.
! With this in mind the code can further be optimized for speed.
! However, for normal sizes it is quite fast.
! 
! Also for testing against Inelastica the use of inner-cell distances is not
! used. Therefore the Hamiltonians can not be numerically compared.
! If this is to be enforced later on. The option is there.
! 
! The use of this module is a straight forward call:
! 
!   call set_HS_matrix(Gamma,ucell,na_u,no_u,no_s,maxnh, &
!       xij,numh,listhptr,listh,H,S, &
!       k,Hk,Sk, &
!       xa,iaorb, &
!       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
! 
! xa, iaorb, RemZConnection, RemUCellDistances, RemNFirstOrbitals, RemNLastOrbitals
! are all optional arguments. xa and iaorb are required if RemZConnection or RemUCellDistances
! are true.
! Notice, that any calls to the optional arguments MUST be with keywords! Otherwise
! the program will end!
! We have added so that when iaorb is required, you could also use lasto. This
! makes it more intuitive as lasto is often more accesible.
! It on the other hand has a small overhead when using lasto (negligeble).
!
! Gamma denotes whether it is a Gamma calculation (if true, it will not
! add k-phases, no matter if k /= \Gamma-point.
! na_u,no_u,no_s,maxnh,xij,numh,listhptr,listh,H,S are all variables
! needed in the definition of the entire H and S matrices in the sparse format.
!
! k is the k-point that will be created for the Hamiltonian.
! Hk and Sk are returned for the user.
!
! The RemZConnection can be used to remove any matrix elements <i|H|j> where i and j
! are connections in the next unit cell in the Z-direction.
! This is usefull if one wishes to see the matrix as it would look while doing 
! transmission calculations.
!
! RemUCellDistances can be set to true so that inner-cell distances are removed.
! Several people on the SIESTA mailing list have "complained/asked questions" about 
! the difference in the Hamiltonians which are not always constructed similarly. 
! However, inner cell differences can be neglected. This is merely to get the same
! matrix as is created in Inelastica, for example.
!
! RemNFirstOrbitals can be used to fully remove that many states from the start of the
! Hamiltonian. This is used when there are buffer regions which should be disregarded.
! RemNLastOrbitals is the same, albeit in the end of the Hamiltonian.
!
! If the sizes of the incoming pointers Hk and Sk do not match the above
! they will be reallocated to the correct size (without notifying the user).
!
! Also there is a routine for obtaining an arbitrary transfer matrix.
! The interface is very much the same as that for the set_HS_matrix.
! However, here it has this interface:
!
!   call set_HS_transfermatrix(Gamma,ucell,na_u,no_u,no_s,maxnh, &
!       xij,numh,listhptr,listh,H,S, &
!       k,transfer_cell,HkT,SkT,xa,iaorb, &
!       RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
! 
! 'transfer_cell' is an integer vector of dimen 3 with each indices denoting the transfer matrix in that
! direction.
! For instance taking:
!   transfer_cell(:) = (/ 1 , 3 , 1 /)
! Will generate the transfer cell to the first neighbouring 'x', third neighbouring 'y' and first
! neighbouring 'z' cell.
! Often it can be necessary to have access to only the transfer matrix in one direction.
! In such cases you can supply:
!   transfer_cell(:) = (/ TRANSFER_ALL , TRANSFER_ALL , 2 /)
! The TRANSFER_ALL is a variable in this module and tells that it should
! use all possible transfer cells in that direction.
!
! In this conjunction the routine
!   call set_HS_available_transfers(Gamma,ucell,na_u,no_u,no_s,maxnh, &
!       xij,numh,listhptr,listh,xa,iaorb,transfer_cell)
! is invaluable. It returns in transfer_cell the allowed transfer cell units for
! the system. Notice that transfer_cell is a 2x3 matrix with the formatting like this:
!
!   transfer_cell(1,1) : minumum transfer cell in the x-direction
!   transfer_cell(2,1) : maximum transfer cell in the x-direction
!   transfer_cell(1,2) : minumum transfer cell in the y-direction
!   transfer_cell(2,2) : maximum transfer cell in the y-direction
!   transfer_cell(1,3) : minumum transfer cell in the z-direction
!   transfer_cell(2,3) : maximum transfer cell in the z-direction
!
! Usually the following holds (for obvious reasons):
!   transfer_cell(1,i) = -transfer_cell(2,i)
! This routine can thus be used to retrieve the allowed transfer cells before using
! set_HS_transfermatrix.
!
! Notice that RemZConnection is NOT available (it doesn't make sense).
! Furthermore, xa and iaorb are needed in order to determine the transfer cell.
! Thus they are no longer optional but mandatory arguments.
!
!
! Furthermore there are the routines:
!   call matrix_rem_left_right(no_tot,Hk,Sk,no_L,no_R)
!   and 
!   call matrix_symmetrize(no_tot,Hk,Sk,Ef)
! Which are used to remove connections of left/right regions
! and used to symmetrize and shift the Hamiltonian Ef, respectively.
! 
! NOTICE that a call to matrix_symmetrize is almost always needed!
! EVEN in the case of the transfer matrix.
  
  use precision, only : dp
 
  implicit none

  private

  interface set_HS_matrix
     module procedure set_HS_matrix_1d
     module procedure set_HS_matrix_2d
  end interface set_HS_matrix

  interface set_HS_transfermatrix
     module procedure set_HS_transfermatrix_1d
     module procedure set_HS_transfermatrix_2d
  end interface set_HS_transfermatrix

  interface matrix_rem_left_right
     module procedure matrix_rem_left_right_1d
     module procedure matrix_rem_left_right_2d
  end interface matrix_rem_left_right

  interface matrix_symmetrize
     module procedure matrix_symmetrize_1d
     module procedure matrix_symmetrize_2d
  end interface matrix_symmetrize

  public :: set_HS_matrix
  public :: set_HS_transfermatrix
  public :: matrix_rem_left_right
  public :: matrix_symmetrize
  public :: set_HS_available_transfers

  integer, parameter, public :: TRANSFER_ALL = -999999

contains

!*****************
! Setting the Hamiltonian for a specific k-point.
! It requires that the Node has the full sparse matrix available
!*****************
  subroutine set_HS_matrix_1d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,H,S, &
       k,Hk,Sk, &
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       xa,iaorb,lasto, &
       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use sys,       only : die 
    use alloc,     only : re_alloc
    use geom_helper, only : ucorb
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)           :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)          :: ucell(3,3) ! The unit cell of system
    integer, intent(in)           :: na_u ! Unit cell atoms
    integer, intent(in)           :: no_u ! Unit cell orbitals
    integer, intent(in)           :: no_s ! Supercell orbitals
    integer, intent(in)           :: maxnh ! Hamiltonian size
    real(dp), intent(in)          :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)           :: numh(no_u),listhptr(no_u)
    integer, intent(in)           :: listh(maxnh)
    real(dp), intent(in)          :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)          :: S(maxnh) ! Overlap
    real(dp), intent(in)          :: k(3) ! k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer          :: Hk(:), Sk(:)

! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    real(dp), intent(in),optional :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in), optional :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
    integer, intent(in), optional :: lasto(0:na_u) ! The number of orbitals on each atom (needed for RemUCellDistances)
    logical, intent(in), optional :: RemZConnection, RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3) ! The reciprocal unit cell
    real(dp) :: xo(3), xc
    real(dp), allocatable :: xuo(:)
    real(dp) :: kxij
    complex(dp) :: cphase
    integer :: no_tot
    integer, allocatable :: liaorb(:)
    integer :: i,j,iu,iuo,juo,iind,ind
    logical :: l_RemZConnection, l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_matrix")

    ! Option collecting
    l_RemZConnection = .false.
    if ( present(RemZConnection) ) &
         l_RemZConnection = RemZConnection
    if (l_RemZConnection .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &the z-connection.")
    if (l_RemZConnection .and. &
         (.not. present(iaorb) .and. .not. present(lasto) ) ) &
         call die("You need iaorb or lasto in set_HS_matrix when removing &
         &the z-connection.")
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances
    if (l_RemUCellDistances .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &unit cell distances.")
    if (l_RemUCellDistances .and. &
         (.not. present(iaorb) .and. .not. present(lasto) ) ) &
         call die("You need iaorb or lasto in set_HS_matrix when removing &
         &unit cell distances.")

    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)


    call re_alloc(Hk,1,no_tot*no_tot,name='Hk',routine='set_HS')
    call re_alloc(Sk,1,no_tot*no_tot,name='Sk',routine='set_HS')
    
    if ( l_RemZConnection ) then
       ! Prepare the cell to calculate the index of the atom
       call reclat(ucell,recell,0) ! Without 2*Pi
       
       ! Find the actual coordinates of the orbitals in the form of the sparse matrices
       ! Notice that this array is without the removed orbitals
       allocate(xuo(no_tot))
       call memory('A','D',no_tot,'set_HS')

       if ( present(iaorb) ) then
          do iuo = 1 , no_tot
             i = iaorb(iuo + l_RemNFirstOrbitals)
             xuo(iuo) = &
               xa(1,i) * recell(1,3) + &
               xa(2,i) * recell(2,3) + &
               xa(3,i) * recell(3,3)
          end do !io in uc
       else if ( present(lasto) ) then
          do j = 1 , na_u
             do i = lasto(j-1) + 1 , lasto(j)
                if ( i <= l_RemNFirstOrbitals ) cycle
                if ( no_tot < i - l_RemNFirstOrbitals ) cycle
                xuo(i-l_RemNFirstOrbitals) = &
                     xa(1,j) * recell(1,3) + &
                     xa(2,j) * recell(2,3) + &
                     xa(3,j) * recell(3,3)
             end do !io in uc
          end do
       end if
          
    end if
    
    ! Create the orb => atom index array
    if ( l_RemUCellDistances ) then
       allocate(liaorb(no_tot))
       call memory('A','I',no_tot,'set_HS')
       
       if ( present(iaorb) ) then
          do iuo = 1 , no_tot
             liaorb(iuo) = iaorb(iuo + l_RemNFirstOrbitals)
          end do !io in uc
       else if ( present(lasto) ) then
          ind = 0
          do j = 1 , na_u
             do i = lasto(j-1) + 1 , lasto(j)
                if ( i <= l_RemNFirstOrbitals ) cycle
                if ( no_tot < i - l_RemNFirstOrbitals ) cycle
                ind = ind + 1
                liaorb(ind) = j
             end do !io in uc
          end do
       end if
       
    end if

!
! Setup H,S for this k-point:
!
    do i = 1,no_tot*no_tot
       Hk(i) = dcmplx(0.d0,0.d0)
       Sk(i) = dcmplx(0.d0,0.d0)
    end do

    xo(:) = 0.0_dp

    setup_HS: if (.not.Gamma ) then

       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = ucorb(listh(ind),no_u) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             ! We also wish to remove the connection in
             ! in the inner cell
             ! I suspect we have a problem here
             ! We may need a check for being in the unitcell
             ! i.e.:
             ! if ( l_RemUCellDistances ) then
             !    if ( listh(ind) > no_u ) then ! As we check in the unit cell!
             !       xo(1) = xa(1,liaorb(juo)) &
             !            -  xa(1,liaorb(iuo))
             !       xo(2) = xa(2,liaorb(juo)) &
             !            -  xa(2,liaorb(iuo))
             !       xo(3) = xa(3,liaorb(juo)) &
             !            -  xa(3,liaorb(iuo))
             !    else
             !       xo = 0.0_dp
             !    end if
             ! end if                
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,liaorb(juo)) &
                     -  xa(1,liaorb(iuo))
                xo(2) = xa(2,liaorb(juo)) &
                     -  xa(2,liaorb(iuo))
                xo(3) = xa(3,liaorb(juo)) &
                     -  xa(3,liaorb(iuo))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,kxij))
             i = iuo+(juo-1)*no_tot
             Hk(i) = Hk(i)+H(ind)*cphase
             Sk(i) = Sk(i)+S(ind)*cphase
          end do
       end do

    else setup_HS
       ! It is not a Gamma calculation, thus we do not have any
       ! neighbouring cells etc.
       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             i = iuo+(juo-1)*no_tot
             Hk(i) = Hk(i)+H(ind)
             Sk(i) = Sk(i)+S(ind)
          end do
       end do

    end if setup_HS

    if ( l_RemUCellDistances ) then
       call memory('D','I',no_tot,'set_HS')
       deallocate(liaorb)
    end if

    if ( l_RemZConnection ) then
       call memory('D','D',no_tot,'set_HS')
       deallocate(xuo)
    end if

  end subroutine set_HS_matrix_1d

!*****************
! Setting the Hamiltonian for a specific k-point.
! It requires that the Node has the full sparse matrix available
!*****************
  subroutine set_HS_matrix_2d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,H,S, &
       k,Hk,Sk, &
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       xa,iaorb,lasto, &
       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use sys,       only : die 
    use alloc,     only : re_alloc
    use geom_helper, only : ucorb
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)           :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)          :: ucell(3,3) ! The unit cell of system
    integer, intent(in)           :: na_u ! Unit cell atoms
    integer, intent(in)           :: no_u ! Unit cell orbitals
    integer, intent(in)           :: no_s ! Total orbitals
    integer, intent(in)           :: maxnh ! Hamiltonian size
    real(dp), intent(in)          :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)           :: numh(no_u),listhptr(no_u)
    integer, intent(in)           :: listh(maxnh)
    real(dp), intent(in)          :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)          :: S(maxnh) ! Overlap
    real(dp), intent(in)          :: k(3) ! k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer          :: Hk(:,:), Sk(:,:)
! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    real(dp), intent(in),optional :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in), optional :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
    integer, intent(in), optional :: lasto(0:na_u) ! The number of orbitals on each atom (needed for RemUCellDistances)
    logical, intent(in), optional :: RemZConnection, RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3)
    real(dp) :: xo(3), xc
    real(dp), allocatable :: xuo(:)
    integer :: no_tot
    real(dp) :: kxij
    complex(dp) :: cphase
    integer, allocatable :: liaorb(:)
    integer :: i,j,iuo,iu,juo,iind,ind
    logical :: l_RemZConnection, l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_matrix")

    ! Option collecting
    l_RemZConnection = .false.
    if ( present(RemZConnection) ) &
         l_RemZConnection = RemZConnection
    if (l_RemZConnection .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &the z-connection.")
    if (l_RemZConnection .and. &
         (.not. present(iaorb) .and. .not. present(lasto) ) ) &
         call die("You need iaorb or lasto in set_HS_matrix when removing &
         &the z-connection.")
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances
    if (l_RemUCellDistances .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &unit cell distances.")
    if (l_RemUCellDistances .and. &
         (.not. present(iaorb) .and. .not. present(lasto) ) ) &
         call die("You need iaorb or lasto in set_HS_matrix when removing &
         &unit cell distances.")


    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)


    call re_alloc(Hk,1,no_tot,1,no_tot,name='Hk',routine='set_HS')
    call re_alloc(Sk,1,no_tot,1,no_tot,name='Sk',routine='set_HS')

    if ( l_RemZConnection ) then
       ! Prepare the cell to calculate the index of the atom
       call reclat(ucell,recell,0) ! Without 2*Pi
       
       ! Find the actual coordinates of the orbitals in the form of the sparse matrices
       ! Notice that this array is without the removed orbitals
       allocate(xuo(no_tot))
       call memory('A','D',no_tot,'set_HS')

       if ( present(iaorb) ) then
          do iuo = 1 , no_tot
             i = iaorb(iuo + l_RemNFirstOrbitals)
             xuo(iuo) = &
               xa(1,i) * recell(1,3) + &
               xa(2,i) * recell(2,3) + &
               xa(3,i) * recell(3,3)
          end do !io in uc
       else if ( present(lasto) ) then
          do j = 1 , na_u
             do i = lasto(j-1) + 1 , lasto(j)
                if ( i <= l_RemNFirstOrbitals ) cycle
                if ( no_tot < i - l_RemNFirstOrbitals ) cycle
                xuo(i-l_RemNFirstOrbitals) = &
                     xa(1,j) * recell(1,3) + &
                     xa(2,j) * recell(2,3) + &
                     xa(3,j) * recell(3,3)
             end do !io in uc
          end do
       end if

    end if

    if ( l_RemUCellDistances ) then
       allocate(liaorb(no_tot))
       call memory('A','I',no_tot,'set_HS')

       if ( present(iaorb) ) then
          do iuo = 1 , no_tot
             liaorb(iuo) = iaorb(iuo + l_RemNFirstOrbitals)
          end do !io in uc
       else if ( present(lasto) ) then
          ind = 0
          do j = 1 , na_u
             do i = lasto(j-1) + 1 , lasto(j)
                if ( i <= l_RemNFirstOrbitals ) cycle
                if ( no_tot < i - l_RemNFirstOrbitals ) cycle
                ind = ind + 1
                liaorb(ind) = j
             end do !io in uc
          end do
       end if

    end if


!
! Setup H,S for this k-point:
!
    do juo = 1,no_tot
       do iuo = 1,no_tot
          Hk(iuo,juo) = dcmplx(0.d0,0.d0)
          Sk(iuo,juo) = dcmplx(0.d0,0.d0)
       end do
    end do

    xo(:) = 0.0_dp

    setup_HS: if (.not.Gamma ) then

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = ucorb(listh(ind),no_u) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             ! We also wish to remove the connection in
             ! in the inner cell
             ! I suspect we have a problem here
             ! We may need a check for being in the unitcell
             ! i.e.:
             ! if ( l_RemUCellDistances ) then
             !    if ( listh(ind) > no_u ) then ! As we check in the unit cell!
             !       xo(1) = xa(1,liaorb(juo)) &
             !            -  xa(1,liaorb(iuo))
             !       xo(2) = xa(2,liaorb(juo)) &
             !            -  xa(2,liaorb(iuo))
             !       xo(3) = xa(3,liaorb(juo)) &
             !            -  xa(3,liaorb(iuo))
             !    else
             !       xo = 0.0_dp
             !    end if
             ! end if      
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,liaorb(juo)) &
                     -  xa(1,liaorb(iuo))
                xo(2) = xa(2,liaorb(juo)) &
                     -  xa(2,liaorb(iuo))
                xo(3) = xa(3,liaorb(juo)) &
                     -  xa(3,liaorb(iuo))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,kxij))
             Hk(iuo,juo) = Hk(iuo,juo)+H(ind)*cphase
             Sk(iuo,juo) = Sk(iuo,juo)+S(ind)*cphase
          end do
       end do

    else setup_HS

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             Hk(iuo,juo) = Hk(iuo,juo)+H(ind)
             Sk(iuo,juo) = Sk(iuo,juo)+S(ind)
          end do
       end do

    end if setup_HS

    if ( l_RemUCellDistances ) then
       call memory('D','I',no_tot,'set_HS')
       deallocate(liaorb)
    end if

    if ( l_RemZConnection ) then
       call memory('D','D',no_tot,'set_HS')
       deallocate(xuo)
    end if

  end subroutine set_HS_matrix_2d

!*****************
! Setting the Hamiltonian for a specific k-point as a specific transfer
! This can be used to obtain the transfer matrix for a certain Hamilton and S.
! It requires that the Node has the full sparse matrix available
!*****************
  subroutine set_HS_transfermatrix_1d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,H,S, &
       k,transfer_cell,HkT,SkT,xa,iaorb, &
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use sys,       only : die 
    use alloc,     only : re_alloc
    use geom_helper, only : ucorb
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)           :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)          :: ucell(3,3) ! The unit cell of system
    integer, intent(in)           :: na_u ! Unit cell atoms
    integer, intent(in)           :: no_u ! Unit cell orbitals
    integer, intent(in)           :: no_s ! Supercell orbitals
    integer, intent(in)           :: maxnh ! Hamiltonian size
    real(dp), intent(in)          :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)           :: numh(no_u),listhptr(no_u)
    integer, intent(in)           :: listh(maxnh)
    real(dp), intent(in)          :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)          :: S(maxnh) ! Overlap
    real(dp), intent(in)          :: k(3) ! k-point in [1/Bohr]
    integer, intent(in)           :: transfer_cell(3) ! The transfer cell directions
    real(dp), intent(in)          :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in)           :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer          :: HkT(:), SkT(:)

! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    logical, intent(in), optional :: RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3) ! The reciprocal unit cell
    real(dp) :: xo(3), xijo(3), xc
    real(dp) :: kxij
    complex(dp) :: cphase
    integer :: no_tot
    integer :: i,j,iu,iuo,juo,iind,ind
    logical :: l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_transfermatrix")

    ! Option collecting
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances

    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)

    call re_alloc(HkT,1,no_tot*no_tot,name='HkT',routine='set_HS')
    call re_alloc(SkT,1,no_tot*no_tot,name='SkT',routine='set_HS')

    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,recell,0) ! Without 2*Pi
       
!
! Setup H,S for this transfer k-point:
!
    do i = 1,no_tot*no_tot
       HkT(i) = dcmplx(0.d0,0.d0)
       SkT(i) = dcmplx(0.d0,0.d0)
    end do

    xo(:) = 0.0_dp

    setup_HST: if (.not.Gamma ) then

       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = ucorb(listh(ind),no_u) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! Determine the transfer matrix cell
             xijo(:) = xij(:,ind) - ( &
                  xa(:,iaorb(juo+l_RemNFirstOrbitals)) - xa(:,iaorb(iu)) &
                  )
             xc = sum(xijo(:) * recell(:,1))
             if ( nint(xc) /= transfer_cell(1) .and. transfer_cell(1) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,2))
             if ( nint(xc) /= transfer_cell(2) .and. transfer_cell(2) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,3))
             if ( nint(xc) /= transfer_cell(3) .and. transfer_cell(3) /= TRANSFER_ALL ) cycle

             ! We also wish to remove the connection in
             ! in the inner cell
             ! I suspect we have a problem here
             ! We may need a check for being in the unitcell
             ! i.e.:
             ! if ( l_RemUCellDistances ) then
             !    if ( listh(ind) > no_u ) then ! As we check in the unit cell!
             !       xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(1,iaorb(iu))
             !       xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(2,iaorb(iu))
             !       xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(3,iaorb(iu))
             !    else
             !       xo = 0.0_dp
             !    end if
             ! end if      
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(1,iaorb(iu))
                xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(2,iaorb(iu))
                xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(3,iaorb(iu))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,kxij))
             i = iuo+(juo-1)*no_tot
             HkT(i) = HkT(i)+H(ind)*cphase
             SkT(i) = SkT(i)+S(ind)*cphase
          end do
       end do

    else setup_HST
       ! It is not a Gamma calculation, thus we do not have any
       ! neighbouring cells etc.
       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! Determine the transfer matrix cell
             xijo(:) = xij(:,ind) - ( &
                  xa(:,iaorb(juo+l_RemNFirstOrbitals)) - xa(:,iaorb(iu)) &
                  )
             xc = sum(xijo(:) * recell(:,1))
             if ( nint(xc) /= transfer_cell(1) .and. transfer_cell(1) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,2))
             if ( nint(xc) /= transfer_cell(2) .and. transfer_cell(2) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,3))
             if ( nint(xc) /= transfer_cell(3) .and. transfer_cell(3) /= TRANSFER_ALL ) cycle

             i = iuo+(juo-1)*no_tot
             HkT(i) = HkT(i)+H(ind)
             SkT(i) = SkT(i)+S(ind)
          end do
       end do

    end if setup_HST

  end subroutine set_HS_transfermatrix_1d


  subroutine set_HS_transfermatrix_2d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,H,S, &
       k,transfer_cell,HkT,SkT,xa,iaorb,&
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use sys,       only : die 
    use alloc,     only : re_alloc
    use geom_helper, only : ucorb
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)           :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)          :: ucell(3,3) ! The unit cell of system
    integer, intent(in)           :: na_u ! Unit cell atoms
    integer, intent(in)           :: no_u ! Unit cell orbitals
    integer, intent(in)           :: no_s ! Total orbitals
    integer, intent(in)           :: maxnh ! Hamiltonian size
    real(dp), intent(in)          :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)           :: numh(no_u),listhptr(no_u)
    integer, intent(in)           :: listh(maxnh)
    real(dp), intent(in)          :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)          :: S(maxnh) ! Overlap
    real(dp), intent(in)          :: k(3) ! k-point in [1/Bohr]
    integer, intent(in)           :: transfer_cell(3) ! The transfer cell directions
    real(dp), intent(in)          :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in)           :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer          :: HkT(:,:), SkT(:,:)
! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    logical, intent(in), optional :: RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3)
    real(dp) :: xo(3), xijo(3), xc
    integer :: no_tot
    real(dp) :: kxij
    complex(dp) :: cphase
    integer :: j,iuo,iu,juo,iind,ind
    logical :: l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_transfermatrix")

    ! Option collecting
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances

    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)


    call re_alloc(HkT,1,no_tot,1,no_tot,name='HkT',routine='set_HS')
    call re_alloc(SkT,1,no_tot,1,no_tot,name='SkT',routine='set_HS')


    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,recell,0) ! Without 2*Pi
    
    ! Setup H,S for this k-point:
    do juo = 1,no_tot
       do iuo = 1,no_tot
          HkT(iuo,juo) = dcmplx(0.d0,0.d0)
          SkT(iuo,juo) = dcmplx(0.d0,0.d0)
       end do
    end do

    xo(:) = 0.0_dp

    setup_HST: if (.not.Gamma ) then

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = ucorb(listh(ind),no_u) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! Determine the transfer matrix cell
             xijo(:) = xij(:,ind) - ( &
                  xa(:,iaorb(juo+l_RemNFirstOrbitals)) - xa(:,iaorb(iu)) &
                  )
             xc = sum(xijo(:) * recell(:,1))
             if ( nint(xc) /= transfer_cell(1) .and. transfer_cell(1) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,2))
             if ( nint(xc) /= transfer_cell(2) .and. transfer_cell(2) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,3))
             if ( nint(xc) /= transfer_cell(3) .and. transfer_cell(3) /= TRANSFER_ALL ) cycle

             ! We also wish to remove the connection in
             ! in the inner cell
             ! I suspect we have a problem here
             ! We may need a check for being in the unitcell
             ! i.e.:
             ! if ( l_RemUCellDistances ) then
             !    if ( listh(ind) > no_u ) then ! As we check in the unit cell!
             !       xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(1,iaorb(iu))
             !       xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(2,iaorb(iu))
             !       xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
             !            -  xa(3,iaorb(iu))
             !    else
             !       xo = 0.0_dp
             !    end if
             ! end if      
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(1,iaorb(iu))
                xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(2,iaorb(iu))
                xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(3,iaorb(iu))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,kxij))
             HkT(iuo,juo) = HkT(iuo,juo)+H(ind)*cphase
             SkT(iuo,juo) = SkT(iuo,juo)+S(ind)*cphase
          end do
       end do

    else setup_HST

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! Determine the transfer matrix cell
             xijo(:) = xij(:,ind) - ( &
                  xa(:,iaorb(juo+l_RemNFirstOrbitals)) - xa(:,iaorb(iu)) &
                  )
             xc = sum(xijo(:) * recell(:,1))
             if ( nint(xc) /= transfer_cell(1) .and. transfer_cell(1) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,2))
             if ( nint(xc) /= transfer_cell(2) .and. transfer_cell(2) /= TRANSFER_ALL ) cycle
             xc = sum(xijo(:) * recell(:,3))
             if ( nint(xc) /= transfer_cell(3) .and. transfer_cell(3) /= TRANSFER_ALL ) cycle

             HkT(iuo,juo) = HkT(iuo,juo)+H(ind)
             SkT(iuo,juo) = SkT(iuo,juo)+S(ind)
          end do
       end do

    end if setup_HST
    
  end subroutine set_HS_transfermatrix_2d

  subroutine set_HS_available_transfers(ucell,na_u,xa,lasto,no_u,maxnh, &
       xij,numh,listhptr,listh,transfer_cell)
    use geom_helper, only : ucorb, iaorb
    use cellSubs, only : reclat

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3) ! The unit cell of system
    integer, intent(in)  :: na_u ! Unit cell atoms
    real(dp), intent(in) :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in)  :: lasto(0:na_u) ! last orbital number of equivalent atom
    integer, intent(in)  :: no_u ! Unit cell orbitals
    integer, intent(in)  :: maxnh ! Hamiltonian size
    real(dp), intent(in) :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)  :: numh(no_u),listhptr(no_u)
    integer, intent(in)  :: listh(maxnh)
! ***********************
! * OUTPUT variables    *
! ***********************
    integer, intent(out) :: transfer_cell(2,3)

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3)
    real(dp) :: xijo(3), xc
    integer :: ia, ja
    integer :: i,j,iuo,juo,ind

    ! Initialize the transfer cell to:
    transfer_cell(:,:) = 0

    ! Prepare the cell to calculate the index of the atom
    call reclat(ucell,recell,0) ! Without 2*Pi
    
    do iuo = 1 , no_u
       ia = iaorb(iuo,lasto)
       do j = 1 , numh(iuo)
          ind = listhptr(iuo) + j
          juo = ucorb(listh(ind),no_u)
          ja = iaorb(juo,lasto)
          xijo(:) = xij(:,ind) - ( xa(:,ja)-xa(:,ia) )
          ! Loop over directions
          do i = 1 , 3 
             ! recell is already without 2*Pi
             xc = sum(xijo(:) * recell(:,i))
             transfer_cell(1,i) = min(transfer_cell(1,i),nint(xc))
             transfer_cell(2,i) = max(transfer_cell(2,i),nint(xc))
          end do
       end do
    end do
    
  end subroutine set_HS_available_transfers

  ! Routine for removing left right overlaps of certain regions.
  ! Is used to fully remove the connection between left and right
  ! states in the Hamiltonian
  subroutine matrix_rem_left_right_2d(no_tot,Hk,Sk,no_L,no_R)

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in)        :: no_tot, no_L, no_R

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot,no_tot), Sk(no_tot,no_tot)

! **************************
! * LOCAL variables        *
! **************************
    integer :: i,j

    ! If nothing is to be removed return immidiately...
    if ( no_L == 0 .or. no_R == 0 ) return

    do j = no_tot - no_R + 1 , no_tot
       do i = 1 , no_L
          Hk(i,j) = dcmplx(0.d0,0.d0)
          Sk(i,j) = dcmplx(0.d0,0.d0)
          Hk(j,i) = dcmplx(0.d0,0.d0)
          Sk(j,i) = dcmplx(0.d0,0.d0)
       end do
    end do

  end subroutine matrix_rem_left_right_2d

  subroutine matrix_rem_left_right_1d(no_tot,Hk,Sk,no_L,no_R)
    integer, intent(in)        :: no_tot, no_L, no_R
    complex(dp), intent(inout) :: Hk(no_tot*no_tot), Sk(no_tot*no_tot)
    call matrix_rem_left_right_2d(no_tot,Hk,Sk,no_L,no_R)
  end subroutine matrix_rem_left_right_1d


  ! Routine for symmetrizing and shifting the matrix Ef
  subroutine matrix_symmetrize_2d(no_tot,Hk,Sk,Ef)

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in)        :: no_tot
    real(dp), intent(in)       :: Ef

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot,no_tot), Sk(no_tot,no_tot)
    
! **************************
! * LOCAL variables        *
! **************************
    integer :: iuo,juo

    do iuo = 1,no_tot
       do juo = 1,iuo-1
          
          Sk(juo,iuo) = 0.5d0*( Sk(juo,iuo) + dconjg(Sk(iuo,juo)) )
          Sk(iuo,juo) =  dconjg(Sk(juo,iuo))
          
          Hk(juo,iuo) = 0.5d0*( Hk(juo,iuo) + dconjg(Hk(iuo,juo)) ) &
               - Ef*Sk(juo,iuo)
          Hk(iuo,juo) =  dconjg(Hk(juo,iuo))
          
       end do
       
       Sk(iuo,iuo)=Sk(iuo,iuo) - dcmplx(0d0,dimag(Sk(iuo,iuo)) )
       
       Hk(iuo,iuo)=Hk(iuo,iuo) - dcmplx(0d0,dimag(Hk(iuo,iuo)) ) &
            - Ef*Sk(iuo,iuo) 
    end do

  end subroutine matrix_symmetrize_2d

  subroutine matrix_symmetrize_1d(no_tot,Hk,Sk,Ef)
    integer, intent(in)        :: no_tot
    real(dp), intent(in)       :: Ef
    complex(dp), intent(inout) :: Hk(no_tot*no_tot), Sk(no_tot*no_tot)
    call matrix_symmetrize_2d(no_tot,Hk(1),Sk(1),Ef)
  end subroutine matrix_symmetrize_1d

end module m_hs_matrix
  
