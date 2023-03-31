! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication

module m_tbt_proj

  ! We here implement a feature to calculate a projection
  ! The equations are the following:

  ! \lambda_i, v_i = diag(H\lambda = S)
  ! , where \lambda_i, v_i is the i'th eigen-value, eigen-vector pair.

  ! To calculate the projection states we do:
  ! p_i = S^{1/2} * v_i / sqrt(S^{1/2} * v_i,(S^{1/2} * v_i)^*)
  ! ( note: ^* is conjugate, not dagger)
  ! The projector is then:
  !   P_i = p_i * p_i^\dagger 
  ! i.e. we get a full matrix by doing the outer product of these vectors

  ! In the following we designate a projection by the name
  ! of "molecule".

  use precision, only : dp
  use m_region
  use m_ts_electype
  use m_tbt_save, only : tNodeE, save_parallel

#ifdef NCDF_4
#ifdef MPI
  use m_tbt_save, only : save_attach_buffer
#endif
  use m_tbt_save, only : local_save_DOS
  use m_tbt_save, only : tbt_cdf_precision
  use netcdf_ncdf, only : NF90_MAX_NAME
#endif

  implicit none

  private

#ifdef NCDF_4
  ! Default to not do any kind of projections

  ! An example input format for projection.
  ! %block TBT.Proj
  !   mol-1
  !   mol-2
  !  ...
  !   etc
  ! %endblock TBT.Proj
    
  ! Then we create all those projection files
  ! %block TBT.Proj.mol-1
  !   # Several lists of atoms "concatenates" the regions.
  !   atom from <> to <>
  !   atom from <> to <>
  !   proj <Name>
  !     # We note the Fermi-level as the 0th level. 
  !     #   (for different biases, the level alignment might change :( )
  !     # Several lines will concatenate.
  !     # So the below lines will only be the HOMO-1, HOMO, LUMO and LUMO+1
  !     level from -2 to -1 # HOMO-1 and HOMO
  !     level from 1 to 2   # LUMO and LUMO+1
  !     level from -2 to 2  # HOMO-1, HOMO, LUMO and LUMO+1
  !     level 1             # LUMO 
  !   end [proj|]
  !   ... more proj <Name> blocks, if needed
  ! %endblock TBT.Proj.mol-1

  integer, parameter :: PROJ_NAME_LEN = 32

  type :: tProjMol
    ! Name of the "molecule"
    character(len=PROJ_NAME_LEN) :: name
    ! Whether this projection is k-resolved
    logical :: Gamma = .true., DOS = .true.
    ! Region containing the atoms, orbitals
    type(tRgn) :: atom, orb, pvt
    ! An array of used projections for this "molecule"
    ! The region contains the levels that constitute a single projection
    type(tRgn), allocatable :: proj(:)
    ! The actual projection state for a single level
    ! Instead of saving N times (orb%n,orb%n) we need only saving N times orb%n vectors
    type(tRgn) :: lvls
    complex(dp), allocatable :: p(:,:)
  end type tProjMol
  type(tProjMol), allocatable :: mols(:)
  integer, save :: N_mol = 0
  public :: N_mol, mols

  ! Contains the Gamma projections of an electrode
  type :: tProjMolEl
    ! The molecule that needs projection.
    type(tProjMol), pointer :: mol => null()
    ! The electrode that will be projected onto
    type(Elec), pointer :: El => null()
    ! All <|Gamma|> projections
    ! this contains all <i|Gamma|j>:
    complex(dp), allocatable :: bGk(:,:)
  end type tProjMolEl
  integer, save :: N_proj_ME
  type(tProjMolEl), allocatable, target :: proj_ME(:)
  public :: N_proj_ME, proj_ME

  ! When ever we do a projection we will also calculate
  ! to the full projection
  type :: tLvlMolEl
    type(tProjMolEl), pointer :: ME => null()
    ! An index below zero corresponds to that electrode
    ! Is used, i.e. no projection.
    ! idx = -2 => Elecs(2)
    integer :: idx = 0
  end type tLvlMolEl
  public :: tLvlMolEl

  type :: tProjT
    ! The projector "left" electrode
    type(tLvlMolEl) :: L
    ! The projector "right" electrode(s)
    type(tLvlMolEl), allocatable :: R(:)
  end type tProjT
  integer, save :: N_proj_T
  type(tProjT), allocatable :: proj_T(:)
  public :: N_proj_T, proj_T

  public :: init_proj
  public :: init_proj_T
  public :: open_cdf_proj
  public :: proj_print
  public :: proj_LME_assoc
  public :: init_proj_save
  public :: proj_update
  public :: proj_Mt_mix, proj_bMtk
  public :: proj_cdf_save
  public :: proj_cdf_save_sp_dev
  public :: proj_cdf_save_bgammak
  public :: proj_cdf_save_S_D
  public :: proj_cdf2ascii

#endif

  public :: read_proj_options
  public :: print_proj_options

contains

#ifdef NCDF_4

  subroutine proj_LME_assoc(lhs,rhs)
    type(tLvlMolEl), pointer :: lhs
    type(tLvlMolEl), intent(inout), target :: rhs
    lhs => rhs
  end subroutine proj_LME_assoc
  
  subroutine init_proj( na_u , lasto , a_Dev , o_Dev, save_DATA )
    
    use fdf
    use fdf_extra

    use dictionary

    integer, intent(in) :: na_u
    integer, intent(in) :: lasto(0:na_u)
    type(tRgn), intent(in) :: a_Dev, o_Dev
    type(dictionary_t), intent(inout) :: save_DATA

    ! Local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    integer :: N_proj, N, n_orb
    integer :: im, ip
    character(len=PROJ_NAME_LEN) :: char, name
    type(tRgn) :: r_tmp, r_tmp2

    ! This routine assumes the projection names are already read in.
    if ( N_mol == 0 ) return

    ! If the user has requested only to calculate the 
    ! self-energies, we should not read in the projections.
    if ( ('Sigma-only'.in.save_DATA) ) return

    ! Read in the data
    do im = 1 , N_mol

      name = mols(im)%name
      
      ! Open the block of the molecule
      if ( .not. fdf_block('TBT.Proj.'//trim(name),bfdf) ) then
        call die('Projection TBT.Proj.'//trim(name)//' &
            &has not been defined.')
      end if
      
      ! Reset number of projections for this molecule.
      N_proj = 0
      N      = 0

      ! We read it line by line and count the number of projections
      do while ( fdf_bline(bfdf,pline) ) 
        if ( fdf_bnnames(pline) == 0 ) cycle
        char = fdf_bnames(pline,1)
        if ( leqi(char,'proj') ) then
          ! We have a projection
          N_proj = N_proj + 1
        else if ( leqi(char,'end') ) then
          ! easy way to check that they are also finished
          N = N - 1
          ! If they are not aligned
          ! then they must be malplaced
          if ( - N /= N_proj ) then
            call die('Error in projection block, &
                &a projection has prematurely ended, correct input.')
          end if

        else if ( leqi(char,'atom') .or. leqi(char,'position') ) then

          ! Add the atom to the list
          call fdf_brange(pline,r_tmp, 1, na_u)
          call rgn_copy(mols(im)%atom,r_tmp2)
          call rgn_union(r_tmp2,r_tmp,mols(im)%atom)

        else if ( leqi(char,'Gamma') ) then

          ! Get whether the molecule should be treated in a
          ! Gamma-consideration
          mols(im)%Gamma = fdf_bboolean(pline,1,after=1)

        else if ( leqi(char,'DOS') .and. N_proj == N ) then

          ! The N_proj == N assures that we are not 
          ! inside a proj block

          mols(im)%DOS = fdf_bboolean(pline,1,after=1)

        end if
      end do
      call rgn_delete(r_tmp,r_tmp2)
      if ( N_proj == 0 ) then
        call die('Error in projection block, &
            &you have not specified any projections.')
      end if
      if ( - N /= N_proj ) then
        call die('Error in projection block, &
            &all projections has not been ended correctly.')
      end if
      if ( mols(im)%atom%n == 0 ) then
        call die('Number of atoms in projection is zero, this is not &
            &allowed.')
      end if
      
      ! Convert the molecule to orbitals
      call rgn_sort(mols(im)%atom)
      mols(im)%atom%name = 'Atoms'
      ! We check that the projection is completely contained
      ! in the device region, if not then a projection of 
      ! eigenstates does not really make sense...
      call rgn_union(a_Dev,mols(im)%atom,r_tmp)
      if ( r_tmp%n /= a_Dev%n ) then
        call die('Projection not fully contained in the device &
            &region, please correct input')
      end if
      call rgn_delete(r_tmp)

      call rgn_Atom2Orb(mols(im)%atom,na_u,lasto,mols(im)%orb)
      mols(im)%orb%sorted = .false.
      mols(im)%orb%name = '[O] projection: '//trim(mols(im)%name)

      ! Retain the number of orbitals
      n_orb = mols(im)%orb%n

      ! Re-read, and then we read in the projections... :)
      call fdf_brewind(bfdf)

      ! Allocate all projections and read them in...
      allocate(mols(im)%proj(N_proj))

      ! Initialize levels 
      call rgn_delete(mols(im)%lvls)

      ip = 0
      do while ( fdf_bline(bfdf,pline) ) 
        if ( fdf_bnnames(pline) == 0 ) cycle
        char = fdf_bnames(pline,1)
        if ( leqi(char,'proj') ) then
          ! We have a projection
          if ( fdf_bnnames(pline) < 2 ) then
            call die('Name for projection has not been provided.')
          end if

          ip = ip + 1

          ! Save the name
          name = fdf_bnames(pline,2)
          if ( index(name,'.') > 0 ) then
            call die('Projections cannot be named with .!')
          end if

          ! Read in the projection
          do while ( fdf_bline(bfdf,pline) ) 

            if ( fdf_bnnames(pline) == 0 ) cycle
            char = fdf_bnames(pline,1)
            if ( leqi(char,'end') ) then
              ! We have ended this level, exit so that we can 
              ! read the next one
              exit

            else if ( leqi(char,'level') .or. &
                leqi(char,'lvl') ) then

              ! We do not know whether all levels
              ! are below/above the Fermi level.
              call fdf_brange(pline,r_tmp,-n_orb,n_orb)
              call rgn_copy(mols(im)%proj(ip),r_tmp2)
              call rgn_union(r_tmp2,r_tmp,mols(im)%proj(ip))

            else if ( leqi(char,'clear') .or. &
                leqi(char,'clear-level') .or. &
                leqi(char,'clear-lvl') ) then

              ! We do not know whether all levels
              ! are below/above the Fermi level.
              call fdf_brange(pline,r_tmp,-n_orb,n_orb)
              call rgn_copy(mols(im)%proj(ip),r_tmp2)
              call rgn_complement(r_tmp,r_tmp2, mols(im)%proj(ip))

            end if

          end do

          ! Remove 0 (Ef) from the levels
          call rgn_init(r_tmp,1,val=0)
          call rgn_complement(r_tmp,mols(im)%proj(ip), &
              mols(im)%proj(ip))
          call rgn_delete(r_tmp,r_tmp2)

          ! A projection has to have at least one projection state
          if ( mols(im)%proj(ip)%n < 1 ) then
            call die('A projection went wrong, it MUST have at least &
                &one state assigned, the Fermi-level (0) is not a &
                &well-defined state.')
          end if

          ! Sort the levels
          call rgn_sort(mols(im)%proj(ip))

          ! Save the name to the levels
          mols(im)%proj(ip)%name = name

          ! Create union of all levels 
          call rgn_union(mols(im)%lvls,mols(im)%proj(ip),r_tmp)
          call rgn_copy(r_tmp,mols(im)%lvls)

#ifdef TBT_PHONON
          if ( any(mols(im)%lvls%r < 1) ) then
            write(*,*)'Error in specifying phonon eigenstates.'
            write(*,*)'There is no such thing as a negative indexed &
                &phonon eigenstate.'
            call die('Phonon molecular levels *must* start from 1')
          end if
#endif
          
        end if
        
      end do ! loop on reading in projections
      
      ! We have now completed reading in all projections.
      
      ! Allocate |> for all unique levels
      allocate(mols(im)%p(n_orb,mols(im)%lvls%n))
      ! Sort the levels
      call rgn_sort(mols(im)%lvls)

      ! Take all projections and make them index based
      do ip = 1 , N_proj

        call rgn_copy(mols(im)%proj(ip),r_tmp)

        mols(im)%proj(ip)%r(:) = rgn_pivot(mols(im)%lvls,r_tmp%r(:))
        ! * NOTE * 
        ! Projections are sorted, lvls are sorted, hence
        ! the pivoting scheme *must* is intrinsically sorted.

      end do

    end do
    
    call rgn_delete(r_tmp,r_tmp2)
    
  end subroutine init_proj
  

  subroutine open_cdf_proj(fname, ncdf)
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif

    character(len=*), intent(in) :: fname
    type(hNCDF), intent(inout) :: ncdf

    ! quick return
    if ( N_proj_ME == 0 ) return

#ifdef MPI
    ! Open the netcdf file
    if ( save_parallel ) then
      call ncdf_open(ncdf,fname, mode=ior(NF90_WRITE,NF90_MPIIO), &
          comm = MPI_COMM_WORLD )
    else
      call ncdf_open(ncdf,fname, mode=NF90_WRITE)
    end if
#else
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#endif
    
  end subroutine open_cdf_proj

  subroutine proj_print( N_Elec, Elecs )

    use parallel, only : Node
    use fdf, only : leqi
    
    integer, intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)
    
    integer :: im, ip
    integer :: it, iE
    integer :: ipt, idx
    logical :: first_P
    type(Elec), pointer :: El_L, El
    type(tProjMol), pointer :: mol
    type(tRgn) :: r_tmp
    logical, allocatable :: checked(:,:)
    
    if ( Node /= 0 ) return
    if ( N_mol == 0 ) return

    ! Lets first print out the projections

    write(*,'(/,a)')'tbt: Projection regions:'

    do im = 1 , N_mol
      
      ! A molecule MUST not be named the same 
      ! as an electrode!
      do iE = 1 , N_Elec
        if ( leqi(mols(im)%name,Elecs(iE)%name) ) then
          call die('Projections and electrodes &
              &must NOT be named the same. Please differ names!')
        end if
      end do

      write(*,'(2a)')' - Projection ',trim(mols(im)%name)
      if ( mols(im)%Gamma ) then
        write(*,'(a)')'   Gamma projection: True'
      else
        write(*,'(a)')'   Gamma projection: False'
      end if
      ! Currently the Gf DOS projection is not implemented,
      ! hence, we do not write it out...
      !if ( mols(im)%DOS ) then
      !   write(*,'(a)')'   DOS projection: True'
      !else
      !   write(*,'(a)')'   DOS projection: False'
      !end if
      call rgn_print(mols(im)%atom, seq_max = 8 , indent = 3)
      write(*,'(a)') '  * Different projections:'
      do ip = 1 , size(mols(im)%proj)

        call rgn_copy(mols(im)%proj(ip),r_tmp)
        do iE = 1 , r_tmp%n
          r_tmp%r(iE) = mols(im)%lvls%r(r_tmp%r(iE))
        end do

        call rgn_print(r_tmp, seq_max = 12 , indent = 5)
        call rgn_delete(r_tmp)

        do iE = 1 , N_Elec
          if ( leqi(mols(im)%proj(ip)%name,Elecs(iE)%name) ) then
            call die('Projections and electrodes &
                &must NOT be named the same. Please differ names!')
          end if
        end do

        ! Figure out if the permutation exists...
        do it = 1 , N_proj_T
          do iE = 1 , N_Elec
            if ( ProjMolEl_same(proj_T(it)%L,Elecs(iE),mols(im),ip) ) then
              write(*,'(tr6,2a,''.T -> ['')',advance='no') '- ',trim(Elecs(iE)%name)

              ! We loop the RHS here
              if ( .not. allocated(proj_T(it)%R) ) then
                call die('Error in setting up RHS projections')
              end if
              do ipt = 1 , size(proj_T(it)%R)
                if ( ipt > 1 ) write(*,'(a)',advance='no') ','
                if ( mod(ipt,4) == 0 ) write(*,'(/,tr10)',advance='no')
                if ( proj_T(it)%R(ipt)%idx > 0 ) then
                  El => proj_T(it)%R(ipt)%ME%El
                  mol => proj_T(it)%R(ipt)%ME%mol
                  idx = proj_T(it)%R(ipt)%idx 
                  write(*,'(tr1,2(a,''.''),a)',advance='no') trim(El%name), &
                      trim(mol%name), trim(mol%proj(idx)%name)
                else
                  write(*,'(tr1,a)',advance='no') &
                      trim(Elecs(-proj_T(it)%R(ipt)%idx)%name)
                end if
                
              end do
              write(*,'(a)') ']'
            end if
          end do
        end do
      end do
      
      write(*,*) ! new-line
      
    end do
    
    ! Print all pristine, to right projections
    first_P = .true.
    do it = 1 , N_proj_T
      if ( proj_T(it)%L%idx < 0 ) then
        El_L => Elecs(-proj_T(it)%L%idx)
        if ( first_P ) then
          write(*,'(a)') '* Non-projected incoming scattering state projections:'
          first_P = .false.
        end if
          
        ! We have a single L-projection
        write(*,'(tr2,2a,''.T -> ['')',advance='no') '- ', trim(El_L%name)

        ! We loop the RHS here
        if ( size(proj_T(it)%R) > 0 ) then
          do ipt = 1 , size(proj_T(it)%R)
            El => proj_T(it)%R(ipt)%ME%El
            if ( ipt > 1 ) write(*,'(a)',advance='no') ','
            if ( mod(ipt,4) == 0 ) write(*,'(/,tr7)',advance='no')
            if ( proj_T(it)%R(ipt)%idx > 0 ) then
              mol => proj_T(it)%R(ipt)%ME%mol
              idx = proj_T(it)%R(ipt)%idx 
              write(*,'(tr1,a)',advance='no') &
                  trim(proj_ME_name(proj_T(it)%R(ipt)))
            else
              write(*,'(2a)')'Pure transport from: '//trim(El_L%name)//' to ',&
                  trim(El%name)
              call die('Erroneous projection, this is NOT a projection')
            end if
          end do
          write(*,'(a)') ']'
        else
          write(*,'(a)')'Pure transport from: '//trim(El_L%name)//' to nothing.'
          call die('Erroneous projection, this is NOT a projection')
        end if
      end if
    end do

    write(*,*) ! New-line

    ! Here we check that all projections match the electrode
    ! they are projected on
    ! A transport projection on molecular states
    ! ONLY has meaning if the scattering states
    ! lives fully on the molecular states.
    ! Hence the union of the molecule and scattering state
    ! MUST equal that of the molecule
    allocate(checked(N_mol,N_Elec))
    checked(:,:) = .false.
    do it = 1 , N_proj_T
      do iE = 1 , N_Elec
        if ( all(checked(:,iE)) ) cycle
        do im = 1 , N_mol
          if ( checked(im,iE) ) cycle
          ! We can only check a projection which
          ! actually projects...
          if ( proj_T(it)%L%idx > 0 ) then
            ! if the electrode is the same
            if ( ProjMolEl_same(proj_T(it)%L, &
                Elecs(iE),mols(im),proj_T(it)%L%idx) ) then
              ! The Left is the same
              call rgn_union(mols(im)%orb,Elecs(iE)%o_inD,r_tmp)
              if ( r_tmp%n /= mols(im)%orb%n ) then
                ! The overlap region does not fully co-incide
                ! with the molecule region.
                ! Hence the projection is erroneuos
                call print_proj(proj_T(it)%L)
                call rgn_print(mols(im)%orb)
                call rgn_print(Elecs(iE)%o_inD)
                write(*, '(a)') 'The selected projection region does not encapsulate &
                    &the electrodes device region.'
                write(*, '(a)') 'Please select a TBT.Atoms.Device region such that &
                    &the second region is fully encapsulated in the first region.'
                call die('The scattering states are not fully &
                    &encapsulated on a LEFT projection, please change &
                    &TBT.Atoms.Device accordingly. This is not allowed.')
              end if
              checked(im,iE) = .true.
            end if
          end if

          do ip = 1 , size(proj_T(it)%R)
            if ( proj_T(it)%R(ip)%idx <= 0 ) cycle
            ! if the electrode is the same
            if ( ProjMolEl_same(proj_T(it)%R(ip), &
                Elecs(iE),mols(im),proj_T(it)%R(ip)%idx) ) then
              call rgn_union(mols(im)%orb,Elecs(iE)%o_inD,r_tmp)
              if ( r_tmp%n /= mols(im)%orb%n ) then
                ! The overlap region does not fully co-incide
                ! with the molecule region.
                ! Hence the projection is erroneuos
                call print_proj(proj_T(it)%R(ip))
                call rgn_print(mols(im)%orb)
                call rgn_print(Elecs(iE)%o_inD)
                write(*, '(a)') 'The selected projection region does not encapsulate &
                    &the electrodes device region.'
                write(*, '(a)') 'Please select a TBT.Atoms.Device region such that &
                    &the second region is fully encapsulated in the first region.'
                call die('The scattering states are not fully &
                    &encapsulated on a RIGHT projection, please change &
                    &TBT.Atoms.Device accordingly. This is not allowed.')
              end if
            end if
            checked(im,iE) = .true.
          end do
        end do
      end do
    end do

    call rgn_delete(r_tmp)
    
  end subroutine proj_print

  subroutine init_proj_T( N_Elec, Elecs , save_DATA )

    use fdf
    use parallel, only : Node
    use dictionary

    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(dictionary_t), intent(inout) :: save_DATA

    character(len=100) :: char
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    integer :: ipt, ip, iE_p, im_p, ip_p, iE_c, im_c, ip_c, it
    integer :: i, j, lp1, lp2

    ! Number of different possible molecule projections
    ! in reality the maximum should be (N_mol * N_Elec)
    integer, parameter :: N_ME = 1000
    type(tProjMolEl), target :: tmp_ME(N_ME)

    ! Local reading block to accomodate several permutations
    ! of the projections
    ! Currently we "only" allow 5000 permutations
    ! of projections.
    integer, parameter :: N_LME = 5000
    type(tLvlMolEl) :: tmp_LME(N_LME)
    integer :: ilme, ilme_o
    integer :: itmp
    logical :: all_T, reflect, any_skipped
    logical :: checked

    ! number of different LHS projections
    N_proj_T = 0
    ! Number of different molecule-electrode projections
    N_proj_ME = 0
    ! if there are no molecules, we return immediately
    if ( N_mol == 0 ) return

    ! If the projection block exists
    ! it means the user is requesting projections.
    if ( .not. fdf_defined('TBT.Projs.T') ) return

    if ( .not. fdf_block('TBT.Projs.T',bfdf) ) then
      call die('TBT.Projs.T is not a block, please correct')
    end if

    ! Whether we should calculate all T?
    all_T = ('proj-T-all'.in.save_DATA)
    reflect = ('proj-T-sum-out'.in.save_DATA)
    any_skipped = .false.
    ! Now all_T determines which projections we will utilize
    ! The regular electrode loop runs:
    ! do iEl = 1 , N_Elec
    !   do jEl = 1 , N_Elec
    !     if ( .not. all_T .and. jEl <= iEl ) cycle
    !   end do
    ! end do
    ! Hence for cases where not all are calculated
    ! Wheter we should reject any electrodes in the 

    ! First read in number of T-projections
    ! This reading is necessary to pre-allocate
    ! the actual required size of the projections used.
    ! I.e. we calculate the number of different LHS projections
    do while ( fdf_bline(bfdf,pline) )
      
      ! skip empty line
      if ( fdf_bnnames(pline) == 0 ) cycle
      char = fdf_bnames(pline,1)
      if ( .not. leqi(char,'from') ) cycle
      if ( fdf_bnnames(pline) < 2 ) then
        call die('Error in TBT.Projs.T block, from <projection> needed.')
      end if

      ! Count number of LHS projections
      char = fdf_bnames(pline,2)

      ! Calculate number of projections that is requested
      ! This corresponds to all projections on Gamma_L
      call parse_T(N_Elec,Elecs,N_mol,mols,char,iE_p,im_p,ip_p)

      ! In case we do not calculate all projections.
      ! Then we should not take the last electrode
      if ( (.not. all_T) .and. iE_p == N_Elec ) then
        if ( Node == 0 ) then
          if ( .not. any_skipped ) write(*,*) ! newline
          any_skipped = .true.
          write(*,'(a)')'tbt: Projection "from '//trim(char)//'" has &
              &been silently rejected.'
        end if
        cycle ! all these projections will be skipped
      end if

      ! Check that the -> are all calculated
      itmp = 0 ! itmp tracks whether there are any possible Gamma_R projections
      do while ( fdf_bline(bfdf,pline) )

        ! Skip empty lines
        if ( fdf_bnnames(pline) == 0 ) cycle
        char = fdf_bnames(pline,1)
        if ( leqi(char,'end') ) exit

        ! 2 This corresponds to all the projections on Gamma_R
        call parse_T(N_Elec,Elecs,N_mol,mols,char,iE_c,im_c,ip_c)

        ! we do not allow pure projections (i.e. no projectinos)
        if ( ip_p == 0 .and. ip_c == 0 ) cycle

        ! check that we should actually calculate this projection
        checked = .false.
        if ( all_T ) then
          if ( iE_c /= iE_p ) checked = .true.
        else if ( iE_p < iE_c ) then
          checked = .true.
        end if

        if ( reflect .and. iE_c == iE_p ) checked = .true.

        if ( checked ) then
          itmp = 1 ! denote that we have a possible RHS projection

          ! Ensure that the electrode projection on this molecule
          ! exists.
          if ( ip_c /= 0 ) then
            ! Increment the possible molecule-electrode projections
            if ( N_proj_ME == 0 ) then
              N_proj_ME = N_proj_ME + 1
              call init_ME(tmp_ME(N_proj_ME),Elecs(iE_c),mols(im_c))
            else if ( ME_idx(N_proj_ME,tmp_ME,Elecs(iE_c),mols(im_c)) == 0 ) then
              N_proj_ME = N_proj_ME + 1
              call init_ME(tmp_ME(N_proj_ME),Elecs(iE_c),mols(im_c))
            end if
          end if

        end if
        
      end do
      
      if ( ip_p == 0 .and. itmp == 0 ) then
        
        call init_ME(tmp_ME(1),Elecs(iE_p),mols(im_p))
        call init_ME(tmp_ME(2),Elecs(iE_c),mols(im_c))
        tmp_LME(1)%ME => tmp_ME(1)
        tmp_LME(1)%idx = ip_p
        tmp_LME(2)%ME => tmp_ME(2)
        tmp_LME(2)%idx = ip_c
        call print_proj(tmp_LME(1))
        call print_proj(tmp_LME(2))

        call die('Error in input, a non-projected from has no &
            &viable to projections. Please ensure your input is &
            &correct.')

      end if
      
      ! This means there where no reasonable RHS projections.
      if ( itmp == 0 ) cycle

      ! Add all the left molecule electrode projections,
      ! if they do not exist already.
      ! We first do it here to capture the user doing 
      ! something wrong... (i.e. create Left -> Right, i.e. no projection)
      if ( ip_p /= 0 ) then
        if ( N_proj_ME == 0 ) then
          N_proj_ME = N_proj_ME + 1
          call init_ME(tmp_ME(N_proj_ME),Elecs(iE_p),mols(im_p))
        else if ( ME_idx(N_proj_ME,tmp_ME,Elecs(iE_p),mols(im_p)) == 0 ) then
          N_proj_ME = N_proj_ME + 1
          call init_ME(tmp_ME(N_proj_ME),Elecs(iE_p),mols(im_p))
        end if
      end if

      ! number of different LHS projections
      ! NOTE that this loop have not checked whether there are any
      ! dublicate LHS projections.
      ! This will be checked in the next one.

      ! Count the number of unique (and not already added,
      ! LHS projections)
      if ( ip_p < 0 ) then
        lp1 = -ip_p
        lp2 = -ip_p
      else
        lp1 = 1
        lp2 = ip_p
      end if

      ! Count number of non-counted LHS projections.
      if ( ip_p == 0 ) then
        ! This is a pure LHS projection
        if ( perm_exist(N_proj_T,tmp_LME,Elecs(iE_p),mols(im_p),-iE_p) == 0 ) then
          N_proj_T = N_proj_T + 1
          ! Create the level projection
          tmp_LME(N_proj_T)%idx = -iE_p
        end if
      else
        do ip = lp1, lp2
          if ( perm_exist(N_proj_T,tmp_LME,Elecs(iE_p),mols(im_p),ip) == 0 ) then
            N_proj_T = N_proj_T + 1
            ! Create the level projection
            i = ME_idx(N_proj_ME,tmp_ME,Elecs(iE_p),mols(im_p))
            if ( i == 0 ) &
                call die('Error in programming: proj_T [1]')
            tmp_LME(N_proj_T)%ME  => tmp_ME(i)
            tmp_LME(N_proj_T)%idx =  ip
          end if
        end do
      end if

    end do

    call fdf_brewind(bfdf)

    if ( N_proj_T == 0 .or. N_proj_ME == 0 ) then
      call die('The projection block was ill-formatted &
          &or all projections has been rejected. Check input.')
    end if

    ! Actually read them in...
    allocate(proj_T(N_proj_T),proj_ME(N_proj_ME))

    ! Copy over the different projections saved in tmp_pE
    do it = 1 , N_proj_ME
      
      ! copy over the molecule-electrode projections
      proj_ME(it)%mol => tmp_ME(it)%mol
      proj_ME(it)%El  => tmp_ME(it)%El

      ! We project all molecule levels on to the
      ! electrode Gamma, no matter whether we use them.
      ! This will at least provide an analysis tool! :)
      i = proj_ME(it)%mol%lvls%n
      allocate(proj_ME(it)%bGk(i,i))

    end do

    ipt = 0 ! LHS projection counter
    ! Start reading
    do while ( fdf_bline(bfdf,pline) )
      ! Skip empty lines
      if ( fdf_bnnames(pline) == 0 ) cycle

      char = fdf_bnames(pline,1)
      if ( leqi(char,'from') ) then
        ! We have started a new T-projection 

        char = fdf_bnames(pline,2)

        ! Parse the parent designation (hence _p)
        call parse_T(N_Elec,Elecs,N_mol,mols,char,iE_p,im_p,ip_p)

        if ( .not. all_T .and. iE_p == N_Elec ) then
          ! quick skip till the end
          do while ( fdf_bline(bfdf,pline) ) 
            if ( fdf_bnnames(pline) == 0 ) cycle
            char = fdf_bnames(pline,1)
            if ( leqi(char,'end') ) exit
          end do
          cycle
        end if

        ! Now we need to figure out the number of attached permutations
        ! This turns out to be a little bit difficult due to the 
        ! possibility of individual lines
        ! The best thing to do is to read in this block and incrementally 
        ! attach them individually to each of the previous ones.
        ilme = 0 ! number of RHS projections for this projection
        ! Parse the child designation (hence _c)
        do while ( fdf_bline(bfdf,pline) )

          ! Skip empty lines
          if ( fdf_bnnames(pline) == 0 ) cycle
          char = fdf_bnames(pline,1)
          if ( leqi(char,'end') ) exit

          ! We must have a RHS projection
          ! Parse the parent designation (hence _p)
          call parse_T(N_Elec,Elecs,N_mol,mols,char,iE_c,im_c,ip_c)

          ! If we do not calculate the transmission
          ! to this one, then skip the projection
          if ( ((.not. all_T) .and. iE_c < iE_p) .or. &
              ((.not. reflect) .and. iE_c == iE_p) ) then
            if ( Node == 0 ) then
              if ( .not. any_skipped ) write(*,*) ! newline
              any_skipped = .true.
              write(*,'(a)')'tbt: Projection "from '//&
                  trim(Elecs(iE_p)%name)// &
                  ' to '//trim(char)//'" has been silently rejected.'
            end if
            cycle
          end if

          ! Add all permutations of this one to the list
          if ( ip_c < 0 ) then
            lp1 = -ip_c
            lp2 = -ip_c
          else
            lp1 = 1
            lp2 = ip_c
          end if

          ! Count number of non-counted RHS projections.
          if ( ip_c == 0 ) then
            ! This is a pure LHS projection
            if ( perm_exist(ilme,tmp_LME,Elecs(iE_c),mols(im_c),-iE_c) == 0 ) then
              ilme = ilme + 1
              ! Create the level projection
              tmp_LME(ilme)%idx = -iE_c
            end if
          else
            do ip = lp1, lp2
              if ( perm_exist(ilme,tmp_LME,Elecs(iE_c),mols(im_c),ip) == 0 ) then
                ilme = ilme + 1
                ! Create the level projection
                i = ME_idx(N_proj_ME,proj_ME,Elecs(iE_c),mols(im_c))
                if ( i == 0 ) &
                    call die('Error in programming: proj_T [2]')
                tmp_LME(ilme)%ME  => proj_ME(i)
                tmp_LME(ilme)%idx =  ip
              end if
            end do
          end if

        end do

        ! In case we have found no RHS permutations
        ! Then immediately skip it...
        if ( ilme == 0 ) cycle

        ! Create the outer parent loop of creation index
        if ( ip_p < 0 ) then
          lp1 = -ip_p
          lp2 = -ip_p
        else if ( ip_p == 0 ) then
          ! Signal a pure projection
          lp1 = 0
          lp2 = 0
        else
          lp1 = 1
          lp2 = ip_p
        end if

        ! Loop over all LHS projections
        do it = lp1 , lp2
          ip = it
          if ( it == 0 ) ip = -iE_p

          ! Figure out the current index...
          ! in case it exists we append...
          i = perm_exist(ipt,proj_T(:)%L,Elecs(iE_p),mols(im_p),ip)
          if ( i == 0 ) then
            ! New LHS side
            ipt = ipt + 1
            if ( ipt > N_proj_T ) &
                call die('Error in programming, proj_T [3]')
            proj_T(ipt)%L%idx = ip
            ! attach molecule electrode part
            if ( it /= 0 ) then
              i = ME_idx(N_proj_ME,proj_ME,Elecs(iE_p),mols(im_p))
              proj_T(ipt)%L%ME => proj_ME(i)
            end if
            i = ipt
            !else
            ! it already exists
          end if

          ! In case the RHS is already allocated
          ! we simply move all projections to the tmp_LME
          ! array and remove dublicates, then move them back...
          ilme_o = ilme ! for several projections, they need not
          ! have a shared RHS projections, hence
          ! we need to "delete" the added ones
          if ( allocated(proj_T(i)%R) ) then
            ! Move all to the end of tmp_LME
            do j = 1 , size(proj_T(i)%R)
              checked = .false.
              if ( proj_T(i)%R(j)%idx > 0 ) then
                if ( perm_exist(ilme,tmp_LME, &
                    proj_T(i)%R(j)%ME%El, proj_T(i)%R(j)%ME%mol, &
                    proj_T(i)%R(j)%idx) == 0 ) then
                  ! it does not exist
                  checked = .true.
                end if
              else
                if ( perm_exist(ilme,tmp_LME, Elecs(1),mols(1), &
                    proj_T(i)%R(j)%idx) == 0 ) then
                  checked = .true.
                end if
              end if
              if ( checked ) then
                ! add it
                ilme = ilme + 1
                tmp_LME(ilme)%ME => proj_T(i)%R(j)%ME
                tmp_LME(ilme)%idx = proj_T(i)%R(j)%idx
              end if
            end do
            deallocate(proj_T(i)%R)
          end if

          ! Add the RHS projection permutations
          allocate(proj_T(i)%R(ilme))
          do j = 1 , ilme
            proj_T(i)%R(j)%idx = tmp_LME(j)%idx
            if ( tmp_LME(j)%idx > 0 ) then
              proj_T(i)%R(j)%ME => tmp_LME(j)%ME
            end if
          end do

          ! Re-instantiate the RHS projections in this block
          ilme = ilme_o

        end do

      end if

    end do

    if ( Node == 0 .and. any_skipped ) then
      write(*,'(a)') 'tbt: Certain projections have been skipped'
      write(*,'(a)') 'tbt: Check manual for allowing all projections'
    end if

  contains

    ! This subroutine calculates the number of projections
    ! that is attached to the projection
    ! Say:
    !   M1.P1
    !   M1.P2
    !  'E1.M1' will result in 2 projections
    !  'E1.M1.P1' will result in 1 projection.
    subroutine parse_T(N_Elec,Elecs,N_mol,mols,EMP,iE,im,N_p)
      integer, intent(in) :: N_Elec
      type(Elec), intent(in) :: Elecs(N_Elec)
      integer, intent(in) :: N_mol
      type(tProjMol), intent(in) :: mols(N_mol)
      character(len=*), intent(in) :: EMP
      integer, intent(out) :: iE, im, N_p
      character(len=100) :: E, M, P
      
      integer :: ip, idot, i

      ! We only have a single designation, 
      ! which is a non-projected one
      ! The index indicates that it is a "full" identity(1) projection
      im  = 1 ! just to not fuck up wrong molecule pointers
      N_p = 0
      
      idot = index(EMP,'.') ! we know it must exist
      if ( idot > 0 ) then
        E = EMP(1:idot-1)
      else
        E = trim(EMP)
      end if
      ! Find electrode index
      do iE = 1 , N_Elec
        if ( leqi(E,Elecs(iE)%name) ) exit
      end do
      if ( iE > N_Elec ) then
        write(*,*)'tbt: Could not recognize electrode designation in TBT.Proj.T block'
        write(*,*)'tbt: The electrode named '//trim(E)//' could not be found.'
        call die('Error in input')
      end if
      if ( idot < 1 ) return
      
      im = 0
      M = EMP(idot+1:)
      idot = index(M,'.')
      ! if idot is 0 we use all projections of molecule M
      ! else, we have a finer designation
      if ( idot > 0 ) then
        P = M(idot+1:)
        if ( len_trim(P) == 0 ) then
          call die('TBT.Projs.T projection not designated in: '//trim(EMP))
        end if
        M = M(1:idot-1)
        ! Ensure that the projection exists, else die
        do im = 1 , N_mol
          if ( leqi(M,mols(im)%name) ) then
            i = size(mols(im)%proj)
            N_p = 0
            do ip = 1 , i
              if ( leqi(P,mols(im)%proj(ip)%name) ) then
                ! We have a match (we return the index)
                N_p = - ip
                return
              end if
            end do
          end if
        end do
      else
        do im = 1 , N_mol
          if ( leqi(M,mols(im)%name) ) then
            N_p = size(mols(im)%proj)
            return
          end if
        end do
      end if

      call die('Could not parse input: '//trim(EMP)//' some &
          &projections does not exist.')
      
    end subroutine parse_T

    function perm_exist(N_LME,LME,El,mol,ip) result(i)
      integer, intent(in) :: N_LME
      type(tLvlMolEl), intent(in) :: LME(N_LME)
      type(Elec), intent(in) :: El
      type(tProjMol), intent(in) :: mol
      integer, intent(in) :: ip
      integer :: i
      
      i = 0
      if ( N_LME == 0 ) return

      do i = 1 , N_LME
        if ( LME(i)%idx == ip ) then
          if ( ip < 0 ) then 
            ! this is a pure projection
            ! hence no associated molecule or electrode
            return
          end if
          if ( LME(i)%ME%El == El .and. LME(i)%ME%mol%name == mol%name ) then
            ! The index, electrode and molecule is the same
            return
          end if
        end if
      end do
      i = 0

    end function perm_exist

    subroutine init_ME(ME,El,mol)
      type(tProjMolEl), intent(inout) :: ME
      type(Elec), target :: El
      type(tProjMol), target :: mol
      ME%El => El
      ME%mol => mol
    end subroutine init_ME
    function ME_idx(N,ME,El,mol) result(i)
      integer, intent(in) :: N
      type(tProjMolEl), intent(in) :: ME(N)
      type(Elec), intent(in) :: El
      type(tProjMol), intent(in) :: mol
      integer :: i
      logical :: exist

      do i = 1 , N
        exist = (ME(i)%El == El) .and. (mol%name == ME(i)%mol%name)
        if ( exist ) return
      end do
      i = 0
      
    end function ME_idx

  end subroutine init_proj_T

  subroutine print_proj(LME)
    type(tLvlMolEl), intent(in) :: LME
    write(*,'(a)') trim(proj_ME_name(LME))
  end subroutine print_proj

  function proj_ME_name(LME) result(name)
    type(tLvlMolEl), intent(in) :: LME
    character(len=NF90_MAX_NAME) :: name
    if ( LME%idx > 0 ) then
      name = trim(LME%ME%El%name) // '.'
      name = trim(name) // trim(LME%ME%mol%name) // '.'
      name = trim(name) // trim(LME%ME%mol%proj(LME%idx)%name)
    else
      name = ' '
    end if
  end function proj_ME_name


  function ProjMolEl_same(pLME,El,mol,ip) result(same)
    type(tLvlMolEl), intent(in) :: pLME
    type(Elec), intent(in) :: El
    type(tProjMol), intent(in) :: mol
    integer, intent(in) :: ip
    logical :: same
    
    if ( ip < 0 ) then
      ! it is an electrode comparison
      same = .true.
      if ( ip == pLME%idx ) return
    end if
    if ( pLME%idx < 0 ) then
      same = .false.
      return
    end if

    same = ( pLME%ME%El == El ) 
    if ( same ) then
      ! We can check whether the molecule is the same
      same = ( pLME%ME%mol%name == mol%name )
    end if
    if ( same ) then
      ! We can check whether the projection exists
      same = ( pLME%idx == ip )
    end if
    
  end function ProjMolEl_same

  ! Initialize the TBT.Proj.nc file
  subroutine init_proj_save( fname, TSHS , r, btd, ispin, &
      N_Elec, Elecs, raEl, roElpd, btd_El, &
      nkpt, kpt, wkpt, NE , Eta, a_Dev, a_Buf, sp_dev_sc, save_DATA )

    use parallel, only : Node, Nodes, IONode
    use units, only: eV
    use fdf, only : fdf_get

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use intrinsic_missing, only : VNORM
    use m_os, only : file_exist

    use dictionary, assign_int => assign
    use netcdf_ncdf, ncdf_parallel => parallel
    use m_ncdf_io, only : cdf_w_Sp
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only: MPI_Bcast, MPI_Logical, MPI_Comm_World
    use mpi_siesta, only: MPI_Send, MPI_Recv, MPI_Status_Size
    use mpi_siesta, only: MPI_Double_Precision, MPI_Double_Complex
    use mpi_siesta, only: MPI_Comm_Self, MPI_Integer
    use mpi_siesta, only: MPI_Barrier
#endif
    use m_tbt_hs, only : tTSHS
    use m_tbt_diag

    character(len=*), intent(in) :: fname
    type(tTSHS), intent(inout) :: TSHS
    type(tRgn), intent(in) :: r, btd
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(tRgn), intent(in) :: raEl(N_Elec), roElpd(N_Elec), btd_El(N_Elec)
    integer, intent(in) :: nkpt, NE
    real(dp), intent(in) :: Eta
    real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt)
    type(tRgn), intent(in) :: a_Dev
    type(tRgn), intent(in) :: a_Buf
    type(Sparsity), intent(inout) :: sp_dev_sc
    type(dictionary_t), intent(inout) :: save_DATA

    type(hNCDF) :: ncdf, grp, grp2, grp3
    type(tRgn) :: r_tmp
    integer :: cmp_lvl
    integer :: no, im, Np, ip, i, ik, iN, iE
    integer :: it, ipt
    logical :: exist, is_same, isGamma, save_state
    logical :: debug_state, sme
    type(dictionary_t) :: dic
    character(len=NF90_MAX_NAME) :: tmp
    ! Create allocatables, they are easier to maintain
    integer :: iLUMO, mol_nkpt
    real(dp), allocatable :: eig(:)
    real(dp), allocatable :: rv(:,:), rS_sq(:,:)
    real(dp) :: dn, mem
    character(len=2) :: unit
    complex(dp), allocatable :: zv(:,:), zS_sq(:,:)
    complex(dp) :: zn
    integer :: prec_DOS, prec_T, prec_Teig, prec_J, prec_COOP, prec_DM
    integer :: nnzs_dev, N_eigen, no_e
    type(OrbitalDistribution) :: fdit
#ifdef TBT_PHONON
    character(len=*), parameter :: T_unit = 'g0'
    character(len=*), parameter :: COHP_unit = 'Ry'
#else
    character(len=*), parameter :: T_unit = 'G0'
    character(len=*), parameter :: COHP_unit = 'Ry/Ry'
#endif
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    ! There is nothing to initialize
    if ( N_mol == 0 ) return

    exist = file_exist(fname, Bcast = .true. )

    ! In case the user thinks the double precision
    ! is too much
    call tbt_cdf_precision('DOS','single',prec_DOS)
    call tbt_cdf_precision('T','single',prec_T)
    call tbt_cdf_precision('T.Eig','single',prec_Teig)
    call tbt_cdf_precision('Current','single',prec_J)
    call tbt_cdf_precision('COOP','single',prec_COOP)
    call tbt_cdf_precision('DM','single',prec_DM)
    
    if ( 'T-eig' .in. save_DATA ) then
      call assign_int(N_eigen,save_DATA,'T-eig')
    else
      N_eigen = 0
    end if
    
    ! Whether we should print the debug statements...
    debug_state = fdf_get('TBT.Projs.Debug',.false.)

    ! This allows to easily change between algorithms.
    isGamma = all(TSHS%nsc(:) == 1)
    if ( isGamma ) then
      ! We cannot create a k-point resolved projection
      ! using a Gamma-point calculation
      mols(:)%Gamma = .true.
    end if

    ! First we re-create the pivoting table for the orbitals
    do im = 1 , N_mol
      
      ! Save number of orbitals in molecule
      no = mols(im)%orb%n

      ! Sort the orbitals according to the device region
      call rgn_copy(mols(im)%orb,r_tmp)
      call rgn_sort(r_tmp)
      ! initialize the pivoting array
      call rgn_init(mols(im)%pvt,no)
      ip = 0
      do i = 1 , r%n
        if ( in_rgn(r_tmp,r%r(i)) ) then

          ip = ip + 1
          mols(im)%orb%r(ip) = r%r(i)
          mols(im)%pvt%r(ip) = i

        end if
      end do
      ! ensure it is not sorted
      mols(im)%orb%sorted = .false.
      if ( ip /= no ) call die('Error in orbitals sorting')

    end do
    call rgn_delete(r_tmp)

    if ( IONode ) write(*,*) ! newline

    if ( IONode .and. .not. all(mols(:)%Gamma) ) then
      write(*,'(a)')'tbt: *********'
      write(*,'(a)')'tbt: k-resolved projections only work if dispersion'
      write(*,'(a)')'tbt: does not create band-crossings.'
      write(*,'(a)')'tbt: IT IS YOUR RESPONSIBILITY TO ENSURE THIS!'
      write(*,'(a)')'tbt: Do specific k-points seperately at band-crossings.'
      write(*,'(a/)')'tbt: *********'

    end if

    if ( exist ) then

      ! We just make sure the indices are correct
      ! We do not assure the projections
      ! This will actually allow the user to replace
      ! their own projections if they want.

      call ncdf_open(ncdf,fname,mode=NF90_NOWRITE)

      dic = ('no_u'.kv. TSHS%no_u)//('na_u'.kv. TSHS%na_u )
      dic = dic //('no_d'.kv.r%n) // ('na_d'.kv.a_Dev%n)
      dic = dic //('n_btd'.kv.btd%n)
      if ( a_Buf%n > 0 ) then
        dic = dic // ('na_b'.kv.a_Buf%n)
      end if
      call ncdf_assert(ncdf,is_same,dims=dic)
      call check(dic,is_same,'Dimensions in the PROJ.nc file does not conform &
          &to the current simulation.')

      ! Clean up
      call delete(dic)

      ! Check the variables
      dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) )
      dic = dic // ('xa'.kvp. TSHS%xa) // ('cell'.kvp.TSHS%cell)
      call rgn_copy(a_Dev, r_tmp)
      call rgn_sort(r_tmp)
      dic = dic // ('pivot'.kvp.r%r)//('a_dev'.kvp.r_tmp%r)
      dic = dic // ('nsc'.kvp. TSHS%nsc) // ('btd'.kvp.btd%r)
      if ( a_Buf%n > 0 ) then
        dic = dic // ('a_buf'.kvp.a_Buf%r)
      end if
      call ncdf_assert(ncdf,is_same,vars=dic, d_EPS = 1.e-4_dp )
      call check(dic,is_same,'lasto, xa or cell in the PROJ.nc file does &
          &not conform to the current simulation.',.false.)
      call delete(dic, dealloc=.false.)
      call rgn_delete(r_tmp)

      ! Check the k-points
      allocate(rv(3,nkpt))
      do i = 1 , nkpt
        call kpoint_convert(TSHS%cell,kpt(:,i),rv(:,i),1)
      end do
      dic = ('kpt'.kvp.rv) // ('wkpt'.kvp. wkpt)
      call ncdf_assert(ncdf,is_same,vars=dic, d_EPS = 1.e-7_dp )
      if ( .not. is_same ) then
        call die('k-points or k-weights are not the same')
      end if
      call delete(dic,dealloc = .false. )
      deallocate(rv)

      ! We check each molecule in the projection file
      do im = 1 , N_mol

        call ncdf_open_grp(ncdf,mols(im)%name,grp)

        ! Check number of atoms in projection,
        ! and confirm number of projections.
        ! *** Currently we do not allow extending
        !     the projection file. ***
        dic = ('na'.kv.mols(im)%atom%n) // ('no'.kv.mols(im)%orb%n)
        call ncdf_assert(grp,is_same,dims = dic )
        call check(dic,is_same,'Projection atoms, orbitals and/or number of &
            &projections does not conform to the current simulation.')
        call delete(dic)

        ! Check that the projection atoms and levels are the is_same
        ! we have another group for each lvl
        dic = ('atom'.kvp.mols(im)%atom%r)//('orb'.kvp.mols(im)%orb%r)
        call ncdf_assert(grp,is_same,vars = dic )
        call check(dic,is_same,'Projection atom list &
            &does not conform to the current simulation.',.false.)
        call delete(dic, dealloc=.false.)

        ! The variable 'eig' have to be present,
        ! We do not save the wave-function. 
        ! If we had 1000 k-points and a huge molecule, the file
        ! would be immense. We only save the state in the projection
        dic = ('eig'.kv.1)
        call ncdf_assert(grp,is_same,has_vars=dic)
        call check(dic,is_same,'Projection eigen values &
            &does not exist in the current simulation.')
        call delete(dic)

        ! Loop all groups
        Np = size(mols(im)%proj)
        do ip = 1 , Np

          call ncdf_open_grp(grp,mols(im)%proj(ip)%name,grp2)

          dic = ('nlvl'.kv.mols(im)%proj(ip)%n)
          call ncdf_assert(grp2,is_same,dims = dic )
          call check(dic,is_same,'Projection levels does not &
              &conform to the current simulation.')
          call delete(dic)

          call rgn_copy(mols(im)%proj(ip),r_tmp)
          do iE = 1 , r_tmp%n
            r_tmp%r(iE) = mols(im)%lvls%r(r_tmp%r(iE))
          end do

          ! Check the levels
          dic = ('lvl'.kvp.r_tmp%r)
          call ncdf_assert(grp2,is_same,vars = dic )
          call check(dic,is_same,'Projection level list &
              &does not conform to the current simulation.',.false.)
          call delete(dic, dealloc=.false.)

        end do

      end do

      ! The number of energy points must
      ! not have been set
      call ncdf_inq_dim(ncdf,'ne',len=iE)
#ifdef MPI
      call MPI_Bcast(iE,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

      call ncdf_close(ncdf)

      if ( Node == 0 ) then
        write(*,'(a)') 'tbt: Projection tables are re-used from the &
            &TBT.Proj.nc file.'
      end if

      if ( iE /= 0 ) then
        print *,iE
        call die('Currently the re-usage of the projection files only works &
            &if it has not been set before. &
            &I.e. no energy-points must be saved.')
      end if

      return

    end if
    
    ! The projection file does not exist.
    ! We need to create it.
    if ( Node == 0 ) then
      write(*,'(2a)')'tbt: Initializing projection data file: ',trim(fname)
    end if

    mem = 0._dp
    
    call timer('proj_init',1)

    ! For easiness we do not parallelize this
    ! Typically molecules are also much smaller
    ! and needs only one diagonalization
    ! If needed, we could diagonalize one molecule per processor,
    ! then collect. However, more than often this will probably
    ! be restricted to one molecule.
    cmp_lvl = fdf_get('CDF.Compress',0)
    cmp_lvl = fdf_get('TBT.CDF.Compress',cmp_lvl)
    cmp_lvl = fdf_get('TBT.CDF.Proj.Compress',cmp_lvl)
    if ( cmp_lvl < 0 ) cmp_lvl = 0
    if ( 9 < cmp_lvl ) cmp_lvl = 9

    call ncdf_create(ncdf,fname,mode = NF90_NETCDF4 )

    ! Save the current system size
    call ncdf_def_dim(ncdf,'no_u',TSHS%no_u)
    call ncdf_def_dim(ncdf,'na_u',TSHS%na_u)
    call ncdf_def_dim(ncdf,'nkpt',nkpt)
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)
    call ncdf_def_dim(ncdf,'na_d',a_Dev%n)
    call ncdf_def_dim(ncdf,'no_d',r%n)
    call ncdf_def_dim(ncdf,'ne',NE)
    call ncdf_def_dim(ncdf,'n_s',product(TSHS%nsc))
    call ncdf_def_dim(ncdf,'n_btd',btd%n)

    ! Create eigenvalue dimension, if needed
    if ( N_eigen > 0 ) then
      call ncdf_def_dim(ncdf,'neig',N_eigen)
    end if
    if ( a_Buf%n > 0 ) then
      call ncdf_def_dim(ncdf,'na_b',a_Buf%n)
    end if

#ifdef TBT_PHONON
    dic = ('source'.kv.'PHtrans-Proj')
#else
    dic = ('source'.kv.'TBtrans-Proj')
#endif

    tmp = datestring()
    dic = dic//('date'.kv.tmp(1:10))
#ifndef TBT_PHONON
    dic = dic//('info'.kv.'State levels are wrt. HOMO=-1,Ef=0,LUMO=1')
#endif
    if ( all(TSHS%nsc(:) == 1) ) then
      dic = dic//('Gamma'.kv.'true')
    else
      dic = dic//('Gamma'.kv.'false')
    end if
    if ( TSHS%nspin > 1 ) then
      if ( ispin == 1 ) then
        dic = dic//('spin'.kv.'UP')
      else
        dic = dic//('spin'.kv.'DOWN')
      end if
    end if
    call ncdf_put_gatt(ncdf, atts = dic )
    call delete(dic)

    ! Create all the variables needed to save the states
    dic = ('info'.kv.'Last orbitals of the equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, TSHS%na_u)
    
    dic = dic//('info'.kv.'Unit cell')//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, 3, 3)
    
    dic = dic//('info'.kv.'Atomic coordinates')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, 3, TSHS%na_u)
    call delete(dic)

    dic = ('info'.kv.'Supercell offsets')
    call ncdf_def_var(ncdf,'isc_off',NF90_INT,(/'xyz', 'n_s'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, 3, product(TSHS%nsc))
    
    dic = dic//('info'.kv.'Number of supercells in each direction')
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, 3)

    dic = dic//('info'.kv.'Device region orbital pivot table')
    call ncdf_def_var(ncdf,'pivot',NF90_INT,(/'no_d'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, r%n)

    dic = dic//('info'.kv.'Blocks in BTD for the pivot table')
    call ncdf_def_var(ncdf,'btd',NF90_INT,(/'n_btd'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, btd%n)

    dic = dic//('info'.kv.'Index of device atoms')
    call ncdf_def_var(ncdf,'a_dev',NF90_INT,(/'na_d'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, a_Dev%n)

    if ( a_Buf%n > 0 ) then
      dic = dic//('info'.kv.'Index of buffer atoms')
      call ncdf_def_var(ncdf,'a_buf',NF90_INT,(/'na_b'/), &
          atts = dic)
      mem = mem + calc_mem(NF90_INT, a_Buf%n)
    end if

    dic = dic//('info'.kv.'k point weights')
    call ncdf_def_var(ncdf,'wkpt',NF90_DOUBLE,(/'nkpt'/), &
        atts = dic)
    dic = dic//('info'.kv.'k point')//('unit'.kv.'b')
    call ncdf_def_var(ncdf,'kpt',NF90_DOUBLE,(/'xyz ','nkpt'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, 4, nkpt) ! kpt and wkpt

#ifdef TBT_PHONON
    dic = dic//('info'.kv.'Frequency')//('unit'.kv.'Ry')
#else
    dic = dic//('info'.kv.'Energy')//('unit'.kv.'Ry')
#endif
    call ncdf_def_var(ncdf,'E',NF90_DOUBLE,(/'ne'/), atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, NE)

    dic = dic//('info'.kv.'Imaginary part for device')
#ifdef TBT_PHONON
    dic = dic//('unit'.kv.'Ry**2')
#endif
    call ncdf_def_var(ncdf,'eta',NF90_DOUBLE,(/'one'/), atts = dic)

    call delete(dic)

    call ncdf_put_var(ncdf,'nsc',TSHS%nsc)
    call ncdf_put_var(ncdf,'isc_off',TSHS%isc_off)
    call ncdf_put_var(ncdf,'pivot',r%r)
    call ncdf_put_var(ncdf,'cell',TSHS%cell)
    call ncdf_put_var(ncdf,'xa',TSHS%xa)
    call ncdf_put_var(ncdf,'lasto',TSHS%lasto(1:TSHS%na_u))
    call rgn_copy(a_Dev, r_tmp)
    call rgn_sort(r_tmp)
    call ncdf_put_var(ncdf,'a_dev',r_tmp%r)
    call rgn_delete(r_tmp)
    call ncdf_put_var(ncdf,'btd',btd%r)
    if ( a_Buf%n > 0 ) then
      call ncdf_put_var(ncdf,'a_buf',a_Buf%r)
    end if

    ! Save all k-points
    ! Even though they are in an unlimited dimension,
    ! we save them instantly.
    ! This ensures that a continuation can check for 
    ! the same k-points in the same go.
    allocate(rv(3,nkpt))
    do i = 1 , nkpt
      call kpoint_convert(TSHS%cell,kpt(:,i),rv(:,i),1)
    end do
    call ncdf_put_var(ncdf,'kpt',rv)
    call ncdf_put_var(ncdf,'wkpt',wkpt)
    deallocate(rv)

    call ncdf_put_var(ncdf,'eta',Eta)

    sme = 'proj-orb-current' .in. save_DATA
    sme = sme .or. ('proj-COOP-A' .in. save_DATA)
    sme = sme .or. ('proj-COHP-A' .in. save_DATA)
    sme = sme .or. ('proj-DM-A' .in. save_DATA)
    if ( sme ) then

      ! In case we need to save the device sparsity pattern
      ! Create dimensions
      nnzs_dev = nnzs(sp_dev_sc)
      call ncdf_def_dim(ncdf,'nnzs',nnzs_dev)

      call delete(dic)

      dic = ('info'.kv.'Number of non-zero elements per row')
      call ncdf_def_var(ncdf,'n_col',NF90_INT,(/'no_u'/), &
          atts=dic)
       mem = mem + calc_mem(NF90_INT, TSHS%no_u)

      dic = dic//('info'.kv. &
          'Supercell column indices in the sparse format ')
      call ncdf_def_var(ncdf,'list_col',NF90_INT,(/'nnzs'/), &
          compress_lvl=cmp_lvl,atts=dic )
       mem = mem + calc_mem(NF90_INT, nnzs_dev)

#ifdef MPI
      call newDistribution(TSHS%no_u,MPI_Comm_Self,fdit,name='TBT-fake dist')
#else
      call newDistribution(TSHS%no_u,-1           ,fdit,name='TBT-fake dist')
#endif

      call cdf_w_Sp(ncdf,fdit,sp_dev_sc)
      call delete(fdit)

    end if

    call delete(dic)

    do iE = 1 , N_Elec
      
      call ncdf_def_grp(ncdf,trim(Elecs(iE)%name),grp)
      
      ! Define atoms etc.
      i = TotUsedAtoms(Elecs(iE))
      call ncdf_def_dim(grp,'na',i)
      
      dic = dic//('info'.kv.'Electrode atoms')
      call rgn_range(r_tmp, ELecs(iE)%idx_a, ELecs(iE)%idx_a + i - 1)
      call ncdf_def_var(grp,'a',NF90_INT,(/'na'/), atts = dic)
      call ncdf_put_var(grp,'a',r_tmp%r)
      mem = mem + calc_mem(NF90_INT, r_tmp%n)
      call rgn_delete(r_tmp)
      
      call ncdf_def_dim(grp,'na_down',raEl(iE)%n)
      dic = dic//('info'.kv.'Electrode + downfolding atoms')
      call ncdf_def_var(grp,'a_down',NF90_INT,(/'na_down'/), atts = dic)
      call ncdf_put_var(grp,'a_down',raEl(iE)%r)
      mem = mem + calc_mem(NF90_INT, raEl(iE)%n)
      
      ! Save generic information about electrode
      dic = dic//('info'.kv.'Bloch expansion')
      call ncdf_def_var(grp,'bloch',NF90_INT,(/'xyz'/), atts = dic)
      call ncdf_put_var(grp,'bloch',Elecs(iE)%Bloch%B)
      mem = mem + calc_mem(NF90_INT, 3)

      call ncdf_def_dim(grp,'no_down',roElpd(iE)%n)

      dic = dic//('info'.kv.'Downfolding region orbital pivot table')
      call ncdf_def_var(grp,'pivot_down',NF90_INT,(/'no_down'/), atts = dic)
      call ncdf_put_var(grp,'pivot_down',roElpd(iE)%r)
      mem = mem + calc_mem(NF90_INT, roElpd(iE)%n)
      
      call ncdf_def_dim(grp,'n_btd',btd_El(iE)%n)

      dic = dic//('info'.kv.'Blocks in BTD downfolding for the pivot_down table')
      call ncdf_def_var(grp,'btd',NF90_INT,(/'n_btd'/), atts = dic)
      call ncdf_put_var(grp,'btd',btd_El(iE)%r)
      mem = mem + calc_mem(NF90_INT, btd_El(iE)%n)

      no_e = Elecs(iE)%o_inD%n
      call ncdf_def_dim(grp,'no_e',no_e)

      dic = ('info'.kv.'Orbital pivot table for self-energy')
      call ncdf_def_var(grp,'pivot',NF90_INT,(/'no_e'/), atts = dic)
      call ncdf_put_var(grp,'pivot',Elecs(iE)%o_inD%r)
      mem = mem + calc_mem(NF90_INT, no_e)

      dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
      call ncdf_def_var(grp,'mu',NF90_DOUBLE,(/'one'/), atts = dic)
      call ncdf_put_var(grp,'mu',Elecs(iE)%mu%mu)

#ifdef TBT_PHONON
      dic = dic//('info'.kv.'Phonon temperature')
#else
      dic = dic//('info'.kv.'Electronic temperature')
#endif
      call ncdf_def_var(grp,'kT',NF90_DOUBLE,(/'one'/), atts = dic)
      call ncdf_put_var(grp,'kT',Elecs(iE)%mu%kT)

      dic = dic//('info'.kv.'Imaginary part of self-energy')
#ifdef TBT_PHONON
      dic = dic//('unit'.kv.'Ry**2')
#endif
      call ncdf_def_var(grp,'eta',NF90_DOUBLE,(/'one'/), atts = dic)
      call ncdf_put_var(grp,'eta',Elecs(iE)%Eta)

      dic = dic//('info'.kv.'Accuracy of the self-energy')//('unit'.kv.'Ry')
      call ncdf_def_var(grp,'Accuracy',NF90_DOUBLE,(/'one'/), atts = dic)
      call ncdf_put_var(grp,'Accuracy',Elecs(iE)%accu)
      call delete(dic)

      mem = mem + calc_mem(NF90_DOUBLE, 4) ! mu, kT, eta, Accuracy

    end do

    call delete(dic)

    do im = 1 , N_mol

      ! Whether this molecule is Gamma-calculated
      isGamma = mols(im)%Gamma
      if ( isGamma ) then
        mol_nkpt = 1
      else
        mol_nkpt = nkpt
      end if
      ! Correct for parallel execution
      if ( mod(mol_nkpt,Nodes) > 0 ) then
        mol_nkpt = mol_nkpt + Nodes - mod(mol_nkpt,Nodes)
      end if

      if ( Node == 0 ) then
        write(*,'(2a)') 'tbt: Calculating eigenvalues and |> of projection: ', &
            trim(mols(im)%name)
      end if

      ! Whether the full states are to be saved
      save_state = fdf_get('TBT.Proj.'//trim(mols(im)%name)//'.States',.false.)

      ! # of orbitals for this molecule
      no = mols(im)%orb%n
      call ncdf_def_grp(ncdf,trim(mols(im)%name),grp)

      ! Define the molecule
      call ncdf_def_dim(grp,'na',mols(im)%atom%n)
      call ncdf_def_dim(grp,'no',no)
      ! the different number of projection levels 
      ! must equal the # of states
      call ncdf_def_dim(grp,'nlvl',mols(im)%lvls%n)

      ! A list the used projections
      dic = ('info'.kv.'Used projections indexed with respect to E_F')
      call ncdf_def_var(grp,'lvl',NF90_INT,(/'nlvl'/),atts=dic)
      mem = mem + calc_mem(NF90_INT, mols(im)%lvls%n)
      call ncdf_put_var(grp,'lvl',mols(im)%lvls%r)

      dic = dic//('info'.kv.'|i> = S^(1/2)|v_i> for unique projections')
      if ( isGamma ) then
        call ncdf_def_var(grp,'state',NF90_DOUBLE,(/'no  ','nlvl'/),atts=dic , &
            compress_lvl = cmp_lvl , chunks = (/no,1/) )
        mem = mem + calc_mem(NF90_DOUBLE, mols(im)%lvls%n, no)
      else
        call ncdf_def_var(grp,'state',NF90_DOUBLE_COMPLEX, &
            (/'no  ','nlvl','nkpt'/),atts=dic, &
            compress_lvl = cmp_lvl , chunks = (/no,1,1/) )
        mem = mem + calc_mem(NF90_DOUBLE, mols(im)%lvls%n, no, nkpt, 2)
      end if

      ! Define variables to contain the molecule
      dic = dic//('info'.kv.'Projection atoms')
      call ncdf_def_var(grp,'atom',NF90_INT,(/'na'/),atts=dic)
      call ncdf_put_var(grp,'atom',mols(im)%atom%r)
      mem = mem + calc_mem(NF90_INT, mols(im)%atom%n)

      dic = dic//('info'.kv.'Projection orbitals')
      call ncdf_def_var(grp,'orb',NF90_INT,(/'no'/),atts=dic)
      call ncdf_put_var(grp,'orb',mols(im)%orb%r)
      mem = mem + calc_mem(NF90_INT, no)

      if ( save_state ) then
        dic = dic//('info'.kv.'State |i> = |v_i> for all i')
        if ( isGamma ) then
          call ncdf_def_var(grp,'states',NF90_DOUBLE,(/'no','no'/),atts=dic, &
              compress_lvl = cmp_lvl , chunks = (/no,1/) )
          mem = mem + calc_mem(NF90_DOUBLE, no, no)
        else
          call ncdf_def_var(grp,'states',NF90_DOUBLE_COMPLEX, &
              (/'no  ','no  ','nkpt'/),atts=dic, &
              compress_lvl = cmp_lvl , chunks = (/no,1,1/) )
          mem = mem + calc_mem(NF90_DOUBLE, no, no, nkpt, 2)
        end if
      end if
#ifdef TBT_PHONON
      dic = dic//('info'.kv.'Eigen frequency')//('unit'.kv.'Ry')
#else
      dic = dic//('info'.kv.'Eigenstate energy')//('unit'.kv.'Ry')
#endif
      if ( isGamma ) then
        call ncdf_def_var(grp,'eig',NF90_DOUBLE,(/'no'/),atts=dic, &
            compress_lvl = cmp_lvl )
          mem = mem + calc_mem(NF90_DOUBLE, no)
      else
        call ncdf_def_var(grp,'eig',NF90_DOUBLE,(/'no  ','nkpt'/),atts=dic, &
            compress_lvl = cmp_lvl , chunks = (/no,1/) )
        mem = mem + calc_mem(NF90_DOUBLE, no, nkpt)
      end if
      call delete(dic)

      ! If the projection should save the DOS
      ! add the DOS
      if ( mols(im)%DOS ) then

        ! Save the diagonal |><| overlap variable
        !dic = dic//('info'.kv.'Single projected diagonal overlap matrix')
        !call ncdf_def_var(grp,'kb_SD',NF90_DOUBLE_COMPLEX,(/'no  ','nlvl','nkpt'/), atts = dic)
        !mem = mem + calc_mem(NF90_DOUBLE, no, nlvl, nkpt, 2)

        !dic = dic//('info'.kv.'<|Gf|> / Pi')
        !call ncdf_def_var(grp,'bGfk',NF90_DOUBLE_COMPLEX,(/'nlvl','ne  ','nkpt'/), atts = dic)
        !mem = mem + calc_mem(NF90_DOUBLE, nlvl, NE, nkpt, 2)

        ! Create the DOS variable (actually just kb_SD * bGfk)
        !dic = dic//('info'.kv.'Single projected density of states <|Gf|>kb_SD')
        !call ncdf_def_var(grp,'DOS',prec_DOS,(/'no  ','nlvl','ne  ','nkpt'/), atts = dic)
        !mem = mem + calc_mem(prec_DOS, nlvl, NE, nkpt)

      end if

      ! Number of projections on this molecule...
      Np = size(mols(im)%proj)

      ! Create the projections
      do ip = 1 , Np

        call delete(dic)
        
        call ncdf_def_grp(grp,mols(im)%proj(ip)%name,grp2)
        call ncdf_def_dim(grp2,'nlvl',mols(im)%proj(ip)%n)

        dic = ('info'.kv.'State levels associated with this projection')
        call ncdf_def_var(grp2,'lvl',NF90_INT,(/'nlvl'/), atts = dic)
        mem = mem + calc_mem(NF90_INT, mols(im)%proj(ip)%n)

        ! Create the correct indices
        call rgn_copy(mols(im)%proj(ip),r_tmp)
        do iE = 1 , r_tmp%n
          r_tmp%r(iE) = mols(im)%lvls%r(r_tmp%r(iE))
        end do
        call ncdf_put_var(grp2,'lvl',r_tmp%r)

        ! Save the diagonal |><| overlap variable
        !dic = dic//('info'.kv.'Projected diagonal overlap matrix')
        !call ncdf_def_var(grp2,'kb_SD',NF90_DOUBLE_COMPLEX,(/'no  ','nkpt'/), atts = dic)

        !dic = dic//('info'.kv.'<|Gf|>')
        !call ncdf_def_var(grp2,'bGfk',NF90_DOUBLE_COMPLEX,(/'ne  ','nkpt'/), atts = dic)

        ! Create the DOS variable (actually just kb_SD * bGfk)
        !dic = dic//('info'.kv.'Projected density of states <|Gf|>kb_SD')
        !call ncdf_def_var(grp2,'DOS',prec_DOS,(/'no  ','ne  ','nkpt'/), atts = dic)

        ! Loop all transport stuff and create the variables needed
        Elec_p: do iE = 1 , N_Elec
          do it = 1 , N_proj_T
            call delete(dic)
            
            if ( ProjMolEl_same(proj_T(it)%L, &
                Elecs(iE),mols(im),ip) ) then
              ! we have transport from this electrode molecular projection
              call ncdf_def_grp(grp2,trim(Elecs(iE)%name),grp3)

              dic = dic//('unit'.kv.'1/Ry')
              if ( 'proj-DOS-A' .in. save_DATA ) then
                dic = dic//('info'.kv.'Spectral function density of states')
                call ncdf_def_var(grp3,'ADOS',prec_DOS,(/'no_d','ne  ','nkpt'/), &
                    atts = dic, &
                    compress_lvl = cmp_lvl, chunks = (/r%n,1,1/))
                mem = mem + calc_mem(prec_DOS, r%n, NE, nkpt)
              end if

              if ( 'proj-DM-A' .in. save_DATA ) then
                dic = dic//('info'.kv.'Spectral function density matrix')
                call ncdf_def_var(grp3,'DM',prec_DM,(/'nnzs','ne  ','nkpt'/), &
                    atts = dic , chunks = (/nnzs_dev/), compress_lvl=cmp_lvl)
                mem = mem + calc_mem(prec_DM, nnzs_dev, NE, nkpt)
              end if

              if ( 'proj-COOP-A' .in. save_DATA ) then
                dic = dic//('info'.kv.'Crystal orbital overlap population')
                call ncdf_def_var(grp3,'COOP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
                    atts = dic , chunks = (/nnzs_dev/), compress_lvl=cmp_lvl)
                mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
              end if

              if ( 'proj-COHP-A' .in. save_DATA ) then
                dic = dic//('info'.kv.'Crystal orbital Hamilton population')//('unit'.kv.COHP_unit)
                call ncdf_def_var(grp3,'COHP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
                    atts = dic , chunks = (/nnzs_dev/), compress_lvl=cmp_lvl)
                mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
              end if

              ! Prepare for orb-current
              call delete(dic)
              dic = ('unit'.kv.T_unit)
              
              if ( 'proj-orb-current' .in. save_DATA ) then
                dic = dic//('info'.kv.'Orbital current')
                call ncdf_def_var(grp3,'J',prec_J,(/'nnzs', 'ne  ', 'nkpt'/), &
                    atts = dic, chunks = (/nnzs_dev/), compress_lvl = cmp_lvl)
                mem = mem + calc_mem(prec_J, nnzs_dev, NE, nkpt)
              end if

              ! Now we create all transport related quantities
              do ipt = 1 , size(proj_T(it)%R)

                dic = dic//('info'.kv.'Transmission')

                i = proj_T(it)%R(ipt)%idx
                if ( i < 0 ) then
                  tmp = trim(Elecs(-i)%name)
                  if ( i == -iE ) then
                    dic = dic//('info'.kv.'Gf transmission')
                    call ncdf_def_var(grp3,trim(tmp)//'.T',prec_T, (/'ne  ','nkpt'/), &
                        atts = dic, chunks =(/NE, 1/) )
                    dic = dic//('info'.kv.'Out transmission correction')
                    mem = mem + calc_mem(prec_T, NE, nkpt)
                    tmp = trim(tmp)//'.C'
                  else
                    tmp = trim(tmp)//'.T'
                  end if
                else
                  tmp = proj_ME_name(proj_T(it)%R(ipt))
                  if ( proj_T(it)%R(ipt)%ME%El == Elecs(iE) ) then
                    dic = dic//('info'.kv.'Gf transmission')
                    call ncdf_def_var(grp3,trim(tmp)//'.T',prec_T, (/'ne  ','nkpt'/), &
                        atts = dic, chunks =(/NE, 1/) )
                    dic = dic//('info'.kv.'Out transmission correction')
                    mem = mem + calc_mem(prec_T, NE, nkpt)
                    tmp = trim(tmp)//'.C'
                  else
                    tmp = trim(tmp)//'.T'
                  end if
                end if

                call ncdf_def_var(grp3,tmp,prec_T, (/'ne  ','nkpt'/), &
                    atts = dic, chunks = (/NE, 1/) )
                mem = mem + calc_mem(prec_T, NE, nkpt)

                if ( N_eigen > 0 ) then
                  dic = dic//('info'.kv.'Transmission eigenvalues')
                  call ncdf_def_var(grp3,trim(tmp)//'.Eig',prec_Teig, &
                      (/'neig','ne  ','nkpt'/), &
                      atts = dic, chunks =(/N_eigen, NE, 1/) )
                  mem = mem + calc_mem(prec_Teig, N_eigen, NE, nkpt)
                end if

              end do

              ! this electrode will only occur once
              cycle Elec_p

            end if
          end do
        end do Elec_p

      end do
      call delete(dic)
      call rgn_delete(r_tmp)

      ! Allocate space for calculating the eigen-values of the MPSH
      allocate(eig(no))
      if ( isGamma ) then 
        allocate(rS_sq(no,no),rv(no,no))
      else
        allocate(zS_sq(no,no),zv(no,no))
      end if

      ! We need to pre-calculate the Fermi-level
      ! to get a common Fermi-level for all k
      if ( isGamma ) then
        call calc_Eig(TSHS%H_2D,TSHS%S_1D,mols(im)%orb,eig,rv)
      else
        call calc_Eig(TSHS%H_2D,TSHS%S_1D,mols(im)%orb,eig)
      end if

#ifdef TBT_PHONON
      ! Calculate the frequency
      do i = 1 , no 
        if ( eig(i) > 0._dp ) then
          eig(i) =  sqrt( eig(i) )
        else
          ! Signal instability
          eig(i) = -sqrt(abs(eig(i)))
        end if
      end do

      ! In phonon transport there is no "unoccupied" levels
      ! However, it makes no sense to count from -size()
      ! Hence we count from 1
      iLUMO = 1

#else
      ! figure out the LUMO level
      iLUMO = no
      do i = 1 , no 
        if ( eig(i) > 0._dp ) then
          iLUMO = i
          exit
        end if
      end do
      ! iLUMO now contains the index of the LUMO lvl

      ! Create attribute to contain the index of the HOMO level
      dic = ('HOMO_index'.kv.iLUMO-1)
      call ncdf_put_gatt(grp,atts=dic)
      call delete(dic)

#endif

      ikpt: do ik = 1 + Node , mol_nkpt , Nodes

        if ( debug_state ) then
          if ( isGamma ) then
            write(*,'(/a)')'tbt-proj-DEBUG: Gamma-point'
          else
            write(*,'(/a,3(f10.5))')'tbt-proj-DEBUG: k-point',kpt(:,ik)
          end if
        end if

        if ( isGamma ) then
          ! We have already calculated the eigen-values for
          ! the Gamma point
        else if ( ik <= nkpt ) then
          call calc_Eig(TSHS%H_2D,TSHS%S_1D,product(TSHS%nsc), &
              TSHS%sc_off,mols(im)%orb,eig,kpt(:,ik),zv)
#ifdef TBT_PHONON
          ! Rescale eigenvalues to frequency
          do i = 1 , no 
            if ( eig(i) > 0._dp ) then
              eig(i) =  sqrt( eig(i) )
            else
              ! Signal instability
              eig(i) = -sqrt(abs(eig(i)))
            end if
          end do
#endif
        end if

        ! Save the eigen-values
        if ( isGamma ) then
          call ncdf_put_var(grp,'eig',eig)
        else

          call ncdf_put_var(grp,'eig',eig,start=(/1,ik/))
#ifdef MPI
          if ( Node == 0 ) then
            ! The eigen-values
            do iN = 1 , Nodes - 1
              if ( ik + iN > nkpt ) exit
              call MPI_Recv(eig,no,MPI_Double_Precision,iN,iN, &
                  MPI_Comm_World,status,MPIerror)
              call ncdf_put_var(grp,'eig',eig,start=(/1,ik+iN/))
            end do
          else if ( ik <= nkpt ) then
            call MPI_Send(eig,no,MPI_Double_Precision,0,Node, &
                MPI_Comm_World,MPIerror)
          end if
#endif
        end if

        ! Save the states
        if ( save_state ) then
          if ( isGamma ) then
            call ncdf_put_var(grp,'states',rv)
          else
            call ncdf_put_var(grp,'states',zv,start=(/1,1,ik/))
#ifdef MPI
            if ( Node == 0 ) then
              do iN = 1 , Nodes - 1
                if ( ik + iN > nkpt ) exit
                call MPI_Recv(zS_sq(1,1),no*no,MPI_Double_Complex,iN,iN, &
                    MPI_Comm_World,status,MPIerror)
                call ncdf_put_var(grp,'states',zS_sq,start=(/1,1,ik+iN/))
              end do
            else if ( ik <= nkpt ) then
              call MPI_Send(zv(1,1),no*no,MPI_Double_Complex,0,Node, &
                  MPI_Comm_World,MPIerror)
            end if
#endif
          end if
        end if

        if ( isGamma ) then
          call calc_sqrt_S(TSHS%S_1D,mols(im)%orb,rS_sq)
          call norm_Eigenstate(mols(im)%orb%n,rv,rS_sq)
        else if ( ik <= nkpt ) then
          call calc_sqrt_S(TSHS%S_1D,product(TSHS%nsc),TSHS%sc_off, &
              mols(im)%orb,zS_sq,kpt(:,ik))
          call norm_Eigenstate(mols(im)%orb%n,zv,zS_sq)
        end if

        if ( debug_state ) then

          write(*,'(a)')' <j|S^1/2 S^1/2|i> = \delta_ij'

          ! Ensure the orthogonality
          if ( isGamma ) then
            call dgemm('T','N',no,no,no,1._dp, &
                rv,no,rv,no, &
                0._dp, rS_sq, no)
            ! Print the norm and the diagonal element
            do i = 1 , no
              dn = VNORM(rS_sq(:,i))
              write(*,'(tr3,i4,2(a,e10.5))') i,' <:|i> = ',dn, &
                  ' <i|i> = ',rS_sq(i,i)
            end do
          else
            call zgemm('C','N',no,no,no,dcmplx(1._dp,0._dp), &
                zv,no,zv,no, &
                dcmplx(0._dp,0._dp), zS_sq, no)
            ! Print the norm and the diagonal element
            do i = 1 , no
              zn = VNORM(zS_sq(:,i))
              write(*,'(tr3,i4,2(a,2(tr1,e10.5)))') i,' <:|i> =',zn, &
                  ' <i|i> =',zS_sq(i,i)
            end do
          end if

          write(*,'(a)')' \sum S^1/2|i><i|S^1/2 = I'

          ! Ensure that \sum |i><i> = I
          if ( isGamma ) then
            rS_sq = 0._dp
            ! Print the norm and the diagonal element
            do i = 1 , no
              call dger(no,no,1._dp,rv(1,i),1,rv(1,i),1,rS_sq(1,1),no)
            end do
            do i = 1 , no
              dn = VNORM(rS_sq(:,i))
              write(*,'(tr3,i4,2(a,e10.5))') i,' \sum_:i = ',dn, &
                  ' \sum_ii = ',rS_sq(i,i)
            end do
          else
            zS_sq = 0._dp
            ! Print the norm and the diagonal element
            do i = 1 , no
              call zgerc(no,no,dcmplx(1._dp,0._dp),zv(1,i),1,zv(1,i),1, &
                  zS_sq(1,1),no)
            end do
            do i = 1 , no
              zn = VNORM(zS_sq(:,i))
              write(*,'(tr3,i4,2(a,2(tr1,e10.5)))') i,' \sum_:i =',zn, &
                  ' \sum_ii =',zS_sq(i,i)
            end do
          end if

        end if

        ! Copy over state levels
        Np = mols(im)%lvls%n
        do ip = 1 , Np

          ! Calculate the index for the state
          i = mols(im)%lvls%r(ip)
          if ( i > 0 ) then
            ! We are asking for LU<>+i
            i = i + iLUMO - 1
          else
            i = i + iLUMO
          end if
          ! Check that the state actually exists
          if ( i < 1 .or. no < i ) then
            write(*,'(a)')'tbt: Projection eigenvalues [eV]:'
            do i = 1 , no
              if ( i < iLUMO ) then
                write(*,'(a,tr1,i4,tr1,e12.5)')'Eigenvalue: ', &
                    i-iLUMO,eig(i)/eV
              else
                write(*,'(a,tr1,i4,tr1,e12.5)')'Eigenvalue: ', &
                    i-iLUMO+1,eig(i)/eV
              end if
            end do
            call die('Requested levels for projections does not exist, &
                &please check eigenvalues before proceeding...')
          end if

          if ( isGamma ) then
            rS_sq(:,ip) = rv(:,i)
          else
            zS_sq(:,ip) = zv(:,i)
          end if

        end do

        ! Save states
        if ( isGamma ) then
          call ncdf_put_var(grp,'state',rS_sq(:,1:Np), &
              start=(/1,1/))
        else
          call ncdf_put_var(grp,'state',zS_sq(:,1:Np), &
              start=(/1,1,ik/))
#ifdef MPI
          if ( Node == 0 ) then
            do iN = 1 , Nodes - 1
              if ( ik + iN > nkpt ) exit
              call MPI_Recv(zS_sq(1,1),no*Np,MPI_Double_Complex,iN,iN, &
                  MPI_Comm_World,status,MPIerror)
              call ncdf_put_var(grp,'state',zS_sq(:,1:Np), &
                  start=(/1,1,ik+iN/))
            end do
          else if ( ik <= nkpt ) then
            call MPI_Send(zS_sq(1,1),no*Np,MPI_Double_Complex,0,Node, &
                MPI_Comm_World,MPIerror)
          end if
#endif

        end if

      end do ikpt

      if ( isGamma ) then
        deallocate(rS_sq,rv)
      else
        deallocate(zS_sq,zv)
      end if
      deallocate(eig)

    end do

    ! Loop on all electrode projections
    dic = dic//('info'.kv.'Projected scattering rate <i|Gam|j>')
    dic = dic//('unit'.kv.'Ry')
    do it = 1 , N_proj_ME

      call ncdf_open_grp(ncdf,trim(proj_ME(it)%mol%name),grp)
      ! Append electrode name to create variable
      tmp = trim(proj_ME(it)%El%name)//'.bGk'

      i = proj_ME(it)%mol%lvls%n
      call ncdf_def_var(grp,tmp,NF90_DOUBLE_COMPLEX, &
          (/'nlvl','nlvl','ne  ','nkpt'/), atts = dic, &
          compress_lvl = cmp_lvl , chunks = (/i,i,1,1/) )
      mem = mem + calc_mem(NF90_DOUBLE, i*i, NE, nkpt, 2)

    end do

    ! clean-up
    call delete(dic)

    ! At this point we still need to add the "non-projected"
    ! LHS projections.
    dic = 'unit'.kv.T_unit
    do it = 1 , N_proj_T

      i = proj_T(it)%L%idx
      if ( i > 0 ) cycle

      ! We have a simple projection
      ! Create electrode group
      call ncdf_open_grp(ncdf,Elecs(-i)%name,grp)

      ! Now we create all transport related quantities
      do ipt = 1 , size(proj_T(it)%R)

        ! Create the transport to
        if ( proj_T(it)%R(ipt)%idx < 0 ) then
          call die('This is a pure transport problem... Do NOT do that!')
        end if
        tmp = proj_ME_name(proj_T(it)%R(ipt))
        if ( proj_T(it)%R(ipt)%ME%El == Elecs(-i) ) then
          dic = dic//('info'.kv.'Gf transmission')
          call ncdf_def_var(grp,trim(tmp)//'.T',prec_T, (/'ne  ','nkpt'/), &
              atts = dic)
          dic = dic//('info'.kv.'Out transmission correction')
          mem = mem + calc_mem(prec_T, NE, nkpt)
          tmp = trim(tmp)//'.C'
        else
          dic = dic//('info'.kv.'Transmission')
          tmp = trim(tmp)//'.T'
        end if

        call ncdf_def_var(grp,tmp,prec_T, (/'ne  ','nkpt'/), &
            atts = dic)
        mem = mem + calc_mem(prec_T, NE, nkpt)

        if ( N_eigen > 0 ) then
          dic = dic//('info'.kv.'Transmission eigenvalues')
          call ncdf_def_var(grp,trim(tmp)//'.Eig',prec_Teig, &
              (/'neig','ne  ','nkpt'/), atts = dic )
          mem = mem + calc_mem(prec_Teig, N_eigen, NE, nkpt)

        end if

      end do

    end do

    call delete(dic)

    call ncdf_close(ncdf)

    if ( Node == 0 ) then
      call pretty_memory(mem, unit)
      write(*,'(3a,f8.3,tr1,a/)') 'tbt: Estimated file size of ', trim(fname), ':', &
          mem, unit
      
      write(*,'(a)') 'tbt: Please ensure the projection eigenvalues &
          &are aligned as you suspect.'
      write(*,'(a)') 'tbt: Molecular states hybridize in proximity &
          &and energy levels might shift.'
      write(*,*) 
    end if

#ifdef MPI
    ! This ensures the timing is correct
    ! AND that the below die command will not be 
    ! executed prematurely.
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

    call timer('proj_init',2)

    if ( fdf_get('TBT.Projs.Init',.false.) ) then
      call die('You have requested to only initialize the &
          &projection tables. We die by your request.')
    end if

  contains

    subroutine check(dic,same,msg,dealloc)
      type(dictionary_t), intent(inout) :: dic
      logical, intent(inout) :: same
      character(len=*), intent(in) :: msg
      logical, intent(in), optional :: dealloc
      integer :: i
      call delete(dic,dealloc=dealloc)
#ifdef MPI
      call MPI_Bcast(same,1,MPI_Logical,0, &
          MPI_Comm_World, i)
#endif
      if ( .not. same ) then
        call die(msg)
      end if
    end subroutine check

    pure function calc_mem(prec_nf90, n1, n2, n3, n4) result(kb)
      use precision, only: dp
      integer, intent(in) :: prec_nf90, n1
      integer, intent(in), optional :: n2, n3, n4
      real(dp) :: kb

      kb = real(n1, dp) / 1024._dp
      if ( present(n2) ) kb = kb * real(n2, dp)
      if ( present(n3) ) kb = kb * real(n3, dp)
      if ( present(n4) ) kb = kb * real(n4, dp)
      
      select case ( prec_nf90 )
      case ( NF90_INT, NF90_FLOAT )
        kb = kb * 4
      case ( NF90_DOUBLE )
        kb = kb * 8
      end select
      
    end function calc_mem

    pure subroutine pretty_memory(mem, unit)
      use precision, only: dp
      real(dp), intent(inout) :: mem
      character(len=2), intent(out) :: unit

      unit = 'KB'
      if ( mem > 1024._dp ) then
        mem = mem / 1024._dp
        unit = 'MB'
        if ( mem > 1024._dp ) then
          mem = mem / 1024._dp
          unit = 'GB'
          if ( mem > 1024._dp ) then
            mem = mem / 1024._dp
            unit = 'TB'
          end if
        end if
      end if

    end subroutine pretty_memory

  end subroutine init_proj_save

  subroutine proj_cdf_save(ncdf, N_Elec, Elecs, &
      ikpt, nE, N_proj_T, proj_T, pDOS, T, &
      N_eigen, Teig, &
      save_DATA)
    
    use parallel, only : Node, Nodes

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_proj_T
    type(tProjT), intent(in) :: proj_T(N_proj_T)
    real(dp), intent(in) :: pDOS(:,:,:)
    real(dp), intent(in) :: T(:,:)
    integer, intent(in) :: N_eigen
    real(dp), intent(in) :: Teig(:,:,:)
    type(dictionary_t), intent(in) :: save_DATA

    type(hNCDF) :: gmol, gproj, gEl
    integer :: ipt, ip, iE, i, iN
    character(len=NF90_MAX_NAME) :: cmol, cproj, ctmp
    logical :: same_E
    integer :: NDOS, NT
    integer :: cnt(2), idx(2)
#ifdef MPI
    real(dp), allocatable, target :: thisDOS(:)
    real(dp), allocatable :: rT(:,:)
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    ! Grab size of arrays
    NDOS = size(pDOS,dim=1)
    NT = size(T,dim=1)
    
#ifdef MPI
    if ( Node == 0 .and. .not. save_parallel ) then
      ! Allocate the buffer array
      if ( N_eigen > NDOS ) then
        allocate(thisDOS(N_eigen))
      else
        allocate(thisDOS(NDOS))
      end if
      call save_attach_buffer(thisDOS)
      
      allocate(rT(NT,Nodes-1))
      
    end if
#endif
    
    cmol  = ' '
    cproj = ' '
    do ipt = 1 , N_proj_T

      ! We this is the same as a non-projected
      ! spectral function, hence we skip it...
      iE = proj_T(ipt)%L%idx
      if ( iE < 0 ) then

        ! Pristine left hand side scattering states
        call ncdf_open_grp(ncdf,Elecs(-iE)%name,gEl)

      else

        if ( cmol /= proj_T(ipt)%L%ME%mol%name ) then
          cmol = proj_T(ipt)%L%ME%mol%name
          call ncdf_open_grp(ncdf,cmol,gmol)
          cproj = ' ' ! forces the projection to be read in
        end if
        if ( cproj /= proj_T(ipt)%L%ME%mol%proj(iE)%name ) then
          cproj = proj_T(ipt)%L%ME%mol%proj(iE)%name
          call ncdf_open_grp(gmol,cproj,gproj)
        end if

        call ncdf_open_grp(gproj,proj_T(ipt)%L%ME%El%name,gEl)

        ! We save the DOS calculated from the spectral function
        if ( 'proj-DOS-A' .in. save_DATA ) then

          call local_save_DOS(gEl,'ADOS',ikpt,nE,NDOS,pDOS(:,2,ipt))

        end if

      end if

#ifdef MPI
      if ( .not. save_parallel ) then
        if ( Node == 0 ) then
          do iN = 1 , Nodes - 1
            if ( nE%iE(iN) > 0 ) then
              call MPI_Recv(rT(1,iN),NT,Mpi_double_precision, &
                  iN, iN, Mpi_comm_world,status,MPIerror)
            end if
          end do
        else if ( nE%iE(Node) > 0 ) then
          call MPI_Send(T(1,ipt),NT,Mpi_double_precision, &
              0, Node, Mpi_comm_world,MPIerror)
        end if
      end if
#endif

      ! Loop different -> terminal transports
      do ip = 1 , size(proj_T(ipt)%R)

        ! Create name
        i = proj_T(ipt)%R(ip)%idx
        same_E = .false.
        if ( i < 0 ) then
          ctmp = trim(Elecs(-i)%name)
          ! We do not allow pure to pure, hence this will
          ! always succeed
          if ( Elecs(-i) == proj_T(ipt)%L%ME%El ) then
            same_E = .true.
          end if
        else
          ctmp = proj_ME_name(proj_T(ipt)%R(ip))
          if ( iE < 0 ) then
            same_E = (Elecs(-iE) == proj_T(ipt)%R(ip)%ME%El)
          else if ( proj_T(ipt)%R(ip)%ME%El == proj_T(ipt)%L%ME%El ) then
            same_E = .true.
          end if
        end if
        if ( same_E ) then
          ctmp = trim(ctmp) // '.C'
        else
          ctmp = trim(ctmp) // '.T'
        end if

        idx = (/nE%iE(Node),ikpt/)
        cnt(:) = 1
        if ( idx(1) <= 0 ) then
          idx(1) = 1
          ! Denote no storage 
          cnt = 0
        end if

        ! Store data
        call ncdf_put_var(gEl,ctmp,T(ip,ipt), start=idx, count=cnt)

#ifdef MPI
        if ( Node == 0 .and. .not. save_parallel ) then
          do iN = 1 , Nodes - 1
            if ( nE%iE(iN) > 0 ) then
              call ncdf_put_var(gEl,ctmp,rT(ip,iN), &
                  start = (/nE%iE(iN),ikpt/) )
            end if
          end do
        end if
#endif

        if ( N_eigen > 0 ) then
          call local_save_DOS(gEl,trim(ctmp)//'.Eig',ikpt,nE,&
              N_eigen,Teig(:,ip,ipt))
        end if

      end do

    end do
    
#ifdef MPI
    if ( allocated(thisDOS) ) deallocate(thisDOS)
    if ( allocated(rT) ) deallocate(rT)
#endif
    
  end subroutine proj_cdf_save

  subroutine proj_cdf_save_sp_dev(ncdf, ikpt, nE, &
      var_name, LME, sp_dev)

    use parallel, only : Node, Nodes
    use class_dSpData1D

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
    use mpi_siesta, only : MPI_Send, MPI_Recv
    use mpi_siesta, only : MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    character(len=*), intent(in) :: var_name
    type(tLvlMolEl), intent(inout) :: LME
    type(dSpData1D), intent(inout) :: sp_dev
    
    type(hNCDF) :: gmol, gproj, gEl
    integer :: nnzs_dev, idx(3), cnt(3)
    real(dp), pointer :: D(:)
#ifdef MPI
    integer :: iN
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif
    
    D => val(sp_dev)
    nnzs_dev = size(D)

    ! We save the orbital current
    !  1. open molecule
    call ncdf_open_grp(ncdf,trim(LME%ME%mol%name),gmol)
    !  2. open projection
    call ncdf_open_grp(gmol,trim(LME%ME%mol%proj(LME%idx)%name),gproj)
    !  3. open electrode (spectral function)
    call ncdf_open_grp(gproj,trim(LME%ME%El%name),gEl)

    if ( save_parallel ) then

      idx(1) = 1
      idx(3) = ikpt
      cnt(1) = nnzs_dev
      cnt(2) = 1
      cnt(3) = 1

      if ( nE%iE(Node) > 0 ) then
        idx(2) = nE%iE(Node)
      else
        idx(2) = 1
        cnt = 0
      end if
      call ncdf_put_var(gEl,var_name,D,start = idx, count=cnt)

    else
       
      ! Save the current
      if ( nE%iE(Node) > 0 ) then
        call ncdf_put_var(gEl,var_name,D,start = (/1,nE%iE(Node),ikpt/) )
      end if

#ifdef MPI
      if ( Node == 0 ) then
        do iN = 1 , Nodes - 1
          if ( nE%iE(iN) > 0 ) then
            call MPI_Recv(D,nnzs_dev,Mpi_double_precision, &
                iN, iN, Mpi_comm_world,status,MPIerror)
            call ncdf_put_var(gEl,var_name,D,start = (/1,nE%iE(iN),ikpt/) )
          end if
        end do
      else if ( nE%iE(Node) > 0 ) then
        call MPI_Send(D(1),nnzs_dev,Mpi_double_precision, &
            0, Node, Mpi_comm_world,MPIerror)
      end if
#endif
    end if
    
  end subroutine proj_cdf_save_sp_dev
  
  ! Returns the projection state for the designated
  ! molecule projection.
  subroutine proj_update(ncdf,N_mol,mols,ikpt)

    use netcdf_ncdf
#ifdef MPI
    use mpi_siesta, only: MPI_Bcast
    use mpi_siesta, only: MPI_Double_Complex, MPI_Comm_World
#endif

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: N_mol
    type(tProjMol), intent(inout) :: mols(N_mol)
    integer, intent(in) :: ikpt

    ! Local variables
    type(hNCDF) :: grp
    integer :: no, im
    real(dp), allocatable :: rp(:,:)
#ifdef MPI
    integer :: MPIerror
#endif
    
    do im = 1 , N_mol

      ! We quickly skip this molecule if it is 
      ! a Gamma-projection and the k-point
      ! is higher than 0, in that case the
      ! molecule already has the correct projectors
      if ( mols(im)%Gamma .and. ikpt > 1 ) cycle

      ! Open the corresponding projection group
      call ncdf_open_grp(ncdf,mols(im)%name,grp)

      ! Get size of projection
      no = mols(im)%orb%n

      if ( mols(im)%Gamma ) then
        ! Read in from Gamma file
        allocate(rp(no,mols(im)%lvls%n))
        call ncdf_get_var(grp,'state',rp)
        mols(im)%p(:,:) = rp(:,:)
        deallocate(rp)
      else
        call ncdf_get_var(grp,'state',mols(im)%p,start=(/1,1,ikpt/))
      end if

#ifdef MPI
      if ( .not. save_parallel ) then
        ! B-cast information
        no = no * mols(im)%lvls%n
        call MPI_Bcast(mols(im)%p(1,1),no, &
            MPI_Double_Complex, 0, MPI_Comm_World, MPIerror)
      end if
#endif

    end do

  end subroutine proj_update

  ! Returns the projection state for the designated
  ! molecule projection.
  subroutine proj_cdf_save_S_D(fname,N_mol,mols,ikpt,S_1D,nwork,zwork)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use parallel, only : Node, Nodes

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only: MPI_Send, MPI_Recv, MPI_Get_Count
    use mpi_siesta, only: MPI_STATUS_SIZE
    use mpi_siesta, only: MPI_Double_Complex, MPI_Comm_World
#endif

    character(len=*), intent(in) :: fname
    integer, intent(in) :: N_mol
    type(tProjMol), intent(inout) :: mols(N_mol)
    integer, intent(in) :: ikpt
    type(zSpData1D), intent(inout) :: S_1D
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: zwork(nwork)

    ! Local variables
    type(hNCDF) :: ncdf, grp, grp2
    integer :: no, Np, Npl, im, poff, Ns, Nsl, soff
    complex(dp), allocatable :: zD(:,:), zP(:,:)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    complex(dp), pointer :: M(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: i, io, lio, ind, j, ip

#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif
    
    ! Assert nwork to be of good size
    no = maxval(mols(:)%orb%n)
    im = 1
    do i = 1 , N_mol
      im = max(im,size(mols(i)%proj) + 1)
    end do
    if ( no ** 2 * (im/Nodes) > nwork ) then
      call die('Work size for projections too small, do you &
          &have an excessive amount of projections?')
    end if

    ! Get data
    dit => dist(S_1D)
    sp  => spar(S_1D)
    M   => val (S_1D)

    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    ! Open the file for writing
    call ncdf_open(ncdf,fname, mode = NF90_WRITE)

    do im = 1 , N_mol

      ! If we are not asked to calculate the DOS
      ! we can easily skip this calculation.
      if ( .not. mols(im)%DOS ) cycle

      ! open the group
      call ncdf_open_grp(ncdf,mols(im)%name,grp)

      ! Size of the current projections
      no = mols(im)%orb%n
      Np = mols(im)%lvls%n
      Ns = size(mols(im)%proj)

#ifdef MPI
      ! The level projections
      Npl = Np / Nodes
      ! The offset of the current node
      poff = Npl * Node
      ! Correct each node
      if ( Npl * Nodes + Node < Np ) then
        Npl = Npl + 1
        ! all previous nodes have also added a calculating unit
        poff = poff + Node
      else if ( Npl * Nodes < Np ) then
        ! the offset of this node (which do not it self
        ! get any more calculations) are the missing ones
        poff = poff + Np - Npl * Nodes
      end if
      ! Do the same for the state projections
      Nsl = Ns / Nodes
      soff = Nsl * Node
      if ( Nsl * Nodes + Node < Ns ) then
        Nsl = Nsl + 1
        soff = soff + Node
      else if ( Nsl * Nodes < Ns ) then
        soff = soff + Ns - Nsl * Nodes
      end if
#else
      ! Not strictly needed, but for clarity in the sequential version
      Npl = Np
      poff = 0
      Nsl = Ns
      soff = 0
#endif

      ! We pre-calculate all state projections that will be calculated
      ! on this node (note these are full matrices ordered
      ! for fastest access
      ! They are ordered:
      !   [1: no*Nsl] = [|1>_1<1|,|2>_1<2|,...]
      !   [no*Nsl+1:] = [|1>_2<1|,|2>_2<2|,...]
      ! etc.
      ! This is how we access them later down...
      do i = 1 , mols(im)%orb%n
        do ip = 1 , Nsl
          j = ((i-1)*Nsl+ip-1) * no + 1
          call proj_state_bra(mols(im),mols(im)%proj(soff+ip), &
              i, zwork(j:j+no-1) )
        end do
      end do
      
      ! Allocate size for the level and state projections
      allocate(zD(no,Npl),zP(no,Nsl))

!$OMP parallel do default(shared), private(ip,i,io,lio,ind,j)
      do i = 1 , mols(im)%orb%n
        
        ! Initialize |><| value
        zD(i,:) = dcmplx(0._dp,0._dp)
        zP(i,:) = dcmplx(0._dp,0._dp)
        io = mols(im)%orb%r(i)
        lio = index_global_to_local(dit,io)

        ! Calculate <| M
        do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)
          j = rgn_pivot(mols(im)%orb,l_col(ind))
          if ( j == 0 ) cycle
          ! this is per level in the system
          do ip = 1 , Npl
            zD(i,ip) = zD(i,ip) + dconjg(mols(im)%p(j,poff+ip)) * M(ind)
          end do
          ! this is per projection calculating |\sum>_i<\sum| M
          j = ((i-1)*Nsl+ip-1)*no+j
          do ip = 1 , Nsl
            zP(i,ip) = zP(i,ip) + zwork(j) * M(ind)
          end do
        end do

        ! Calculate remaining |> and take the trace
        do ip = 1 , Npl
          zD(i,ip) = mols(im)%p(i,poff+ip) * zD(i,ip)
        end do

      end do
!$OMP end parallel do

      ! Save the diagonal of the projection S
      call ncdf_put_var(grp,'kb_SD',zD,start=(/1,1,ikpt/))
#ifdef MPI
      if ( Node == 0 ) then
        ! Offset for master node
        poff = Npl + 1
        do i = 1 , Nodes - 1
          if ( poff > Np ) exit
          ! We know that the master node ALWAYS have more or 
          ! an equal amount of processing states. Hence we can re-use
          ! the array
          call MPI_Recv(zD(1,1),no*Npl,MPI_Double_Complex,i,i, &
              MPI_Comm_World, status, MPIerror)
          call MPI_Get_Count(status,MPI_Double_Complex,ind,MPIerror)
          ! calculate remote Npl
          j = ind / no
          call ncdf_put_var(grp,'kb_SD',zD(:,1:j),start=(/1,poff,ikpt/))
          poff = poff + j
        end do
      else if ( Npl > 0 ) then
        call MPI_Send(zD(1,1),no*Npl,MPI_Double_Complex,0,Node, &
            MPI_Comm_World, MPIerror)
      end if
#endif

      if ( Node == 0 ) then

        ! Save projected overlaps of full states
        i = 0
        soff = 0
        do ip = 1 , Ns

          call ncdf_open_grp(grp,mols(im)%proj(ip)%name,grp2)

          ! Save the diagonal of the projection S
          call ncdf_put_var(grp2,'kb_SD',zP(:,soff+ip),start=(/1,ikpt/))
#ifdef MPI
          ! If there are no more projections on the master
          ! we receive the next batch
          if ( soff + ip == Nsl .and. ip < Ns ) then
            ! Calculate zP offset
            soff = - ip 
            ! step processor
            i = i + 1
            ! Get next segment
            call MPI_Recv(zP(1,1),no*Nsl,MPI_Double_Complex,i,i, &
                MPI_Comm_World, status, MPIerror)
            call MPI_Get_Count(status,MPI_Double_Complex,ind,MPIerror)
            ! calculate remote Nsl
            Nsl = ind / no
          end if
#endif

        end do

#ifdef MPI
      else if ( Nsl > 0 ) then
        call MPI_Send(zP(1,1),no*Nsl,MPI_Double_Complex,0,Node, &
            MPI_Comm_World, MPIerror)
#endif
      end if

      deallocate(zD,zP)

    end do

    call ncdf_close(ncdf)

  end subroutine proj_cdf_save_S_D

  ! Save <|Gamma|>
  subroutine proj_cdf_save_bGammak(ncdf,N_proj_ME,proj_ME,ikpt,nE)
    
    use class_OrbitalDistribution
    use class_Sparsity
    use class_zSpData1D

    use parallel, only : Node, Nodes

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only: MPI_Send, MPI_Recv, MPI_Get_Count
    use mpi_siesta, only: MPI_STATUS_SIZE, MPI_STATUSES_IGNORE
    use mpi_siesta, only: MPI_REQUEST_NULL
    use mpi_siesta, only: MPI_Double_Complex, MPI_Comm_World
#endif

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: N_proj_ME
    type(tProjMolEl), intent(inout) :: proj_ME(N_proj_ME)
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE

    ! Local variables
    type(hNCDF) :: gmol
    integer :: iE, nl
    character(len=NF90_MAX_NAME) :: cmol, ctmp
    integer :: idx(4), cnt(4)
#ifdef MPI
    integer :: reqs(N_proj_ME)
    complex(dp), allocatable :: tmp(:)
    integer :: MPIerror, Status(MPI_STATUS_SIZE), iN
#endif

    cmol  = ' '

    if ( save_parallel ) then
      
      idx(:) = 1
      idx(4) = ikpt


      ! Easy storage
      do iE = 1 , N_proj_ME

        ! Number of levels on this projection
        nl = proj_ME(iE)%mol%lvls%n

        ! Reset count
        cnt(1) = nl
        cnt(2) = nl
        cnt(3) = 1
        cnt(4) = 1

        if ( cmol /= proj_ME(iE)%mol%name ) then
          cmol = proj_ME(iE)%mol%name
          call ncdf_open_grp(ncdf,cmol,gmol)
        end if

        ctmp = trim(proj_ME(iE)%El%name)//'.bGk'

        ! Now we can save the data

        ! Determine index and count of storages
        if ( nE%iE(Node) > 0 ) then
          idx(3) = nE%iE(Node)
        else
          ! first index, we *need* to make sure the
          ! position exists
          idx(3) = 1
          ! Tell to not store any data
          cnt = 0
        end if

        ! In the code bGk is _without_ factor "i".
        ! Hence, we here add factor i
        proj_ME(iE)%bGk = proj_ME(iE)%bGk * dcmplx(0._dp, 1._dp)

        ! ALL nodes _have_ to participate
        call ncdf_put_var(gmol,ctmp,proj_ME(iE)%bGk, &
            start = idx, count=cnt )

        ! and back
        proj_ME(iE)%bGk = proj_ME(iE)%bGk * dcmplx(0._dp, -1._dp)

      end do

    else
      
#ifdef MPI
      ! Find the maximum number of levels
      nl = 0
      do iE = 1 , N_proj_ME
        nl = max(nl,proj_ME(iE)%mol%lvls%n)
      end do
      allocate(tmp(nl*nl))
#endif

      do iE = 1 , N_proj_ME

        ! Number of levels on this projection
        nl = proj_ME(iE)%mol%lvls%n

#ifdef MPI
        if ( Node /= 0 ) then
          reqs(iE) = MPI_REQUEST_NULL
        end if
#endif

        if ( cmol /= proj_ME(iE)%mol%name ) then
          cmol = proj_ME(iE)%mol%name
          call ncdf_open_grp(ncdf,cmol,gmol)
        end if

        ctmp = trim(proj_ME(iE)%El%name)//'.bGk'
        ! Now we can save the data

        if ( nE%iE(Node) > 0 ) then
          ! In the code bGk is _without_ factor "i".
          ! Hence, we here add factor i
          proj_ME(iE)%bGk = proj_ME(iE)%bGk * dcmplx(0._dp, 1._dp)
          call ncdf_put_var(gmol,ctmp,proj_ME(iE)%bGk, &
              start = (/1,1,nE%iE(Node),ikpt/) )
          proj_ME(iE)%bGk = proj_ME(iE)%bGk * dcmplx(0._dp, -1._dp)
        end if
#ifdef MPI
        if ( Node == 0 ) then
          do iN = 1 , Nodes - 1
            if ( nE%iE(iN) <= 0 ) cycle
            call MPI_Recv(tmp,nl*nl,Mpi_double_complex, &
                iN, iN, Mpi_comm_world,status,MPIerror)
            tmp = tmp * dcmplx(0._dp, 1._dp)
            call ncdf_put_var(gmol,ctmp,reshape(tmp,(/nl,nl/)), &
                start = (/1,1,nE%iE(iN),ikpt/) )
          end do
        else if ( nE%iE(Node) > 0 ) then
          call MPI_ISend(proj_ME(iE)%bGk,nl*nl,Mpi_double_complex, &
              0, Node, Mpi_comm_world,reqs(iE),MPIerror)
        end if
#endif

      end do

#ifdef MPI
      if ( Node /= 0 ) then
        call MPI_WaitAll(N_proj_ME,reqs(1),MPI_STATUSES_IGNORE,MPIerror)
      end if
      deallocate(tmp)
#endif

    end if
    
  end subroutine proj_cdf_save_bGammak

  subroutine proj_cdf2ascii(fname,N_proj_T,proj_T,save_DATA)
    
    use units, only : eV
    use variable
    use dictionary

    use m_timestamp, only : datestring
    use netcdf_ncdf
    use m_ts_electype

    character(len=*), intent(in) :: fname
    integer, intent(in) :: N_proj_T
    type(tProjT), intent(in) :: proj_T(N_proj_T)
    type(dictionary_t), intent(in) :: save_DATA

  end subroutine proj_cdf2ascii

  subroutine proj_Mt_mix(mol,ip,Mt,bGk)
    type(tProjMol), intent(in) :: mol
    integer, intent(in) :: ip ! projection index
    complex(dp), intent(out) :: Mt(mol%orb%n,mol%orb%n)
    complex(dp), intent(in) :: bGk(mol%lvls%n,mol%lvls%n)
    complex(dp) :: p(mol%orb%n), tmp(mol%orb%n)
    integer :: i, j, gi, gj

    if ( ip == 0 ) then
      call die('Error in programming, proj_Mt_mix')
    end if

    ! loop over number of levels associated with
    ! this projection
    do j = 1 , mol%proj(ip)%n
      gj = mol%proj(ip)%r(j)

      p(:) = dcmplx(0._dp,0._dp)
      do i = 1 , mol%proj(ip)%n
        gi = mol%proj(ip)%r(i)
        ! Create summation |i> . <i|Gam|j>
        p(:) = p(:) + mol%p(:,gi) * bGk(gi,gj)
      end do

      ! Do last product |i> . <i|Gam|j> . <j|
      ! and take the transpose
      tmp(:) = dconjg(mol%p(:,gj))
      if ( j == 1 ) then
        do i = 1 , mol%orb%n
          Mt(:,i) = p(i) * tmp(:)
        end do
      else
        do i = 1 , mol%orb%n
          Mt(:,i) = Mt(:,i) + p(i) * tmp(:)
        end do
      end if

    end do
    
  end subroutine proj_Mt_mix

  ! Projects the projection state onto a transposed matrix
  subroutine proj_bMtk(mol,orb,Mt,bMk,nwork,work)
    type(tProjMol), intent(in) :: mol
    type(tRgn), intent(in) :: orb ! The orbitals of the current matrix
    complex(dp), intent(in) :: Mt(orb%n,orb%n)
    complex(dp), intent(out) :: bMk(mol%lvls%n,mol%lvls%n)
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)
    complex(dp), pointer :: tmp(:), pl(:)
    integer :: i, j
    complex(dp), external :: zdotc
    
    if ( nwork < orb%n*(1+mol%lvls%n) ) then
      call die('Projection proj_bMtk, not enough work space.')
    end if

    do j = 1 , mol%lvls%n
      ! point to the |j>
      pl => work((j-1)*orb%n+1:j*orb%n)
      call proj_sort(orb,mol,j,pl)
    end do
    j = mol%lvls%n + 1
    tmp => work((j-1)*orb%n+1:j*orb%n)

    ! |j>
    do j = 1 , mol%lvls%n

      pl => work((j-1)*orb%n+1:j*orb%n)

      ! Note that Mt is a transposed matrix, hence we need to 
      ! transpose back
      call zgemv('T',orb%n,orb%n,dcmplx(1._dp,0._dp),Mt(1,1),orb%n, &
          pl,1,dcmplx(0._dp,0._dp),tmp,1)

      ! <i|
      do i = 1 , mol%lvls%n
        pl => work((i-1)*orb%n+1:i*orb%n)
        bMk(i,j) = zdotc(orb%n,pl,1,tmp,1)
      end do

    end do
    
  end subroutine proj_bMtk

  ! Calculates the projection vector for
  ! the projection levels
  ! |1>_j * <1| + |2>_j * <2|  for all states in the projector
  subroutine proj_state_bra(mol,proj,j,p)
    type(tProjMol), intent(in) :: mol
    type(tRgn), intent(in) :: proj
    integer, intent(in) :: j ! column j (corresponds to the index of p)
    complex(dp), intent(out) :: p(mol%orb%n)
    
    integer :: ip, i

    p(:) = dcmplx(0._dp,0._dp)
    do ip = 1 , proj%n
      i = proj%r(ip)
      !         add |ip>_j <ip|
      p(:) = p(:) + mol%p(j,i) * dconjg(mol%p(:,i))
    end do
    
  end subroutine proj_state_bra

  ! Takes a consecutive projection and sorts
  ! it to a region:
  ! Consider a projection on the orbitals:
  !   [1,2,3,4]
  ! consider now a region made up of orbitals in this order:
  !   [2,4,1]
  ! then this will return the projection:
  !   [2,4,1]
  ! The region MUST be a subset of the molecule orbitals
  subroutine proj_sort(r,mol,il,psort)
    ! The region on which we wish to create
    ! an aligned projector.
    type(tRgn), intent(in) :: r
    ! The molecule
    type(tProjMol), intent(in) :: mol
    ! The level that we want to create the ket of
    integer, intent(in) :: il
    ! The full projector aligned to the requested projector
    complex(dp), intent(out) :: psort(r%n)

    ! Local variables
    integer :: i

    do i = 1 , r%n
      ! Copy over the projector to the assigned index
      psort(i) = mol%p(rgn_pivot(mol%orb,r%r(i)),il)
    end do
    
  end subroutine proj_sort

  subroutine read_proj_options( save_DATA )
    
    use fdf

    use dictionary

    type(dictionary_t), intent(inout) :: save_DATA

    ! Local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    logical :: ltmp
    integer :: i

    ! Reset the projection molecules
    N_mol = 0
    
    ! If the user has requested only to calculate the 
    ! self-energies, we should not read in the projections.
    if ( ('Sigma-only'.in.save_DATA) ) return

    ! If the projection block exists
    ! it means the user is requesting projections.
    if ( .not. fdf_defined('TBT.Projs') ) return

    if ( .not. fdf_block('TBT.Projs',bfdf) ) then
      call die('TBT.Projs is not a block, please correct')
    end if

    ! First read in the molecules
    do while ( fdf_bline(bfdf,pline) )
      ! skip empty line
      if ( fdf_bnnames(pline) == 0 ) cycle
      N_mol = N_mol + 1
    end do
    allocate(mols(N_mol))
    
    ! rewind to read again
    call fdf_brewind(bfdf)

    ! Retrieve the names
    i = 0
    do while ( fdf_bline(bfdf,pline) )
      ! empty line
      if ( fdf_bnnames(pline) == 0 ) cycle
      i = i + 1
      mols(i)%name = fdf_bnames(pline,1)
      if ( index(mols(i)%name,'.') > 0 ) then
        call die('Projections cannot be named with .!')
      end if
    end do


    ! Whether we should assert and calculate
    ! all transmission amplitudes for the projections
    ltmp = fdf_get('TBT.Projs.T.Elecs.All', ('T-all'.in.save_DATA) )
    ltmp = fdf_get('TBT.Projs.T.All', ltmp )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-T-all'.kv.1)
    end if
    ltmp = fdf_get('TBT.Projs.T.Out', ('T-sum-out'.in.save_DATA) )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-T-sum-out'.kv.1)
    end if
    ltmp = fdf_get('TBT.Projs.Only', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-only'.kv.1)
    end if
    
    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.Projs.DOS.A', ('DOS-A'.in.save_DATA))
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-DOS-A'.kv.1)
    end if

    ltmp = fdf_get('TBT.Projs.Current.Orb', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-DOS-A'.kv.1)
      save_DATA = save_DATA // ('proj-orb-current'.kv.1)
    end if

    ltmp = fdf_get('TBT.Projs.DM.A', .false.)
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-DOS-A'.kv.1)
      save_DATA = save_DATA // ('proj-DM-A'.kv.1)
    end if

    ltmp = fdf_get('TBT.Projs.COOP.A', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-DOS-A'.kv.1)
      save_DATA = save_DATA // ('proj-COOP-A'.kv.1)
    end if

    ltmp = fdf_get('TBT.Projs.COHP.A', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('proj-DOS-A'.kv.1)
      save_DATA = save_DATA // ('proj-COHP-A'.kv.1)
    end if

  end subroutine read_proj_options
  
  subroutine print_proj_options( save_DATA )
    
    use parallel, only: IONode
    use dictionary
    
    type(dictionary_t), intent(inout) :: save_DATA

    character(len=*), parameter :: f1 ='(''tbt-proj: '',a,t53,''='',tr4,l1)'

    if ( .not. IONode ) return
    if ( N_mol == 0 ) return

    write(*,f1) 'Calc. T between all electrodes',('proj-T-all'.in.save_DATA)
    write(*,f1) 'Calc. total T out of electrodes',('proj-T-sum-out'.in.save_DATA)
    write(*,f1) 'Saving DOS from spectral functions',('proj-DOS-A' .in. save_DATA)
    write(*,f1) 'Saving bond currents (orb-orb)',('proj-orb-current'.in.save_DATA)
    write(*,f1) 'Saving DM from spectral functions',('proj-DM-A'.in.save_DATA)
    write(*,f1) 'Saving COOP from spectral functions',('proj-COOP-A'.in.save_DATA)
    write(*,f1) 'Saving COHP from spectral functions',('proj-COHP-A'.in.save_DATA)

  end subroutine print_proj_options
#else
  
  subroutine read_proj_options( save_DATA )
    use dictionary
    type(dictionary_t), intent(inout) :: save_DATA
  end subroutine read_proj_options
  
  subroutine print_proj_options( save_DATA )
    use dictionary
    type(dictionary_t), intent(inout) :: save_DATA
  end subroutine print_proj_options
#endif

end module m_tbt_proj
