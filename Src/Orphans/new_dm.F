! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_new_dm

!     Prepares a starting density matrix for a new geometry iteration
!     This DM can be:
!     1. Synthesized directly from atomic occupations (not idempotent)
!     2. Read from file
!     3. Extrapolated from (two) previous steps
!     3.a The DM of the previous iteration
!
!     In cases 2 and 3, a check is done to guarantee that the structure
!     of the read or extrapolated DM conforms to the current sparsity.
!     If it does not, the information is re-arranged.
!
!     Special cases:
!            Harris: The matrix is always initialized
!            Force calculation: The DM should be written to disk
!                               at the time of the "no displacement"
!                               calculation and read from file at
!                               every subsequent step.
!            Variable-cell calculation:
!              If the auxiliary cell changes, the DM is forced to be
!              initialized (conceivably one could rescue some important
!              information from an old DM, but it is too much trouble
!              for now). NOTE that this is a change in policy with respect
!              to previous versions of the program, in which a (blind?)
!              re-use was allowed, except if 'ReInitialiseDM' was 'true'.
!              Now 'ReInitialiseDM' is 'true' by default. Setting it to
!              'false' is not recommended. (This fdf variable maps to the
!              'initdmaux' module variable)
!
!              In all other cases (including "server operation"), the
!              default is to allow DM re-use (with possible extrapolation)
!              from previous geometry steps.
!              The fdf variables 'DM.AllowReuse' and 'DM.AllowExtrapolation'
!              (mapped to 'allow_dm_reuse' and 'allow_dm_extrapolation', and
!              both 'true' by default) can be used to change this behavior.
!
!              There is no re-use of the DM for "Forces", and "Phonon"
!              dynamics types (i.e., the DM is re-initialized)
!
!              For "CG" calculations, the default is not to extrapolate the
!              DM (unless requested by setting 'DM.AllowExtrapolation' to
!              "true"). The previous step's DM is reused.
!
!      For the purposes of extrapolation, this module keeps DMsaved, which
!      holds the DM from a previous geometry step (used, together with the
!      current DM to extrapolate).
!      The module also keeps the sparsity at the time of initialization,
!      so that it can be compared to the current one. If restructuring is
!      needed, the saved sparsity is updated.
!
!     Alberto Garcia, September 2007
!
      use m_sparse, only: same_sparsity, change_sparsity
      use sys, only: die
      use precision, only: dp
      use alloc, only: re_alloc, de_alloc
      use parallel,  only: IOnode

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      character(len=*),parameter:: modName = 'm_new_dm '

      ! Private variables to implement DM extrapolation and reuse

      ! If extrapolation is allowed, DMsaved stores the DM of
      ! previous geometry iteration
      real(dp), pointer, private :: DMsaved(:,:) => null()
      !
      ! Even if extrapolation is not allowed, these variables
      ! store the sparsity pattern of the previous geometry
      ! iteration.
      integer, pointer, private  :: listold(:) => null()
      integer, pointer, private  :: listptrold(:) => null()
      integer, pointer, private  :: numold(:) => null()

      private
      public :: new_dm

      CONTAINS

      subroutine  new_dm( auxchanged, numh, listhptr, listh, Dscf)

      USE siesta_options
      use siesta_geom
      use atomlist, only: datm, no_s, iaorb, lasto, no_u, no_l
      use m_spin, only: nspin, h_spin_dim, e_spin_dim
      use m_steps, only: istp
      USE m_normalize_dm, only: normalize_dm
      use sparse_matrices, only: Escf     !! To be moved
#ifdef TRANSIESTA
      use m_ts_global_vars, only : TSrun
#endif /* TRANSIESTA */

      implicit none

      logical, intent(in) :: auxchanged ! Has auxiliary supercell changed?
      integer, intent(in) :: numh(:)
      integer, intent(in) :: listhptr(:), listh(:)
      real(dp), pointer   :: Dscf(:,:)

!     Local variables

      logical :: dminit     ! Initialize density matrix?
      logical :: try_to_read_from_file
      integer :: nnz

      if (IOnode) then
         write(6,"(a,i5)") "New_DM. Step: ", istp
      endif


!     In principle we allow the re-use of the DM (i.e, we do not initialize it)
!
      dminit = .false.
      try_to_read_from_file = usesavedm     ! As per defaults
!
!     Except if there are explicit instructions
!
      if (.not. allow_dm_reuse) then
         dminit = .true.
         try_to_read_from_file = .false.  ! In case the user has a fossil DM.UseSaveDM
         if (IOnode) then
            write(6,"(a)") "DM re-use not allowed. Resetting always"
            if (usesavedm) then
               write(6,"(a)") "DM.UseSaveDM  overriden !!"
            endif
         endif
      endif
!
!     or using Harris...
!
      if (harrisfun) dminit = .true.
!
!     or we are in the first step, or performing force-constant calculations
!
      if (istp .eq. 1) then
         dminit = .true.
      else
         if ((idyn .eq. 6) ! Force Constants
     $         .and. usesavedm .and. writedm)  dminit = .true.
         if ((idyn .eq. 7) ! Phonon series (writedm??)
     $         .and. usesavedm)  dminit = .true.
      endif

!
!     ... or if the auxiliary cell has changed
!     (in this case we have to  avoid reading back saved copy from file)
!
      if (initdmaux.and.auxchanged) then
        dminit = .true.
        try_to_read_from_file = .false.
        if (IOnode) then
           write(6,"(a)") "DM reset as supercell changed."
        endif
      endif


      if (dminit) then
         if (IOnode) then
            write(6,"(a)") "Initializing Density Matrix..."
         endif

         ! This is the point of resizing of the DM
         ! to the bounds appropriate to the new sparsity
         nnz=sum(numh)
         call re_alloc(Dscf,1,nnz,1,h_spin_dim,name="Dscf",
     $                 routine="sparseMat",copy=.false.,shrink=.true.)

         ! Reset DMsaved to null
         call de_alloc(DMsaved,name="DMsaved",routine="new_dm")

         call initdm(Datm, Dscf, lasto, na_s,
     .               no_s, no_l, h_spin_dim, na_u, no_l, 
     .               numh, listhptr, listh, iaorb, inspn,
     .               try_to_read_from_file)
         !
         ! Update saved sparsity pattern.
         ! no_l should not have changed, but in the future
         ! one might conceivably implement a domain decomposition,
         ! and that would involve re-checking this issue
         call re_alloc(numold,1,no_l)
         call re_alloc(listptrold,1,no_l)
         ! NOTE possible shrinking, as the test for changed
         ! sparsity might involve the length of the listh array
         call re_alloc(listold,1,size(listh),shrink=.true.)
         numold = numh
         listold = listh
         listptrold = listhptr


      else    ! not initializing the DM

         if (IOnode) then
            write(6,"(a)") "Re-using DM from previous geometry..."
         endif

        ! Extrapolate density matrix between steps
         if (idyn == 0) then
            ! do not extrapolate for CG relaxation, in case
            ! there are big jumps
            allow_dm_extrapolation = fdf_get(
     $           'DM.AllowExtrapolation', .FALSE.)
            if (allow_dm_extrapolation) then
               ! The user insisted...
               if (ionode) then
                  write(6,"(a)")
     $                 "NOTE: DM.AllowExtrapolation set " //
     $                 "with CG dynamics option"
                  write(0,"(a)")
     $                 "NOTE: DM.AllowExtrapolation set " //
     $                 "with CG dynamics option"
               endif
            endif
         endif

         call new_extrapol(no_l, numh, listhptr, listh, Dscf)

      endif

      ! Initialize energy-density matrix to zero for first call to overfsm
#ifdef TRANSIESTA
      ! Only part of Escf is updated in TS, so if it is put as zero here
      ! a continuation run gives bad forces.
      if(.not.TSrun) Escf(:,:) = 0.0_dp
#else
      Escf(:,:) = 0.0_dp
#endif /* TRANSIESTA */

!     Normalize density matrix to exact charge
#ifdef TRANSIESTA      
      if (.not.TSrun ) call normalize_dm( first= .true. )
#else
      call normalize_dm(first=.true.)
#endif  /* TRANSIESTA */

      END subroutine new_dm

      subroutine new_extrapol(nrows,num, listptr,list, DM)

      USE siesta_options

      ! This routine extrapolates the density matrix DM
      ! If it is called for the first time after initialization
      ! it just passes DM through, possibly re-structured.

      ! On INPUT, DM conforms to the previous geometry iteration
      ! sparsity pattern (stored in module variables)
      ! On OUTPUT, DM conforms to the current geometry iteration
      ! sparsity pattern

      integer, intent(in) :: nrows
      integer, intent(in) :: num(nrows), listptr(nrows)
      integer, intent(in) :: list(:)
      real(dp), pointer   :: DM(:,:)


      real(dp), pointer   :: DMtmp(:,:) => null()
      logical  :: sparsity_has_changed
      real(dp) :: msave
      integer  :: i, j, ispin, h_spin_dim, nnz, new_nnz, is

#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif

      sparsity_has_changed =
     $  (.not.  same_sparsity (nrows,numold,listold, num,list))

#ifdef MPI
      call MPI_AllReduce(sparsity_has_changed,lbuffer,1,MPI_logical,
     $                   MPI_lor,
     $                   MPI_Comm_World,MPIerror)
      sparsity_has_changed = lbuffer
#endif

      ! Here is an interesting thought:
      ! Can we do the re-structuring node by node, only on those
      ! in which it is necessary? I very much think so, so the
      ! above block is superfluous.

      ! In future it might be possible that no_l changes during
      ! the calculation (i.e., because of the implementation of
      ! orbital distribution by domain decomposition).
      ! In that case we might simply request an initialization of
      ! the density matrix, or write a super extrapolator.

      if (sparsity_has_changed) then

         ! Fix structure of DM
         ! This section now uses a temporary matrix and a safer
         ! restructuring routine.

         h_spin_dim = size(DM,dim=2)
         nnz = size(DM,dim=1)      ! This should be sum(numold)
         call re_alloc(DMtmp,1,nnz,1,h_spin_dim,
     $                 name="DMtmp",routine="new_dm", shrink=.true.)
         DMtmp(:,:) = DM(:,:)
         new_nnz = sum(num(:))
         call re_alloc(DM,1,new_nnz,1,h_spin_dim,
     $                 name="DM",routine="new_dm",shrink=.true.)

         do is = 1, h_spin_dim
            call change_sparsity(nrows,
     $                           numold,listptrold,listold,
     $                           num, listptr, list,
     $                           DMtmp(:,is), DM(:,is) )
         enddo
         call de_alloc(DMtmp,name="DMtmp",routine="new_dm")

      endif

      ! At this point DM conforms already to the current sparsity pattern.

      if (associated(DMsaved)) then

         if (IOnode) then
            write(6,"(a)") "Extrapolating Density Matrix..."
         endif

         if (sparsity_has_changed) then
            !Fix structure of DMsaved
            h_spin_dim = size(DMsaved,dim=2)
            nnz = size(DMsaved,dim=1)
            call re_alloc(DMtmp,1,nnz,1,h_spin_dim,
     $           name="DMtmp",routine="new_dm", shrink=.true.)
            DMtmp(:,:) = DMsaved(:,:)
            new_nnz = sum(num(:))
            call re_alloc(DMsaved,1,new_nnz,1,h_spin_dim,
     $           name="DMsaved",routine="new_dm",shrink=.true.)

            do is = 1, h_spin_dim
               call change_sparsity(nrows,
     $                              numold,listptrold,listold,
     $                              num, listptr, list,
     $                              DMtmp(:,is), DMsaved(:,is) )
            enddo
            call de_alloc(DMtmp,name="DMtmp",routine="new_dm")

         endif

         ! Extrapolate DM and update DMsaved
         ! This is a simple linear extrapolation

         h_spin_dim = size(DM,dim=2)
         nnz = size(DM,dim=1)

         do ispin = 1, h_spin_dim
            do i = 1, nnz
               msave = DM(i,ispin)
               DM(i,ispin) = 2.0_dp * DM(i,ispin) - DMsaved(i,ispin)
               DMsaved(i,ispin) = msave
            enddo
         enddo

      else    ! not associated(DMsaved)

         if (IOnode) then
            if (allow_dm_extrapolation) then
               write(6,"(a)") "Re-using DM without extrapolation " //
     $                        "for lack of history"
            else
               write(6,"(a)") "Re-using DM without extrapolation " //
     $                        "(not allowed)"
            endif
         endif

         if (allow_dm_extrapolation) then
            ! Copy DM to DMsaved for next iteration
            h_spin_dim = size(DM,dim=2)
            nnz = size(DM,dim=1)
            call re_alloc(DMsaved,1,nnz,1,h_spin_dim,name="DMsaved",
     $                    routine="new_dm",copy=.false.,shrink=.true.)
            DMsaved(1:nnz,1:h_spin_dim) = DM(1:nnz,1:h_spin_dim)
         endif
      endif

      if (sparsity_has_changed) then

         if (IOnode) then
            write(6,"(a)") "Density Matrix sparsity pattern changed."
         endif

         ! Store sparsity
         call re_alloc(numold,1,size(num),shrink=.true.)
         call re_alloc(listptrold,1,size(listptr),shrink=.true.)
         call re_alloc(listold,1,size(list),shrink=.true.)
         numold = num
         listold = list
         listptrold = listptr

      endif

      end subroutine new_extrapol
      subroutine initdm(Datm, Dscf, lasto, maxa,
     .                  maxo, maxuo, nspin, nua, no_l, 
     .                  numh, listhptr, listh, iaorb, inspn,
     .                  try_dm_from_file)

c *******************************************************************
c Density matrix initialization
c
c    If Try_Dm_From_File is true, it is read from file if present.
c    Otherwise it is generated assuming atomic charging
c      (filling up atomic orbitals). The DM originated that way is
c      not a good DM due to overlaps, but the SCF cycling corrects
c      that for the next cycle.
c    Spin polarized calculations starting from atoms:
c      Default: All atoms with maximum polarization compatible with
c               atomic configuration. In Ferromagnetic ordering (up).
c      If DM.InitSpinAF is true, as default but in Antiferro order:
c               even atoms have spin down, odd up.
c      If fdf %block DM.InitSpin is present it overwrites previous
c         schemes: magnetic moments are explicitly given for some atoms.
c         Atoms not mentioned in the block are initialized non polarized.
c
c Written by E. Artacho. December 1997. Taken from the original piece
c of siesta.f written by P. Ordejon.
c Non-collinear spin added by J.M.Soler, May 1998.
c ********* INPUT ***************************************************
c logical try_dm_from_file     : whether DM has to be read from files or not
c logical found         : whether DM was found in files
c logical inspn         : true : AF ordering according to atom ordering
c                                if no DM files, no DM.InitSpin, ispin=2
c                         false: Ferro ordering  (fdf DM.InitSpinAF)
c integer nua           : Number of atoms in the unit cell
c integer no_l          : Number of orbitals in the unit cell
c integer nspin         : Number of spin components in density matrix
c integer maxa          : Max num. atoms for dimension
c integer maxo          : Max. number of orbitals (globally)
c integer maxuo         : Max. number of orbitals (locally)
c integer lasto(0:maxa) : List with last orbital of each atom
c integer numh(:)       : Dscf matrix sparse information
c integer listhptr(:)   : Pointer in CSR sparsity matrix
c integer listh(:)      : Column index of equivalent orbital
c integer iaorb(maxo)   : List saying to what atom an orbital belongs
c double Datm(no)       : Occupations of basis orbitals in free atom
c ********* OUTPUT **************************************************
c double Dscf(:,h_spin_dim=nspin) : Density matrix in sparse form
c *******************************************************************

C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes, IOnode
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use fdf
      use parsing
      use sys,          only : die
      use m_iodm,       only : read_dynamic_dm
      use m_sparse,     only : same_sparsity, change_sparsity
      use alloc,        only : re_alloc, de_alloc
#ifdef TRANSIESTA
      use sparse_matrices, only: Escf
      use m_ts_iodm, only    : read_dynamic_ts_dm, ts_init_dm
      use m_energies, only: ef  ! Transiesta uses the EF obtained in a initial SIESTA run
                                ! to place the electrodes and scattering region energy
                                ! levels at the appropriate relative position, so it is
                                ! stored in the TSDE file.
      use m_ts_options,   only : TSmode
#endif /* TRANSIESTA */

#ifdef MPI
      use mpi_siesta
#endif
      use units, only : pi
      use sparse_matrices, only: S

      implicit          none

      logical           found, inspn, try_dm_from_file
      integer           no_l, nua, maxo, maxuo, h_spin_dim, maxa,
     .                  nspin
      integer           lasto(0:maxa), numh(maxuo), iaorb(maxo)
      integer       ::  listhptr(:), listh(:)
      real(dp), pointer :: Dscf(:,:)
      real(dp)          Datm(maxo)
      external          memory

c ---------------------------------------------------------------------

C Internal variables and arrays
      character(len=*),parameter:: myName = 'initdm'
      character         updo*1, msg*80
      logical           noncol, peratm, badsyntax
      integer           nh, ni, nn, nr, nv, iat, nat, ia, iu,
     .                  i1, i2, in, ind, ispin, jo, io,
     .                  iio, maxatnew

      integer, save ::  maxat

      integer           integs(4), lastc, lc(0:3)

      integer, pointer, save ::  atom(:)
#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif
      real(dp)          aspin, cosph, costh, epsilon,
     .                  qio, rate, reals(4),
     .                  sinph, sinth, spinat, spio, values(4)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp), pointer, save :: phi(:), spin(:), theta(:)

      integer, pointer :: list_tmp(:), num_tmp(:), listptr_tmp(:)
      real(dp), pointer :: D_tmp(:,:), work(:)
#ifdef TRANSIESTA
      logical :: ts_found
      real(dp), pointer :: EDM_tmp(:,:)
#endif /* TRANSIESTA */

      integer :: nh_tmp, h_spin_dim_read, is
      logical :: same_sp

      data maxat / 1000 /
      data epsilon / 1.d-8 /

      h_spin_dim = size(Dscf,dim=2)

C Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------
      nullify(list_tmp,D_tmp,num_tmp,listptr_tmp)
#ifdef TRANSIESTA
      nullify(EDM_tmp)
#endif /* TRANSIESTA */

      if (try_dm_from_file) then

         call re_alloc(num_tmp,1,no_l)
         call re_alloc(listptr_tmp,1,no_l)

#ifndef TRANSIESTA
         call read_dynamic_dm(nh_tmp, no_l, h_spin_dim_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, found )
#else
         found = .false.   ! initialize
         if(TSmode) then ! Transiesta solution method
C First try to read the .TSDE file
           call read_dynamic_ts_dm(nh_tmp, no_l, h_spin_dim_read, 
     .             num_tmp, listptr_tmp, list_tmp, D_tmp, EDM_tmp, 
     .             ef, ts_found )

C Print information and do the some TS initializations
           call ts_init_dm(ts_found)

C If it does not find a .TSDE file, try to read a .DM file
           if(.not. ts_found) 
     .        call read_dynamic_dm(nh_tmp, no_l, h_spin_dim_read, 
     .             num_tmp, listptr_tmp, list_tmp, D_tmp, found )     
C Either .TSDE or .DM may have been found
           found = ts_found .or. found
         else ! Not Transiesta solution mathod
           ts_found = .false.
           call read_dynamic_dm(nh_tmp, no_l, h_spin_dim_read, num_tmp,
     $             listptr_tmp, listd=list_tmp, dm=D_tmp, found=found )
         end if ! Transiesta solution method
#endif /* TRANSIESTA */

      else
        found = .false.
      endif

C If found, check and update, otherwise initialize with neutral atoms

      if (found) then
        ! Various degrees of sanity checks
        if (h_spin_dim_read /= h_spin_dim) then
           if (Node.eq.0) then
              write(6,"(a,i6,/,a)")
     $        "WARNING: Wrong h_spin_dim in DM file: ",
     $             h_spin_dim_read,
     $        "WARNING: Falling back to atomic initialization of DM."
              write(0,"(a,i6,/,a)")
     $        "WARNING: Wrong h_spin_dim in DM file: ",
     $             h_spin_dim_read,
     $        "WARNING: Falling back to atomic initialization of DM."
           endif
           found = .false.
        endif
      endif
      
      if (found) then

         same_sp = same_sparsity(no_l,num_tmp,list_tmp,
     $                                 numh,listh)
#ifdef MPI
        call MPI_AllReduce(same_sp,lbuffer,1,MPI_logical,MPI_land,
     .    MPI_Comm_World,MPIerror)
        same_sp = lbuffer
#endif
        if (same_sp) then
           Dscf(:,:) = D_tmp(:,:)
#ifdef TRANSIESTA
           if(ts_found) Escf(:,:) = EDM_tmp(:,:)
#endif /* TRANSIESTA */
        else
           if (IOnode) then
              write(6,*) "Read DM has different structure. Fixing..."
           endif
           nullify(work)
           do is= 1, h_spin_dim
              call change_sparsity (no_l,num_tmp,listptr_tmp,list_tmp,
     $                                  numh,listhptr,listh,
     $                                  D_tmp(:,is), Dscf(:,is),
     $                                  work=work)
#ifdef TRANSIESTA
! change_sparsity allocates work, but should work, since it does not 
! deallocate inside the spin loop
              if(ts_found.and.is<=e_spin_dim) 
     .          call change_sparsity (no_l,num_tmp,listptr_tmp,list_tmp,
     $                                  numh,listhptr,listh,
     $                                  EDM_tmp(:,is), Escf(:,is),
     $                            work)
               
#endif /* TRANSIESTA */
           enddo
           call de_alloc(work)  ! Allocated by change_sparsity
        endif


      else  ! Initialize with neutral atoms

C See whether specific initial spins are given in a DM.InitSpin block
C and read them in a loop on atoms where lines are read and parsed
C   integer nat       : how many atoms to polarize
C   integer atom(nat) : which atoms
C   double  spin(nat) : what polarization -----------------------------

        noncol = .false.
        peratm = fdf_block('DM.InitSpin',bfdf)
        if (Node.eq.0) then
          if (peratm .and. h_spin_dim.lt.2) write(6,'(/,a)')
     .    'initdm: WARNING: DM.InitSpin not used because h_spin_dim < 2'
        endif

        if (peratm .and. h_spin_dim.ge.2) then

C Allocate local memory
          nullify(atom,phi,spin,theta)
          call re_alloc( atom, 1, maxat, 'atom', 'initdm' )
          call re_alloc( phi, 1, maxat, 'phi', 'initdm' )
          call re_alloc( spin, 1, maxat, 'spin', 'initdm' )
          call re_alloc( theta, 1, maxat, 'theta', 'initdm' )

          nat = 0
          badsyntax = .FALSE.
          do while(fdf_bline(bfdf,pline) .and. (nat .lt. nua) .and.
     .            (.not. badsyntax))

            nn = fdf_bnnames(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

            if (ni .eq. 1) then
              if (nat .eq. maxat) then
                maxatnew = nat + nint(0.1*nat)
C
                call re_alloc( atom, 1, maxatnew, 'atom', 'initdm',
     &                         copy=.true. )
C
                call re_alloc( phi, 1, maxatnew, 'phi', 'initdm',
     &                         copy=.true. )
                call re_alloc( spin, 1, maxatnew, 'spin', 'initdm',
     &                         copy=.true. )
                call re_alloc( theta, 1, maxatnew, 'theta', 'initdm',
     &                         copy=.true. )
C
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = fdf_bintegers(pline,1)

              if (nn .eq. 0) then
C Read value of spin
                if (nr .eq. 3) then
C Read spin value and direction
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = fdf_breals(pline,2) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,3) * pi/180.0d0
                elseif (nr .eq. 1) then
C Read spin value. Default direction.
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              elseif (nn .eq. 1) then
C Read spin as + or - (maximun value)
                updo = fdf_bnames(pline,1)
                if (updo .eq. '+') then
                  spin(nat) =  100.d0
                elseif (updo .eq. '-') then
                  spin(nat) = -100.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
                if (nr .eq. 2) then
                  theta(nat) = fdf_breals(pline,1) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,2) * pi/180.0d0
                elseif (nr .eq. 0) then
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              else
C Print bad-syntax error and stop
                badsyntax = .TRUE.
              endif

              if ((atom(nat) .lt. 1) .or. (atom(nat) .gt. nua)) then
                write(msg,'(a,a,i4)') 'intdm: ERROR: Bad atom ' //
     .            'index in DM.InitSpin, line', nat+1
                call die(TRIM(msg))
              endif
              if (abs(theta(nat)) .gt. 1.d-12) noncol = .true.
            else
C Print bad-syntax error and stop
              badsyntax = .TRUE.
            endif
          enddo

          if (badsyntax) then
            write(msg,'(a,i4)')
     .        'initdm: ERROR: bad syntax in DM.InitSpin, line', nat+1
            call die(msg)
          endif

          if (noncol .and. h_spin_dim.lt.4) then
            if (Node.eq.0) then
            write(6,'(/,2a)') 'initdm: WARNING: noncollinear spins ',
     .                 'in DM.InitSpin not used because h_spin_dim < 4'
            endif
            noncol = .false.
          endif

C Initialize to 0
          Dscf(:,:) = 0.0d0

C Initialize all paramagnetic

          do ia = 1, nua
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    Dscf(ind,1) = 0.5d0 * Datm(io)
                    Dscf(ind,2) = Dscf(ind,1)
                  endif
                enddo
              endif
            enddo
          enddo

C Loop on atoms with spin

          do iat = 1, nat
            ia = atom(iat)

C Find maximum atomic moment that the atoms involved can carry

            spinat = 0.d0
            do io = lasto(ia-1) + 1, lasto(ia)
              spinat = spinat + min( Datm(io), 2.d0 - Datm(io) )
            enddo
            if (spinat.lt.epsilon .and. Node.eq.0) print'(a,i6,a)',
     .        'initdm: WARNING: atom ', atom(iat),
     .        ' has a closed-shell and cannot be polarized'

C If given spin is larger than possible, make it to max atomic

            aspin = abs(spin(iat))
            if ((aspin .gt. spinat) .and. (aspin .gt. epsilon))
     .         spin(iat) = spinat*spin(iat)/aspin

C Initialize orbitals with same rate as atom

            rate = spin(iat) / (spinat+epsilon)
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                qio = Datm(io)
                spio = rate * min( Datm(io), 2.d0 - Datm(io) )
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    if ( noncol ) then
C Store non-collinear-spin density matrix as
C   ispin=1 => D11, ispin=2 => D22;
C   ispin=3 => Real(D12); ispin=4 => Imag(D12)
                      costh = cos(theta(iat))
                      sinth = sin(theta(iat))
                      cosph = cos(phi(iat))
                      sinph = sin(phi(iat))
                      Dscf(ind,1) = (qio + spio * costh) / 2
                      Dscf(ind,2) = (qio - spio * costh) / 2
                      Dscf(ind,3) =   spio * sinth * cosph / 2
                      Dscf(ind,4) = - spio * sinth * sinph / 2
                      if ( h_spin_dim == 8 ) then
                        Dscf(ind,5) = 0.0_dp
                        Dscf(ind,6) = 0.0_dp
                        Dscf(ind,7) = Dscf(ind,3)
                        Dscf(ind,8) = Dscf(ind,4)
                      end if
                    else
                      Dscf(ind,1) = (qio + spio) / 2
                      Dscf(ind,2) = (qio - spio) / 2
                    end if
                  end if
                end do
              end if
            end do

          end do

C Deallocate local memory
          call de_alloc( atom, 'atom', 'initdm' )
          call de_alloc( phi, 'phi', 'initdm' )
          call de_alloc( spin, 'spin', 'initdm' )
          call de_alloc( theta, 'theta', 'initdm' )

c ---------------------------------------------------------------------

        else

C Initialize to 0
          Dscf(:,:) = 0.0d0

C Automatic, for non magnetic (h_spin_dim=1) or for Ferro or Antiferro -----
          do io = 1, no_l
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
              ind = listhptr(io)+in
              jo = listh(ind)
              if (iio .eq. jo) then
                if (h_spin_dim .eq. 1) then

C No spin polarization

                  Dscf(ind,1) = Datm(iio)
                else

C Spin polarization

                  i1 = 1
                  i2 = 2

C Ferro or antiferro according to DM.InitSpinAF (inspn)

                  if (inspn) then
                    if (mod(iaorb(iio),2).eq.0) then
                      i1 = 2
                      i2 = 1
                    endif
                  endif
                  Dscf(ind,i1) = min( Datm(iio), 1.d0 )
                  Dscf(ind,i2) = Datm(iio) - Dscf(ind,i1)
                endif
              endif
            enddo
          enddo

        endif

      ! We have initialized with atomic information. Correct in case we
      ! are using such a small cell that there are direct interactions
      ! of orbitals with their own images, and we insist on using the
      ! Gamma-point only. Otherwise S(diagonal) is always 1.0 and the
      ! simple atomic-orbital superposition works as intended.
        

        do io = 1, no_l
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
               ind = listhptr(io)+in
               jo = listh(ind)
               if (iio .eq. jo) then  ! diagonal element
                  Dscf(ind,:) = Dscf(ind,:) / S(ind)
               endif
            enddo
         enddo

      endif

      if ( associated(num_tmp) )
     &    call de_alloc(num_tmp,name='num_tmp',routine=myName)
      if ( associated(listptr_tmp) )
     &    call de_alloc(listptr_tmp,name='listptr_tmp',routine=myName)
      if ( associated(list_tmp) )
     &    call de_alloc(list_tmp,name='listd')
      if ( associated(D_tmp) )
     &    call de_alloc(D_tmp,name='dm')
#ifdef TRANSIESTA
      if ( associated(EDM_tmp) )
     &    call de_alloc(EDM_tmp)
#endif /* TRANSIESTA */

      call print_initial_spin()


      CONTAINS

      subroutine print_initial_spin()
      use m_mpi_utils, only: Globalize_sum
      use sparse_matrices, only: S
      use m_spin, only: nspin => spinor_dim
      
      real(dp) :: qspin(nspin)
      integer  :: io, j, ispin, ind
#ifdef MPI
      real(dp) :: qtmp(nspin)
#endif
      
      ! We do nothing if there is only one spin-component
      if ( nspin == 1 ) return

      ! Print spin polarization
      do ispin = 1 , nspin
         qspin(ispin) = 0.0_dp
         do io = 1 , no_l
            do j = 1 , numh(io)
               ind = listhptr(io)+j
               qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            end do
         end do
      end do

#ifdef MPI
      ! Global reduction of spin components
      call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
      qspin(1:nspin) = qtmp(1:nspin)
#endif

      ! Only let the io-node print information
      if ( .not. IONode ) return
      
      write(6,'(/,a,f12.6)')
     .     'initdm: Initial spin polarization (Qup-Qdown) =', 
     .     qspin(1) - qspin(2)
      
      end subroutine print_initial_spin

      end subroutine initdm

      end module m_new_dm
