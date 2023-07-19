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
      use m_fdf_global, only: fdf_global_get

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

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
      use m_spin,   only: nspin
      use m_steps, only: istp
      use sparse_matrices, only: Escf     !! To be moved

#ifdef TRANSIESTA
      use M_ts_global_vars, only : TSrun
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

#ifdef TRANSIESTA
      integer :: ts_alloc_stat
#endif /* TRANSIESTA */

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
         call re_alloc(Dscf,1,nnz,1,nspin,name="Dscf",
     $                 routine="new_dm",copy=.false.,shrink=.true.)

         ! Reset DMsaved to null
         call de_alloc(DMsaved,name="DMsaved",routine="new_dm")
         nullify(DMsaved)
         !
         call initdm(Datm, Dscf, lasto, na_s,
     .               no_s, no_l, nspin, na_u, no_l, nspin,
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
            call fdf_global_get(allow_dm_extrapolation,
     $                          'DM.AllowExtrapolation'  , .FALSE.)
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

      END subroutine new_dm

!---------------------------------------------------

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
      logical   ::  sparsity_has_changed
      real(dp)  ::  msave
      integer    ::  i, j, ispin, nspin, nnz, new_nnz, is

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

         nspin = size(DM,dim=2)
         nnz = size(DM,dim=1)      ! This should be sum(numold)
         call re_alloc(DMtmp,1,nnz,1,nspin,
     $                 name="DMtmp",routine="new_dm", shrink=.true.)
         DMtmp(:,:) = DM(:,:)
         new_nnz = sum(num(:))
         call re_alloc(DM,1,new_nnz,1,nspin,
     $                 name="DM",routine="new_dm",shrink=.true.)
 
         do is = 1, nspin 
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
            nspin = size(DMsaved,dim=2)
            nnz = size(DMsaved,dim=1)
            call re_alloc(DMtmp,1,nnz,1,nspin,
     $           name="DMtmp",routine="new_dm", shrink=.true.)
            DMtmp(:,:) = DMsaved(:,:)
            new_nnz = sum(num(:))
            call re_alloc(DMsaved,1,new_nnz,1,nspin,
     $           name="DMsaved",routine="new_dm",shrink=.true.)
 
            do is = 1, nspin 
               call change_sparsity(nrows,
     $                              numold,listptrold,listold,
     $                              num, listptr, list,
     $                              DMtmp(:,is), DMsaved(:,is) )
            enddo
            call de_alloc(DMtmp,name="DMtmp",routine="new_dm")

         endif

         ! Extrapolate DM and update DMsaved
         ! This is a simple linear extrapolation

         nspin = size(DM,dim=2)
         nnz = size(DM,dim=1)

         do ispin = 1, nspin
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
            nspin = size(DM,dim=2)
            nnz = size(DM,dim=1)
            call re_alloc(DMsaved,1,nnz,1,nspin,name="DMsaved",
     $                    routine="new_dm",copy=.false.,shrink=.true.)
            DMsaved(1:nnz,1:nspin) = DM(1:nnz,1:nspin)
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

!
      subroutine initdm(Datm, Dscf, lasto, maxa,
     .                  maxo, maxuo, maxspn, nua, no_l, nspin, 
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
c integer no_l           : Number of orbitals in the unit cell
c integer nspin         : Number of spin components
c integer maxa          : Max num. atoms for dimension
c integer maxo          : Max. number of orbitals (globally)
c integer maxuo         : Max. number of orbitals (locally)
c integer maxspn        : Max number of spin components in dimensions
c integer lasto(0:maxa) : List with last orbital of each atom
c integer numh(:)       : Dscf matrix sparse information
c integer listhptr(:)   :  "
c integer listh(:)      :  "
c integer iaorb(maxo)   : List saying to what atom an orbital belongs 
c double Datm(no)       : Occupations of basis orbitals in free atom
c ********* OUTPUT **************************************************
c double Dscf(:,maxspn) : Density matrix in sparse form
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
      use siesta_options, only : isolve
#endif /* TRANSIESTA */

#ifdef MPI
      use mpi_siesta
#endif
      use units, only : pi

      implicit          none

      logical           found, inspn, try_dm_from_file
      integer           no_l, nua, maxo, maxuo, nspin, maxa, 
     .                  maxspn
      integer           lasto(0:maxa), numh(maxuo), iaorb(maxo)
      integer       ::  listhptr(:), listh(:)
      real(dp), pointer :: Dscf(:,:)
      real(dp)          Datm(maxo)
      external          memory

c ---------------------------------------------------------------------

C Internal variables and arrays
 
      character         updo*1, line*130, names*80
      logical           noncol, peratm
      integer           nh, ni, nn, nr, nv, iat, nat, ia, iu,
     .                  i1, i2, in, ind, ispin, jo, io, 
     .                  iio, maxatnew

      integer, save ::  maxat

      integer           integs(4), lastc, lc(0:3)

      integer, dimension(:), allocatable, save ::
     .                  atom, ibuffer
#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif
      real(dp)          aspin, cosph, costh, epsilon, 
     .                  qio, rate, reals(4),
     .                  sinph, sinth, spinat, spio, values(4)

      real(dp), dimension(:), allocatable, save ::
     .                  phi, spin, theta, dpbuffer1

      integer, pointer :: list_tmp(:), num_tmp(:), listptr_tmp(:)
      real(dp), pointer :: D_tmp(:,:), work(:)
#ifdef TRANSIESTA
      logical :: ts_found
      real(dp), pointer :: EDM_tmp(:,:)
#endif /* TRANSIESTA */

      integer :: nh_tmp, nspin_read, is
      logical :: same_sp

      data maxat / 1000 /
      data epsilon / 1.d-8 /

C Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------

      if (try_dm_from_file) then

         nullify(list_tmp,D_tmp,num_tmp,listptr_tmp)
         call re_alloc(num_tmp,1,no_l)
         call re_alloc(listptr_tmp,1,no_l)
#ifdef TRANSIESTA
         nullify(EDM_tmp)
#endif /* TRANSIESTA */


#ifndef TRANSIESTA
         call read_dynamic_dm(nh_tmp, no_l, nspin_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, found )
#else
         found = .false.   ! initialize
         if(isolve == 2) then ! Transiesta solution method
C First try to read the .TSDE file
           call read_dynamic_ts_dm(nh_tmp, no_l, nspin_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, EDM_tmp, ef, ts_found )

C Print information and do the some TS initializations
           call ts_init_dm(ts_found)

C If it does not find a .TSDE file, try to read a .DM file
           if(.not. ts_found) 
     .        call read_dynamic_dm(nh_tmp, no_l, nspin_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, found )           

C Either .TSDE or .DM may have been found
           found = ts_found .or. found
         else ! Not Transiesta solution mathod
           ts_found = .false.
           call read_dynamic_dm(nh_tmp, no_l, nspin_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, found )
         end if ! Transiesta solution method
#endif /* TRANSIESTA */

      else
        found = .false.
      endif

C If found, check and update, otherwise initialize with neutral atoms

      if (found) then
        ! Various degrees of sanity checks
        if (nspin_read /= nspin) then
           if (Node.eq.0) then
              write(6,"(a,i6,/,a)")
     $        "WARNING: Wrong nspin in DM file: ",
     $             nspin_read,
     $        "WARNING: Falling back to atomic initialization of DM."
              write(0,"(a,i6,/,a)")
     $        "WARNING: Wrong nspin in DM file: ",
     $             nspin_read,
     $        "WARNING: Falling back to atomic initialization of DM."
           endif
           found = .false.
           call de_alloc(num_tmp)
           call de_alloc(listptr_tmp)
           call de_alloc(list_tmp)
           call de_alloc(D_tmp)
#ifdef TRANSIESTA
           if(ts_found) call de_alloc(EDM_tmp)
#endif /* TRANSIESTA */
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
           do is= 1, nspin
              call change_sparsity (no_l,num_tmp,listptr_tmp,list_tmp,
     $                                  numh,listhptr,listh,
     $                                  D_tmp(:,is), Dscf(:,is),
     $                            work)     
#ifdef TRANSIESTA
! change_sparsity allocates work, but should work, since it does not 
! deallocate inside the spin loop
              if(ts_found) 
     .          call change_sparsity (no_l,num_tmp,listptr_tmp,list_tmp,
     $                                  numh,listhptr,listh,
     $                                  EDM_tmp(:,is), Escf(:,is),
     $                            work)
               
#endif /* TRANSIESTA */
           enddo
           call de_alloc(num_tmp)
           call de_alloc(listptr_tmp)
           call de_alloc(list_tmp)
           call de_alloc(D_tmp)
           call de_alloc(work)  ! Allocated by change_sparsity
#ifdef TRANSIESTA
           if(ts_found) call de_alloc(EDM_tmp)
#endif /* TRANSIESTA */
        endif


      else

C See whether specific initial spins are given in a DM.InitSpin block
C and read them in a loop on atoms where lines are read and parsed
C   integer nat       : how many atoms to polarize
C   integer atom(nat) : which atoms
C   double  spin(nat) : what polarization -----------------------------

        noncol = .false.
        if (Node.eq.0) then
          peratm = fdf_block('DM.InitSpin',iu)
          if (peratm .and. nspin.lt.2) write(6,'(/,a)')
     .    'initdm: WARNING: DM.InitSpin not used because nspin < 2'
        endif
#ifdef MPI
        call MPI_Bcast(peratm,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

        if (peratm .and. nspin.ge.2) then

C Allocate local memory
          allocate(atom(maxat))
          call memory('A','I',maxat,'initdm')
          allocate(phi(maxat))
          call memory('A','D',maxat,'initdm')
          allocate(spin(maxat))
          call memory('A','D',maxat,'initdm')
          allocate(theta(maxat))
          call memory('A','D',maxat,'initdm')

          nat = 0
          do iat = 1, nua+1
C Read and parse a line of the data block
            if (Node.eq.0) then
              read(iu,'(a)', end=50) line
            endif
#ifdef MPI
            call MPI_Bcast(line,130,MPI_character,0,MPI_Comm_World,
     .        MPIerror)
#endif
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (nn.ge.1 .and. names(lc(0)+1:lc(1)).eq.'%endblock') then
C End data reading
              goto 50
            elseif (ni .eq. 1) then
              if (nat.eq.maxat) then
                maxatnew = nat + nint(0.1*nat)
C
                allocate(ibuffer(maxat))
                call memory('A','I',maxat,'initdm')
                ibuffer(1:maxat) = atom(1:maxat)
                call memory('D','I',size(atom),'initdm')
                deallocate(atom)
                allocate(atom(maxatnew))
                call memory('A','I',maxatnew,'initdm')
                atom(1:maxat) = ibuffer(1:maxat)
                call memory('D','I',size(ibuffer),'initdm')
                deallocate(ibuffer)
C
                allocate(dpbuffer1(maxat))
                call memory('A','D',maxat,'initdm')
                dpbuffer1(1:maxat) = phi(1:maxat)
                call memory('D','D',size(phi),'initdm')
                deallocate(phi)
                allocate(phi(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                phi(1:maxat) = dpbuffer1(1:maxat)
                dpbuffer1(1:maxat) = spin(1:maxat)
                call memory('D','D',size(spin),'initdm')
                deallocate(spin)
                allocate(spin(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                spin(1:maxat) = dpbuffer1(1:maxat)
                dpbuffer1(1:maxat) = theta(1:maxat)
                call memory('D','D',size(theta),'initdm')
                deallocate(theta)
                allocate(theta(maxatnew))
                call memory('A','D',maxatnew,'initdm')
                theta(1:maxat) = dpbuffer1(1:maxat)
                call memory('D','D',size(dpbuffer1),'initdm')
                deallocate(dpbuffer1)
C
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = integs(1)
            else
C Print bad-syntax error and stop
              goto 40
            endif
            if (nn .eq. 0) then
C Read value of spin
              if (nr .eq. 3) then
C Read spin value and direction
                spin(nat)  = reals(1)
                theta(nat) = reals(2) * pi/180.0d0
                phi(nat)   = reals(3) * pi/180.0d0
              elseif (nr .eq. 1) then
C Read spin value. Default direction.
                spin(nat)  = reals(1)
                theta(nat) = 0.d0
                phi(nat)   = 0.d0
              else
                goto 40
              endif
            else if (nn .eq. 1) then
C Read spin as + or - (maximun value)
              updo = names(lc(0)+1:lc(1))
              if (updo .eq. '+') then
                spin(nat) =  100.d0
              elseif (updo .eq. '-') then
                spin(nat) = -100.d0
              else
                goto 40
              endif
              if (nr .eq. 2) then
                theta(nat) = reals(1) * pi/180.0d0
                phi(nat)   = reals(2) * pi/180.0d0
              elseif (nr .eq. 0) then
                theta(nat) = 0.d0
                phi(nat)   = 0.d0
              else
                goto 40
              endif
            else
              goto 40
            endif
            if (atom(nat).lt.1 .or. atom(nat).gt.nua) then
              if (Node.eq.0) then
                write(6,'(/,a)')
     .          'intdm: ERROR: Bad atom index in DM.InitSpin, line', iat
              endif
              call die('intdm: ERROR: Bad atom index in DM.InitSpin')
            endif
            if (abs(theta(nat)).gt.1.d-12) noncol = .true.
          enddo

          if (Node.eq.0) then
            write(6,'(a)') 
     .         'initdm: ERROR: Too many atom entries in DM.InitSpin'
           call die('initdm: ERROR: Too many' //
     $              ' atom entries in DM.InitSpin')
          else
            call die()
          endif

   40     continue
          if (Node.eq.0) then
            write(6,*)
     .       'initdm: ERROR: bad syntax in DM.InitSpin, line', iat
            call die('initdm: ERROR: bad syntax in DM.InitSpin')
          else
            call die()
          endif

   50     continue

          if (noncol .and. nspin.lt.4) then
            if (Node.eq.0) then
            write(6,'(/,2a)') 'initdm: WARNING: noncolinear spins ',
     .                 'in DM.InitSpin not used because nspin < 4'
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
                    if (noncol) then
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
                    else
                      Dscf(ind,1) = (qio + spio) / 2
                      Dscf(ind,2) = (qio - spio) / 2
                    endif
                  endif
                enddo
              endif
            enddo

          enddo

C Deallocate local memory
          call memory('D','I',size(atom),'initdm')
          deallocate(atom)
          call memory('D','D',size(phi),'initdm')
          deallocate(phi)
          call memory('D','D',size(spin),'initdm')
          deallocate(spin)
          call memory('D','D',size(theta),'initdm')
          deallocate(theta)

c ---------------------------------------------------------------------

        else

C Initialize to 0

          Dscf(:,:) = 0.0d0

C Automatic, for non magnetic (nspin=1) or for Ferro or Antiferro -----

          do io = 1, no_l
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
              ind = listhptr(io)+in
              jo = listh(ind)
              if (iio .eq. jo) then
                if (nspin .eq. 1) then

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

      endif

      call print_initial_spin()


      CONTAINS

      subroutine print_initial_spin()
      use m_mpi_utils, only: Globalize_sum
      use sparse_matrices, only: S

      real(dp) :: qspin(nspin)
      integer  :: io, j, ispin, ind
#ifdef MPI
      real(dp) :: qtmp(nspin)
#endif

! Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo

#ifdef MPI
! Global reduction of spin components
        call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
        qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
           if (IOnode) then
              write(6,'(/,a,f12.6)')
     .             'initdm: Initial spin polarization (Qup-Qdown) =', 
     .             qspin(1) - qspin(2)
           endif
        endif
      endif
      end subroutine print_initial_spin

      end subroutine initdm

      END MODULE m_new_dm




