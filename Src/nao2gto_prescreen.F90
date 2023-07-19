      module prescreen

! prescreen ERIs using schwarz inequality to get list_ij
! list_uv and list_mn
! xmqin, October 2013.

      use kinds,          only : dp, int_8
      use atm_types,      only : species, species_info, &
                                 nspecies,l_max, nco, nso 
      use hfx_types,      only : D2Sindx, subshell, pair_list_type
      use hfx_types,      only : hfx_input_parameter, eri_prescreen
      use hfx_types,      only : hfx_screen_coeff_type,&
                                 log_zero,&
                                 powell_min_log
      use hfx_types,      only : pair_dist_radii_pgf, &
                                 sfc_pgf,sfc_shell
      use parallel,       only : Node, Nodes
      use atmfuncs,       only : lofio, mofio, rcut
      use atomlist,       only : indxuo, rmaxo, lasto
      use atomlist,       only : rco, rmaxkb, lastkb
      use listsc_module,  only : listsc
      use extended_index, only : indexsc
      use cell_pbc,       only : pbc, trans_pbc
      use alloc,          only : re_alloc, de_alloc
      use sorting
      use m_recipes,      only : sort
      use neighbour,      only : jna=>jan, xij, r2ij, maxna=>maxnna
      use neighbour,      only : mneighb
      use contract_eri,         only : calc_contract_eri2
      use primitive_eri,        only : calc_primitive_screen                                   
      use libint_wrapper_types, only : lib_int
      use powell,               only : opt_state_type,&
                                       powell_optimize
      use nao2gto_util,         only : exp_radius_very_extended

      
      implicit none
      public ::  calc_prescreen_eri, build_pair_list, &
                 calc_pair_dist_radii, calc_screening_functions
      private

      CONTAINS

      subroutine calc_prescreen_eri(lib, norb, iaorb, iphorb, nuo, nuotot, &
                                    na, isa, xa, indxua, cell, rcell, &
                                    maxnh, numh, listhptr, listh, &
                                    hfx_parameter, max_eri )
! Input variables and arrays:
       type(lib_int)         :: lib
       type(hfx_input_parameter) :: hfx_parameter
       integer, intent(in)   ::  &
        maxnh, na, norb, nuo, nuotot, &
        iaorb(norb), indxua(na), iphorb(norb), isa(na), &
        listh(maxnh),listhptr(nuo), numh(nuo)

       real(dp), intent(in)  :: xa(3,na), cell(3,3), rcell(3,3)
! Out 
!       real(dp), intent(out) ::  eri_prescreen(maxnh)
       real(dp), intent(out) ::  max_eri
! Internal variables and arrays:
       integer ia, ioa, io, iu, is, i, j,       &
               ja, joa, jo, ju, js, ind,        &
               l_i, l_j, m_i, m_j, ncoi, ncoj,   &
               npgfi, npgfj

       real(dp) ri(3), rj(3), r_temp(3), r_pbc(3)

       type(species_info), pointer  :: ispp, jspp
       real(dp), dimension(:,:,:,:), allocatable  :: eri
    
!       external trans_pbc, reclat
        do io   =  1, nuotot
           ia   =  iaorb(io)
           is   =  isa(ia)
           ispp => species(is)
           ioa  =  iphorb(io)
           l_i  =  lofio(is, ioa)
           m_i  =  mofio(is, ioa)
           if(m_i .ne. -l_i) cycle
           ri(:) = xa(:, ia)
           npgfi =  ispp%orbnl_contract(ispp%orb_index(ioa)) 
           ncoi  =  nco(l_i)*npgfi
          
           do j   =  1, numh(io)
              ind =  listhptr(io) + j
              jo  =  listh(ind)
              ja  =  iaorb(jo)
              js  =  isa(ja)
              jspp => species(js)
              joa =  iphorb(jo)
              l_j =  lofio(js, joa)
              m_j =  mofio(js, joa)
              if(m_j .ne. -l_j) cycle

              rj(:)  = xa(:, ja)
              r_temp = rj - ri
           !   r_pbc = pbc(r_temp, cell, rcell )
              call trans_pbc(r_temp,cell,rcell,r_pbc)
              rj = ri + r_pbc

              npgfj = jspp%orbnl_contract(jspp%orb_index(joa)) 
              ncoj = nco(l_j)*npgfj
            
             allocate( eri(nso(l_i),nso(l_j),nso(l_i),nso(l_j)) )
             eri=0.0d0

             call calc_contract_eri2( lib, cell, rcell, ri, rj, ri, rj, &
                                     npgfi, npgfj, npgfi, npgfj, &
                                     l_i, l_j, l_i, l_j,         &
                                     ncoi, ncoj, ncoi, ncoj,     &
                 ispp%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)),  &
                 jspp%orbnl_zeta(1:npgfj, jspp%orb_index(joa)),  & 
                 ispp%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)),  &
                 jspp%orbnl_zeta(1:npgfj, jspp%orb_index(joa)),  &
                 ispp%sphi(1:ncoi,ioa:ioa+nso(l_i)-1),           &  
                 jspp%sphi(1:ncoj,joa:joa+nso(l_j)-1),           &  
                 ispp%sphi(1:ncoi,ioa:ioa+nso(l_i)-1),           &
                 jspp%sphi(1:ncoj,joa:joa+nso(l_j)-1),           &
                 hfx_parameter, eri )

             eri_prescreen(ind) = 2.d0*max(maxval(eri),-minval(eri))

             deallocate(eri) 
           enddo
        enddo
       
       !----------end prescreen use schwarz inequality----------
 
       max_eri=max(maxval(eri_prescreen),-minval(eri_prescreen))
       max_eri=dsqrt(max_eri)

       if(Node.eq.0) write(6,*) 'max_eri_prescreen',  max_eri
       
      end subroutine calc_prescreen_eri

! ---------------------------------------------------------------------------
!  New subroutine to find shell pair list (uv) and (mn)
!
!  Coded by xmqin, 12. 06, 2016,  different from shanghui's method.
!    
!  (1) Find neighbours for all atoms in supercell using rmax , 
!      Number of neighbours  : nnia    max : maxnna = 200
!      Neighbours' index: jna(nnia)
!      PBC vectors to neighbours : xij(jna) 
!      Squared distances to neighbours : r2ij(jna), 
!
!      Actually, I call the original subritine "neighbour()" of SISTEA
!      program to find atomic pair_list without additional definition 
!      for pair_atomic_list '(type)'.  
!      Loop over atoms to give xij and rij in batch !
!      
!  (2) Find all shell pairs over orbital loop
!     (a) Hermite and trans symmetry : ishell_unit > jshell_unit
!     (b) NAO overlap criteria: rcut(io)+rcut(jo) > rij, io overlap with jo,
!     (c) Schwarz inequality included denity matrix
!---------------------------------------------------------------------------
!
! Useful information:
!
! (1)  To deal with periodic boundary conditions, the unit cell is 'extended'
!      on each side. 
!      There are two kinds of index for images( unit cell, atoms and obitals ) :
!      one is in the PBC supercell, the other is in the normal (or auxiliary) supercell.
!  
!  PBC unit cell : Wigner-Seitz cell A(i) = -1/2a(i), 1/2a(i) , a, Ai : lattice vector.
!  PBC supercell : To extend Wigner cell on each side of direction: -nsc(i)*A/2, nsc(i)*A/2
!  Normal supercell : Direct supercell : A(i) = 1, nsc(i)*a !
! 
!  The index of images should be shift between pbc supercell and normal supercell 
!  because we usually use normal index. 
!
!  (1) In order to take into account the 8-fold symmetry under pbc :
!
!   (u0v[R]|m[R']n[R"])        =   (v0u[-R]|m[R'-R]n[R"-R]) 
! = (u0v[R]|n[R"]m[R'])        =   (v0u[-R]|n[R"-R]m[R'-R])
! = (m0n[R"-R']|u[-R']v[R-R']) =   (n0m[R'-R"]|u[-R"]v[R-R"]) 
! = (m0n[R"-R']|v[R-R']u[-R']) =   (n0m[R'-R"]|v[R-R"]u[-R"])
!
!  The first index must be always in unit cell, so we should transfer R to 0.
!
!   pair_list by iterating in the following way.
!   we should transfer R to 0 
!
!   DO iatom=1,natom               
!      iatom_unit = mod(iatom-1, unit_atoms)+1      
!
!      DO jatom=iatom, negbours(iatom)                      
!         jatom_pbc = ind1 (ind2(jatom)-ind2(iatom)+ind2(iatom_unit))
!         jatom_unit_pbc = mod(jatom-1, unit_atoms)+1
!         IF(jatom_unit<iatom_unit) CYCLE        
!
!         atom_ij = ncells*(iatom_unit)*(iatom_unit-1)/2
!                  +mod(jatom_pbc,unit_atoms)*iatom_unit+jatom_unit_pbc
!     
!  atom_ij is the uv or mn pair list index gives u[0]v[R] or
!  m[0]n[R],  it is not equal to number of loop.
! 
!  Here "_pbc" denotes images' index in pbc supercell.  
!  v[0] means that we choose v[R] as reference unitcell.
!  then u, m, n must be shift to extended supercell related to v[R].
!
!!
!! IF(katom+latom<=iatom+jatom)  THEN
!!   IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
! or 
!----------------------------------------------------------------------------
! integer nia: Atom index of orbital (u)s into the loop .
!              To build (uv), nia = nua  for u
!                       (mn), nia = na   for m

! integer norb         : Number of orbitals in supercell
! real*8  scell(3,3)   : Supercell vectors SCELL(IXYZ,IVECT)
! integer nsc(3)       : Num. of unit cells in each supercell direction
! real*8  xa(3,na)     : Atomic positions in cartesian coordinates
! integer lasto(0:na)  : Last orbital of each atom in array iphorb
! integer iphorb(norb) : Orbital index of each orbital in its atom,
!                       where no=lasto(na)
! integer iphorb(norb) 
! integer isa(na)      : Species index of each atom

!---------------------------------------------------------------------------  
      subroutine build_pair_list( nia, na, nua, isa, xa, indxua, nuotot, norb,   &
                                  iaorb, iphorb, maxnh, numh, listhptr, listh,  &
                                  cell, rcell, hfx_parameter, &
                                  max_eri, list_orb )

! Input variables and arrays:
      implicit none
      integer,  intent(in)    :: nia, na, nua, isa(na), indxua(na), norb, maxnh
      integer,  intent(in)    :: iphorb(norb), iaorb(norb)
      integer, intent(in)     :: &
        listh(maxnh), listhptr(nuotot), &
        numh(nuotot)

      real(dp), intent(in)    :: cell(3,3), rcell(3,3), xa(3,na)
      real(dp), intent(in)    :: max_eri
      type(hfx_input_parameter) :: hfx_parameter

      type(pair_list_type)      :: list_orb
!------------------------------------------------------------------
       integer ncells, nuotot, n_shells
       integer ia, ind, io, ioa, is, i, j,       &
               ja, jn, jo, joa, js, jx, nnia,      &
               iuo, juo, jos, juos, iua, jua, josc
       
        
       integer l_i, l_j, m_i, m_j, index_ij, ishell, jshell, &
               ishell_unit, jshell_unit,  unit_shells
       real(dp) ri(3), rj(3), r_temp(3), r_pbc(3)
       real(dp) rci, rcj, rij, r2 
! For order vertor :
       integer, dimension(:),  pointer :: index => null()
       real(dp) mcell(3,3), rmcell(3,3)
       integer :: eri_prescreen_ind
!       logical :: find_in_Dscf
       real(dp) ::  eri_prescreen_value

! ......................
!  rmaxo = max rcut of all orbitals
!      if(node.eq.0) then
!       write(6,'(a)') 'To find the neighbours of an atom in supercell.'
!       write(6,'(a,f12.6,a)') "using rmaxo= ", rmaxo+rmaxkb, " bohr." 
!      endif

! Initialize neighb subroutine 

! Allocate local memory
!      nullify(list_orb)

      unit_shells = subshell(nuotot)
      ncells = norb / nuotot

!      call mneighb( cell, 2.0_dp*rmaxo, na, xa, 0, 0, nnia )
!
!      write(6,'(a,2i8)') 'Number of cells and unit oritals', ncells, nuotot
!      call re_alloc(index,1, maxna, name="index")

      list_orb%nelement=0

      do ia = 1, nia ! the last cell 5*5*5 scell
         iua = mod(ia-1, nua)+1
         ri(1:3)=xa(1:3, ia)

!         call mneighb( cell, 2._dp*rmaxo, na, xa, ia, 0, nnia)

!       sort by by distance
!         call sort( nnia, r2ij, index )
!         call iorder( jna, 1, nnia, index )
!         call order(  r2ij, 1, nnia, index )
!         call order(  xij, 3, nnia, index )

         do io  = lasto(ia-1)+1,lasto(ia)
            ioa = iphorb(io)
            is  = isa(iaorb(io))
            l_i = lofio(is,ioa)
            m_i = mofio(is,ioa)
            if(m_i .ne.-l_i) cycle
            ishell = subshell(io)
            iuo = mod(io-1,nuotot)+1
            ishell_unit = subshell(iuo)

!  Find neighbor atoms using rmaxo = max(rcut(natom))
!            do jn  = 1, nnia
!               ja  = jna(jn)
             do ja = 1, na
               !rj(1:3) = xa(1:3,ia) + xij(1:3,jn)
               !rij = dsqrt(r2ij(jn))
               r_temp(1:3)=xa(1:3,ja)-xa(1:3,ia)
               call trans_pbc(r_temp, cell, rcell, r_pbc)  ! xmqin for debug
               r2=r_pbc(1)**2+r_pbc(2)**2+r_pbc(3)**2
               rj = ri + r_pbc
               rij = dsqrt(r2) 
               
               do jo  = lasto(ja-1)+1,lasto(ja)
                  joa = iphorb(jo)
                  js  = isa(iaorb(jo))
                  l_j = lofio(js,joa)
                  m_j = mofio(js,joa)
                  if(m_j .ne.-l_j) cycle
                  jos = indexsc(io, iuo, jo)
                  juo = mod(jos-1,nuotot)+1                 
                  jshell =  subshell(jos)
                  jshell_unit = subshell(juo)
                
                  if(jshell_unit .le. ishell_unit) then
                     !write(99,*) "iii", io, jo , rcut(is,ioa)+rcut(js,joa) 
                     !write(99,*) rij
  
                     if(rcut(is,ioa)+rcut(js,joa) .le. rij) cycle

                     eri_prescreen_value = eri_prescreen(D2Sindx(io,jo))
                     
                     if( dsqrt(dabs(eri_prescreen_value))*max_eri.ge. &
                         hfx_parameter%eps_pairlist ) then

                         n_shells= ncells*(ishell_unit)*(ishell_unit-1)/2 &
                                   +((jshell-1)/unit_shells)*ishell_unit+jshell_unit

                         list_orb%nelement=list_orb%nelement+1
                         i=list_orb%nelement
                         list_orb%element(i)%pair(1) = io
                         list_orb%element(i)%pair(2) = jo
                         list_orb%element(i)%nl_index = n_shells
                         list_orb%element(i)%r1(1:3) = ri(1:3)
                         list_orb%element(i)%r2(1:3) = rj(1:3)
                         list_orb%element(i)%dist2 = r2

                     endif
                   endif

               enddo
             enddo

        enddo
      enddo
!      call de_alloc(index,name="index")
! Finish timer

      end subroutine build_pair_list



    SUBROUTINE calc_pair_dist_radii(hfx_parameter, radii_pgf)
!    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hfx_input_parameter)                 :: hfx_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:,:,:,:,:,:), POINTER   :: radii_pgf

    INTEGER :: handle, i, ikind, ipgf, isave, iset, jkind, jpgf, jset, la, &
      lb, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, nsgfa, &
                                                nsgfb

    REAL(dp) :: cutoff, DATA(2,0:100), ff, max_contraction_a, &
      max_contraction_b, prefactor, R1, R_max, ra(3), rab(3), rab2, radius, &
      rap(3), rb(3), rp(3), x(2), zetp
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_a, set_radius_b
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    integer ia, ioa, io, iu, is,        &
            ja, joa, jo, ju, js, ind,        &
            l_i, l_j, m_i, m_j, ncoi, ncoj,   &
            npgfi, npgfj, ishell, jshell
    type(species_info), pointer  :: ispp, jspp
    real(dp) :: eps_new

! ----------------------------------------------------------------
    if(hfx_parameter%eps_schwarz>1.0E-6_dp ) then
      eps_new = 1.0E-6_dp 
    else
      eps_new = hfx_parameter%eps_schwarz
    endif

    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp
    DO is=1,nspecies
       ispp => species(is)
      DO js = 1,nspecies
         jspp => species(js)

        ishell = 0
        DO io = 1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if (m_i .ne. -l_i) cycle     ! not a normal orbital
           ishell=ishell+1
           npgfi = ispp%orbnl_contract(ishell)
           ncoi  = nco(l_i)*npgfi
!           if(ishell == 2) cycle

          jshell = 0
          DO jo = 1,jspp%norbs
             l_j = jspp%orb_l(jo)
             m_j = jspp%orb_m(jo)
             if (m_j .ne. -l_j) cycle     ! not a normal orbital
             jshell=jshell+1
             npgfj = jspp%orbnl_contract(jshell)
             ncoj  = nco(l_j)*npgfj

!            if(jshell == 2) cycle
            DO ipgf = 1,npgfi
              DO jpgf = 1,npgfj
                 radius=ispp%pgf_radius(ipgf,ishell)+jspp%pgf_radius(jpgf,jshell)
             !   write(18,*)"ra",  ipgf,ishell, ispp%pgf_radius(ipgf,ishell)
             !   write(18,*) "rb", jpgf,jshell, jspp%pgf_radius(jpgf,jshell)
             !   write(18,*)  radius
                DO i=0,100
                  rb(1) = 0.0_dp + 0.01_dp * radius*i
                  R_max = 0.0_dp
                      zetp = ispp%orbnl_zeta(ipgf,ishell) + jspp%orbnl_zeta(jpgf,jshell)
                      ff = jspp%orbnl_zeta(jpgf,jshell)/zetp
                      rab = 0.0_dp
                      rab(1) = rb(1)
                      rab2 = rb(1)**2
                      prefactor = EXP(-ispp%orbnl_zeta(ipgf,ishell)*ff*rab2)
                      rap(:) =ff*rab(:)
                      rp(:) = ra(:) + rap(:)
                      rb(:) = ra(:) + rab(:)
                      cutoff = 1.0_dp

                      R1=exp_radius_very_extended(l_i,l_i,l_j,l_j,ra=ra,rb=rb,rp=rp, &
                                       zetp=zetp,eps=eps_new, prefactor=prefactor, &
                                       cutoff=cutoff,epsabs=1.0E-12_dp)
                      R_max = MAX(R_max,R1)

                  DATA(1,i) = rb(1)
                  DATA(2,i) = R_max
!                  write(17,*) jpgf,ipgf,jshell,ishell, i
!                  write(17,*) R1,  DATA(1,i), DATA(2,i)
                END DO

                CALL optimize_it(DATA,x,0.0_dp)
                radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x
!                write(21,*) jpgf,ipgf,jshell,ishell,js,is
!                write(21,*) x
              END DO !jpgf
            END DO !ipgf
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE calc_pair_dist_radii

  SUBROUTINE calc_screening_functions(lib, hfx_parameter, cell, rcell,  &
                                      coeffs_pgf, coeffs_set, radii_pgf)
    TYPE(lib_int)                            :: lib
    TYPE(hfx_input_parameter)                :: hfx_parameter
    real(dp), intent(in)  :: cell(3,3), rcell(3,3)

!    TYPE(hfx_screen_coeff_type), &
!      DIMENSION(:, :), POINTER               :: coeffs_kind
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
     DIMENSION(:, :, :, :, :, :), POINTER   :: radii_pgf, coeffs_pgf

! -----------------------------------------------------------------------------
    INTEGER                                  :: handle, i, ikind, ipgf, iset, &
                                                jkind, jpgf, jset, la, lb, &
                                                ncoa, ncob, nkind, nseta, &
                                                nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    REAL(dp) :: DATA(2,0:100), kind_radius_a, kind_radius_b, &
      max_contraction_a, max_contraction_b, max_val, max_val_temp, R1, ra(3), &
      radius, rb(3), x(2), max_val2, max_val_temp2
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_a, set_radius_b
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: tmp_R_1
    integer ia, ioa, io, iu, is,        &
            ja, joa, jo, ju, js, ind,        &
            l_i, l_j, m_i, m_j, ncoi, ncoj,   &
            npgfi, npgfj, ishell, jshell
    real(dp) zeta_i, zeta_j
    type(species_info), pointer  :: ispp, jspp
       real(dp), dimension(:,:,:,:), allocatable  :: eri
! ------------------------------------------------------------------------------

    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp

    DO is=1, nspecies
       ispp => species(is)
      DO js = 1, nspecies
         jspp => species(js)

        ishell = 0
        DO io = 1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if (m_i .ne. -l_i) cycle     ! not a normal orbital
           ishell=ishell+1
           npgfi = ispp%orbnl_contract(ishell)
           ncoi  = nco(l_i)*npgfi
           max_contraction_a = maxval((/(sum(abs(ispp%sphi(1:ncoi, &
                         i))),i=io,io+nso(l_i)-1)/))

         jshell = 0
          DO jo = 1,jspp%norbs
             l_j = jspp%orb_l(jo)
             m_j = jspp%orb_m(jo)
             if (m_j .ne. -l_j) cycle     ! not a normal orbital
             jshell=jshell+1
            ! if((jshell .ne. 1) .or.(jshell .ne. 3)) cycle
             npgfj = jspp%orbnl_contract(jshell)
             ncoj  = nco(l_j)*npgfj
             max_contraction_b =  maxval((/(sum(abs(jspp%sphi(1:ncoj, &
                         i))),i=jo,jo+nso(l_j)-1)/))
            !write(41,*) ishell,jshell, max_contraction_a, max_contraction_b
            DO ipgf = 1, npgfi
               zeta_i = ispp%orbnl_zeta(ipgf,ishell)
              DO jpgf = 1, npgfj
                 zeta_j = jspp%orbnl_zeta(jpgf,jshell)
                 radius=ispp%pgf_radius(ipgf,ishell)+jspp%pgf_radius(jpgf,jshell)
           !       write(51,*) is, js, ishell,jshell
           !       write(51,*) ipgf,jpgf,ispp%pgf_radius(ipgf,ishell),jspp%pgf_radius(jpgf,jshell)
           !       write(51,*) ZETA_I,ZETA_J
           !       write(101,"(4I4)") ishell,jshell, ipgf, jpgf
           !       write(301,*) "###", ishell,jshell, ipgf, jpgf

                DO i=0,100
                  rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
                  max_val = 0.0_dp
                  max_val_temp = 0.0_dp

                  R1 = MAX(0.0_dp, radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(1)*rb(1)**2 &
                       + radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(2))

                  CALL calc_primitive_screen(lib, ra, rb, ra, rb, &
                                             zeta_i, zeta_j, zeta_i, zeta_j, &
                                             l_i, l_j, l_i ,l_j, &
                                             max_val_temp, hfx_parameter, R1, R1)
!                 write(43,*) R1, max_val_temp
             !     write(501,*) zeta_i,zeta_j,l_i,l_j
                  max_val = MAX(max_val, max_val_temp)
 !                max_val = 2.0_dp*MAX(max_val, max_val_temp)
             !     write(43,*) max_val_temp, max_val
             !     write(301,*)  rb(1),  max_val
                  max_val = SQRT(max_val)

                  max_val = max_val * max_contraction_a * max_contraction_b
                  DATA(1,i) = rb(1)
                  IF(max_val == 0.0_dp ) THEN
                    DATA(2,i) = powell_min_log
                  ELSE
                    DATA(2,i) = LOG10((max_val))
                  END IF
             !     write(101,*) i, DATA(1, i), DATA(2, i)
                END DO
!                  write(101,*) "#xxxxxxxxxxxxxxxxxxx"
!                  write(301,*) "#xxxxxxxxxxxxxxxxxxx"

                CALL optimize_it(DATA,x,powell_min_log)
                coeffs_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x
!                write(22,*) jpgf,ipgf,jshell,ishell,js,is
!                write(22,*) x


!            write(201,"(6I4)") is, js, ishell, jshell, ipgf,jpgf
!             write(201,*)
!             write(201,*) ishell, jshell,ipgf,jpgf
!            DO i=0,100
!              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
!              max_val2 = x(1)*rb(1)**2 + x(2)
!              write(201,*) rb(1), DATA(2, i),  max_val2
!            ENDDO
!             write(201,*) "xxxxxxxxxxxxxxxxxxx"

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp

    DO is=1,nspecies
       ispp => species(is)
      DO js = 1,nspecies
         jspp => species(js)

       
        ishell = 0
        DO io = 1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if (m_i .ne. -l_i) cycle     ! not a normal orbital
           ishell=ishell+1
           npgfi = ispp%orbnl_contract(ishell)
           ncoi  = nco(l_i)*npgfi

           max_contraction_a =  maxval((/(sum(abs(ispp%sphi(1:ncoi, &
                          i))),i=io,io+nso(l_i)-1)/))

          jshell = 0
          DO jo = 1,jspp%norbs
             l_j = jspp%orb_l(jo)
             m_j = jspp%orb_m(jo)
             if (m_j .ne. -l_j) cycle     ! not a normal orbital
             jshell=jshell+1
             npgfj = jspp%orbnl_contract(jshell)
             ncoj  = nco(l_j)*npgfj
             max_contraction_b =  maxval((/(sum(abs(jspp%sphi(1:ncoj, &
                         i))),i=jo,jo+nso(l_j)-1)/))
             radius=ispp%shell_radius(ishell)+jspp%shell_radius(jshell)

             tmp_R_1 => radii_pgf(:,:,jshell,ishell,js,is)

!            write(203,*)is, js, ishell, jshell
!            write(103,*)is, js, ishell, jshell

            DO i=0,100
              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              max_val_temp = 0.0_dp
!              max_val2 = 0.0_dp
!              max_val_temp2 = 0.0_dp

             allocate( eri(nso(l_i),nso(l_j),nso(l_i),nso(l_j)) )
             eri=0.0d0
             call calc_contract_eri2( lib, cell, rcell, ra, rb, ra, rb, &
                                     npgfi, npgfj, npgfi, npgfj, &
                                     l_i, l_j, l_i, l_j,         &
                                     ncoi, ncoj, ncoi, ncoj,     &
                 ispp%orbnl_zeta(1:npgfi, ishell),  &
                 jspp%orbnl_zeta(1:npgfj, jshell),  &
                 ispp%orbnl_zeta(1:npgfi, ishell),  &
                 jspp%orbnl_zeta(1:npgfj, jshell),  &
                 ispp%sphi(1:ncoi,io:io+nso(l_i)-1),           &
                 jspp%sphi(1:ncoj,jo:jo+nso(l_j)-1),           &
                 ispp%sphi(1:ncoi,io:io+nso(l_i)-1),           &
                 jspp%sphi(1:ncoj,jo:jo+nso(l_j)-1),           &
                 hfx_parameter, eri )

!              write(1111,*) eri
              max_val = MAX(max_val, maxval(dabs(eri)))
              max_val = 2.0_dp*max_val
!              write(102,*) i, rb(1), max_val
!              write(102,*)
              max_val = DSQRT(max_val)


!             DO ipgf = 1,npgfi
!                 zeta_i = ispp%orbnl_zeta(ipgf,ishell)
!                 DO jpgf = 1,npgfj
!                    zeta_j = jspp%orbnl_zeta(jpgf,jshell)

!                    R1 = MAX(0.0_dp, tmp_R_1(jpgf,ipgf)%x(1)*rb(1)**2  &
!                             + tmp_R_1(jpgf,ipgf)%x(2))

!                    CALL calc_primitive_screen(lib, ra, rb, ra, rb,    &
!                                       zeta_i, zeta_j, zeta_i, zeta_j, &
!                                       l_i, l_j, l_i ,l_j,             &
!                                       max_val_temp, hfx_parameter, R1, R1)
!
!                    max_val = MAX(max_val, max_val_temp)
!                 END DO !jpgf
!              END DO !ipgf

!              max_val = SQRT(max_val)
!              max_val = max_val * max_contraction_a * max_contraction_b
             
              DATA(1,i) = rb(1)
              IF(max_val == 0.0_dp ) THEN
                DATA(2,i) = powell_min_log
              ELSE
                DATA(2,i) = LOG10((max_val))
              END IF
             deallocate(eri)

           END DO
            CALL optimize_it(DATA,x,powell_min_log)
            coeffs_set(jshell,ishell,js,is)%x = x
!            write(23,*) is, js, ishell, jshell
!            write(23,*) x
!            write(23,*) "xxxxxxxxxxxxxxxxxxxxxxxx"

!            write(204,*) is, js, ishell, jshell
!            DO i=0,100
!              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
!              max_val2 = coeffs_set(jshell,ishell,js,is)%x(1)*rb(1)**2 +   &
!                         coeffs_set(jshell,ishell,js,is)%x(2)
!              write(204,*) rb(1), DATA(2,i), max_val2
!            ENDDO
!            write(204,*) "xxxxxxxxxxxxxxxxxxxxxxxx"


          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE calc_screening_functions

!
!
! little drive routine for the powell minimizer
! data is the data to fit, x is of the form (x(1)*DATA(1)**2+x(2))
! only values of DATA(2) larger than fmin are taken into account
! it constructs an approximate upper bound of the fitted function
!
!
  SUBROUTINE optimize_it(DATA,x,fmin)

    REAL(KIND=dp), INTENT(IN)                :: DATA(2,0:100)
    REAL(KIND=dp), INTENT(OUT)               :: x(2)
    REAL(KIND=dp), INTENT(IN)                :: fmin

    INTEGER                                  :: i, k
    REAL(KIND=dp)                            :: f, large_weight, &
                                                small_weight, weight
    TYPE(opt_state_type)                     :: opt_state

! initial values

     x(1) = 0.0_dp
     x(2) = 0.0_dp

     ! we go in two steps, first we do the symmetric weight to get a good,
     ! unique initial guess
     ! we restart afterwards for the assym.
     ! the assym function appears to have several local minima, depending on the
     ! data to fit
     ! the loop over k can make the switch gradual, but there is not much need,
     ! seemingly
     DO k=0,4,4

        small_weight=1.0_dp
        large_weight=small_weight*(10.0_dp**k)

        ! init opt run
        opt_state%state = 0
        opt_state%nvar = 2
        opt_state%iprint = 3
        opt_state%unit = 6
        opt_state%maxfun = 100000
        opt_state%rhobeg = 0.1_dp
        opt_state%rhoend = 0.000001_dp

        DO

          ! compute function
          IF ( opt_state%state == 2 ) THEN
            opt_state%f = 0.0_dp
            DO i=0,100
              f = x(1)*DATA(1,i)**2 +  x(2)
              IF( f > DATA(2,i) ) THEN
                weight = small_weight
              ELSE
                weight = large_weight
              END IF
              IF( DATA(2,i) > fmin ) opt_state%f = opt_state%f + weight * (f-DATA(2,i))**2
            END DO
          END IF

          IF ( opt_state%state == -1 ) EXIT
          CALL powell_optimize (opt_state%nvar, x, opt_state)
        END DO

        ! dealloc mem
        opt_state%state = 8
        CALL powell_optimize (opt_state%nvar, x, opt_state)

     ENDDO

  END SUBROUTINE optimize_it

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i,j 2d index
!> \param N matrix size
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  PURE FUNCTION get_1D_idx(i,j,N)
    INTEGER, INTENT(IN)                      :: i, j
    INTEGER(int_8), INTENT(IN)               :: N
    INTEGER(int_8)                           :: get_1D_idx

    INTEGER(int_8)                           :: min_ij

    min_ij = MIN(i,j)
    get_1D_idx = min_ij*N + MAX(i,j) - (min_ij-1)*min_ij/2 - N

  END FUNCTION get_1D_idx

      end module prescreen 


