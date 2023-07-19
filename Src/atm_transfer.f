! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine atm_transfer

      use atm_types, only: maxnorbs, nspecies
      use atm_types, only: species, species_info

      use radial
      use atmparams, only:NTBMAX
!----------------------------------------------------------------
      use old_atmfuncs, only: nsmax
!
!     old_atmfuncs arrays
!
      use old_atmfuncs, only: tabpol, table, tab2
      use old_atmfuncs, only: coretab, tab2pol
      use old_atmfuncs, only: qtb, slfe
      use old_atmfuncs, only: chloctab, vlocaltab
      use old_atmfuncs, only: lmxosave, npolorbsave
      use old_atmfuncs, only: nzetasave, nsemicsave, nkblsave
!
!     old_atmfuncs procedures
!
      use old_atmfuncs, only: labelfis, izofis, zvalfis
      use old_atmfuncs, only: massfis, lomaxfis, nofis
      use old_atmfuncs, only: cnfigfio, lofio, mofio
      use old_atmfuncs, only: atmpopfio, epskb, rcut
      use old_atmfuncs, only: lmxkbfis, nkbfis

!----------------------------------------------------------------

      use periodic_table, only: symbol
      use sys,            only: die

      implicit none

      type(species_info), pointer        :: spp
      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp

      integer is, io, i , n, ntot, l, m, max_npjnl
      integer max_norbnl, nsm, izeta, j
      integer norb, indx, ipol, num_normal, num_pol

      integer, dimension(maxnorbs) :: index_normal, z_normal,
     $     nsm_normal, index_pol, z_pol, nsm_pol
      integer, dimension(maxnorbs) :: mark  ! overkill

!     Allocate main structures for new atmfuncs -----------

      nspecies = nsmax           ! From old_atmfuncs

      allocate(species(nspecies))
!-----------------------------------------------------------

      max_npjnl = 0
      max_norbnl = 0
      do is=1,nspecies
         spp => species(is)

         spp%read_from_file = .false.

         spp%label = labelfis(is)
         spp%z     = izofis(is)
         if (spp%z.eq.-100) then
            spp%symbol = 'BS'
         else
            spp%symbol = symbol(abs(spp%z))
         endif
         spp%zval  = zvalfis(is)
         spp%mass  = massfis(is)
         spp%self_energy  = slfe(is)


         spp%lmax_basis = lomaxfis(is)
         spp%norbs = nofis(is)

         do io = 1,  spp%norbs
            spp%orb_n(io) = cnfigfio(is,io)  !! Not sure about this
            spp%orb_l(io) = lofio(is,io)
            spp%orb_m(io) = mofio(is,io)
            spp%orb_pop(io) = atmpopfio(is,io)
         enddo
!
!           Double check
!

         index_normal = 0
         index_pol = 0
         do io = 1,  spp%norbs
            norb=0 
            indx=0
            do  l=0,lmxosave(is)
               do nsm=1,nsemicsave(l,is)+1
                  do izeta=1,nzetasave(l,nsm,is)
                     norb=norb+(2*l+1)
                     indx=indx+1
                     if(norb.ge.io) then
                        index_normal(io) = indx
                        z_normal(io) = izeta
                        nsm_normal(io) = nsm
                        goto 20
                     endif
                  enddo 
               enddo 
            enddo 

            indx=0
            do  l=0,lmxosave(is)
               do nsm=1,nsemicsave(l,is)+1
                  do ipol=1, npolorbsave(l,nsm,is)
                     norb=norb+(2*(l+1)+1)
                     indx=indx+1
                     if(norb.ge.io) then
                        index_pol(io) = indx
                        z_pol(io) = ipol
                        nsm_pol(io) = nsm
                        goto 20
                     endif
                  enddo 
               enddo 
            enddo  

 20         continue

         enddo

         num_normal = maxval(index_normal)
         num_pol = maxval(index_pol)

         spp%n_orbnl = num_normal + num_pol


         allocate(spp%orbnl(spp%n_orbnl))
!
         mark = 0
         do io = 1, spp%norbs
            i = index_normal(io)
            if (i .eq. 0) cycle     ! not a normal orbital

            spp%orb_index(io) = i

            if (mark(i) .ne. 0) cycle   ! nl orb already set up
            mark(i) = 1

            l = spp%orb_l(io)
            spp%orbnl_l(i) = l
            spp%orbnl_n(i) = spp%orb_n(io)
!!            spp%orbnl_n(i) = nsm_pol(io)
            spp%orbnl_z(i) = z_normal(io)
            spp%orbnl_pop(i) = qtb(io,is) * (2*l+1)
            spp%orbnl_ispol(i) = .false.

            op => spp%orbnl(i)

            call rad_alloc(op,NTBMAX)
            op%delta  =   table(1,i,is)
            op%cutoff =   op%delta *  (NTBMAX - 1)
            op%f(1:)    = table(3:,i,is)
            op%d2(1:)   = tab2(1:,i,is)

         enddo
!
!        Polarization orbitals
!
         do io = 1, spp%norbs
            j = index_pol(io)
            if (j .eq. 0) cycle     ! not a polarization orbital
            i = num_normal + j      ! pol orbs come after normal ones
            spp%orb_index(io) = i

            if (mark(i) .ne. 0) cycle   ! nl orb already set up
            mark(i) = 1

            l = spp%orb_l(io)
            spp%orbnl_l(i) = l
            spp%orbnl_n(i) = spp%orb_n(io)
!!            spp%orbnl_n(i) = nsm_pol(io)
            spp%orbnl_z(i) = z_pol(io)
            spp%orbnl_pop(i) = qtb(io,is) * (2*l+1)
            spp%orbnl_ispol(i) = .true.

            op => spp%orbnl(i)
            call rad_alloc(op,NTBMAX)
            op%delta  =   tabpol(1,j,is)
            op%cutoff =   op%delta *  (NTBMAX - 1)
            op%f(1:)    = tabpol(3:,j,is)
            op%d2(1:)   = tab2pol(1:,j,is)
         enddo

!
!        KB projectors (relatively easy...)
!
         spp%nprojs = nkbfis(is)
         spp%lmax_projs = lmxkbfis(is)

         do i = 1, spp%nprojs
            io = - i                        !! Old convention
!!!!        spp%pj_n(i)  ???????????? useful??
            spp%pj_l(i) = lofio(is,io)
            spp%pj_m(i) = mofio(is,io)
         enddo

!
!        This piece of code assumes that the projectors are ordered
!        in their usual manner
!
         n = 0
         ntot = 0
         do l = 0, spp%lmax_projs
            do i = 1, nkblsave(l,is)
               n = n + 1
               spp%pjnl_n(n) = i
               spp%pjnl_l(n) = l
               do m = 1, 2*l+1
                  ntot = ntot + 1
                  spp%pj_index(ntot) = n
               enddo
            enddo
         enddo
         spp%n_pjnl = n
         if (ntot .ne. spp%nprojs) call die('KB indexing...')


         allocate(spp%pjnl(spp%n_pjnl))

         do i = 1, spp%n_pjnl
            pp => spp%pjnl(i)
            call rad_alloc(pp,NTBMAX)
            pp%delta  =   table(1,-i,is)
            pp%f(1:)    = table(3:,-i,is)
            pp%d2(1:)   = tab2(1:,-i,is)
         enddo
!
!        Fill in the KB energy array and the cutoffs
!        A bit redundant

         do i = 1, spp%nprojs
            io = -i
            indx = spp%pj_index(i)
            spp%pjnl_ekb(indx) = epskb(is,io)
            pp => spp%pjnl(indx)
            pp%cutoff = rcut(is,io)
         enddo



         call rad_alloc(spp%vna,NTBMAX)
         spp%vna%f(1:)       = table(3:,0,is)
         spp%vna%cutoff      = table(2,0,is)
         spp%vna%delta       = table(1,0,is)
         spp%vna%d2(1:)      = tab2(1:,0,is)

         call rad_alloc(spp%chlocal,NTBMAX)
         spp%chlocal%delta      = chloctab(1,1,is)
         spp%chlocal%cutoff     = chloctab(1,1,is)*(NTBMAX-1)
         spp%chlocal%f(1:)      = chloctab(2:,1,is)
         spp%chlocal%d2(1:)     = chloctab(2:,2,is)

         call rad_alloc(spp%reduced_vlocal,NTBMAX)
         spp%reduced_vlocal%delta      = vlocaltab(1,1,is)
         spp%reduced_vlocal%cutoff     = vlocaltab(1,1,is)*(NTBMAX-1)
         spp%reduced_vlocal%f(1:)      = vlocaltab(2:,1,is)
         spp%reduced_vlocal%d2(1:)     = vlocaltab(2:,2,is)

         spp%there_is_core      = (coretab(1,2,is) .eq. 1)

         call rad_alloc(spp%core,NTBMAX)
         spp%core%delta         = coretab(1,1,is)
         spp%core%cutoff        = coretab(1,1,is)*(NTBMAX-1)
         spp%core%f(1:)         = coretab(2:,1,is)
         spp%core%d2(1:)        = coretab(2:,2,is)

      enddo

      end subroutine atm_transfer







