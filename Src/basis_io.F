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
      module basis_io
!
!     Support for dumping and reading PAO and KB information from
!     ASCII or NetCDF files.
!
!     Alberto Garcia, 2000, 2001
!
      use chemical
      use sys, only: die
      use precision
      use atom_options, only: write_ion_plot_files
#ifdef CDF
      use atmparams, only:NTBMAX
#endif
      use atm_types
      use basis_types, only: write_basis_specs, basis_parameters
      use pseudopotential
      use radial
      use fdf
!      use flib_wxml, only: str
      use xml, only: xml_dump_attribute, xml_dump_element, str

      implicit none

      public :: dump_basis_ascii, read_basis_ascii
      public :: dump_basis_xml
      public :: dump_basis_netcdf, read_basis_netcdf
      public :: read_ion_ascii

      type(species_info), pointer        :: spp

!!!!!!!      private :: xml_element_content

      private

      CONTAINS

      subroutine read_basis_netcdf(ns)

#ifndef CDF
      integer, intent(out) :: ns

      call die(
     $ '*** You need netCDF to read the new user-defined basis files...'
     $ )
      ns = 0
      end subroutine  read_basis_netcdf
#else 

      use netcdf

      integer, intent(out) :: ns

      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp

      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $     pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $     kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $     cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vna_id, chlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l, nrp_tables, core_flag, nor, nk, m

      character(len=40) filename
      character(len=40) dummy
                                !
      call read_chemical_types
      nspecies = number_of_species()
      ns = nspecies

      allocate(species(nspecies))

      do is = 1, nspecies
         spp => species(is)
         spp%label = species_label(is)
         spp%read_from_file = .true.
         write(filename,'(a,a)') trim(spp%label), ".ion.nc"
         iret = nf90_open(trim(filename),NF90_NOWRITE,ncid)

         iret = nf90_inq_dimid(ncid,'norbs',norbs_id)
         iret = nf90_inquire_dimension(ncid,norbs_id,len=norbs)
         if (norbs .gt. maxnorbs)
     $  call die("read_user_basis: Increase maxnorbs in atm_types.f")

         spp%n_orbnl = norbs

         iret = nf90_inq_dimid(ncid,'nkbs',nkbs_id)
         iret = nf90_inquire_dimension(ncid,nkbs_id,len=nkbs)
         spp%n_pjnl = nkbs

!
!        For now, it is assumed that *all* the radial arrays have
!        the same length.
!
         iret = nf90_inq_dimid(ncid,'ntb',ntb_id)
         iret = nf90_inquire_dimension(ncid,ntb_id,len=nrp_tables)
         if (nrp_tables .ne. NTBMAX) call die("NTBMAX mismatch")

         allocate(spp%orbnl(norbs))
         allocate(spp%pjnl(nkbs))


         iret = nf90_get_att(ncid,nf90_global,'Element',spp%symbol)
!!!        iret = nf90_get_att(ncid,nf90_global,'Label',dummy)
                                !! Sanity check here??

         iret = nf90_get_att(ncid,nf90_global,'Atomic_number',spp%z)
         if (atomic_number(is) .ne. spp%z)
     $                            call die("Atomic number mismatch")

         iret = nf90_get_att(ncid,nf90_global,'Valence_charge',spp%zval)
         iret = nf90_get_att(ncid,nf90_global,'Mass',spp%mass)
         iret = nf90_get_att(ncid,nf90_global,'Self_energy',
     $                                         spp%self_energy)
         iret = nf90_get_att(ncid,nf90_global,
     $        'Number_of_orbitals',spp%norbs)
         iret = nf90_get_att(ncid,nf90_global,
     $        'L_max_basis',spp%lmax_basis)
         iret = nf90_get_att(ncid,nf90_global,
     $        'Number_of_projectors',spp%nprojs)
         iret = nf90_get_att(ncid,nf90_global,
     $        'L_max_projs',spp%lmax_projs)

                                !! Orbitals

         iret = nf90_inq_varid(ncid,'orbnl_l',orbnl_l_id)
         iret = nf90_inq_varid(ncid,'orbnl_n',orbnl_n_id)
         iret = nf90_inq_varid(ncid,'orbnl_z',orbnl_z_id)
         iret = nf90_inq_varid(ncid,'orbnl_ispol',orbnl_ispol_id)
         iret = nf90_inq_varid(ncid,'orbnl_pop',orbnl_pop_id)

         iret = nf90_inq_varid(ncid,'cutoff',cutoff_id)
         iret = nf90_inq_varid(ncid,'delta',delta_id)

                                !!      Projectors

         iret = nf90_inq_varid(ncid,'pjnl_l',pjnl_l_id)
         call check(iret)
         iret = nf90_inq_varid(ncid,'pjnl_n',pjnl_n_id)
         iret = nf90_inq_varid(ncid,'pjnl_ekb',pjnl_ekb_id)
         iret = nf90_inq_varid(ncid,'kbcutoff',kbcutoff_id)
         iret = nf90_inq_varid(ncid,'kbdelta',kbdelta_id)
         call check(iret)

         iret = nf90_inq_varid(ncid,'orb',orb_id)
         call check(iret)
                                !
                                !       Local potential
                                !
         iret = nf90_inq_varid(ncid,'vna',vna_id)
         iret = nf90_get_att(ncid,vna_id,
     $        'Vna_cutoff',spp%vna%cutoff)
         iret = nf90_get_att(ncid,vna_id,
     $        'Vna_delta',spp%vna%delta)
                                !
                                !       Local potential charge density
                                !
         iret = nf90_inq_varid(ncid,'chlocal',chlocal_id)
         iret = nf90_get_att(ncid,chlocal_id,
     $        'Chlocal_cutoff',spp%chlocal%cutoff)
         iret = nf90_get_att(ncid,chlocal_id,
     $        'Chlocal_delta',spp%chlocal%delta)
                                !
                                !       Core charge
                                !
         iret = nf90_get_att(ncid,nf90_global,'Core_flag',core_flag)
         spp%there_is_core = (core_flag .eq. 1)

         if (spp%there_is_core) then
            iret = nf90_inq_varid(ncid,'core',core_id)
            iret = nf90_get_att(ncid,core_id,
     $           'Core_cutoff',spp%core%cutoff)
            iret = nf90_get_att(ncid,core_id,
     $           'Core_delta',spp%core%delta)
         else
            call rad_zero(spp%core)
         endif

         call check(iret)
         iret = nf90_inq_varid(ncid,'proj',proj_id)
         call check(iret)

         iret = nf90_get_var(ncid,pjnl_l_id,spp%pjnl_l,count=(/nkbs/))
         call check(iret)
         iret = nf90_get_var(ncid,pjnl_n_id,spp%pjnl_n,count=(/nkbs/))
         call check(iret)
         iret = nf90_get_var(ncid,pjnl_ekb_id,spp%pjnl_ekb,
     $        count=(/nkbs/))
         call check(iret)

         iret=nf90_get_var(ncid,orbnl_l_id,spp%orbnl_l,count=(/norbs/))
         iret=nf90_get_var(ncid,orbnl_n_id,spp%orbnl_n,count=(/norbs/))
         iret=nf90_get_var(ncid,orbnl_z_id,spp%orbnl_z,count=(/norbs/))

         iret = nf90_get_var(ncid,orbnl_ispol_id,aux,count=(/norbs/))
         do i = 1, norbs
            spp%orbnl_ispol(i) =  aux(i) .eq. 1
         enddo
         call check(iret)
         iret = nf90_get_var(ncid,orbnl_pop_id,
     $        spp%orbnl_pop,count=(/norbs/))
         call check(iret)

c     
         nk = 0
         do i = 1, nkbs
            pp => spp%pjnl(i)
            call rad_alloc(pp,NTBMAX)
            iret = nf90_get_var(ncid,proj_id,pp%f(1:),
     $           start=(/1,i/),count=(/NTBMAX,1/))
            call check(iret)
            iret = nf90_get_var(ncid,kbcutoff_id,pp%cutoff,
     $           start=(/i/))
            call check(iret)
            iret = nf90_get_var(ncid,kbdelta_id,pp%delta,
     $           start=(/i/))
            call check(iret)
            call rad_setup_d2(pp,yp1=0.0_dp,ypn=huge(1.0_dp))
            l = spp%pjnl_l(i)
            do m = -l,l
               nk = nk+1
               spp%pj_n(nk) = spp%pjnl_n(i)
               spp%pj_l(nk) = spp%pjnl_l(i)
               spp%pj_m(nk) = m
               spp%pj_index(nk) = i
            enddo
         enddo
         spp%nprojs = nk
                                !
                                !       Local potential
                                !
         call rad_alloc(spp%vna,NTBMAX)
         iret = nf90_get_var(ncid,vna_id,spp%vna%f(1:),
     $        start=(/1/),count=(/NTBMAX/))
         call check(iret)
         call rad_setup_d2(spp%vna,yp1=0.0_dp,ypn=huge(1.0_dp))

                                !
                                !       Local potential charge density
                                !
         call rad_alloc(spp%chlocal,NTBMAX)
         iret = nf90_get_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $        start=(/1/),count=(/NTBMAX/))
         call check(iret)
         call rad_setup_d2(spp%chlocal,yp1=0.0_dp,ypn=huge(1.0_dp))

         if (spp%there_is_core) then
            call rad_alloc(spp%core,NTBMAX)
            iret = nf90_get_var(ncid,core_id,spp%core%f(1:),
     $           start=(/1/),count=(/NTBMAX/))
            call check(iret)
            call rad_setup_d2(spp%core,yp1=0.0_dp,ypn=huge(1.0_dp))
         endif

         nor = 0
         do i = 1, norbs
            op => spp%orbnl(i)
            call rad_alloc(op,NTBMAX)
            iret = nf90_get_var(ncid,orb_id,op%f(1:),
     $           start=(/1,i/),count=(/NTBMAX,1/))
            call check(iret)
            iret = nf90_get_var(ncid,cutoff_id,op%cutoff,
     $           start=(/i/))
            call check(iret)
            iret = nf90_get_var(ncid,delta_id,op%delta,
     $           start=(/i/))
            call check(iret)
            call rad_setup_d2(op,yp1=0.0_dp,ypn=huge(1.0_dp))
            l = spp%orbnl_l(i)
            do m = -l,l
               nor = nor+1
               spp%orb_n(nor) = spp%orbnl_n(i)
               spp%orb_l(nor) = spp%orbnl_l(i)
               spp%orb_m(nor) = m
               spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
               spp%orb_index(nor) = i
            enddo
         enddo
         spp%norbs = nor
         iret = nf90_close(ncid)
         call check(iret)

      enddo
!
      CONTAINS

      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         call die()
      endif
      end subroutine check

      end subroutine read_basis_netcdf
#endif

!----------------------------------------------------------------
#ifdef CDF
#endif
!=======================================================================

      subroutine read_basis_ascii(ns)
      integer, intent(out) :: ns

      integer is

      call read_chemical_types
      nspecies = number_of_species()
      ns = nspecies

      allocate(species(nspecies))

      do is = 1, nspecies
         spp => species(is)
         spp%label = species_label(is)
         spp%read_from_file = .true.
         call read_ion_ascii(spp)
      enddo

      end subroutine read_basis_ascii
!
!----------------------
      subroutine read_ion_ascii(spp)
      type(species_info), pointer  :: spp

      character(len=20) filename
      integer i, l, m, lun, nor, nk, ispol

      write(filename,'(a,a)') trim(spp%label), ".ion"
      call io_assign(lun)
      open(lun,file=filename,status='old',form='formatted')
      rewind(lun)

      call read_header(spp,lun)
      read(lun,*) 
         do i=1,spp%n_orbnl
            read(lun,*)
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i)
            spp%orbnl_ispol(i) =  ispol.eq.1
            call radial_read_ascii(spp%orbnl(i),lun,
     $                             yp1=0.0_dp,ypn=huge(1.0_dp))
         enddo
!
!        Update indexes
!
         nor = 0
         do i = 1, spp%n_orbnl
            l = spp%orbnl_l(i)
            do m = -l,l
               nor = nor+1
               spp%orb_n(nor) = spp%orbnl_n(i)
               spp%orb_l(nor) = spp%orbnl_l(i)
               spp%orb_m(nor) = m
               spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
               spp%orb_index(nor) = i
            enddo
         enddo
         spp%norbs = nor

! KBs
         read(lun,*)
         do i=1,spp%n_pjnl
            read(lun,*)
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i)
            call radial_read_ascii(spp%pjnl(i),lun,
     $                             yp1=0.0_dp,ypn=huge(1.0_dp))
         enddo
!
!        Update indexes
!

         nk = 0
         do i = 1, spp%n_pjnl
            l = spp%pjnl_l(i)
            do m = -l,l
               nk = nk+1
               spp%pj_n(nk) = spp%pjnl_n(i)
               spp%pj_l(nk) = spp%pjnl_l(i)
               spp%pj_m(nk) = m
               spp%pj_index(nk) = i
            enddo
         enddo
         spp%nprojs = nk
!
!Vna
         read(lun,*)
         call radial_read_ascii(spp%vna,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))

!
!Chlocal
         read(lun,*)
         call radial_read_ascii(spp%chlocal,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))
!
!Core
         read(lun,*,end=9999)
         call radial_read_ascii(spp%core,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))

 9999    continue
         call io_close(lun)

      CONTAINS

      subroutine read_header(p,unit)

      type(species_info), pointer :: p
      integer, intent(in)         :: unit
      
      character(len=78) line

      read(unit,'(a)') line
      if (trim(line) .eq. '<preamble>') then
 1       continue
         read(unit,'(a)') line
         if (trim(line) .ne. '</preamble>') goto 1
      endif
         
      read(unit,'(a2)') p%symbol
      read(unit,'(a20)') p%label
      read(unit,*) p%z
      read(unit,*) p%zval
      read(unit,*) p%mass
      read(unit,*) p%self_energy
      read(unit,*) p%lmax_basis, p%n_orbnl
      read(unit,*) p%lmax_projs, p%n_pjnl

      allocate(p%orbnl(p%n_orbnl))
      allocate(p%pjnl(p%n_pjnl))

      end subroutine read_header

      end subroutine read_ion_ascii


#ifndef CDF
      subroutine dump_basis_netcdf
!     Do nothing
      end subroutine dump_basis_netcdf
#else
      subroutine dump_basis_netcdf

      use netcdf

      type(rad_func), pointer            :: pp
      type(rad_func), pointer            :: op

      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $        pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $        kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $        cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vna_id, chlocal_id, reduced_vlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l
      character*40 filename

      do is = 1, nspecies
         spp => species(is)
         if (spp%read_from_file) cycle   !! Do not dump

         write(filename,'(a,a)') trim(spp%label), ".ion.nc"
         write(6,'(2a)') 'Dumping basis to NetCDF file ',
     $                   trim(filename)

         iret = nf90_create(trim(filename),NF90_CLOBBER,ncid)

        nkbs =  spp%n_pjnl
        norbs = spp%n_orbnl

        iret = nf90_def_dim(ncid,'norbs',norbs,norbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'nkbs',nkbs,nkbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'ntb',NTBMAX,ntb_id)
        call check(iret)


!! Orbitals

        iret = nf90_put_att(ncid,nf90_global,'Element',spp%symbol)
        iret = nf90_put_att(ncid,nf90_global,'Label',spp%label)
        iret = nf90_put_att(ncid,nf90_global,'Atomic_number',spp%z)
        iret = nf90_put_att(ncid,nf90_global,'Valence_charge',spp%zval)
        iret = nf90_put_att(ncid,nf90_global,'Mass',spp%mass)
        iret = nf90_put_att(ncid,nf90_global,'Self_energy',
     $                                       spp%self_energy)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'Number_of_orbitals',spp%norbs)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'L_max_basis',spp%lmax_basis)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'Number_of_projectors',spp%nprojs)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'L_max_projs',spp%lmax_projs)

        iret = nf90_def_var(ncid,'orbnl_l',nf90_int,norbs_id,orbnl_l_id)
        iret = nf90_def_var(ncid,'orbnl_n',nf90_int,norbs_id,orbnl_n_id)
        iret = nf90_def_var(ncid,'orbnl_z',nf90_int,norbs_id,orbnl_z_id)
        iret = nf90_def_var(ncid,'orbnl_ispol',nf90_int,
     $                            norbs_id,orbnl_ispol_id)
        iret = nf90_def_var(ncid,'orbnl_pop',nf90_double,
     $                            norbs_id,orbnl_pop_id)

        iret = nf90_def_var(ncid,'cutoff',nf90_double,
     $                           norbs_id,cutoff_id)
        iret = nf90_def_var(ncid,'delta',nf90_double,
     $                           norbs_id,delta_id)

!!      Projectors

        iret = nf90_def_var(ncid,'pjnl_l',nf90_int,nkbs_id,pjnl_l_id)
        call check(iret)
        iret = nf90_def_var(ncid,'pjnl_n',nf90_int,nkbs_id,pjnl_n_id)
        iret = nf90_def_var(ncid,'pjnl_ekb',nf90_double,
     $                           nkbs_id,pjnl_ekb_id)
        iret = nf90_def_var(ncid,'kbcutoff',nf90_double,
     $                           nkbs_id,kbcutoff_id)
        iret = nf90_def_var(ncid,'kbdelta',nf90_double,
     $                           nkbs_id,kbdelta_id)
        call check(iret)

        iret = nf90_def_var(ncid,'orb',nf90_double,
     $                      (/ntb_id,norbs_id/),orb_id)
        call check(iret)
!
!       Local potential
!
        iret = nf90_def_var(ncid,'vna',nf90_double,
     $                      (/ntb_id/),vna_id)
        iret = nf90_put_att(ncid,vna_id,
     $                      'Vna_cutoff',spp%vna%cutoff)
        iret = nf90_put_att(ncid,vna_id,
     $                      'Vna_delta',spp%vna%delta)
!
!       Local potential charge density
!
        iret = nf90_def_var(ncid,'chlocal',nf90_double,
     $                      (/ntb_id/),chlocal_id)
        iret = nf90_put_att(ncid,chlocal_id,
     $                      'Chlocal_cutoff',spp%chlocal%cutoff)
        iret = nf90_put_att(ncid,chlocal_id,
     $                      'Chlocal_delta',spp%chlocal%delta)
!
!       Reduced Local potential (rV+2*Zval)
!
        iret = nf90_def_var(ncid,'reduced_vlocal',nf90_double,
     $                      (/ntb_id/),reduced_vlocal_id)
        iret = nf90_put_att(ncid,reduced_vlocal_id,
     $         'Reduced_vlocal_cutoff',spp%reduced_vlocal%cutoff)
        iret = nf90_put_att(ncid,reduced_vlocal_id,
     $         'Reduced_vlocal_delta',spp%reduced_vlocal%delta)
!
!       Core charge
!
        if (spp%there_is_core) then
           iret = nf90_put_att(ncid,nf90_global,
     $                         'Core_flag',1)
           iret = nf90_def_var(ncid,'core',nf90_double,
     $                      (/ntb_id/),core_id)
           iret = nf90_put_att(ncid,core_id,
     $                      'Core_cutoff',spp%core%cutoff)
           iret = nf90_put_att(ncid,core_id,
     $                      'Core_delta',spp%core%delta)
        else
           iret = nf90_put_att(ncid,nf90_global,
     $                         'Core_flag',0)
        endif

        call check(iret)
        iret = nf90_def_var(ncid,'proj',nf90_double,
     $                      (/ntb_id,nkbs_id/),proj_id)
        call check(iret)

!!!!!!!!
        iret = nf90_enddef(ncid)
        call check(iret)

        iret = nf90_put_var(ncid,pjnl_l_id,spp%pjnl_l,count=(/nkbs/))
        call check(iret)
        iret = nf90_put_var(ncid,pjnl_n_id,spp%pjnl_n,count=(/nkbs/))
        call check(iret)
        iret = nf90_put_var(ncid,pjnl_ekb_id,spp%pjnl_ekb,
     $                           count=(/nkbs/))
        call check(iret)

        iret = nf90_put_var(ncid,orbnl_l_id,spp%orbnl_l,count=(/norbs/))
        iret = nf90_put_var(ncid,orbnl_n_id,spp%orbnl_n,count=(/norbs/))
        iret = nf90_put_var(ncid,orbnl_z_id,spp%orbnl_z,count=(/norbs/))

        if (norbs .gt. maxnorbs)
     $        call die("dump_basis_netcdf: Increase maxnorbs")
        aux(1:norbs) = 0
        do i = 1, norbs
           if (spp%orbnl_ispol(i)) aux(i)=1
        enddo
        iret = nf90_put_var(ncid,orbnl_ispol_id,aux,count=(/norbs/))
        call check(iret)
        iret = nf90_put_var(ncid,orbnl_pop_id,
     $                           spp%orbnl_pop,count=(/norbs/))
        call check(iret)

c
        do i = 1, nkbs
           pp => spp%pjnl(i)
           iret = nf90_put_var(ncid,proj_id,pp%f(1:),
     $                      start=(/1,i/),count=(/NTBMAX,1/))
           call check(iret)
           iret = nf90_put_var(ncid,kbcutoff_id,pp%cutoff,
     $                      start=(/i/))
           call check(iret)
           iret = nf90_put_var(ncid,kbdelta_id,pp%delta,
     $                      start=(/i/))
           call check(iret)
        enddo
!
!       Local potential
!
        iret = nf90_put_var(ncid,vna_id,spp%vna%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
        call check(iret)

!
!       Local potential charge density
!
        iret = nf90_put_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
        call check(iret)

!
!       Reduced Local potential
!
        iret = nf90_put_var(ncid,reduced_vlocal_id,
     $       spp%reduced_vlocal%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
        call check(iret)

        if (spp%there_is_core) then
           iret = nf90_put_var(ncid,core_id,spp%core%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
           call check(iret)
        endif

        do i = 1, norbs
           op => spp%orbnl(i)
           iret = nf90_put_var(ncid,orb_id,op%f(1:),
     $                      start=(/1,i/),count=(/NTBMAX,1/))
           call check(iret)
           iret = nf90_put_var(ncid,cutoff_id,op%cutoff,
     $                      start=(/i/))
           call check(iret)
           iret = nf90_put_var(ncid,delta_id,op%delta,
     $                      start=(/i/))
           call check(iret)
        enddo
        iret = nf90_close(ncid)
        call check(iret)

      enddo

      contains
      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         call die()
      endif
      end subroutine check

      end subroutine dump_basis_netcdf

#endif


      subroutine dump_basis_ascii

      integer is

      do is = 1, nspecies
         call dump_ion_ascii(is)
      enddo
      end subroutine dump_basis_ascii
!---
      subroutine dump_ion_ascii(is)
      integer, intent(in)      :: is

      type(species_info), pointer  ::   spp
      character*40 filename
      character*30 fileid
      integer i, l, ispol, lun, lun2, zeta, n_series

      spp => species(is)
      if (spp%read_from_file) return    !! Do not dump

      write(filename,'(a,a)') trim(spp%label), ".ion"
      call io_assign(lun)
      open(lun,file=filename,status='replace',form='formatted')

      write(lun,'(a)') '<preamble>'
      call write_basis_specs(lun,is)
      call pseudo_header_print(lun,
     $                   basis_parameters(is)%pseudopotential)
      write(lun,'(a)') '</preamble>'
      call write_header(spp,lun)
      write(lun,'(a)') "# PAOs:__________________________"
         n_series = 0
         do i=1,spp%n_orbnl
            ispol = 0
            if (spp%orbnl_ispol(i)) ispol = 1
            zeta = spp%orbnl_z(i)
            if (zeta == 1) n_series = n_series + 1

            write(lun,'(4i3,f10.6,2x,a)')
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i),
     $           " #orbital l, n, z, is_polarized, population"
            call radial_dump_ascii(spp%orbnl(i),lun)

           if (write_ion_plot_files) then
            write(fileid,'(a,i1,a,i1,a,a)') "ORB.S", n_series, ".",
     $            zeta, ".", trim(spp%label)
            call io_assign(lun2)
            open(unit=lun2,file=fileid,status='replace',
     $                     form='formatted')
            write(lun2,'(2a,4i3,f10.4)') "# ", trim(spp%label),
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i)
            write(lun2,'(a)')
     $       "#(species label, l, n, z, is_polarized, popul)"
            call radial_dump_ascii(spp%orbnl(i),lun2,header=.false.)
            call io_close(lun2)
           endif ! (write_ion_plot_files)

         enddo
         
         write(lun,'(a)') "# KBs:__________________________"
         n_series = 0
         do i=1,spp%n_pjnl
            zeta = spp%pjnl_n(i)
            l = spp%pjnl_l(i)
            write(lun,'(2i3,f22.16,2x,a)')
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i),
     $           " #kb l, n (sequence number), Reference energy"
            call radial_dump_ascii(spp%pjnl(i),lun)

           if (write_ion_plot_files) then
            write(fileid,'(a,i1,a,i1,a,a)') "KB.L", l , ".",
     $            zeta, ".", trim(spp%label)
            call io_assign(lun2)
            open(unit=lun2,file=fileid,status='replace',
     $                     form='formatted')
            write(lun2,'(2a,2i3,f12.5)') "# ",trim(spp%label), 
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i)
            write(lun2,'(a)')
     $           "#kb l, n (sequence number), Reference energy"
            call radial_dump_ascii(spp%pjnl(i),lun2,header=.false.)
            call io_close(lun2)
           endif ! (write_ion_plot_files)


         enddo

         write(lun,'(a)') "# Vna:__________________________"
         call radial_dump_ascii(spp%vna,lun)

        if (write_ion_plot_files) then
         write(fileid,'(a,a)') "VNA.", trim(spp%label)
         call io_assign(lun2)
         open(unit=lun2,file=fileid,status='replace',
     $        form='formatted')
         write(lun2,'(3a)') "# ",trim(spp%label), " Vna"
         call radial_dump_ascii(spp%vna,lun2,header=.false.)
         call io_close(lun2)
        endif ! (write_ion_plot_files)

         write(lun,'(a)') "# Chlocal:__________________________"
         call radial_dump_ascii(spp%chlocal,lun)

        if (write_ion_plot_files) then
         write(fileid,'(a,a)') "CHLOCAL.", trim(spp%label)
         call io_assign(lun2)
         open(unit=lun2,file=fileid,status='replace',
     $        form='formatted')
         write(lun2,'(3a)') "# ",trim(spp%label), " ChLocal"
         call radial_dump_ascii(spp%chlocal,lun2,header=.false.)
         call io_close(lun2)
!
!        Vlocal does not go to the .ion file
!
         write(fileid,'(a,a)') "RED_VLOCAL.", trim(spp%label)
         call io_assign(lun2)
         open(unit=lun2,file=fileid,status='replace',
     $        form='formatted')
         write(lun2,'(3a)') "# ",trim(spp%label), " Red_Vlocal"
         call
     $    radial_dump_ascii(spp%reduced_vlocal,lun2,header=.false.)
         call io_close(lun2)
        endif ! (write_ion_plot_files)

!
         if (spp%there_is_core) then
            write(lun,'(a)') "# Core:__________________________"
            call radial_dump_ascii(spp%core,lun)

           if (write_ion_plot_files) then
            write(fileid,'(a,a)') "CHCORE.", trim(spp%label)
            call io_assign(lun2)
            open(unit=lun2,file=fileid,status='replace',
     $           form='formatted')
            write(lun2,'(3a)') "# ",trim(spp%label), " ChCore"
            call radial_dump_ascii(spp%core,lun2,header=.false.)
            call io_close(lun2)
           endif ! (write_ion_plot_files)

         endif
         call io_close(lun)

         CONTAINS

         subroutine write_header(p,unit)
         
         type(species_info), pointer :: p
         integer, intent(in)         :: unit

         write(unit,'(a2,28x,a)') p%symbol, "# Symbol"
         write(unit,'(a20,10x,a)') p%label, "# Label"
         write(unit,'(i5,25x,a)') p%z, "# Atomic number"
         write(unit,'(g22.12,25x,a)') p%zval, "# Valence charge"
         write(unit,'(g22.12,4x,a)') p%mass, "# Mass "
         write(unit,'(g22.12,4x,a)') p%self_energy, "# Self energy "
         write(unit,'(2i4,22x,a)') p%lmax_basis, p%n_orbnl,
     $        "# Lmax for basis, no. of nl orbitals "
         write(unit,'(2i4,22x,a)') p%lmax_projs, p%n_pjnl,
     $        "# Lmax for projectors, no. of nl KB projectors "

         end subroutine write_header

      end subroutine dump_ion_ascii
!
!-----------------------------------------------------------------
      subroutine dump_basis_xml

      integer is

      do is = 1, nspecies
         call dump_ion_xml(is)
      enddo
      end subroutine dump_basis_xml
!---
      subroutine dump_ion_xml(is)

      integer, intent(in)      :: is

      type(species_info), pointer  ::   spp
      character*40 filename
      integer i, ispol, lun

      spp => species(is)
      if (spp%read_from_file) return    !! Do not dump

      write(filename,'(a,a)') trim(spp%label), ".ion.xml"
      call io_assign(lun)
      open(lun,file=filename,status='replace',form='formatted')

      write(lun,'(a)') '<ion version="0.1">'
      call xml_dump_element(lun,'symbol',str(spp%symbol))
      call xml_dump_element(lun,'label',str(spp%label))
      call xml_dump_element(lun,'z',str(spp%z))
      call xml_dump_element(lun,'valence',str(spp%zval))
      call xml_dump_element(lun,'mass',str(spp%mass))
      call xml_dump_element(lun,'self_energy',str(spp%self_energy))
      call xml_dump_element(lun,'lmax_basis',str(spp%lmax_basis))
      call xml_dump_element(lun,'norbs_nl',str(spp%n_orbnl))
      call xml_dump_element(lun,'lmax_projs',str(spp%lmax_projs))
      call xml_dump_element(lun,'nprojs_nl',str(spp%n_pjnl))

      write(lun,'(a)') '<preamble>'
      call write_basis_specs(lun,is)
      call pseudo_header_print(lun,
     $                   basis_parameters(is)%pseudopotential)
      write(lun,'(a)') '</preamble>'

      write(lun,'(a)') "<paos>"
         do i=1,spp%n_orbnl
            write(lun,'(a)') "<orbital "
            ispol = 0
            if (spp%orbnl_ispol(i)) ispol = 1
            call xml_dump_attribute(lun,'l',str(spp%orbnl_l(i)))
            call xml_dump_attribute(lun,'n',str(spp%orbnl_n(i)))
            call xml_dump_attribute(lun,'z',str(spp%orbnl_z(i)))
            call xml_dump_attribute(lun,'ispol',str(ispol))
            call xml_dump_attribute(lun,'population',
     $                            str(spp%orbnl_pop(i)))
            write(lun,'(a)') " >"
            call radial_dump_xml(spp%orbnl(i),lun)
            write(lun,'(a)') "</orbital>"
         enddo
      write(lun,'(a)') "</paos>"
         
      write(lun,'(a)') "<kbs>"
      do i=1,spp%n_pjnl
         write(lun,'(a)') "<projector "
         call xml_dump_attribute(lun,'l',str(spp%pjnl_l(i)))
         call xml_dump_attribute(lun,'n',str(spp%pjnl_n(i)))
         call xml_dump_attribute(lun,'ref_energy',
     $                                   str(spp%pjnl_ekb(i)))
         write(lun,'(a)') " >"
         call radial_dump_xml(spp%pjnl(i),lun)
         write(lun,'(a)') "</projector>"
       enddo
       write(lun,'(a)') "</kbs>"

       write(lun,'(a)') "<vna>"
       call radial_dump_xml(spp%vna,lun)
       write(lun,'(a)') "</vna>"

       write(lun,'(a)') "<chlocal>"
       call radial_dump_xml(spp%chlocal,lun)
       write(lun,'(a)') "</chlocal>"

       write(lun,'(a)') "<reduced_vlocal>"
       call radial_dump_xml(spp%reduced_vlocal,lun)
       write(lun,'(a)') "</reduced_vlocal>"

       if (spp%there_is_core) then
          write(lun,'(a)') "<core>"
          call radial_dump_xml(spp%core,lun)
          write(lun,'(a)') "</core>"
       endif

       write(lun,'(a)') "</ion>"

       call io_close(lun)

      end subroutine dump_ion_xml
!
!
!---------------------------------------------------------------------
c$$$      subroutine xml_element_content(xf,name,content)
c$$$      use flib_wxml
c$$$      type(xmlf_t), intent(inout) :: xf            ! For new XML output
c$$$      character(len=*), intent(in)  :: name
c$$$      character(len=*), intent(in)  :: content
c$$$
c$$$      call xml_NewElement(xf,name)
c$$$      call xml_AddPcdata(xf,content)
c$$$      call xml_EndElement(xf,name)
c$$$      end subroutine xml_element_content
c$$$!---------------------------------------------------------------------

      end module basis_io












