
! Calculation of the spherical orbital transformation matrices.
! Literature
!     H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)
!
! Shanghui modified in atm_transfer following cp2k code, 2010
!
! Modified by xmqin for input from fdf format, December, 2013 
! output NAO2GTO for user's check
!
!  We call dgemm()  to calculate c2s (Cartesian to Spherical )

      subroutine nao2gto_transfer

      use precision, only: dp
      use atm_types, only: maxnorbs, nspecies
      use atm_types, only: species, species_info
      
      use atm_types, only: ncon_max,l_max,nco,nso,ncosum,co,coset,indco, &
                           orbtramat, maxn_orbnl,maxn_contract
      use atmfuncs,  only: lofio,mofio,labelfis, rcut
      use chemical,  only: species_label
      use mathconstants
      use sys,       only: die 
      use fdf
      use parsing 
      use nao2gto_util, only: exp_radius, exp_radius_very_extended
      use radial,    only: rad_get, rad_setup_d2
      use radial,    only: rad_func, rad_alloc
      use alloc,     only: re_alloc, de_alloc
      use hfx_types, only: hfx_input_parameter, read_hfx_info
      use nao2gto_nonlin
     
!----------------------------------------------------------------
!
      implicit none

      type(species_info), pointer        :: spp

      integer is, io, i ,l, m, n, j, iunit
      integer max_norbnl, izeta
      integer, dimension(maxnorbs) :: index_normal, z_normal,&
               nsm_normal, index_pol, z_pol, nsm_pol

      integer  ipgf,jpgf,lx,ly,lz,nco_max,nso_max,ico
      real(dp) expzet,fnorm,zeta,zetb,prefac,gcca,gccb
      real(dp) radius_pgf, radius_shell, radius_kind 
      integer  nao2gto_contract(maxn_orbnl,nspecies)
      real(dp) nao2gto_zeta(maxn_contract,maxn_orbnl,nspecies), &
               nao2gto_coefficient(maxn_contract,maxn_orbnl,nspecies)
      integer i_cphi,num_cphi(maxnorbs) ! maxnorb=100,so cphi is always less than 100.
      character*50  relabel
      integer  num_gaus,  ren, reorb , rel, rezeta
      logical  found
     
      integer ::  inlz, jnlz, irad, errno
      real(dp) :: dummy
      real(dp), dimension(:,:),   pointer :: fit_trial => null()
      real(dp), dimension(:),     pointer :: rad_pts => null()
      real(dp), dimension(:),     pointer :: orb_pts => null()
      real(dp), dimension(:),     pointer :: fval => null()
      real(dp), dimension(:),     pointer :: orb_dr => null()
  
      type(rad_func), pointer      :: func => null()
      type(hfx_input_parameter) :: hfx_opts

!e-----------------------for nco,nso,co,coset,orbtramat etc------------------
        ncosum(-1:l_max)=0
        DO l=0,l_max
          nco(l) = (l + 1)*(l + 2)/2
          nso(l) = 2*l+1
          ncosum(l) = ncosum(l-1) + nco(l)
        END DO
     
        DO lx=0,l_max
          DO ly=0,l_max
            DO lz=0,l_max
              l = lx + ly + lz
              IF (l > l_max) CYCLE
              co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
              coset(lx,ly,lz) = ncosum(l-1) + co(lx,ly,lz)
            END DO
          END DO
        END DO
     
         indco(:,:) = 0
     
        DO l=0,l_max
          DO lx=0,l
            DO ly=0,l-lx
              lz = l - lx - ly
              indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
            END DO
          END DO
        END DO

        allocate (orbtramat(0:l_max))
        do l=0,l_max
         nco_max = nco(l)
         nso_max = nso(l)
        allocate (orbtramat(l)%c2s(nso_max,nco_max))
        enddo
        call  calc_c2s_matrix(orbtramat,l_max,co)

!---------------------end for nco,nso,co,coset,orbtramat etc-----------------

      hfx_opts%dump_fit_data = fdf_boolean("HFX.DumpFitData",.true.)
      hfx_opts%npts_fit  = fdf_integer( "HFX.FitDataPoints",500)
      hfx_opts%max_num_gaus_s = fdf_integer("HFX.MaximumNumberGaussians-S",6)
      hfx_opts%max_num_gaus_p = fdf_integer("HFX.MaximumNumberGaussians-P",6)
      hfx_opts%max_num_gaus_d = fdf_integer("HFX.MaximumNumberGaussians-D",6)
      hfx_opts%tolerance_gaus = fdf_double( "HFX.ToleranceFit", 1.d-3)
      hfx_opts%threshold_exp_gaus = fdf_double("HFX.SeparationExponents", 1.4d0)
      hfx_opts%gto_eps = fdf_double("HFX.GaussianEPS", 1.0d-4)

      
      nao2gto_contract = 0
      nao2gto_zeta = 0.0_dp
      nao2gto_coefficient  = 0.0_dp
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

      found = fdf_block('NAO2GTO',iunit)
      if (found) then
        write(6,*) "Read NAO2GTO from input fdf file"
!      endif
        do is=1,nspecies
           spp => species(is)
           read(iunit,*) relabel
         do i=1,spp%n_orbnl
           read(iunit,*)   ren, rel ,rezeta, nao2gto_contract(i,is)
           
!           write(6,*) is ,

           do j=1,nao2gto_contract(i,is)
         
             read(iunit,*)  nao2gto_zeta(j,i,is), &
                      nao2gto_coefficient(j,i,is)
         !    write(6,*) nao2gto_zeta(j,i,is),nao2gto_coefficient(j,i,is)
           enddo
         enddo
        enddo

      else

!       call read_hfx_info(hfx_opts)

      call re_alloc(rad_pts, 1,  hfx_opts%npts_fit, 'rad_pts', &
        'nao2gto_transfer')
      call re_alloc(orb_pts, 1,  hfx_opts%npts_fit, 'orb_pts', &
        'nao2gto_transfer')
!     Allocate the temporal variable that will contain the
!     coefficients of the Gaussians that fit a NAO
      call re_alloc(fit_trial, 1, 2, 1,  hfx_opts%max_num_gaus_s, 'nao2gto_transfer')


      do is=1,nspecies
           spp => species(is)
           spp%label = species_label(is)

           inlz = 0

          do io = 1, spp%norbs
!       Identify the angular momentum quantum number
             l = spp%orb_l(io)
             m = spp%orb_m(io)

             if ( m .ne. -l ) cycle           

             inlz = inlz + 1
             
             do irad=1, hfx_opts%npts_fit
                rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
     &                 (real(hfx_opts%npts_fit, dp) - 1)
                call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)
             enddo
              orb_pts(hfx_opts%npts_fit) = 0.0_dp
              fit_trial(:,:) = 0.0_dp
              if( l == 0)  num_gaus = hfx_opts%max_num_gaus_s
              if( l == 1)  num_gaus = hfx_opts%max_num_gaus_p
              if( l == 2)  num_gaus = hfx_opts%max_num_gaus_d
!              if (spp%orbnl_ispol(inlz)) hfx_opts%max_num_gaus = 3
              write(6,*) inlz, l, spp%orbnl_ispol(inlz), num_gaus 

              call nao2gto_gaussfit( num_gaus , rad_pts, orb_pts, &
     &                               hfx_opts, fit_trial, errno )
              ipgf = 0
          do jnlz = 1, num_gaus
!!           For debugging
            write(6,'(a,3i5,3f12.5,i5,l5)') &
 &            'nao2gto_transfer: species, orbital, gauss, exp, coef, epsilon, error = ', &
 &            is, inlz, jnlz, fit_trial(1,jnlz), fit_trial(2,jnlz), &
 &            epsilon(1.0_dp), errno, abs(fit_trial(1,jnlz)) < epsilon(1.0_dp)
!!           End debugging
            if ( abs(fit_trial(1,jnlz)) < epsilon(1.0_dp) ) exit
            ipgf = jnlz
          end do
          if ( ipgf == 0 ) errno = 10

          if ( errno == 0 ) then
              nao2gto_contract(inlz,is) = ipgf
              do jnlz = 1, ipgf
                nao2gto_zeta(jnlz,inlz,is)        = fit_trial(1,jnlz)
                nao2gto_coefficient(jnlz,inlz,is) = fit_trial(2,jnlz)
              end do
          else
            write(*, fmt='("GAUSSIAN FITTING ERROR ",I4.4)') errno
          end if

         enddo
       enddo

      call de_alloc( rad_pts,   'rad_pts',   'nao2gto_transfer' )
      call de_alloc( orb_pts,   'orb_pts',   'nao2gto_transfer' )
      call de_alloc( fit_trial, 'fit_trial', 'nao2gto_transfer' )


      endif


!        do is=1,nspecies
!           spp => species(is)
!           i = 0
!           do io = 1, spp%norbs
!              l=spp%orb_l(io)
!              m=spp%orb_m(io)
!              if (m .ne. -l) cycle     ! not a normal orbital
!              i=i+1
!              spp%orbnl_contract(i)= nao2gto_contract(i,is)
!              expzet = 0.25_dp*REAL(2*l+3, dp)
!              prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
!              DO n = 1, spp%orbnl_contract(i)
!                 gcca = nao2gto_coefficient(n,i,is)
!                 zeta = nao2gto_zeta(n,i,is)
!                 spp%orbnl_zeta(n,i)= zeta
!                 spp%orbnl_coefficient(n,i) =prefac*zeta**expzet*gcca
!               END DO
!            spp%orbnl_adjoined_zeta(i) = &
!                     minval(spp%orbnl_zeta(1:spp%orbnl_contract(i),i))
!             END DO
!           END DO


      do is=1,nspecies
         spp => species(is)
!-------------- GTOs for normal orbitals    ---------------
          i=0
         do io = 1, spp%norbs
             l=spp%orb_l(io)
             m=spp%orb_m(io)
             if (m .ne. -l) cycle     ! not a normal orbital
            i=i+1
            spp%orbnl_contract(i)= nao2gto_contract(i,is)
            do ipgf=1,spp%orbnl_contract(i)
               spp%orbnl_zeta(ipgf,i)= nao2gto_zeta(ipgf,i,is)
               spp%orbnl_coefficient(ipgf,i)= nao2gto_coefficient(ipgf,i,is)
            enddo
            spp%orbnl_adjoined_zeta(i) = &
                     minval(spp%orbnl_zeta(1:spp%orbnl_contract(i),i))
         enddo
! xmqin add cutoff radius for primitive gto and contracted gto, actually
! PAO
! PAO itself has cutoff radius in siesta, I set shell cutoff just for
! comparison.
        spp%pgf_radius(1:maxn_contract,1:maxn_orbnl)=0.0d0
        spp%shell_radius(1:maxn_orbnl)=0.0d0
         i = 0
         do io = 1, spp%norbs
             l=spp%orb_l(io)
             m=spp%orb_m(io)
             if (m .ne. -l) cycle     ! not a normal orbital
            i=i+1
            spp%orbnl_contract(i)= nao2gto_contract(i,is)
            do ipgf=1,spp%orbnl_contract(i)
               gcca =  nao2gto_coefficient(ipgf,i,is)
               zeta =  nao2gto_zeta(ipgf,i,is)
!               if(zeta .le.0.1_dp) then
!               spp%pgf_radius(ipgf,i)= & 
!                        exp_radius(l,zeta, 1.0d-2, 1.0_dp)
!               else
               spp%pgf_radius(ipgf,i)= &
                        exp_radius(l,zeta, hfx_opts%gto_eps, 1.0_dp)
!               endif

          write(6,'(a,3i5,f12.5)')                                           &
 &          'nao2gto_transfer: isp, inlz, ipgf, pgf_radius = ',              &
 &            is, i, ipgf, spp%pgf_radius(ipgf,i)
!!         End debugging

            enddo
            spp%shell_radius(i) = &
               maxval(spp%pgf_radius(1:spp%orbnl_contract(i),i))
        write(6,'(a,2i5,2f12.5)')                                              &
 &        'nao2gto_transfer: is, inlz, fittedcut, origincut = ',                     &
 &          is, i, spp%shell_radius(i), rcut(is,io)

         enddo
         spp%kind_radius=maxval(spp%shell_radius(1:maxn_orbnl))

!----------------------------this is for spp%orb_*_cphi----------------------

         i_cphi=1
         do io = 1, spp%norbs
           if(spp%orb_m(io).lt.2) then
                num_cphi(i_cphi)=io
                i_cphi=i_cphi+1
           else if(spp%orb_m(io).eq.2) then
                num_cphi(i_cphi)=io
                num_cphi(i_cphi+1)=io
                i_cphi=i_cphi+2
           else if(spp%orb_m(io).eq.3) then
                num_cphi(i_cphi)=io
                num_cphi(i_cphi+1)=io
                num_cphi(i_cphi+2)=io
                num_cphi(i_cphi+3)=io
                i_cphi=i_cphi+4
          endif
         enddo
         i_cphi=i_cphi-1 ! total number of cphi   

!----------------------------this is for spp%norbscphi----------------------
! Here, we only consider the case that lmax <= 3 (s, p, d, and f)
         if( spp%lmax_basis.lt.2) then
         spp%norbs_cphi=spp%norbs   !s,p
         else 
         spp%norbs_cphi=i_cphi !d,f
         endif


         do i_cphi = 1, spp%norbs_cphi
         io=num_cphi(i_cphi)
         spp%orb_n_cphi(i_cphi)=spp%orb_n(io)
         spp%orb_l_cphi(i_cphi)=spp%orb_l(io)
         spp%orb_index_cphi(i_cphi)=spp%orb_index(io)
         enddo



        io=0
        do i=1,spp%n_orbnl
           l=spp%orbnl_l(i)
     
           if(i.eq.1) then
           spp%orbnl_index_cphi(i)=1
           spp%orbnl_index_sphi(i)=1
           else
           spp%orbnl_index_cphi(i)=spp%orbnl_index_cphi(i-1) &
                                  +nco(spp%orbnl_l(i-1))
           spp%orbnl_index_sphi(i)=spp%orbnl_index_sphi(i-1) &
                                  +nso(spp%orbnl_l(i-1))
           endif
     
           do ico=ncosum(l-1)+1,ncosum(l)
           io=io+1
           spp%orb_cphi_lx(io)=indco(1, ico) 
           spp%orb_cphi_ly(io)=indco(2, ico)
           spp%orb_cphi_lz(io)=indco(3, ico)
           enddo
        enddo

!-----------------------this is for spp%norm_cphi------------------------------
         do io = 1, spp%norbs_cphi
          l = spp%orb_l_cphi(io) 

          expzet = 0.5_dp*REAL(2*l + 3,dp)

          fnorm = 0.0_dp
          i = spp%orb_index_cphi(io)
            DO ipgf=1,spp%orbnl_contract(i)
             gcca = spp%orbnl_coefficient(ipgf,i)
             zeta = spp%orbnl_zeta(ipgf,i)
              DO jpgf=1,spp%orbnl_contract(i)
                gccb = spp%orbnl_coefficient(jpgf,i)
                zetb = spp%orbnl_zeta(jpgf,i)
                fnorm = fnorm + gcca*gccb/(zeta + zetb)**expzet
              END DO
             END DO

          fnorm = 0.5_dp**l*pi**1.5_dp*fnorm

             lx = spp%orb_cphi_lx(io)
             ly = spp%orb_cphi_ly(io)
             lz = spp%orb_cphi_lz(io)
             prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
             spp%norm_cphi(io) = 1.0_dp/SQRT(prefac*fnorm)
         enddo

!---------------------------------calc cphi and sphi----------------------------
        spp%cphi(1:ncon_max,1:spp%norbs_cphi)=0.0d0
        spp%sphi(1:ncon_max,1:spp%norbs)=0.0d0
   
           !write(6,*) 'lx,ly,lz,co(lx,ly,lz),io'
        do io=1,spp%norbs_cphi 
           i=spp%orb_index_cphi(io) 
           lx=spp%orb_cphi_lx(io)
           ly=spp%orb_cphi_ly(io)
           lz=spp%orb_cphi_lz(io)
           !write(6,*) lx,ly,lz,co(lx,ly,lz),io
           do j=1,spp%orbnl_contract(i)
          spp%cphi((j-1)*nco(spp%orb_l_cphi(io))+co(lx,ly,lz),io) =  &
                 spp%orbnl_coefficient(j,i) * spp%norm_cphi(io)
           enddo
       enddo
    
       call  cphi2sphi(ncon_max,spp%norbs_cphi,spp%norbs,l_max,&
                      spp%n_orbnl,spp%orbnl_l,nco,nso,&
                      spp%orbnl_index_cphi,spp%orbnl_index_sphi,&
                      spp%cphi,spp%sphi,orbtramat)

          i=0
         do io = 1, spp%norbs
             l=spp%orb_l(io)
             m=spp%orb_m(io)
             if (m .ne. -l) cycle     ! not a normal orbital
            i=i+1
        !   write(72,*) spp%sphi(1:nco(l)*spp%orbnl_contract(i),io)
            spp%orbnl_contraction_coeff(i)=  &
                  maxval((/(sum(abs(spp%sphi(1:nco(l)*spp%orbnl_contract(i), &
                          j))),j=io,io+nso(l)-1)/))
!            spp%orbnl_contraction_coeff(i)=  &
!                  maxval((/(sum(abs(spp%sphi(1:nco(l)*spp%orbnl_contract(i), &
!                          io))))/))
         enddo


!---------------------------------end calc cphi and sphi-------------------
!
!
       enddo  ! for is =1, nspecies


         write(6,'(A)') "%block NAO2GTO"                
!         write(6,*) 'the order of index norbs'
         do is=1,nspecies
            spp=>species(is)
            spp%label = species_label(is)
!             write(6,'(a)')  &
!                "#(species label, l, n, z, is_polarized, popul)"
!             write(6,*)'Number of GTOs'
!             write(6,*) '        zeta          coefficent'
             write(6,'(A,1X,I3)') trim(spp%label), spp%n_orbnl
          do io = 1,  spp%n_orbnl
!             write(6,'(a)')  &
!                "#(species label, l, n, z, is_polarized, popul)"
!             write(6,*) io
!             write(6,'(2a,4i3,f10.4)') "# ", trim(spp%label), &
!                spp%orbnl_l(io), spp%orbnl_n(io), spp%orbnl_z(io),&
!                spp%orbnl_contract(io)
!             i=spp%orb_index(io)
             write(6,'(I1,3(1X,I2))') spp%orbnl_n(io), spp%orbnl_l(io), &
                     spp%orbnl_z(io), spp%orbnl_contract(io)
          
!             write(6,*)  spp%orbnl_contract(io)
           
           do n=1,spp%orbnl_contract(io)
              write(6,'(2x,2f16.9)') spp%orbnl_zeta(n,io),&
                              spp%orbnl_coefficient(n,io)
           enddo
          enddo
         enddo
           write(6,'(A)') "%endblock NAO2GTO"



    hfx_opts%is_fitted_nao = fdf_boolean("HFX.UseFittedNAOs", .true.)

    if(hfx_opts%is_fitted_nao) then
    write(6,*) "Use fitted NAOs as new basis"
    call re_alloc(rad_pts, 1, hfx_opts%npts_fit, 'rad_pts', &
                  'nao2gto_itransfer')

    do is=1,nspecies
      spp => species(is)

      do inlz = 1,spp%n_orbnl

        l = spp%orbnl_l(inlz)
         expzet = 0.5_dp*REAL(2*l+3, dp)
         prefac = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
        fnorm = 0.0_dp
        do ipgf =1,spp%orbnl_contract(inlz)
             gcca = spp%orbnl_coefficient(ipgf,inlz)
             zeta = spp%orbnl_zeta(ipgf,inlz)
             do jpgf=1, spp%orbnl_contract(inlz)
                gccb = spp%orbnl_coefficient(jpgf,inlz)
                zetb = spp%orbnl_zeta(jpgf,inlz)
                fnorm = fnorm+gcca*gccb*prefac/(zeta+zetb)**expzet
             enddo
        enddo
        fnorm = 1.0_dp/SQRT(fnorm)

        func => spp%orbnl(inlz)
!         call reset_rad_func(func)
!         spp%orbnl(inlz)%cutoff = gtos(isp)%shell_radius(inlz)
!         spp%orbnl(inlz)%delta =
!         spp%orbnl(inlz)%cutoff/(dble(hfx_opts%npts_fit-1)+1.0d-20)
        func%cutoff = spp%shell_radius(inlz)
        func%n = hfx_opts%npts_fit
        func%delta =func%cutoff/(dble(func%n-1)+1.0d-20)
!        write(6,*) "cutoff 1 :: ", func%cutoff, func%n, func%delta
!        write(6,*) "cutoff 1 :: ",
!        gtos(isp)%shell_radius(inlz),spp%orbnl(inlz)%cutoff
!        write(6,*) "n, delta 1 :: ",  spp%orbnl(inlz)%n, spp%orbnl(inlz)%delta

!        orb_pts(hfx_opts%npts_fit) = 0.0_dp

        do irad=1,func%n
           rad_pts(irad) = func%cutoff * real(irad-1, dp) / &
&                 (real(func%n, dp) - 1)
!          call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)

        enddo

        func%f(:) = 0.0_dp

        do jnlz=1,spp%orbnl_contract(inlz)
           func%f(:) = func%f(:) + spp%orbnl_coefficient(jnlz,inlz) *  &
 &                    exp(-rad_pts(:)**2*spp%orbnl_zeta(jnlz,inlz))
        enddo
        func%f(func%n) = 0.0_dp
        func%f = func%f*fnorm

        call rad_setup_d2(func,yp1=0.0_dp,ypn=huge(1.0_dp))

       enddo
     enddo
      call de_alloc( rad_pts,   'rad_pts',   'nao2gto_transfer' )
   endif

      end subroutine nao2gto_transfer

      subroutine calc_c2s_matrix(orbtramat,l_max,co)
        USE kinds,                           ONLY: dp
        USE mathconstants
     
        implicit none
        !------------input and output-------------------------
        integer ns,nc,l_max,co(0:l_max,0:l_max,0:l_max)
          !real(kind=dp)  c2s(ns,nc,0:l_max)
        TYPE orbtramat_type
        REAL(KIND = dp), DIMENSION(:,:), POINTER :: c2s
        END TYPE orbtramat_type
     
        TYPE(orbtramat_type), DIMENSION(0:l_max) :: orbtramat
     
        !-----------local variables---------------------------
        INTEGER   expo,i,ic,is,j,k,l,lx,ly,lz,m,ma
        REAL(KIND=dp)  s, s1, s2, binomial
     
     
        DO l=0,l_max
     
!     *** Build the orbital transformation matrix for the     ***
!     *** transformation from Cartesian to spherical orbitals ***
!     *** (c2s, formula 15)                                   ***

         DO lx=0,l
           DO ly=0,l-lx
             lz = l - lx - ly
             ic = co(lx,ly,lz)
             DO m=-l,l
               is = l + m + 1
               ma = ABS(m)
               j = lx + ly - ma
               IF ((j >= 0).AND.(MODULO(j,2) == 0)) THEN
                 j = j/2
                 s1 = 0.0_dp
                 DO i=0,(l-ma)/2
                   s2 = 0.0_dp
                   DO k=0,j
                     IF (((m < 0).AND.(MODULO(ABS(ma-lx),2) == 1)).OR. &
                     ((m > 0).AND.(MODULO(ABS(ma-lx),2) == 0))) THEN
                       expo = (ma - lx + 2*k)/2
                       s = (-1.0_dp)**expo*SQRT(2.0_dp)
                     ELSE IF ((m == 0).AND.(MODULO(lx,2) == 0)) THEN
                       expo = k - lx/2
                       s = (-1.0_dp)**expo
                     ELSE
                       s = 0.0_dp
                     END IF
                     s2 = s2 + binomial(j,k)*binomial(ma,lx-2*k)*s
                   END DO
                   s1 = s1 + binomial(l,i)*binomial(i,j)* &
                        (-1.0_dp)**i*fac(2*l-2*i)/fac(l-ma-2*i)*s2
                 END DO
                 orbtramat(l)%c2s(is,ic) = &
             SQRT((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l-ma))/ &
                  (fac(lx)*fac(ly)*fac(lz)*fac(2*l)*fac(l+ma)))*s1/ &
                  (2.0_dp**l*fac(l))*(-1)**ma
               ELSE
                   orbtramat(l)%c2s(is,ic) = 0.0_dp
               END IF
             END DO
           END DO
         END DO
       END DO
      end subroutine calc_c2s_matrix
    
   
      subroutine cphi2sphi(ncon_max,norb_cphi,norb_sphi,l_max,norb_nl,&
                           orbnl_l,nco,nso,index_cphi,index_sphi,   &
                           cphi,sphi, orbtramat)
    
       USE kinds,                           ONLY: dp
       implicit none
       !----------input variables---------------------------------
       integer  ncon_max,norb_cphi,norb_sphi,l_max,norb_nl
       integer  orbnl_l(norb_nl),nco(0:l_max),nso(0:l_max),&
                index_cphi(norb_nl),index_sphi(norb_nl)
       real(dp) cphi(ncon_max,norb_cphi),sphi(ncon_max,norb_sphi)
                !c2s(nso(l_max),nco(l_max),0:l_max)    
    
       TYPE orbtramat_type
       REAL(KIND = dp), DIMENSION(:,:), POINTER :: c2s
       END TYPE orbtramat_type
    
       TYPE(orbtramat_type), DIMENSION(0:l_max) :: orbtramat
       !-----------local variables-------------------------------- 
       integer i,l,ncgf,nsgf,first_cgf,first_sgf     
       
!       write(6,*) "testsphi",norb_sphi,norb_cphi
    
       do i=1,norb_nl
          l=orbnl_l(i)
          ncgf=nco(l)
          nsgf=nso(l)   
          first_cgf=index_cphi(i)
          first_sgf=index_sphi(i) 
          CALL dgemm("N","T",ncon_max,nsgf,ncgf,&
                    1.0_dp,cphi(1,first_cgf),ncon_max,&
                    orbtramat(l)%c2s(1,1),nsgf,&
                    0.0_dp,sphi(1,first_sgf),ncon_max)
    
       enddo
    
      end subroutine cphi2sphi 
    
    
      FUNCTION binomial(n,k) RESULT(n_over_k)
    
       USE kinds,                           ONLY: dp
       USE mathconstants
       implicit none
       
       INTEGER, INTENT(IN)                      :: n, k
       REAL(KIND=dp)                            :: n_over_k
    
       IF ((k >= 0).AND.(k <= n)) THEN
          n_over_k = fac(n)/(fac(n-k)*fac(k))
       ELSE
          n_over_k = 0.0_dp
       END IF
    
      END FUNCTION binomial
