program test6

! Writes the xc energy density eps_xc for all available functionals
! J.M.Soler, May.2014

  use siestaXC, only: dp      ! double precision kind
  use siestaXC, only: ldaxc   ! LDA xc
  use siestaXC, only: ggaxc   ! GGA xc

! Set mesh parameters
  implicit none
  integer, parameter::   &
    irel = 0,            &! relativistic exchange? 0=>no, 1=>yes
    nspin = 2,           &! number of spin polarization components (1|2)
    nz  = 2,             &! number of spin polarization values
    nsx = 20,            &! number of s mesh values for Fx(s)
    nkf = 20,            &! number of kf mesh values for epsc(kf,kg)
    nkg = 30              ! number of kg mesh values for epsc(kf,kg)
  real(dp),parameter::   &
    smin = 0.0_dp,       &! min. value of s in Fx(s) interpol.
    smax = 1.e2_dp,      &! max. value of s in Fx(s) interpol.
    kfmin = 1.e-3_dp,    &! min. value of kf in epsc(kf,kg) interpol.
    kfmax = 10.0_dp,     &! max. value of kf in epsc(kf,kg) interpol.
    kgmin =  0.0_dp,     &! min. value of kg in epsc(kf,kg) interpol.
    kgmax = 10.0_dp,     &! max. value of kg in epsc(kf,kg) interpol.
    dsnds1 = 1.e3_dp,    &! (s(n)-s(n-1))/(s(2)-s(1))
    dkfndkf1 = 1.e3_dp,  &! (kf(n)-kf(n-1))/(kf(2)-kf(1))
    dkgndkg1 = 2.0_dp     ! (kg(n)-kg(n-1))/(kg(2)-kg(1))

  ! List of functionals to be tested
  integer, parameter:: nfTot = 18      ! total number of functionals
  integer, parameter:: nf = nfTot-2    ! GGA functionals
  integer:: indexf(nf) = (/3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

  ! All functionals available
  !                  1,           2,          3,           4,   
  !                  5,           6,          7,           8, 
  !                  9,          10,         11,          12,
  !                 13,          14,         15,          16,
  !                 17,          18
  character(len=3):: &
    func(nfTot)=(/'LDA',       'LDA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'GGA',       'GGA', &
                  'GGA',       'GGA'       /)
  character(len=10):: &
    auth(nfTot)=(/'PZ        ','PW92      ','PW91      ','PBE       ', &
                  'RPBE      ','revPBE    ','LYP       ','WC        ', &
                  'PBEJsJrLO ','PBEJsJrHEG','PBEGcGxLO ','PBEGcGxHEG', &
                  'PBESOL    ','AM05      ','PW86      ','B88       ', &
                  'C09       ','BH        '/) 

  integer:: &
    ifunc, ikf, ikg, isx, iz, jf, mspin
  real(dp):: &
    a, b, dExdD(nspin), dEcdD(nspin), dExdGD(3,nspin), dEcdGD(3,nspin), &
    dVcdD(nspin,nspin), dVxdD(nspin,nspin), &
    epsc(nkf,nkg,nz,nf), epsc0, epsc1, epsx(nkf,nkg,nz,nf), epsx0, epsx1, &
    fx(nsx,nf), fxc(nkf,nkg,nz,nf), grho(3,nspin), &
    kf, kfmesh(nkf), kg, kgmesh(nkg), &
    pi, rho(nspin), rhoTot, s, smesh(nsx), zeta, zmesh(nz)

! Set mesh values
  do iz = 1,nz
    zmesh(iz) = real(iz-1,dp)/(nz-1)
  enddo
  a = log(dsnds1) / (nsx-1)
  b = smax / (exp(a*(nsx-1)) - 1)
  do isx = 1,nsx
    smesh(isx) = b*exp(a*(isx-1)) - b
  enddo
  a = log(dkfndkf1) / (nkf-1)
  b = (kfmax-kfmin) / (exp(a*(nkf-1)) - 1)
  do ikf = 1,nkf
    kfmesh(ikf) = kfmin + b*exp(a*(ikf-1)) - b
  enddo
  a = log(dkgndkg1) / (nkg-1)
  b = (kgmax-kgmin) / (exp(a*(nkg-1)) - 1)
  do ikg = 1,nkg
    kgmesh(ikg) = b*exp(a*(ikg-1)) - b
  enddo

! Loop on functionals
  do jf = 1,nf
    ifunc = indexf(jf)

    ! Find exchange enhancement factors at s mesh values
    pi = acos(-1.0_dp)
    kf = 1                                ! Fermi wavevector
    rho(1) = kf**3/3/pi**2                ! electron density
    epsx0 = -3*kf/4/pi                    ! LDA exchange energy
    do isx = 1,nsx
      s = smesh(isx)                      ! |grad(rho)|/(2*rho*kf)
      kg = 2*kf*s                         ! |grad(rho)|/rho
      grho = 0                            ! grad(rho)
      grho(3,1) = rho(1)*kg
      mspin = 1                           ! unpolarized to call ggaxc
      call ggaxc( auth(ifunc), irel, mspin, rho, grho, &
                  epsx1, epsc1, dEXdD, dECdD, dEXdGD, dECdGD )
      fx(isx,jf) = epsx1/epsx0            ! exchange enhancement factor
    enddo

    ! Find exchange-correlation energies at (kf,kg) mesh
    do ikf = 1,nkf
      do ikg = 1,nkg
        do iz = 1,nz
          kf = kfmesh(ikf)
          kg = kgmesh(ikg)
          zeta = zmesh(iz)
          rhoTot = kf**3/3/pi**2
          rho(1) = rhoTot*(1+zeta)/2
          rho(2) = rhoTot*(1-zeta)/2
          grho = 0
          grho(3,1:nspin) = rho(1:nspin)*kg
          call ldaxc( 'PW92', irel, nspin, rho, &
                      epsx0, epsc0, dEXdD, dECdD, dVxdD, dVcdD )
          call ggaxc( auth(ifunc), irel, nspin, rho, grho, &
                      epsx1, epsc1, dEXdD, dECdD, dEXdGD, dECdGD )
          epsx(ikf,ikg,iz,jf) = epsx1
          epsc(ikf,ikg,iz,jf) = epsc1
          fxc(ikf,ikg,iz,jf) = (epsx1+epsc1)/epsx0
        enddo
      enddo
    enddo

  enddo ! jf

! Write results
  open(2,file='funcs.out')
  write(2,'(20a14)') (trim(auth(indexf(jf))),jf=1,nf)
  close(2)

  open(2,file='fx.out')
!  write(2,'(a14,2x,20a14)') 's',(trim(auth(indexf(jf))),jf=1,nf)
  do isx = 1,nsx
    write(2,'(e14.6,2x,20e14.6)') smesh(isx),(fx(isx,jf),jf=1,nf)
  enddo
  close(2)

  open(2,file='epsx.out')
!  write(2,'(3a14,2x,20a14)') 'kf','kg','zeta',(trim(auth(indexf(jf))),jf=1,nf)
  do iz = 1,nz
    do ikg = 1,nkg
      do ikf = 1,nkf
        write(2,'(3e14.6,2x,20e14.6)') &
          kfmesh(ikf), kgmesh(ikg), zmesh(iz), (epsx(ikf,ikg,iz,jf),jf=1,nf)
      enddo
    enddo
  enddo
  close(2)

  open(2,file='epsc.out')
!  write(2,'(3a14,2x,20a14)') 'kf','kg','zeta',(trim(auth(indexf(jf))),jf=1,nf)
  do iz = 1,nz
    do ikg = 1,nkg
      do ikf = 1,nkf
        write(2,'(3e14.6,2x,20e14.6)') &
          kfmesh(ikf), kgmesh(ikg), zmesh(iz), (epsc(ikf,ikg,iz,jf),jf=1,nf)
      enddo
    enddo
  enddo
  close(2)

  open(2,file='fxc.out')
!  write(2,'(3a14,2x,20a14)') 'kf','kg','zeta',(trim(auth(indexf(jf))),jf=1,nf)
  do iz = 1,nz
    do ikg = 1,nkg
      do ikf = 1,nkf
        write(2,'(3e14.6,2x,20e14.6)') &
          kfmesh(ikf), kgmesh(ikg), zmesh(iz), (fxc(ikf,ikg,iz,jf),jf=1,nf)
      enddo
    enddo
  enddo
  close(2)

end program test6

