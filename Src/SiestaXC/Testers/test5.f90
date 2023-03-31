PROGRAM siestaXCtest5

  ! Finds and writes the integrand of the nonlocal vdW functional
  ! for fixed kf2 and kg2=grad(n2)/n2 in a mesh of kf1, kg1, and r12. 
  ! J.M.Soler. Jul.2012

  ! Used module procedures
  USE siestaXC, only: setXC
  USE m_vdwxc,  only: phiofr
  USE m_vdwxc,  only: vdw_get_qmesh
  USE m_vdwxc,  only: vdw_set_kcut
  USE m_vdwxc,  only: vdw_theta
  USE m_vv_vdwxc, only: vv_vdw_phi_val

  ! Used module parameters
  USE siestaXC, only: dp
  USE siestaXC, only: gp => grid_p

  implicit none

  ! Tester parameters
  integer, parameter:: nspin = 1      ! Number of spin components
  integer, parameter:: nfunc = 1      ! Number of functionals
  integer, parameter:: nr    = 6      ! Number of radial points
  real(dp),parameter:: r(nr) = (/0.1_dp, 1.0_dp, 2.0_dp, 5.0_dp, &
                                 10._dp, 20._dp/) ! |r1-r2| distances
  real(dp),parameter:: kcut  = 15._dp ! Planewave vector cutoff
  real(dp),parameter:: kf2   = 1.0_dp ! Fermi wavevector at point 2
  real(dp),parameter:: kg2   = 2.5_dp ! grad(n)/n at point 2
  real(dp),parameter:: kfmax = 5._dp  ! Max. Fermi wavevector at point 1
  real(dp),parameter:: kgmax = 5._dp  ! Max. grad(n)/n at point 1
  real(dp),parameter:: dkf   = 0.1_dp ! Fermi wavevector mesh interval
  real(dp),parameter:: dkg   = 0.1_dp ! grad(n)/n mesh interval
  character(len=*),parameter:: func = 'VDW'
  character(len=*),parameter:: auth = 'VV'  ! 'DRSLL'|'LMKLL'|'VV'
                                      ! KBM, C09, and BH use DRSLL kernel

  ! Tester variables and arrays
  integer :: ik, ikf, ikg, ir, nk, nkf, nkg
  real(dp):: kf1, kg1, gn1(3,nspin), gn2(3,nspin), n1(nspin), n2(nspin), &
             n1n2phi_exact, n1n2phi_interp, pi, r12, wc(nfunc), wx(nfunc)
  real(dp),allocatable:: &
             dtdn1(:,:), dtdn2(:,:), dtdgn1(:,:,:), dtdgn2(:,:,:), &
             n1phi(:), phi(:,:), theta1(:), theta2(:)

  ! Initialize XC functional
  wx = 1._dp/nfunc
  wc = 1._dp/nfunc
  call setXC( nfunc, func, auth, wx, wc )

  ! Set planewave cutoff
  call vdw_set_kcut( kcut )

  ! Find number of interpolation points and allocate interpolation arrays
  call vdw_get_qmesh( nk )
  allocate( dtdn1(nk,nspin), dtdn2(nk,nspin) )
  allocate( dtdgn1(3,nk,nspin), dtdgn2(3,nk,nspin) )
  allocate( n1phi(nk), phi(nk,nk), theta1(nk), theta2(nk) )

  ! Open file for output
  open(1,file='n1n2phi.table')

  ! Iterate on kf1, kg1, and r12
  nkf = nint(kfmax/dkf)
  nkg = nint(kgmax/dkg)
  do ir = 1,nr
    r12 = r(ir)
    call phiofr( r12, phi )
!    print'(a,/,(6f12.6))','phi(:,13)=',phi(:,13)
!    print'(a,/,(6f12.6))','phi(13,:)=',phi(13,:)
    do ikf = 0,nkf
      do ikg = 0,nkg
        kf1 = ikf*dkf
        kg1 = ikg*dkg
        kf1 = max(kf1,1.e-6)
        pi = acos(-1._dp)
        n1 = 0
        n2 = 0
        gn1 = 0
        gn2 = 0
        n1(1) = kf1**3/(3*pi**2)
        n2(1) = kf2**3/(3*pi**2)
        gn1(1,1) = n1(1)*kg1
        gn2(1,1) = n2(1)*kg2
        call vdw_theta( nspin, n1, gn1, theta1, dtdn1, dtdgn1 )
        call vdw_theta( nspin, n2, gn2, theta2, dtdn2, dtdgn2 )
!        if (ir==1 .and. ikf==0) then
!          print'(a,6f12.6)','r12,kf1,kg1=',r12,kf1,kg1
!          print'(a,/,(6f12.6))','theta1=',theta1
!        endif
        do ik = 1,nk
          n1phi(ik) = sum(theta1(:)*phi(:,ik))
        end do

        ! Assume phi_val returns phi, with which we want to compare
!        n1n2phi_interp = sum(n1phi*theta2)/n1(1)/n2(1)

        ! Assume phi_val returns sqrt(n1*n2)*phi
!        n1n2phi_interp = sum(n1phi*theta2) /sqrt(n1(1)*n2(1))

        ! Assume phi_val returns kf1*kf2*phi
!        n1n2phi_interp = sum(n1phi*theta2) *kf1*kf2/n1(1)/n2(1)

        ! Assume phi_val returns (kf1*kf2)**2*phi
        n1n2phi_interp = sum(n1phi*theta2) *(kf1*kf2)**2/n1(1)/n2(1)

        ! Assume phi_val returns n1*n2*phi
!        n1n2phi_interp = sum(n1phi*theta2)

        n1n2phi_exact = vv_vdw_phi_val( kf1, kf2, kg1, kg2, r12 )
        write(1,'(3f9.3,3e15.6)') r12, kf1, kg1, &
          n1n2phi_exact, n1n2phi_interp, n1n2phi_interp-n1n2phi_exact
      end do
    end do
  end do

  close(1)

END PROGRAM siestaXCtest5

