      subroutine radial_nao2gto(op, n_contract, zeta,
     $                                coeff, Rcut)
      type(rad_func)    :: op
      integer           :: n_contract
      real(dp)          :: zeta(maxn_contract), coeff(maxn_contract)
      real(dp)          :: Rcut

      integer  :: j, m
      real(dp) ::  x, y, dy, s_error, new_delta, new_Rcut
!      logical  :: print_header
!
!     The standard dump is to unit "lun"
!     and includes a header with npts, delta, and cutoff
!

!      if (print_header) then
!          write(lun,'(i4,2g22.12,a)') op%n,
!     $          op%delta, op%cutoff, " # npts, delta, cutoff"
!      endif

       delt=Rcut/(dble(ntbmax-1)+1.0d-20)

!      s_error = 0.d0
      do j=1,op%n
         x = (j-1)*delta
         y =0.0_dp
         do m = 1, n_contract
            y = y+coeff(m)*dexp(-zeta(m)*x**2)
         enddo

         write(6,'(2g22.12)') x, y
       !  write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
      !   s_error = s_error + (op%f(j)-y)**2
      enddo

!      write(lun, *) "ERROR:", s_error

      end subroutine radial_nao2gto_ascii

