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
      module m_cell
        use precision, only: dp
        use siesta_geom, only: ucell
        implicit none
        real(dp), public, save  :: celli(3,3) = 0.0_dp

        public :: cart2frac, frac2cart
        public :: write_canonical_ucell

        private

        CONTAINS

        subroutine cart2frac(cart,frac)
        real(dp), intent(in)  :: cart(3)
        real(dp), intent(out)  :: frac(3)

        frac =  matmul(transpose(celli),cart)
        end subroutine cart2frac

        subroutine frac2cart(frac,cart)
        real(dp), intent(in)  :: frac(3)
        real(dp), intent(out)  :: cart(3)

        cart =  matmul(ucell,frac)
        end subroutine frac2cart

      subroutine write_canonical_ucell(iunit,filename)
      use units, only: Ang
!
!     Writes out unit cell information in fdf-compatible format
!     with LatticeVectors in angstrom.
!
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: iunit

        integer  :: iu, ix, iv
        character(len=90)    :: fname

      if (present(iunit)) then
         iu = iunit
      else
         if (present(filename)) then
            fname = filename
         else
            fname = "OUT.UCELL"
         endif
         call io_assign( iu )
         open(iu, file=trim(fname), form='formatted',
     $        position='rewind', status='unknown')
      endif

      write(iu,"(a)") "LatticeConstant 1.0 Ang"
      write(iu,"(a)") "%block LatticeVectors"
      write(iu,'(3x,3f18.9)')
     .     ((ucell(ix,iv)/Ang,ix=1,3),iv=1,3)
      write(iu,"(a)") "%endblock LatticeVectors"

      if (.not. present(iunit)) then
         call io_close(iu)
      endif

      end subroutine write_canonical_ucell

      end module m_cell
!---------------------------------------------------
