! *** Module: nao2gto_transform ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief I/O module for Gaussian-based Hartree-Fock exchange
!!
!!  This module bridges the input file of SIESTA with the NAO2GTO routines,
!!  which calculate the Hartree-Fock exchange interaction using Gaussians.
!!
!! \note
!!      This file currently works with a version of LIBINT configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2017 Reviewed for inclusion in SIESTA [Xinming Qin]
!!      - 01.2018 Brought together from separate files [Yann Pouillon]
!!      - 01.2022 modified  l_max = 5 and Y(lmn)=(-1)^l for SIESTA v4.1.5 [Xinming Qin]
!!  Now the spher_harm.f is unchanged, otherwise positive coeff is required when l=odd.
!!        
! *****************************************************************************
module nao2gto_transform

  use precision, only: dp

  implicit none

  private

  public :: &
&   calc_c2s_matrix, &
&   cphi2sphi, &
&   orbtramat_type

    !> \brief Spherical harmonics and orbital transformation matrices
    !!
    !! This type stores spherical harmonics and the corresponding orbital
    !! transformation matrices.
    !!
    !! \par Literature
    !!      H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)
    type :: orbtramat_type
      real(dp), dimension(:,:), pointer :: c2s
    end type orbtramat_type

contains

! *****************************************************************************
!> \brief Computes the cartesian-to-spherical transformation matrix
!!
!! This routine builds the orbital transformation matrix for the
!! conversion from cartesian gaussians to real spherical harmonic
!! gaussians (c2s). 
!! This is based on Eq. (15) given in \cite Schlegel-95
!!
!! Let us define the Cartesian Gaussian function as
!!
!! \f{eqnarray*}{
!! g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} ) =
!!        N(\alpha, l_x, l_y, l_z ) x^{l_{x}} y^{l_{y}} z^{l_{z}} 
!!        e^{-\alpha r^{2}}
!! \f}
!!
!! where \f$ N \f$ is a normalization factor, and for a given angular momentum
!! \f$ l = l_{x} + l_{y} + l_{z} \f$. 
!! 
!! The electron repulsion integrals (ERIs) are calculated over this set
!! of Cartesian Gaussian functions within LIBINT.
!! But SIESTA uses a basis set of real spherical harmonics.
!!
!! In this context, we define a basis set of spherical harmonic Gaussians as
!!
!! \f{eqnarray*}{
!! \tilde{g}_{\rm Spher}^{\rm complex} (\alpha, l, m, n, \vec{r} ) =
!!        \tilde{N}(\alpha, n ) r^{n} Y_{lm}^{\rm complex} 
!!        e^{-\alpha r^{2}}
!! \f}
!!
!! where \f$ \tilde{N} \f$ is a normalization factor, and \f$ Y_{lm}^{\rm complex} \f$ 
!! is a complex spherical harmonic.
!!
!! The pure spherical harmonic Gaussians can be constructed from a linear
!! combination of the appropriate cartesian Gaussians
!!
!! \f{eqnarray*}{
!! r^{l-n} \tilde{g}_{\rm Spher}^{\rm complex} (\alpha, l, m, n, \vec{r} ) =
!!        \sum_{l_x + l_y + l_z = l}
!!        c(l, m, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} )
!! \f}
!! This corresponds with Eq.(3) of Ref. \cite Schlegel-95. 
!! According to the previous Equations, we might conclude that
!!
!! \f{eqnarray*}{
!!        \tilde{N}(\alpha, n ) r^{l} Y_{lm}^{\rm complex} 
!!        e^{-\alpha r^{2}} = 
!!        \sum_{l_x + l_y + l_z = l}
!!        c(l, m, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} )
!! \f}
!!
!! The transformation coefficients are given in Eq. (15) of 
!! Ref. \cite Schlegel-95.
!!
!! \f{eqnarray*}{
!! c(l, m, l_x, l_y, l_z ) = 
!!     \sqrt{\frac{(2l_x)! (2l_y)! (2l_z)! l! (l - \vert m \vert)! } 
!!                {(2l)! l_x! l_y! l_z! (l + \vert m \vert)! } }
!!     \frac{1}{2^l l!} \sum_{i=0}^{(l - \vert m \vert)/2} 
!!     \left( \begin{array}{c}
!!          l \\\
!!          i \end{array} \right)
!!     \left( \begin{array}{c}
!!          i \\\
!!          j \end{array} \right)
!!     \frac{(-1)^{i} (2l-2i)!}{(l - \vert m \vert -2i)!}
!!     \sum_{k=0}^{j} 
!!     \left( \begin{array}{c}
!!          j \\\
!!          k \end{array} \right)
!!     \left( \begin{array}{c}
!!          \vert m \vert \\\
!!          l_x-2k \end{array} \right)
!!          (-1)^{\pm(\vert m \vert - l_x + 2k)/2},
!! \f}
!! where \f$ j = ( l_x + l_y - \vert m \vert) / 2 \f$, and is an integer;
!! \f$ c(l, m, n, l_x, l_y, l_z ) =   0 \f$ if \f$ j \f$ is a half-integer.
!! The binomial coefficients \f$ \left( \begin{array}{c}
!!          j \\\
!!          k \end{array} \right) \f$ are zero for \f$ q < 0\f$ and \f$ q>p\f$.
!!  
!! In SIESTA, we use real spherical harmonics, so an extra transformation
!! has to be considered. 
!! 
!! The \f$ m = 0 \f$ spherical harmonics are already real.
!!
!! For \f$ m \ne 0 \f$ pairs, pure spherical harmonics can be combined into 
!! two real spherical harmonics functions,
!! \f$ (Y_{lm}^{\rm complex} + Y_{l,-m}^{\rm complex}) / \sqrt{2} \f$ and
!! \f$ (Y_{lm}^{\rm complex} - Y_{l,-m}^{\rm complex}) / \sqrt{-2} \f$ and.
!! Since \f$ m \f$ always appear in absolute value in the equations for the
!! coefficients, then
!!
!! \f{eqnarray*}{
!!        \tilde{N}(\alpha, n ) r^{l} Y_{lm}^{\rm real} 
!!        e^{-\alpha r^{2}} = 
!!        \sum_{l_x + l_y + l_z = l}
!!        \frac{2}{\sqrt{2}} c(l, m, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} )
!!        = \sum_{l_x + l_y + l_z = l}
!!        \sqrt{2} c(l, m, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} )
!!        = \sum_{l_x + l_y + l_z = l}
!!        c^{\rm real}(l, m, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} )
!! \f}
!! if \f$ m \ne 0 \f$.
!!
!! This subroutine computes the coefficients of the linear transformation
!! from cartesian Gaussian functions to real spherical harmonic Gaussian 
!! functions (\f$ c^{\rm real} \f$), including the factor of square root of 2.
!!
!! Examples:
!!
!! For the real spherical harmonic \f$ p_x \f$ orbital 
!! \f$ \Rightarrow l = 1, m = 1  \f$, and ignoring the normalization functions
!!
!! then 
!!
!! \f{eqnarray*}{
!!        r Y_{l=1,m=1}^{\rm real} 
!!        e^{-\alpha r^{2}} = 
!!        x e^{-\alpha r^{2}} = 
!!        \sum_{l_x + l_y + l_z = 1}
!!        c^{\rm real}(l=1, m=1, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} ) =
!!        c^{\rm real}(l=1, m=1, n, l_x=1, l_y=0, l_z=0 )
!!        x e^{-\alpha r^{2}} + 
!!        c^{\rm real}(l=1, m=1, n, l_x=0, l_y=1, l_z=0 )
!!        y e^{-\alpha r^{2}} + 
!!        c^{\rm real}(l=1, m=1, n, l_x=0, l_y=0, l_z=1 )
!!        z e^{-\alpha r^{2}} ,
!! \f}
!! and therefore
!! \f{eqnarray*}{c^{\rm real}(l=1, m=1, n, l_x=1, l_y=0, l_z=0 ) = 1 \f}
!! \f{eqnarray*}{c^{\rm real}(l=1, m=1, n, l_x=0, l_y=1, l_z=0 ) = 0 \f}
!! \f{eqnarray*}{c^{\rm real}(l=1, m=1, n, l_x=0, l_y=0, l_z=1 ) = 0 \f}
!!
!! For the real spherical harmonic \f$ d_{x^{2}-y^{2}} \f$ orbital 
!! \f$ \Rightarrow l = 2, m = 2  \f$, and ignoring the normalization functions
!!
!! then 
!!
!! \f{eqnarray*}{
!!        r^{2} Y_{l=2,m=2}^{\rm real} 
!!        e^{-\alpha r^{2}} = 
!!        (x^{2}-y^{2})  e^{-\alpha r^{2}} = 
!!        \sum_{l_x + l_y + l_z = 2}
!!        c^{\rm real}(l=2, m=2, n, l_x, l_y, l_z )
!!        g_{\rm Cart} (\alpha, l_x, l_y, l_z, \vec{r} ) =
!!        c^{\rm real}(l=2, m=2, n, l_x=2, l_y=0, l_z=0 )
!!        x^{2} e^{-\alpha r^{2}} + 
!!        c^{\rm real}(l=2, m=2, n, l_x=0, l_y=2, l_z=0 )
!!        y^{2} e^{-\alpha r^{2}} + 
!!        c^{\rm real}(l=2, m=2, n, l_x=0, l_y=0, l_z=2 )
!!        z^{2} e^{-\alpha r^{2}} +
!!        c^{\rm real}(l=2, m=2, n, l_x=1, l_y=1, l_z=0 )
!!        xy e^{-\alpha r^{2}} +
!!        c^{\rm real}(l=2, m=2, n, l_x=1, l_y=0, l_z=1 )
!!        xz e^{-\alpha r^{2}} +
!!        c^{\rm real}(l=2, m=2, n, l_x=0, l_y=1, l_z=1 )
!!        yz e^{-\alpha r^{2}} 
!! \f}
!! and therefore
!! \f{eqnarray*}{c^{\rm real}(l=2, m=2, n, l_x=2, l_y=0, l_z=0 ) = 
!!               - c^{\rm real}(l=2, m=2, n, l_x=0, l_y=2, l_z=0 ) \f}
!! \f{eqnarray*}{c^{\rm real}(l=2, m=2, n, l_x=0, l_y=0, l_z=2 ) = 0 \f}
!! \f{eqnarray*}{c^{\rm real}(l=2, m=2, n, l_x=1, l_y=1, l_z=0 ) = 0 \f}
!! \f{eqnarray*}{c^{\rm real}(l=2, m=2, n, l_x=1, l_y=0, l_z=1 ) = 0 \f}
!! \f{eqnarray*}{c^{\rm real}(l=2, m=2, n, l_x=0, l_y=1, l_z=1 ) = 0 \f}
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] l_max: Maximum angular momentum
!! \param[in] co: matrix that transforms the three exponents on the cartesian
!!                components of a Cartesian Gaussian function (lx, ly, lz)
!!                into a single index that identifies the Cartesian Gaussian
!!                function of a given angular momentum
!!                For a given l, co runs between 1 and (l + 1)*(l + 2)/2
!! \param[out] orbtramat: Coefficients of the linear transformation
!!                from cartesian Gaussian functions
!!                to real spherical harmonic Gaussian functions
! *****************************************************************************
  subroutine calc_c2s_matrix(l_max, co, orbtramat)

    use nao2gto_utils, only: fac     ! Factorial function

    implicit none

    ! Arguments
    integer, intent(in) :: l_max     ! Maximum angular momentum
    integer, intent(in) :: co(0:l_max,0:l_max,0:l_max)
                                     ! Matrix that transforms the three 
                                     !   exponents on the cartesian
                                     !   components of a Cartesian Gaussian 
                                     !   function (lx, ly, lz)
                                     !   into a single index that identifies 
                                     !   the Cartesian Gaussian
                                     !   function of a given angular momentum
                                     !   For a given l, co runs between 
                                     !      1 and (l + 1)*(l + 2)/2
    type(orbtramat_type), dimension(0:l_max), intent(out) :: orbtramat
                                     ! Coefficients of the linear transformation
                                     !    from cartesian Gaussian functions
                                     !    to real spherical harmonic Gaussian
                                     !    functions

    ! Local variables
    integer  :: expo, ic, k
    integer  :: lx             ! Exponent of the x-cartesian coordinate in 
                               !   a cartesian gaussian function
    integer  :: ly             ! Exponent of the y-cartesian coordinate in 
                               !   a cartesian gaussian function
    integer  :: lz             ! Exponent of the z-cartesian coordinate in 
                               !   a cartesian gaussian function
                               ! They must comply with
                               !     lx + ly + lz = l
    integer  :: l              ! Counter for loop on angular quantum numbers
    integer  :: m              ! Magnetic quantum number
    integer  :: ma             ! Absolute value of the magnetic quantum number
    integer  :: i              ! Counter for the first sum of Eq. (15) in the 
                               !   paper by Schlegel and Frisch.
                               !   Same notation as in the paper
    integer  :: j              ! Counter for the second sum of Eq. (15) in the 
                               !   paper by Schlegel and Frisch.
                               !   Same notation as in the paper
    integer  :: isp            ! Index that run between 1 and (2*l + 1) for
                               !   a given l
    real(dp) :: s              ! Last line of Eq. (15) of the 
                               !   paper by Schlegel and Frisch
                               !   (-1)**expo
    real(dp) :: s1             ! First sum appearing in Eq. (15) of the 
                               !   paper by Schlegel and Frisch
    real(dp) :: s2             ! Second sum appearing in Eq. (15) of the 
                               !   paper by Schlegel and Frisch

    ! -------------------------------------------------------------------------

!!   For debugging
!    write(6,'(a,i5)')                        &
! &    'calc_c2s_matrix: l_max = ',  l_max  
!!   End debugging

!   Loop over all possible angular momenta
    do l = 0, l_max 
      do lx = 0, l
        do ly = 0, l - lx

!         In a basis of Gaussian Cartessian functions, for a given angular
!         momentum l, l = lx + ly + lz
!         See right after Eq. (2) of the paper by Schlegel and Frisch
          lz = l - lx - ly

          ic = co(lx,ly,lz)
!         Loop over magnetic quantum number
          do m = -l, l
!           In principle, m runs between -l and +l.
!           We can perform a change so m runs from 1 to 2*l+1.
!           Useful to store information in matrices where the entrances depend
!           on m. 
!           This is done in the following line
            isp = l + m + 1
!           In equation (15) of the paper by Schlegel and Frisch,
!           only the absolute value of the magnetic quantum number appears
            ma = abs(m)

!           Define the upper limit of the second-summatory in 
!           Eq. (15) of the paper by Schlegel and Frisch
            j = lx + ly - ma

!           The coefficient is non zero only if j is an even integral number,
!           as said right below Eq. (15) of the paper by
!           Schlegel and Frisch
            if ( (j >= 0) .and. (modulo(j,2) == 0) ) then
!             The value of j defined right after Eq. (15) of the 
!             previous paper corresponds to the former one divided by 2
              j = j/2
              s1 = 0.0_dp

!             Perform the first of the sums in Eq. (15) of the reference paper
              do i = 0, (l-ma)/2
                s2 = 0.0_dp
!               Perform the second of the sums in Eq. (15) of the reference pape
!               The factor sqrt(2.0_dp) if m \ne 0 is due to the
!               transformations from complex to real spherical harmonics
!               (see the notes in the documentation in the header).
                do k = 0, j
                  if ( ((m < 0) .and. (modulo(abs(ma-lx),2) == 1)) .or. &
&                      ((m > 0) .and. (modulo(abs(ma-lx),2) == 0)) ) then
                    expo = (ma - lx + 2*k)/2
                    s = (-1.0_dp)**expo*sqrt(2.0_dp)
                  else if ( (m == 0) .and. (modulo(lx,2) == 0) ) then
                    expo = k - lx/2
                    s = (-1.0_dp)**expo
                  else
                    s = 0.0_dp
                  end if
                  s2 = s2 + binomial(j,k)*binomial(ma,lx-2*k)*s
                end do ! End loop to compute the second sum of Eq. (15)
                s1 = s1 + binomial(l,i)*binomial(i,j)*             &
&                    (-1.0_dp)**i*fac(2*l-2*i)/fac(l-ma-2*i)*s2
              end do  ! End loop to compute the first sum of Eq. (15)

!             Multiply here all the prefactors in front of the two sums of 
!             Eq. (15)
              orbtramat(l)%c2s(isp,ic) = &
&               sqrt((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l-ma))/ &
&                    (fac(lx)*fac(ly)*fac(lz)*fac(2*l)*fac(l+ma)))*s1/ &
&                    ((2.0_dp**l)*fac(l))*(-1)**ma
            else
              orbtramat(l)%c2s(isp,ic) = 0.0_dp
            end if
!!           For debugging
!            write(6,'(a,7i5,f12.5)')                                       &
! &            'calc_c2s_matrix: l, lx, ly, lz, m, isp, ic, orbtramat  = ', &
! &            l, lx, ly, lz, m, isp, ic, orbtramat(l)%c2s(isp,ic)  
!!           End debugging
          end do    ! End loop over magnetic quantum number
        end do      ! End loop on ly
      end do        ! End loop on lx
    end do          ! End loop over angular momenta

  end subroutine calc_c2s_matrix

! *****************************************************************************
!> \brief Computes cartesian to spherical coordinates transforms
!!
!! This routine converts Gaussian coefficients expressed in a cartesian
!! coordinates system to spherical coordinates.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] ncon_max: Maximum number of Cartesian Gaussian type orbitals
!! \param[in] norb_cphi: Total number of Cartesian Gaussian type orbitals 
!!                       for a given atom
!! \param[in] norb_sphi: Total number of Spherical Harmonics type orbitals
!!                       for a given atom
!! \param[in] l_max: Maximum value of the angular momentum
!! \param[in] norb_nl: Total number of shells of NAO in the basis set
!! \param[in] orbnl_l: Orbital angular momentum of a given shell
!! \param[in] nco: Number of Cartesian Gaussian type orbitals for a given shell
!! \param[in] nso: Number of Spherical Harmonics for a given shell
!! \param[in] index_cphi: Index of the first Cartesian Gaussian type orbital
!!           with the angular momentum l in the list of orbitals of a given atom
!! \param[in] index_sphi: Index of the first Spherical Gaussian type orbital
!!           with the angular momentum l in the list of orbitals of a given atom
!! \param[in] cphi: ...
!! \param[in,out] sphi: ...
!! \param[in] orbtramat: Coefficients of the linear transformation
!!           from cartesian Gaussian functions
!!           to real spherical harmonic Gaussian functions
!!
!! \note The sphi variable must be declared inout because of LAPACK's DGEMM.
! *****************************************************************************
  subroutine cphi2sphi(ncon_max, norb_cphi, norb_sphi, l_max, norb_nl, &
&                orbnl_l, nco, nso, index_cphi, index_sphi, cphi, sphi, &
&                orbtramat)

    implicit none

    ! Arguments
    integer, intent(in)     :: norb_nl     ! Total number of shells of NAO
                                           ! Example: for Si DZP, 
                                           ! norb_nl = 5 
                                           ! Two for the 3 s-shell
                                           ! Two for the 3 p-shell
                                           ! One for the 3 d-shell
    integer, intent(in)     :: ncon_max    ! Maximum number of Gaussian type
                                           !   orbitals for a given atom
    integer, intent(in)     :: norb_cphi   ! Total number of Cartesian Gaussian
                                           !   type orbitals for a given atom
    integer, intent(in)     :: norb_sphi   ! Total number of Spherical Harmonic
                                           !   type orbitals for a given atom
    integer, intent(in)     :: l_max       ! Maximum angular momentum shell
    integer, intent(in)     :: orbnl_l(norb_nl) 
                                           ! Orbital angular momentum of a given
                                           ! shell
    integer, intent(in)     :: nco(0:l_max)
                                           ! Number of Cartesian Gausians type
                                           !   orbitals for a given shell
                                           !   1 for the s
                                           !   3 for the p
                                           !   6 for the d
                                           !  10 for the f
    integer, intent(in)     :: nso(0:l_max)
                                           ! Number of Spherical Harmonics 
                                           !   for a given shell
                                           !   1 for the s
                                           !   3 for the p
                                           !   5 for the d
                                           !   7 for the f
    integer, intent(in)     :: index_cphi(norb_nl), index_sphi(norb_nl)
    real(dp), intent(in)    :: cphi(ncon_max,norb_cphi)
    real(dp), intent(inout) :: sphi(ncon_max,norb_sphi)
    type(orbtramat_type), dimension(0:l_max), intent(in) :: orbtramat

    ! Local variables
    integer :: i         ! Counter for loop on the shells
    integer :: l         ! Orbital angular momentum of a given shell
    integer :: ncgf      ! Number of Cartesian Gaussian function for a 
                         !   particular shell
    integer :: nsgf      ! Number of Spherical Harmonics for a given shell
    integer :: first_cgf ! Index of the first Cartessian Gaussian type orbital
                         !   with a given value of l in the list of CG type 
                         !   orbitals of a given atom
    integer :: first_sgf ! Index of the first Spherical Harmonic type orbital
                         !   with a given value of l in the list of spherical
                         !   harmonic type orbitals of a given atom
    integer :: icounter  ! Counter for loops for debugging
    integer :: jcounter  ! Counter for loops for debugging

! -------------------------------------------------------------------------
!!     For debugging
!      write(6,'(a,i5)')       &
! &      'cphi2sphi: l_max     = ', l_max
!      write(6,'(a,i5)')       &
! &      'cphi2sphi: ncon_max  = ', ncon_max
!      write(6,'(a,i5)')       &
! &      'cphi2sphi: norb_cphi = ', norb_cphi
!      write(6,'(a,i5)')       &
! &      'cphi2sphi: norb_sphi = ', norb_sphi
!       do i = 1, norb_nl
!         write(6,'(a,3i5)')                                   &
! &         'cphi2sphi: iorb_nl, index_cphi, index_sphi = ',   &
! &         i, index_cphi(i), index_sphi(i)
!       enddo
!!     End debugging

!   Loop on all the shells of NAO in the basis of given atom
    do i = 1, norb_nl

!     Identify the angular momentum shell
      l = orbnl_l(i)
!     Identify the number of Cartesian Gaussian functions for shell i
      ncgf = nco(l)
!     Identify the number of Spherical hamonics for shell i
      nsgf = nso(l)
!     Indentify the index of the first Cartesian Gaussian type orbital
!     with the angular momentum l in the list of orbitals of a given atom
      first_cgf = index_cphi(i)
!     Indentify the index of the first Spherical Harmonic type orbital
!     with the angular momentum l in the list of orbitals of a given atom
      first_sgf = index_sphi(i)


!!     For debugging
!      write(6,'(a,7i5)') &
! &      'cphi2sphi: i, l, ncgf, nsgf, first_cgf, first_sgf = ', &
! &                  i, l, ncgf, nsgf, first_cgf, first_sgf 
!      do icounter = 1, ncon_max
!        write(6,'(a,2i5,f12.5)')                      &
! &        'cphi2sphi: icounter, first_cgf, cphi  = ', &
! &                    icounter, first_cgf, cphi(icounter,first_cgf) 
!      enddo 
!      do icounter = 1, nsgf
!        do jcounter = 1, ncgf
!          write(6,'(a,2i5,f12.5)')                           &
! &          'cphi2sphi: icounter, jcounter, orbtramat  = ',  &
! &                      icounter, jcounter, orbtramat(l)%c2s(icounter,jcounter)
!        enddo 
!      enddo 
!!     End debugging

!     DGEMM  performs one of the matrix-matrix operations
!     C := alpha*op( A )*op( B ) + beta*C,
!     where  op( X ) is one of
!     op( X ) = X   or   op( X ) = X**T,
!     alpha and beta are scalars, and A, B and C are matrices, with op( A )
!     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!     Here: alpha = 1.0_dp
!           beta  = 0.0_dp
!           A     = cphi,  that is taken as is (not transposed)
!           B     = orbtramat(l)%c2s(1,1),  that is taken as transposed
!           C     = sphi

!     write(6,*) "nsgf, ncgf", nsgf, ncgf
!     write(6,*) "cphi: ", cphi(1:ncon_max,1:ncgf)


      call dgemm("N", "T", ncon_max, nsgf, ncgf, 1.0_dp, &
&       cphi(1,first_cgf), ncon_max, orbtramat(l)%c2s(1,1), nsgf, 0.0_dp, &
&       sphi(1,first_sgf), ncon_max)

!    write(6,*) "sphi: ",sphi(1:ncon_max,1:nsgf)

    enddo

  end subroutine cphi2sphi

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Returns the number of binomial combinations of two integers
!!
!! This function computes the binomial combinations of two integers k and n
!! as (n!/(k!(n-k)!)).
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] n: sample size
!! \param[in] k: number of sample elements selected
!! \return number of binomial combinations of k elements among n
! *****************************************************************************
  function binomial(n,k) result(n_over_k)

    use nao2gto_utils, only: fac

    implicit none

    ! Arguments
    integer, intent(in) :: n, k

    ! Local variables
    real(dp) :: n_over_k

    ! -------------------------------------------------------------------------

    if ( (k >= 0) .and. (k <= n) ) then
      n_over_k = fac(n)/(fac(n-k)*fac(k))
    else
      n_over_k = 0.0_dp
    endif

  end function binomial

end module nao2gto_transform
