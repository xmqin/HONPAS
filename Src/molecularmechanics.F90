module molecularmechanics
!
! Add additional interactions through pair-potentials
! Implementation: Julian Gale (Curtin, AU)
! Modified by Alberto Garcia  (stress sign fix, plots of V(r))
! Implemented Grimme's cutoff function (A. Garcia, April 2008)
!

 use precision, only: dp

 implicit none

    ! The original implementation by Julian had the wrong sign
    ! for the stress tensor.
    ! Define sign_change as +1 to recover that behavior

 integer, parameter    :: sign_change = -1

 real(dp), parameter   :: d_grimme_default = 20.0_dp    ! 2006 Grimme paper

 private

 integer, parameter    :: maxMMpot = 10
 integer, save         :: nMMpot = 0
 integer, save         :: nMMpotptr(2,maxMMpot)
 integer, save         :: nMMpottype(maxMMpot)
 logical, save         :: PotentialsPresent = .false.
 real(dp), save        :: MMcutoff
 real(dp), save        :: s6_grimme
 real(dp), save        :: d_grimme
 real(dp), save        :: MMpotpar(6,maxMMpot)

 public :: inittwobody, twobody

 CONTAINS

 subroutine inittwobody()

   use fdf
   use units,   only : eV, Ang
   use parsing, only : parse
   use sys,     only : die
   use parallel,only : Node
#  ifdef MPI
   use mpi_siesta
#  endif

   real(dp), parameter   :: s6_grimme_default = 1.66_dp    ! Fit by Roberto Peverati for DZ basis sets

      interface
         function leqi(s1,s2)
         logical leqi
         character(len=*), intent(in)   :: s1, s2
         end function leqi
      end interface

   character(len=130) :: line
   character(len=80)  :: names
   character(len=80)  :: scale
   integer            :: il
   integer            :: integs(4)
   integer            :: iu
   integer            :: lastc
   integer            :: lc(0:3)
   integer            :: maxlin
   integer            :: ni
   integer            :: nn
   integer            :: nr
   integer            :: nv
#  ifdef MPI
   integer            :: MPIerror
#  endif
   real(dp)           :: Dscale
   real(dp)           :: Escale
   real(dp)           :: reals(4)
   real(dp)           :: values(4)

 ! Allocation of arrays formerly done here moved
 ! to top of module, as they are really static for now

   MMpotpar(1:6,1:maxMMpot) = 0.0_dp
 !
 ! Get potential cutoff
 !
#ifdef MPI
  if (Node.eq.0) then
    MMcutoff = fdf_physical('MM.Cutoff',30.0d0,'Bohr')
  endif

  call MPI_Bcast(MMcutoff,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
#else
  MMcutoff = fdf_physical('MM.Cutoff',30.0d0,'Bohr')
#endif
!
! Set MM units of energy for potential parameters
!
  if (Node.eq.0) then
    scale = fdf_string( 'MM.UnitsEnergy','eV' )
    if (scale.eq.'eV') then
      Escale = eV
    else
      Escale = 1.0_dp
    endif
    scale = fdf_string( 'MM.UnitsDistance','Ang' )
    if (scale.eq.'Ang') then
      Dscale = Ang
    else
      Dscale = 1.0_dp
    endif
  endif
!
! Read in data from block
!
  nMMpot = 0
#ifdef MPI
  if (Node.eq.0) then
    PotentialsPresent = fdf_block('MM.Potentials',iu)
  endif
  call MPI_Bcast(PotentialsPresent,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#else
  PotentialsPresent = fdf_block('MM.Potentials',iu)
#endif
  if (PotentialsPresent) then
#ifdef MPI
    if (Node.eq.0) then
#endif
    write(6,"(a)") "Reading two-body potentials"
    maxlin = maxMMpot
    do il = 1,maxlin
!
! Read and parse data line
!
      read(iu,'(a)',end=50) line
      if (index(line,'%endblock').ne.0) exit
      lastc = index(line,'#') - 1
      if (lastc .le. 0) lastc = len(line)
      call parse( line(1:lastc), nn, lc, names, nv, values, ni, integs, nr, reals )
!
      if (nn.gt.0) then
        if (index(names(1:lc(1)),'c6').gt.0.or.index(names(1:lc(1)),'C6').gt.0) then
          nMMpot = nMMpot + 1
          if (nMMpot.gt.maxMMpot) then
            call die('MM: Too many MM potentials - increase maxMMpot!')
          endif
          nMMpottype(nMMpot) = 1
          if (ni.lt.2) then
            call die('MM: Species numbers missing in potential input!')
          endif
          nMMpotptr(1,nMMpot) = integs(1)
          nMMpotptr(2,nMMpot) = integs(2)
          write(6,"(a,i3,a,i3)") "C6 - two-body potential between ", integs(1), " and ", integs(2)
          if (nr.ge.2) then
! C6 : Parameter one is C6 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
! C6 : Parameter two is damping exponent
            MMpotpar(2,nMMpot) = reals(2)/Dscale
          elseif (nr.eq.1) then
! C6 : Parameter one is C6 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
            MMpotpar(2,nMMpot) = 0.0_dp
          endif
        elseif (index(names(1:lc(1)),'c8').gt.0.or.index(names(1:lc(1)),'C8').gt.0) then
          nMMpot = nMMpot + 1
          if (nMMpot.gt.maxMMpot) then
            call die('MM: Too many MM potentials - increase maxMMpot!')
          endif
          nMMpottype(nMMpot) = 2
          if (ni.lt.2) then
            call die('MM: Species numbers missing in potential input!')
          endif
          nMMpotptr(1,nMMpot) = integs(1)
          nMMpotptr(2,nMMpot) = integs(2)
          write(6,"(a,i3,a,i3)") "C8 - two-body potential between ", integs(1), " and ", integs(2)
          if (nr.ge.2) then
! C8 : Parameter one is C8 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
! C8 : Parameter two is damping exponent
            MMpotpar(2,nMMpot) = reals(2)/Dscale
          elseif (nr.eq.1) then
! C8 : Parameter one is C8 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
            MMpotpar(2,nMMpot) = 0.0_dp
          endif
        elseif (index(names(1:lc(1)),'c10').gt.0.or.index(names(1:lc(1)),'C10').gt.0) then
          nMMpot = nMMpot + 1
          if (nMMpot.gt.maxMMpot) then
            call die('MM: Too many MM potentials - increase maxMMpot!')
          endif
          nMMpottype(nMMpot) = 3
          if (ni.lt.2) then
            call die('MM: Species numbers missing in potential input!')
          endif
          nMMpotptr(1,nMMpot) = integs(1)
          nMMpotptr(2,nMMpot) = integs(2)
          write(6,"(a,i3,a,i3)") "C10 - two-body potential between ", integs(1), " and ", integs(2)
          if (nr.ge.2) then
! C10 : Parameter one is C10 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
! C10 : Parameter two is damping exponent
            MMpotpar(2,nMMpot) = reals(2)/Dscale
          elseif (nr.eq.1) then
! C10 : Parameter one is C10 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)
            MMpotpar(2,nMMpot) = 0.0_dp
          endif
        elseif (index(names(1:lc(1)),'harm').gt.0.or.index(names(1:lc(1)),'HARM').gt.0) then
          nMMpot = nMMpot + 1
          if (nMMpot.gt.maxMMpot) then
            call die('MM: Too many MM potentials - increase maxMMpot!')
          endif
          nMMpottype(nMMpot) = 4
          if (ni.lt.2) then
            call die('MM: Species numbers missing in potential input!')
          endif
          nMMpotptr(1,nMMpot) = integs(1)
          nMMpotptr(2,nMMpot) = integs(2)
          write(6,"(a,i3,a,i3)") "Harmonic two-body potential between ", integs(1), " and ", integs(2)
          if (nr.ge.2) then
! Harm : Parameter one is force constant
            MMpotpar(1,nMMpot) = reals(1)*Escale/(Dscale**2)
! Harm : Parameter two is r0
            MMpotpar(2,nMMpot) = reals(2)*Dscale
          elseif (nr.eq.1) then
! Harm : Parameter one is force constant
            MMpotpar(1,nMMpot) = reals(1)*Escale/(Dscale**2)
            MMpotpar(2,nMMpot) = 0.0_dp
          endif
        elseif (leqi(names(1:lc(1)),'Grimme')) then
          nMMpot = nMMpot + 1
          if (nMMpot.gt.maxMMpot) then
            call die('MM: Too many MM potentials - increase maxMMpot!')
          endif
          nMMpottype(nMMpot) = 5
          if (ni.lt.2) then
            call die('MM: Species numbers missing in potential input!')
          endif
          nMMpotptr(1,nMMpot) = integs(1)
          nMMpotptr(2,nMMpot) = integs(2)
          write(6,"(a,i3,a,i3)") "Grimme two-body potential between ", integs(1), " and ", integs(2)
          if (nr.eq.2) then

! C6 : Parameter one is C6 coefficient
            MMpotpar(1,nMMpot) = reals(1)*Escale*(Dscale**6)

! C6 : Parameter two is the sum of the van-der-Waals radii
!      Note 1: This must be already appropriately corrected (i.e., factor of 1.1 ...)
!      Note 2: This is a real length, as opposed to the damping parameters for the Tang-Toenes
!              potentials, so note the correct application of the scale factor.

            MMpotpar(2,nMMpot) = reals(2) * Dscale

          else
            call die('MM: Need both C6 and R0 values in Grimme line!')
          endif
        endif
      endif
    enddo
50  continue

    if (any(nMMpottype .eq. 5)) then
       s6_grimme = fdf_double("MM.Grimme.S6",s6_grimme_default)
       d_grimme = fdf_double("MM.Grimme.D",d_grimme_default)
    endif

#ifdef MPI
    endif      ! Node 0
    call MPI_Bcast(nMMpot,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nMMpottype,nMMpot,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nMMpotptr(1,1),2*nMMpot,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(MMpotpar(1,1),6*nMMpot,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(s6_grimme,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(d_grimme,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
#endif
  endif  

if (Node .eq. 0)  call plot_functions()

end subroutine inittwobody

subroutine twobody(na,xa,isa,cell,emm,ifa,fa,istr,stress)

  use parallel,         only : Node
  use units,            only : kbar
  use alloc,            only : re_alloc, de_alloc

  integer,  intent(in)    :: na
  integer,  intent(in)    :: isa(*)
  integer,  intent(in)    :: ifa   ! compute forces if > 0
  integer,  intent(in)    :: istr  ! compute stress if > 0
  real(dp), intent(in)    :: cell(3,3)
  real(dp), intent(in)    :: xa(3,*)
  real(dp), intent(out)   :: emm
  real(dp), intent(inout) :: fa(3,*)
  real(dp), intent(inout) :: stress(3,3)

  integer                 :: i
  integer                 :: idir
  integer                 :: ii
  integer                 :: imid
  integer                 :: ix
  integer                 :: iy
  integer                 :: iz
  integer                 :: j
  integer                 :: jdir
  integer                 :: jj
  integer                 :: jmid
  integer                 :: kdir
  integer                 :: kk
  integer                 :: kmid
  integer                 :: np
  integer                 :: nvec
  integer                 :: nvecj0
  integer                 :: nveck0
  logical                 :: lallfound1
  logical                 :: lallfound2
  logical                 :: lallfound3
  logical                 :: lanyvalidpot
  logical, pointer        :: lvalidpot(:)
  real(dp)                :: MMcutoff2
  real(dp)                :: arg
  real(dp)                :: a
  real(dp)                :: b
  real(dp)                :: br6
  real(dp)                :: br8
  real(dp)                :: br10
  real(dp)                :: c
  real(dp)                :: alpha
  real(dp)                :: beta
  real(dp)                :: df6
  real(dp)                :: df8
  real(dp)                :: df10
  real(dp)                :: gamma
  real(dp)                :: earg
  real(dp)                :: ebr6
  real(dp)                :: ebr8
  real(dp)                :: ebr10
  real(dp)                :: etrm
  real(dp)                :: etrm1
  real(dp)                :: f6
  real(dp)                :: f8
  real(dp)                :: f10
  real(dp)                :: fg
  real(dp)                :: fgprime
  real(dp)                :: ftrm
  real(dp)                :: factor
  real(dp)                :: proj1
  real(dp)                :: proj2
  real(dp)                :: proj3
  real(dp)                :: r
  real(dp)                :: R0
  real(dp)                :: r2
  real(dp)                :: r2i
  real(dp)                :: r2j
  real(dp)                :: r2k
  real(dp)                :: rcx1
  real(dp)                :: rcy1
  real(dp)                :: rcz1
  real(dp)                :: rcx2
  real(dp)                :: rcy2
  real(dp)                :: rcz2
  real(dp)                :: rcx3
  real(dp)                :: rcy3
  real(dp)                :: rcz3
  real(dp)                :: recipa
  real(dp)                :: recipb
  real(dp)                :: recipc
  real(dp)                :: rnorm
  real(dp)                :: rx
  real(dp)                :: ry
  real(dp)                :: rz
  real(dp)                :: rxi
  real(dp)                :: ryi
  real(dp)                :: rzi
  real(dp)                :: rxj
  real(dp)                :: ryj
  real(dp)                :: rzj
  real(dp)                :: rvol
  real(dp)                :: vol
  real(dp), external      :: volcel
  real(dp)                :: x
  real(dp)                :: y
  real(dp)                :: z

  real(dp)                :: mm_stress(3,3)
  integer                 :: jx
!
! Start timer
!
  call timer('MolMec', 1 )
!
! Allocate workspace arrays
!
  nullify(lvalidpot)
  call re_alloc(lvalidpot,1,nMMpot,name="lvalidpot",routine="twobody")
!
! Initialise energy and mm_stress
!
  emm = 0.0_dp
  mm_stress(1:3,1:3) = 0.0_dp
!
! Find number of cell images required
!
  MMcutoff2 = MMcutoff**2
  call uncell(cell,a,b,c,alpha,beta,gamma,1.0_dp)
  recipa = 1.0_dp/a
  recipb = 1.0_dp/b
  recipc = 1.0_dp/c
!
! Find volume if required
!
  if (istr.ne.0) then
    vol = volcel(cell)
    rvol = 1.0_dp/vol
  endif
!
! Loop over first atom
!
  do i = 1,na
!                 
! Loop over second atom 
!               
    do j = 1,i
      if (i.eq.j) then
        factor = 0.5_dp
      else
        factor = 1.0_dp 
      endif
!
! Find valid potentials
!
      lanyvalidpot = .false.
      do np = 1,nMMpot
        if (nMMpotptr(1,np).eq.isa(i).and.nMMpotptr(2,np).eq.isa(j)) then
          lanyvalidpot = .true.
          lvalidpot(np) = .true.
        elseif (nMMpotptr(1,np).eq.isa(j).and.nMMpotptr(2,np).eq.isa(i)) then
          lanyvalidpot = .true.
          lvalidpot(np) = .true.
        else
          lvalidpot(np) = .false.
        endif
      enddo
      if (lanyvalidpot) then
!
! Find image of j nearest to i 
!
        x = xa(1,j) - xa(1,i)
        y = xa(2,j) - xa(2,i)
        z = xa(3,j) - xa(3,i)
!     
! Find projection of cell vector 3 on to i - j vector
!     
        rnorm = x*x + y*y + z*z
        if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
        proj3 = rnorm*recipc*(x*cell(1,3) + y*cell(2,3) + z*cell(3,3))
        kmid = nint(proj3)
        x = x - kmid*cell(1,3)
        y = y - kmid*cell(2,3) 
        z = z - kmid*cell(3,3)
!
! Find projection of cell vector 2 on to i - j vector
!
        rnorm = x*x + y*y + z*z
        if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
        proj2 = rnorm*recipb*(x*cell(1,2) + y*cell(2,2) + z*cell(3,2))
        jmid = nint(proj2) 
        x = x - jmid*cell(1,2)
        y = y - jmid*cell(2,2)
        z = z - jmid*cell(3,2)
!     
! Find projection of cell vector 1 on to i - j vector
!     
        rnorm = x*x + y*y + z*z 
        if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
        proj1 = rnorm*recipa*(x*cell(1,1) + y*cell(2,1) + z*cell(3,1))
        imid = nint(proj1)
        x = x - imid*cell(1,1)
        y = y - imid*cell(2,1)
        z = z - imid*cell(3,1)
!
! Initialise counter for number of valid vectors
!
        nvec = 0
!
! Outer loop over first cell vector direction 
!
        do idir = 1,-1,-2
!     
! Reinitialise distance squared 
!     
          r2i = 10000.0_dp*MMcutoff2
!
! Loop over first cell vector
!
          lallfound1 = .false.
          if (idir.eq.1) then
            ii = 0
          else
            ii = - 1
          endif
!
! Set initial coordinate vector
!
          rxi = x + dble(ii)*cell(1,1)
          ryi = y + dble(ii)*cell(2,1)
          rzi = z + dble(ii)*cell(3,1)
!
! Set increment vector
!
          rcx1 = dble(idir)*cell(1,1)
          rcy1 = dble(idir)*cell(2,1)
          rcz1 = dble(idir)*cell(3,1)
!
          do while (.not.lallfound1)
!
! Save number of vectors before search over second direction
!
            nvecj0 = nvec
!
! Outer loop over second cell vector direction 
!
            do jdir = 1,-1,-2
!     
! Reinitialise saved distance squared
!     
              r2j = 10000.0_dp*MMcutoff2
!
! Loop over second cell vector
!
              lallfound2 = .false.
              if (jdir.eq.1) then
                jj = 0
              else
                jj = - 1
              endif
!
! Set initial coordinate vector
!
              rxj = rxi + dble(jj)*cell(1,2)
              ryj = ryi + dble(jj)*cell(2,2)
              rzj = rzi + dble(jj)*cell(3,2)
!
! Set increment vector
!
              rcx2 = dble(jdir)*cell(1,2)
              rcy2 = dble(jdir)*cell(2,2)
              rcz2 = dble(jdir)*cell(3,2)
!
              do while (.not.lallfound2)
!
! Save number of vectors before search over third direction
!
                nveck0 = nvec
!
! Outer loop over third cell vector direction 
!
                do kdir = 1,-1,-2
!
! Reinitialise saved distance squared
!           
                  r2k = 10000.0_dp*MMcutoff2
!
! Loop over third!cell vector
!
                  lallfound3 = .false.
                  if (kdir.eq.1) then
                    kk = 0
                  else
                    kk = - 1
                  endif
!
! Set initial coordinate vector
!
                  rx = rxj + dble(kk)*cell(1,3)
                  ry = ryj + dble(kk)*cell(2,3)
                  rz = rzj + dble(kk)*cell(3,3)
!
! Set increment vector
!
                  rcx3 = dble(kdir)*cell(1,3)
                  rcy3 = dble(kdir)*cell(2,3)
                  rcz3 = dble(kdir)*cell(3,3)
!
                  do while (.not.lallfound3)
!
! Calculate square of distance
!
                    r2 = rx*rx + ry*ry + rz*rz
!
! Check distance squared against cutoff squared
!
                    if (r2.le.MMcutoff2) then
                      if (r2.gt.1.0d-10) then
!
! Valid distance, so increment counter
!
                        nvec = nvec + 1
!
! Evaluate potentials for this valid distance
!
                        do np = 1,nMMpot
                          if (lvalidpot(np)) then
                            if (nMMpottype(np).eq.1) then
                              if (MMpotpar(2,np).eq.0.0_dp) then
                                etrm = MMpotpar(1,np)/(r2**3)
                                ftrm = 6.0_dp*factor*etrm/r2
                                etrm = - etrm
                              else
                                br6 = MMpotpar(2,np)*sqrt(r2)
                                ebr6 = exp(-br6)
                                f6 = 1.0_dp + br6*(1.0_dp + 0.5_dp*br6*(1.0_dp + (br6/3.0_dp)*( &
                                     1.0_dp + 0.25_dp*br6*(1.0_dp + 0.2_dp*br6*(1.0_dp + (br6/6.0_dp))))))
                                f6 = 1.0_dp - f6*ebr6
                                etrm1 = MMpotpar(1,np)/(r2**3)
                                df6 = ebr6*(br6*(br6**6))/720.0_dp
                                etrm = - etrm1*f6
                                ftrm = factor*(6.0_dp*etrm1*f6 - etrm1*df6)/r2
                              endif
                            elseif (nMMpottype(np).eq.2) then
                              if (MMpotpar(2,np).eq.0.0_dp) then
                                etrm = MMpotpar(1,np)/(r2**4)
                                ftrm = 8.0_dp*factor*etrm/r2
                                etrm = - etrm
                              else
                                br8 = MMpotpar(2,np)*sqrt(r2)
                                ebr8 = exp(-br8)
                                f8 = 1.0_dp + br8*(1.0_dp + 0.5_dp*br8*(1.0_dp + (br8/3.0_dp)*( &
                                     1.0_dp + 0.25_dp*br8*(1.0_dp + 0.2_dp*br8*(1.0_dp + &
                                     (br8/6.0_dp)*(1.0_dp+(br8/7.0_dp)*(1.0_dp + 0.125_dp*br8)))))))
                                f8 = 1.0_dp - f8*ebr8
                                etrm1 = MMpotpar(1,np)/(r2**4)
                                df8 = ebr8*(br8*(br8**8))/40320.0_dp
                                etrm = - etrm1*f8
                                ftrm = factor*(8.0_dp*etrm1*f8 - etrm1*df8)/r2
                              endif
                            elseif (nMMpottype(np).eq.3) then
                              if (MMpotpar(2,np).eq.0.0_dp) then
                                etrm = MMpotpar(1,np)/(r2**5)
                                ftrm = 10.0_dp*factor*etrm/r2
                                etrm = - etrm
                              else
                                br10 = MMpotpar(2,np)*sqrt(r2)
                                ebr10 = exp(-br10)
                                f10 = 1.0_dp + br10*(1.0_dp + 0.5_dp*br10*(1.0_dp + &
                                      (br10/3.0_dp)*(1.0_dp + 0.25_dp*br10*(1.0_dp + 0.2_dp*br10*( &
                                      1.0_dp + (br10/6.0_dp)*(1.0_dp + (br10/7.0_dp)*(1.0_dp + &
                                      0.125_dp*br10*(1.0_dp + (br10/9.0_dp)*(1.0_dp + 0.1_dp*br10)))))))))
                                f10 = 1.0_dp - f10*ebr10
                                etrm1 = MMpotpar(1,np)/(r2**5)
                                df6 = ebr10*(br10*(br10**10))/3628800.0_dp
                                etrm = - etrm1*f10
                                ftrm = factor*(10.0_dp*etrm1*f10 - etrm1*df10)/r2
                              endif
                            elseif (nMMpottype(np).eq.4) then
                              r = sqrt(r2)
                              etrm = MMpotpar(1,np)*(r - MMpotpar(2,np))
                              ftrm = factor*etrm/r
                              etrm = 0.5_dp*etrm*(r - MMpotpar(2,np))
                            elseif (nMMpottype(np).eq.5) then    ! Grimme
                              r = sqrt(r2)
                              R0 = MMpotpar(2,np)
                              arg = - d_grimme*(r/R0 - 1.0_dp)
                              earg = exp(arg)
                              fg = 1.0_dp / ( 1.0_dp + earg )
                              etrm1 = s6_grimme * MMpotpar(1,np)/(r2**3)
                              fgprime = (d_grimme/R0) * earg / ( 1 + earg )**2
                              etrm = - etrm1*fg
                              ftrm = factor*(6.0_dp*etrm1*fg - etrm1*r*fgprime)/r2
                            endif
                            emm = emm + factor*etrm
                            if (ifa.ne.0) then
!
! Gradients
!
                              fa(1,i) = fa(1,i) + ftrm*rx
                              fa(2,i) = fa(2,i) + ftrm*ry
                              fa(3,i) = fa(3,i) + ftrm*rz
                              fa(1,j) = fa(1,j) - ftrm*rx
                              fa(2,j) = fa(2,j) - ftrm*ry
                              fa(3,j) = fa(3,j) - ftrm*rz
                            endif
                            if (istr.ne.0) then
!
! Stress
!
                              ftrm = ftrm*rvol
                              ! *** Sign_change should be -1 for the correct stress
                              ! *** See parameter definition above
                              ! ***
                   mm_stress(1,1) = mm_stress(1,1) - sign_change*ftrm*rx*rx
                   mm_stress(2,1) = mm_stress(2,1) - sign_change*ftrm*ry*rx
                   mm_stress(3,1) = mm_stress(3,1) - sign_change*ftrm*rz*rx
                   mm_stress(1,2) = mm_stress(1,2) - sign_change*ftrm*rx*ry
                   mm_stress(2,2) = mm_stress(2,2) - sign_change*ftrm*ry*ry
                   mm_stress(3,2) = mm_stress(3,2) - sign_change*ftrm*rz*ry
                   mm_stress(1,3) = mm_stress(1,3) - sign_change*ftrm*rx*rz
                   mm_stress(2,3) = mm_stress(2,3) - sign_change*ftrm*ry*rz
                   mm_stress(3,3) = mm_stress(3,3) - sign_change*ftrm*rz*rz
                            endif
                          endif
                        enddo
                      endif
                    endif
!
! Increment by third vector
!
                    kk = kk + kdir
                    rx = rx + rcx3
                    ry = ry + rcy3
                    rz = rz + rcz3
!                     
! Check to see if this direction is complete
!                     
                    lallfound3 = (r2.gt.r2k.and.r2.gt.MMcutoff2)
                    r2k = r2
                  enddo
                enddo
!
! Increment by second vector
!
                jj = jj + jdir
                rxj = rxj + rcx2
                ryj = ryj + rcy2
                rzj = rzj + rcz2
!
! Check to see if this direction is complete
!
                lallfound2 = (r2.gt.r2j.and.r2.gt.MMcutoff2.and.nvec.eq.nveck0)
                r2j = r2
              enddo
            enddo
!
! Increment by first vector
!
            ii = ii + idir
            rxi = rxi + rcx1
            ryi = ryi + rcy1
            rzi = rzi + rcz1
!
! Check to see if this direction is complete
!
            lallfound1 = (r2.gt.r2i.and.r2.gt.MMcutoff2.and.nvec.eq.nvecj0)
            r2i = r2
          enddo
        enddo
!
! Endif over valid potentials
!
      endif
!
! End loops over atoms
!
    enddo
  enddo

!
! Print and add MM contribution to stress
!
  if (istr.ne.0) then

     if (Node .eq. 0 .and. PotentialsPresent)  then
        write(6,'(/,a,6f12.2)')  'MM-Stress (kbar):',   &
                (mm_stress(jx,jx)/kbar,jx=1,3),       &
                 mm_stress(1,2)/kbar,                 &
                 mm_stress(2,3)/kbar,                 &
                 mm_stress(1,3)/kbar                  
     endif

     stress = stress + mm_stress
  endif

!
! Free workspace arrays
!
  call de_alloc(lvalidpot,name="lvalidpot")
!
! Stop timer
!
  call timer('MolMec', 2 )

end subroutine twobody


subroutine plot_functions()
!
! Writes out V(r) info to files of the form MMpot.NN
! Units: energy: eV, distance: Ang
!
  use units,   only : eV, Ang

integer :: np
character(len=20) :: fname
integer, parameter   :: npts = 1000

real(dp) :: rmin, rmax, delta, range
real(dp) :: etrm, ftrm, r1, r2, factor, br6, ebr6, f6, etrm1
real(dp) :: df6, ebr8, br8, f8, df8, br10, ebr10, f10, df10, r
real(dp) :: fg, fgprime, arg, earg, R0
integer  :: i, iu

factor = 1.0_dp       !! ??

do np = 1,nMMpot

write(fname,"(a,i2.2)") "MMpot.", np
call io_assign(iu)
open(iu,file=trim(fname),form="formatted",status="replace",  &
     position="rewind",action="write")

     rmin = 0.1_dp
     rmax = min(20.0_dp,MMcutoff)
     range = rmax - rmin
     delta = range/npts
     do i = 0, npts
        r1 = rmin + delta * i
        r2 = r1*r1

        if (nMMpottype(np).eq.1) then
           if (MMpotpar(2,np).eq.0.0_dp) then
              etrm = MMpotpar(1,np)/(r2**3)
              ftrm = 6.0_dp*factor*etrm/r2
              etrm = - etrm
           else
              br6 = MMpotpar(2,np)*sqrt(r2)
              ebr6 = exp(-br6)
              f6 = 1.0_dp + br6*(1.0_dp + 0.5_dp*br6*(1.0_dp + (br6/3.0_dp)*( &
                   1.0_dp + 0.25_dp*br6*(1.0_dp + 0.2_dp*br6*(1.0_dp + (br6/6.0_dp))))))
              f6 = 1.0_dp - f6*ebr6
              etrm1 = MMpotpar(1,np)/(r2**3)
              df6 = ebr6*(br6*(br6**6))/720.0_dp
              etrm = - etrm1*f6
              ftrm = factor*(6.0_dp*etrm1*f6 - etrm1*df6)/r2
           endif
        elseif (nMMpottype(np).eq.2) then
           if (MMpotpar(2,np).eq.0.0_dp) then
              etrm = MMpotpar(1,np)/(r2**4)
              ftrm = 8.0_dp*factor*etrm/r2
              etrm = - etrm
           else
              br8 = MMpotpar(2,np)*sqrt(r2)
              ebr8 = exp(-br8)
              f8 = 1.0_dp + br8*(1.0_dp + 0.5_dp*br8*(1.0_dp + (br8/3.0_dp)*( &
                   1.0_dp + 0.25_dp*br8*(1.0_dp + 0.2_dp*br8*(1.0_dp + &
                   (br8/6.0_dp)*(1.0_dp+(br8/7.0_dp)*(1.0_dp + 0.125_dp*br8)))))))
              f8 = 1.0_dp - f8*ebr8
              etrm1 = MMpotpar(1,np)/(r2**4)
              df8 = ebr8*(br8*(br8**8))/40320.0_dp
              etrm = - etrm1*f8
              ftrm = factor*(8.0_dp*etrm1*f8 - etrm1*df8)/r2
           endif
        elseif (nMMpottype(np).eq.3) then
           if (MMpotpar(2,np).eq.0.0_dp) then
              etrm = MMpotpar(1,np)/(r2**5)
              ftrm = 10.0_dp*factor*etrm/r2
              etrm = - etrm
           else
              br10 = MMpotpar(2,np)*sqrt(r2)
              ebr10 = exp(-br10)
              f10 = 1.0_dp + br10*(1.0_dp + 0.5_dp*br10*(1.0_dp + &
                   (br10/3.0_dp)*(1.0_dp + 0.25_dp*br10*(1.0_dp + 0.2_dp*br10*( &
                   1.0_dp + (br10/6.0_dp)*(1.0_dp + (br10/7.0_dp)*(1.0_dp + &
                   0.125_dp*br10*(1.0_dp + (br10/9.0_dp)*(1.0_dp + 0.1_dp*br10)))))))))
              f10 = 1.0_dp - f10*ebr10
              etrm1 = MMpotpar(1,np)/(r2**5)
              df6 = ebr10*(br10*(br10**10))/3628800.0_dp
              etrm = - etrm1*f10
              ftrm = factor*(10.0_dp*etrm1*f10 - etrm1*df10)/r2
           endif
        elseif (nMMpottype(np).eq.4) then
           r = sqrt(r2)
           etrm = MMpotpar(1,np)*(r - MMpotpar(2,np))
           ftrm = factor*etrm/r
           etrm = 0.5_dp*etrm*(r - MMpotpar(2,np))
        elseif (nMMpottype(np).eq.5) then    ! Grimme
           r = sqrt(r2)
           R0 = MMpotpar(2,np)
           arg = - d_grimme*(r/R0 - 1.0_dp)
           earg = exp(arg)
           fg = 1.0_dp / ( 1.0_dp + earg )
           etrm1 = s6_grimme * MMpotpar(1,np)/(r2**3)
           fgprime = (d_grimme/R0) * earg / ( 1 + earg )**2
           etrm = - etrm1*fg
           ftrm = factor*(6.0_dp*etrm1*fg - etrm1*r*fgprime)/r2
        endif

        write(iu,*) r1/Ang, etrm/eV, ftrm*Ang/eV

     enddo  !! points
     call io_close(iu)
   enddo     !! nPots

  end subroutine plot_functions

end module molecularmechanics
