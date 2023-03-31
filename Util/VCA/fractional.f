! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      program fractional
!
!     Program to multiply the Pseudopotential strength by a given
!     fraction to simulate "fractional occupation" in the VCA
!
      use precision,       only: dp
      use pseudopotential, only: pseudopotential_t,
     $                           pseudo_read, pseudo_write_formatted
      use periodic_table,  only: cnfig, qvlofz
      use f2kcli

      implicit none

      type(pseudopotential_t), target :: pot1, pmix
      type(pseudopotential_t), pointer :: p1,  p
      
      real(dp) :: xfraction, z1
      integer  :: i, iostat, nargs, l
      character(len=200) :: name, name1, xmixstr
      integer, dimension(0:3) :: cnfig1, config
      real(dp), dimension(0:3) :: occup, occup1


      nargs = command_argument_count()
      if (nargs /= 2) then
         write(0,*) "Usage: fractional ps1 xfraction (no .psf)"
         stop
      endif

      call get_command_argument(1,value=name1,status=iostat)
      if (iostat == 0) then
         call pseudo_read(trim(name1), pot1)
      else
         call die("Cannot get first argument")
      endif

      call get_command_argument(2,value=xmixstr,status=iostat)
      if (iostat == 0) then
         read(xmixstr,fmt=*,iostat=iostat) xfraction
         if (iostat /= 0) call die("Cannot parse xfraction")
      else
         call die("Cannot get xfraction")
      endif

      p1 => pot1
      p => pmix

      print *, "Generation valence: ", p1%gen_zval

      z1 = nucl_z(p1%name)

      call cnfig(int(z1),cnfig1(0:3))
      print *, "Cnfig1: ", cnfig1

      call qvlofz(int(z1),occup1(0:3))
      print *, "Occup1: ", occup1

      do l = 0, 3
         config(l) = cnfig1(l)
         occup(l) = occup1(l)*xfraction
      enddo
!

      p%name = "SX"
      p%nr = p1%nr
      p%nrval = p1%nrval
      p%zval = xfraction*p1%zval 
      p%gen_zval = xfraction*p1%gen_zval
      p%relativistic = p1%relativistic
      p%correlation = p1%correlation
      p%icorr = p1%icorr
      p%irel = p1%irel
      p%nicore = p1%nicore
      p%a = p1%a            ! Use p1's grid
      p%b = p1%b

      p%method(1) = "FRACT"
      p%method(2:3) = p1%method(1:2)
      p%method(4) = " --- "
      p%method(5:6) = " fraction "

      write(p%text,"(f8.5,a)")
     $             xfraction, " fraction of " // p1%name

      p%npotu = p1%npotu
      p%npotd = p1%npotd
      allocate(p%r(size(p1%r)))
      p%r = p1%r
      allocate(p%chcore(size(p1%r)))
      allocate(p%chval(size(p1%r)))
      p%chcore = xfraction * p1%chcore
      p%chval = xfraction * p1%chval
      allocate(p%vdown(p%npotd,size(p1%r)))
      allocate(p%ldown(p%npotd))
      do i = 1, p%npotd
         p%vdown(i,:) = xfraction*p1%vdown(i,:)
         p%ldown(i) = p1%ldown(i)
      enddo
      if (p%npotu /= 0) then
         allocate(p%vup(p%npotu,size(p1%r)))
         allocate(p%lup(p%npotu))
         do i = 1, p%npotu
            p%vup(i,:) = xfraction*p1%vup(i,:)
            p%lup(i) = p1%lup(i)
         enddo
      endif

      write(name,"(a,a,f7.5,a)")
     $        trim(name1), "-Fraction-", xfraction, ".psf"
      call pseudo_write_formatted(name,p)
      write(name,"(a,a,f7.5,a)")
     $        trim(name1), "-Fraction-", xfraction, ".synth"
      open (unit=4,file=name, form="formatted", status="unknown",
     $      action="write", position="rewind")
      write(4,"(a)") "%block SyntheticAtoms"
      write(4,"(i2)") 1
      write(4,"(4i3)") config(0:3)
      write(4,"(4f12.6)") occup(0:3)
      write(4,"(a)") "%endblock SyntheticAtoms"
      close(4)
!
!     This kludge serves to make life easier for shell scripts.
!     They can pick up the label from this file inmediately after
!     invocation of mixps
!
      open (unit=5,file="FRACTLABEL",
     $     form="formatted", status="unknown",
     $      action="write", position="rewind")
      write(5,"(a,a,f7.5,a)")
     $        trim(name1), "-Fraction-", xfraction
      close(5)

      CONTAINS

      function nucl_z(name)
      character(len=2), intent(in) ::  name
      real(dp)                     ::  nucl_z
c
c      function determines the nuclear charge of an element
c
      real(dp), parameter :: one = 1.0_dp
      real(dp) :: charge

      if (name .eq. 'ZN' .or. name .eq. 'ZR') then
         charge = 0.0*one
      else if (name .eq. 'H ' .or. name .eq. ' H') then
         charge = 1*one
      else if (name .eq. 'He') then
         charge = 2*one
      else if (name .eq. 'Li') then
         charge = 3*one
      else if (name .eq. 'Be') then
         charge = 4*one
      else if (name .eq. 'B ' .or. name .eq. ' B') then
         charge = 5*one
      else if (name .eq. 'C ' .or. name .eq. ' C') then
         charge = 6*one
      else if (name .eq. 'N ' .or. name .eq. ' N') then
         charge = 7*one
      else if (name .eq. 'O ' .or. name .eq. ' O') then
         charge = 8*one
      else if (name .eq. 'F ' .or. name .eq. ' F') then
         charge = 9*one
      else if (name .eq. 'Ne') then
         charge = 10*one
      else if (name .eq. 'Na') then
         charge = 11*one
      else if (name .eq. 'Mg') then
         charge = 12*one
      else if (name .eq. 'Al') then
         charge = 13*one
      else if (name .eq. 'Si') then
         charge = 14*one
      else if (name .eq. 'P ' .or. name .eq. ' P') then
         charge = 15*one
      else if (name .eq. 'S ' .or. name .eq. ' S') then
         charge = 16*one
      else if (name .eq. 'Cl') then
         charge = 17*one
      else if (name .eq. 'Ar') then
         charge = 18*one
      else if (name .eq. 'K ' .or. name .eq. ' K') then
         charge = 19*one
      else if (name .eq. 'Ca') then
         charge = 20*one
      else if (name .eq. 'Sc') then
         charge = 21*one
      else if (name .eq. 'Ti') then
         charge = 22*one
      else if (name .eq. 'V ' .or. name .eq. ' V') then
         charge = 23*one
      else if (name .eq. 'Cr') then
         charge = 24*one
      else if (name .eq. 'Mn') then
         charge = 25*one
      else if (name .eq. 'Fe') then
         charge = 26*one
      else if (name .eq. 'Co') then
         charge = 27*one
      else if (name .eq. 'Ni') then
         charge = 28*one
      else if (name .eq. 'Cu') then
         charge = 29*one
      else if (name .eq. 'Zn') then
         charge = 30*one
      else if (name .eq. 'Ga') then
         charge = 31*one
      else if (name .eq. 'Ge') then
         charge = 32*one
      else if (name .eq. 'As') then
         charge = 33*one
      else if (name .eq. 'Se') then
         charge = 34*one
      else if (name .eq. 'Br') then
         charge = 35*one
      else if (name .eq. 'Kr') then
         charge = 36*one
      else if (name .eq. 'Rb') then
         charge = 37*one
      else if (name .eq. 'Sr') then
         charge = 38*one
      else if (name .eq. 'Y ' .or. name .eq. ' Y') then
         charge = 39*one
      else if (name .eq. 'Zr') then
         charge = 40*one
      else if (name .eq. 'Nb') then
         charge = 41*one
      else if (name .eq. 'Mo') then
         charge = 42*one
      else if (name .eq. 'Tc') then
         charge = 43*one
      else if (name .eq. 'Ru') then
         charge = 44*one
      else if (name .eq. 'Rh') then
         charge = 45*one
      else if (name .eq. 'Pd') then
         charge = 46*one
      else if (name .eq. 'Ag') then
         charge = 47*one
      else if (name .eq. 'Cd') then
         charge = 48*one
      else if (name .eq. 'In') then
         charge = 49*one
      else if (name .eq. 'Sn') then
         charge = 50*one
      else if (name .eq. 'Sb') then
         charge = 51*one
      else if (name .eq. 'Te') then
         charge = 52*one
      else if (name .eq. 'I ' .or. name .eq. ' I') then
         charge = 53*one
      else if (name .eq. 'Xe') then
         charge = 54*one
      else if (name .eq. 'Cs') then
         charge = 55*one
      else if (name .eq. 'Ba') then
         charge = 56*one
      else if (name .eq. 'La') then
         charge = 57*one
      else if (name .eq. 'Ce') then
         charge = 58*one
      else if (name .eq. 'Pr') then
         charge = 59*one
      else if (name .eq. 'Nd') then
         charge = 60*one
      else if (name .eq. 'Pm') then
         charge = 61*one
      else if (name .eq. 'Sm') then
         charge = 62*one
      else if (name .eq. 'Eu') then
         charge = 63*one
      else if (name .eq. 'Gd') then
         charge = 64*one
      else if (name .eq. 'Tb') then
         charge = 65*one
      else if (name .eq. 'Dy') then
         charge = 66*one
      else if (name .eq. 'Ho') then
         charge = 67*one
      else if (name .eq. 'Er') then
         charge = 68*one
      else if (name .eq. 'Tm') then
         charge = 69*one
      else if (name .eq. 'Yb') then
         charge = 70*one
      else if (name .eq. 'Lu') then
         charge = 71*one
      else if (name .eq. 'Hf') then
         charge = 72*one
      else if (name .eq. 'Ta') then
         charge = 73*one
      else if (name .eq. 'W ' .or. name .eq. ' W') then
         charge = 74*one
      else if (name .eq. 'Re') then
         charge = 75*one
      else if (name .eq. 'Os') then
         charge = 76*one
      else if (name .eq. 'Ir') then
         charge = 77*one
      else if (name .eq. 'Pt') then
         charge = 78*one
      else if (name .eq. 'Au') then
         charge = 79*one
      else if (name .eq. 'Hg') then
         charge = 80*one
      else if (name .eq. 'Tl') then
         charge = 81*one
      else if (name .eq. 'Pb') then
         charge = 82*one
      else if (name .eq. 'Bi') then
         charge = 83*one
      else if (name .eq. 'Po') then
         charge = 84*one
      else if (name .eq. 'At') then
         charge = 85*one
      else if (name .eq. 'Rn') then
         charge = 86*one
      else if (name .eq. 'Fr') then
         charge = 87*one
      else if (name .eq. 'Ra') then
         charge = 88*one
      else if (name .eq. 'Ac') then
         charge = 89*one
      else if (name .eq. 'Th') then
         charge = 90*one
      else if (name .eq. 'Pa') then
         charge = 91*one
      else if (name .eq. ' U' .or. name .eq. 'U ') then
         charge = 92*one
      else if (name .eq. 'Np') then
         charge = 93*one
      else if (name .eq. 'Pu') then
         charge = 94*one
      else if (name .eq. 'Am') then
         charge = 95*one
      else if (name .eq. 'Cm') then
         charge = 96*one
      else if (name .eq. 'Bk') then
         charge = 97*one
      else if (name .eq. 'Cf') then
         charge = 98*one
      else if (name .eq. 'Es') then
         charge = 99*one
      else if (name .eq. 'Fm') then
         charge = 100*one
      else if (name .eq. 'Md') then
         charge = 101*one
      else if (name .eq. 'No') then
         charge = 102*one
      else if (name .eq. 'Lr') then
         charge = 103*one
      else
         write(6,9000) name
 9000    format(//'element ',a2,' unknown')
         call die("Unknown element")
      end if

      nucl_z = charge

      end function nucl_z


      end program fractional


