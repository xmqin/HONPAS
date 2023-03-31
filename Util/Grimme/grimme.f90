! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This program has been implemented by:
!  Nick Papior, 2014
! Program for creating an output for fdf with the 
! input of a standard input file
program grimme_program

  use f2kcli
  use fdf
  use periodic_table
  use chemical, only: read_chemical_types, number_of_species, species_label
  use chemical, only: atomic_number, is_floating
  use precision, only : dp

  implicit none 

  type grimme
     integer :: Z = 0 ! atomic number H = 1
     character(len=3) :: atom ! atomic name  H = 'H'
     real(dp) :: r0 ! vdw-radius (in Ang)
     real(dp) :: C6 ! parameter (in eV/Ang^6)
     character(len=100) :: opt ! any option given if one needs more different parameters for the same element
     type(grimme), pointer :: next => null()
  end type grimme

  type(grimme), pointer :: first => null()

  character(len=100) :: filein
  logical :: exists
  ! default variables
  integer :: narg, i
  ! species loop variables
  integer :: ns, is, js

  ! Here we start the routines
  filein = 'none'
  narg = command_argument_count()
  do while ( narg > 0 )
     if (narg == 1) then
        call get_command_argument(narg,filein)
     end if
     narg = narg - 1
  end do
  if ( leqi(filein,'none') ) then
     write(0,'(a)') 'Input file can not be piped into this program, &
          &please supply FDF file on command line...'
     stop
  end if

  ! check whether the file exists
  inquire(file=filein,exist=exists)
  if (.not. exists ) then
     write(0,'(a)') 'Input file can not be piped into this program, &
          &please supply FDF file on command line...'
     stop
  end if





  ! add default parameters
  call add_grimme_parameters()

  ! For Iridium arXiv: 1311.7022v2 gives a hint which actually is
  ! the same as found by Grimme. (so apparantly no need to change it.

  ! in DOI: 10.1039/c3cp51791h we find for Platinum:
  call add_grimme(78,0.175_dp,42.44_dp,'nm','J',opt='10.1039/c3cp51791h')

  ! in arXiv: 1011.1113 we find for Gold
  call add_grimme(79,0.1772_dp,40.62_dp,'nm','J',opt='1011.1113')


  ! **** *** *** *** *** ** *** *** ** ** *** *** **
  ! *** insert new Grimme-parameters here
  ! **** *** *** *** *** ** *** *** ** ** *** *** **


  ! open the fdf file
  call fdf_init(filein,"grimme.log")

  ! create the data for the species list
  call read_chemical_types(silent=.true.)
  ns = number_of_species()

  ! Write to standard output the commands needed in the fdf file
  write(*,'(a)') 'MM.UnitsDistance Ang  # what this program prints out DO NOT CHANGE'
  write(*,'(a)') 'MM.UnitsEnergy    eV  # what this program prints out DO NOT CHANGE'
  write(*,'(a)') 'MM.Grimme.S6     0.75 # Grimme-paper for PBE (correct for your functional)'
  write(*,'(a)') 'MM.Grimme.D     20.   # Grimme-paper (correct for your functional)'
  write(*,'(a)') '%block MM.Potentials'
  
  do is = 1 , ns
     ! we don't add the Grimme parameter for floating 
     ! orbitals
     if ( is_floating(is) ) cycle

     ! loop over everything that it needs to contact with
     do js = is , ns
        if ( is_floating(js) ) cycle

        call write_grimme(first,is,js)
     end do

  end do
  write(*,'(a)') '%endblock MM.Potentials'

contains

  subroutine add_grimme(Z,r0,c6,L,E,opt,atom)
    integer, intent(in) :: Z
    real(dp), intent(in) :: r0, c6
    character(len=*), intent(in) :: L,E
    character(len=*), intent(in), optional :: opt,atom

    type(grimme), pointer :: n => null(), new => null()
    character(len=3) :: latom

    real(dp), parameter :: NAvogadro = 6.022e23_dp
    

    latom = symbol(Z)
    if ( present(atom) ) latom = atom

    allocate(new)
    new%Z = Z
    new%atom = latom
    new%r0 = r0 * fdf_convfac(L,'Ang')
    new%c6 = c6 * fdf_convfac(E,'eV') * fdf_convfac(L,'Ang') ** 6 / NAvogadro
    new%opt = '10.1002/jcc.20495'
    if ( present(opt) ) new%opt = opt

    if ( .not. associated(first) ) then
       first => new
       return
    end if

    n => first
    if ( same(n,new) ) then ! the first atom gets replaced
       new%next => first%next
       deallocate(first)
       first => new
       return
    end if
    do while ( associated(n%next) )
       if ( same(n%next,new) ) then ! the next gets replaced
          new%next => n%next%next
          deallocate(n%next)
          n%next => new
          return
       end if
       n => n%next
    end do
    if ( same(n,new) ) then
       stop 'Grimme functional and parameters are the same'
    end if

    n%next => new

  end subroutine add_grimme

  function same(g1,g2) result(ret)
    type(grimme), intent(in) :: g1, g2
    logical :: ret
    ret = g1%z == g2%z
    if (.not. ret ) return
    ret = g1%atom == g2%atom
    if (.not. ret ) return
  end function same

  subroutine write_grimme(branch,is,js)
    type(grimme), pointer :: branch
    integer, intent(in) :: is, js
    type(grimme), pointer :: igrim, jgrim

    call associate_grimme(branch,is,igrim)
    call associate_grimme(branch,js,jgrim)
    
    if ( .not. associated(igrim) ) then
       stop 'Could not find one of the species.. Weird'
    end if
    if ( .not. associated(jgrim) ) then
       stop 'Could not find one of the species.. Weird'
    end if

    if ( is == js ) then
       write(*,'(2(tr1,i2,tr1),a,f10.2,tr1,f10.3,tr1,''# '',a,'', '',a)') &
            is,js,'Grimme', &
            sqrt(igrim%c6*jgrim%c6), & ! C6
            igrim%r0+jgrim%r0, & ! R0
            trim(igrim%atom),trim(igrim%opt)
    else
       write(*,'(2(tr1,i2,tr1),a,f10.2,tr1,f10.3,tr1,''# '',a,'' / '',a)') &
            is,js,'Grimme', &
            sqrt(igrim%c6*jgrim%c6), & ! C6
            igrim%r0+jgrim%r0,trim(igrim%atom),trim(jgrim%atom)
    end if
    
  end subroutine write_grimme

  subroutine associate_grimme(branch,is,grim)
    type(grimme), pointer :: branch
    integer, intent(in) :: is
    type(grimme), pointer :: grim
    integer :: Z
    type(grimme), pointer :: n
    
    ! nullify
    nullify(grim)

    if ( .not. associated(branch) ) return

    Z = atomic_number(is)

    n => branch
    do 
       if ( n%z == Z ) then
          grim => n
          return
       end if
       if ( .not. associated(n%next) ) exit
       n => n%next
    end do

  end subroutine associate_grimme
  
  subroutine add_grimme_parameters()
    integer :: i
    ! All these values comes from the program dftd3 vs: 3.0.2
    ! All have been corrected by the factor 1.1 so as to correspond with
    ! the values found in the paper by Grimme: 10.1002/jcc.20495

    call add_grimme( 1,0.1001_dp,  0.14_dp,'nm','J')
    call add_grimme( 2,0.1012_dp,  0.08_dp,'nm','J')
    call add_grimme( 3,0.0825_dp,  1.61_dp,'nm','J')
    call add_grimme( 4,0.1408_dp,  1.61_dp,'nm','J')
    call add_grimme( 5,0.1485_dp,  3.13_dp,'nm','J')
    call add_grimme( 6,0.1452_dp,  1.75_dp,'nm','J')
    call add_grimme( 7,0.1397_dp,  1.23_dp,'nm','J')
    call add_grimme( 8,0.1342_dp,  0.70_dp,'nm','J')
    call add_grimme( 9,0.1287_dp,  0.75_dp,'nm','J')
    call add_grimme(10,0.1243_dp,  0.63_dp,'nm','J')
    call add_grimme(11,0.1144_dp,  5.71_dp,'nm','J')
    call add_grimme(12,0.1364_dp,  5.71_dp,'nm','J')
    call add_grimme(13,0.1639_dp, 10.79_dp,'nm','J')
    call add_grimme(14,0.1716_dp,  9.23_dp,'nm','J')
    call add_grimme(15,0.1705_dp,  7.84_dp,'nm','J')
    call add_grimme(16,0.1683_dp,  5.57_dp,'nm','J')
    call add_grimme(17,0.1639_dp,  5.07_dp,'nm','J')
    call add_grimme(18,0.1595_dp,  4.61_dp,'nm','J')
    call add_grimme(19,0.1485_dp, 10.80_dp,'nm','J')
    call add_grimme(20,0.1474_dp, 10.80_dp,'nm','J')
    do i = 21 , 30
       call add_grimme(i,0.1562_dp, 10.80_dp,'nm','J')
    end do
    call add_grimme(31,0.1650_dp, 16.99_dp,'nm','J')
    call add_grimme(32,0.1727_dp, 17.10_dp,'nm','J')
    call add_grimme(33,0.1760_dp, 16.37_dp,'nm','J')
    call add_grimme(34,0.1771_dp, 12.64_dp,'nm','J')
    call add_grimme(35,0.1749_dp, 12.47_dp,'nm','J')
    call add_grimme(36,0.1727_dp, 12.01_dp,'nm','J')
    call add_grimme(37,0.1628_dp, 24.67_dp,'nm','J')
    call add_grimme(38,0.1606_dp, 24.67_dp,'nm','J')
    do i = 39 , 48
       call add_grimme(i,0.1639_dp, 24.67_dp,'nm','J')
    end do
    call add_grimme(49,0.1672_dp, 37.32_dp,'nm','J')
    call add_grimme(50,0.1804_dp, 38.71_dp,'nm','J')
    call add_grimme(51,0.1881_dp, 38.44_dp,'nm','J')
    call add_grimme(52,0.1892_dp, 31.74_dp,'nm','J')
    call add_grimme(53,0.1892_dp, 31.50_dp,'nm','J')
    call add_grimme(54,0.1881_dp, 29.99_dp,'nm','J')
    call add_grimme(55,0.1802_dp,315.27_dp,'nm','J')
    call add_grimme(56,0.1762_dp,226.99_dp,'nm','J')
    call add_grimme(57,0.1720_dp,176.25_dp,'nm','J')
    do i = 58 , 71
       call add_grimme(i,0.1753_dp,140.68_dp,'nm','J')
    end do
    call add_grimme(72,0.1787_dp,105.11_dp,'nm','J')
    do i = 73 , 79
       call add_grimme(i,0.1772_dp, 81.24_dp,'nm','J')
    end do
    call add_grimme(80,0.1758_dp, 57.36_dp,'nm','J')
    call add_grimme(81,0.1986_dp, 57.25_dp,'nm','J')
    call add_grimme(82,0.1944_dp, 63.16_dp,'nm','J')
    call add_grimme(83,0.1898_dp, 63.54_dp,'nm','J')
    call add_grimme(84,0.2005_dp, 55.28_dp,'nm','J')
    call add_grimme(85,0.1991_dp, 57.17_dp,'nm','J')
    call add_grimme(86,0.1924_dp, 56.64_dp,'nm','J')

  end subroutine add_grimme_parameters

end program grimme_program

