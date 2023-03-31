! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!*************************************************************************
!*                                                                    
!*  ORBMOL_PROJ projects the states of the adsorbed system (surface + adsorbed  
!*  molecule) on some specified states of the isolated molecule.          
!*                                                                        
!*                                                                        
!* written by R. Rurali and N. Lorente, June 2004                         
!*                                                                        
! Version with k-points, written by Alberto Garcia re-using code by Carolina
! Tabares and Pablo Ordejon. February 2009
!****************************************************************************

! USAGE:                                                                    
!                                                                           
! Run your SIESTA calculation of the adsorbed system and ask for the wave  
! functions in a sufficiently broad range across the Fermi level. In the    
! same calculation set the flag SaveHS to make the SIESTA write the overlap  
! matrix (you get H too, but it is not needed here).   
!
! For the wave functions, use COOP.Write T and use the "energy window"
! fdf symbols
!
!   WFS.EnergyMin -10.0 eV
!   WFS.EnergyMax 10.0 eV
!
! (The reason for the COOP tag is to get a WFSX file, with weight information,
!  and the k-point set used for self-consistency).
!
!  If you need to process legacy calculations, use the wfs2wfsx converter in
!  the Util/ directory.
!
! Separately, run a calculation of the molecule alone with the same 
! coordinates it has when it is adsorbed, without relaxing it (i.e. take 
! the relaxed coordinates of the adsorbed system, remove the surface and
! run a total energy calculation). Ask for the wave functions of all the
! states on which you want to project (HOMO and LUMO, for instance).
!
! **** This version assumes that the molecule appears first in the geometry
! **** list in the full calculation
! 
! The program needs four input files:
!
! 1) Main input file, read by standard input. A sample of input file is:
!
!    --- begin input file ---
!        C60_on_Si111.WFSX
!        C60.WFSX
!        C60_on_Si111.HS
!        delta_e, sigma
!        orbmol.out
!    --- end input file ---
!
!    where:
!    - The first line is the name of the wavefunctions file 
!      of the adsorbed system
!    - The second line is the name of the wavefunctions file
!      of the molecule
!    - The third is the file where the Hamiltonian and Overlap 
!      matrix have been stored 
!    - The fourth line contains two real numbers which specify 
!      the parameter for the histogram construction (energy interval 
!      width and sigma of the gaussian)
!    - The fifth line is the output file
!       
!***************************************************************************


program orbmol_proj

   implicit none

   integer, parameter  :: sp = selected_real_kind(6,30)
   integer, parameter  :: dp = selected_real_kind(14,100)

!  Local variables

   character( len = 200 )  :: fname_sys, fname_mol, fname_hs, &
                             fname_out, fname_out_up,        &
                             fname_out_down                    ! filenames 

   integer :: i, j, n, m, ih, ind, jj, nnz, im, io, is0, iw, k, iw0
   integer :: nao, nkp, nsp, no_s, ik
   integer :: wfs_u = 1, wfs_m=2, iu_hs=3, stdin=5, iu_out=4,            &
              iu_out_up=11, iu_out_down=12
   logical :: gamma, debug, gamma_wfsx, gamma_wfsx_m

! variables of the whole system

   integer :: n_basis_sys,                                   &
              n_wave_sys, wave_sys,                          &
              n_k_point_sys, k_point_sys,                    &
              n_spin_sys, spin_sys

   integer, allocatable :: indxuo(:)
   real(dp), allocatable :: wk(:), pk(:,:), xij(:,:)
   real (sp), allocatable :: wf(:,:)

! variables of the molecule

   integer :: n_basis_mol,                                   &
              n_wave_mol, wave_mol,                          &
              n_k_point_mol, k_point_mol,                    &
              n_spin_mol, spin_mol                          

   real(sp), allocatable    :: c_mol(:,:,:)

! variables of the HS

   integer :: no_u,  n_spin_hs, n_dummy_hs
   integer, allocatable :: numh(:), listhptr(:), listh(:)
   real(dp), allocatable :: Ssparse(:)

! variables of the PDOS

   integer :: is, e_int, ibin, idummy

   real(dp)    :: e_min, e_max, delta_e, factor, ene, sigma, eigval
   real(dp)    :: factor_re, factor_im, qtot, temp
   real(dp)    :: min_energy, max_energy, phase
   integer     :: number_of_wfns, nwfmx
   real(dp), allocatable :: pdos(:,:,:)
   real(dp), parameter :: zero = 0.0d0
   real(dp), parameter :: pi   = 3.141592654d0



! Start of program

! Read from standard input

  open( unit=stdin )
  read(stdin,'(a30)') fname_sys
  read(stdin,'(a30)') fname_mol
  read(stdin,'(a30)') fname_hs
  read(stdin,*)       delta_e, sigma
  read(stdin,'(a30)') fname_out

! -------------------------------------------------------------------
! 
! Now I read the wave function coefficients of the states of the 
! whole system
!
! -------------------------------------------------------------------

  open(wfs_u,file=fname_sys,status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) nsp
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "Whole-system WFSX read: nkp, nsp, nao: ", nkp, nsp, nao

  nwfmx = 0
  min_energy = huge(1.0_dp)
  max_energy = -huge(1.0_dp)

  do ik=1,nkp
     do is=1,nsp

        read(wfs_u) idummy, pk(1:3,ik), wk(ik)
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_u) is0
        read(wfs_u) number_of_wfns
        nwfmx = max(nwfmx,number_of_wfns)

        do iw=1,number_of_wfns
           read(wfs_u) iw0
           read(wfs_u) eigval
           min_energy = min(min_energy,eigval)
           max_energy = max(max_energy,eigval)
           read(wfs_u)
        enddo
     enddo
  enddo

  print *, "K-points used: ", .not. gamma_wfsx
  print *, " Maximum number of whole-system wfs per k-point: ", nwfmx
  print "(a,2f12.4)", "Min_energy, max_energy on whole-system WFS file: ",  &
       min_energy, max_energy
            
   e_min = min_energy
   e_max = max_energy

   n_k_point_sys = nkp
   n_spin_sys = nsp
   n_basis_sys = nao

   print *, "System: k-points, nspin, num_orbs: ", n_k_point_sys, n_spin_sys, n_basis_sys
   print *, "Energy range in whole system WFS file: ", e_min, e_max



! -------------------------------------------------------------------
!
! In exactly the same way, I now read the wave function coefficients 
! of the molecular states I am interested in
!
! -------------------------------------------------------------------

  open(wfs_m,file=fname_mol,status='old',form='unformatted')
  read(wfs_m) nkp, gamma_wfsx_m
  if (.not. gamma_wfsx_m) STOP "molecule uses k-points!"
  read(wfs_m) nsp
  read(wfs_m) nao
  read(wfs_m)        !! Symbols, etc

   n_k_point_mol = nkp
   n_spin_mol    = nsp
   n_basis_mol   = nao

  if (debug) print *, "Molecule WFSX read: nkp, nsp, nnao: ", nkp, nsp, nao

  nwfmx = 0
  min_energy = huge(1.0_dp)
  max_energy = -huge(1.0_dp)

  do ik=1,n_k_point_mol
     do is=1,n_spin_mol

        read(wfs_m) idummy  ! , pkm(1:3,ik), wkm(ik)
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_m) is0
        read(wfs_m) number_of_wfns
        nwfmx = max(nwfmx,number_of_wfns)

        if ( ( ik == 1 ) .and. ( is == 1 ) ) then
           n_wave_mol = number_of_wfns
           allocate( c_mol( n_basis_sys, number_of_wfns, nsp ) )
           c_mol = zero
        end if
        do iw=1,number_of_wfns
           read(wfs_m) iw0
           read(wfs_m) eigval
           min_energy = min(min_energy,eigval)
           max_energy = max(max_energy,eigval)
           read(wfs_m) c_mol(1:n_basis_mol,iw,is)
        enddo
     enddo
  enddo

  print *, "K-points used: ", .not. gamma_wfsx_m
  print *, " Maximum number of molecule wfs per k-point: ", nwfmx
  print "(a,2f12.4)", "Min_energy, max_energy on molecule WFSX file: ",  &
       min_energy, max_energy
            
  close(wfs_m)

! -------------------------------------------------------------------
! Finally I read the S matrix in sparse form 

   open( iu_hs, file=fname_hs, form='unformatted', status='unknown' )

   read(iu_hs) no_u, no_s, n_spin_hs, n_dummy_hs

   if ( no_u /= n_basis_sys ) then
      write(*,*) '             ERROR! ERROR!                        '
      write(*,*) '    Dimension conflict: the S matrix is not       '
      write(*,*) '     coherent with the system wave function       '
      stop
   end if

   allocate( numh(no_u), listhptr(no_u))

   write(*,*) '*****************************************'
   write(*,*) '*    S MATRIX OF THE ADSORBED SYSTEM    *'
   write(*,*) '*****************************************'

   read(iu_hs) gamma   ! Flag to signal k-point usage
   if (.not.gamma) then
      allocate(indxuo(no_s))
      read(iu_hs) (indxuo(ih),ih=1,no_s) 
   endif

   do i = 1, no_u
      read(iu_hs) numh(i)
   end do
   nnz = sum(numh(1:no_u))

   listhptr(1) = 0
   do ih = 2,no_u
      listhptr(ih) = listhptr(ih-1) + numh(ih-1)
   enddo

   allocate(listh(nnz),Ssparse(nnz),xij(3,nnz))
   do ih = 1,no_u
      do im = 1,numh(ih)
         read(iu_hs) listh(listhptr(ih)+im)
      enddo
   enddo

   do n = 1, n_spin_hs
      do i = 1, no_u
         do j = 1, numh(i)
            read(iu_hs)  !!! Hamiltonian
         end do
      end do
   end do

! S
   do ih = 1,no_u
      do im = 1,numh(ih)
         read(iu_hs) Ssparse(listhptr(ih)+im)
      enddo
   enddo

   read(iu_hs) qtot,temp
   if (.not.gamma) then
!C Read interorbital vectors for K point phasing
       do ih = 1,no_u
           do im = 1,numh(ih)
               read(iu_hs) (xij(k,listhptr(ih)+im),k=1,3)
           enddo
       enddo
   endif

   close(iu_hs)

! -------------------------------------------------------------------
!
! Let's finally calculate the projection over the specified 
! molecular orbitals
!
! -------------------------------------------------------------------

! At first I loop over the different molecular orbitals on which
! I want to project the states

   write(*,*)
   write(*,*) 'The calculation of the PDOS has started'


     if (gamma) then
        allocate(wf(1,1:no_u))
     else
        allocate(wf(2,1:no_u))
     endif

   e_int = ( e_max - e_min ) / delta_e + 1
   allocate( pdos( e_int, n_wave_mol, n_spin_sys ) )

   pdos = zero

   debug = .false.
        rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        if (debug) print *, "Number of k-points, spins: ", nkp, nsp
        do ik=1,n_k_point_sys
           if (debug) print *, "k-point: ", ik
           do is=1,n_spin_sys
              read(wfs_u)
              read(wfs_u)
              read(wfs_u)  number_of_wfns
              if (debug) print *, "  Number of wfns: ", number_of_wfns
              do iw=1,number_of_wfns
                 if (debug) print *, "     wfn: ", iw
                 read(wfs_u) 
                 read(wfs_u) eigval
                 ! Early termination of iteration
                 ! Note that we keep a few more states on the sides, due to
                 ! the smearing
                 !!if (eigval < low_e .or. eigval > high_e) then
                 !!   read(wfs_u)   ! Still need to read this
                 !!   CYCLE
                 !!endif

                 read(wfs_u) (wf(:,io), io=1,no_u)

                 do n = 1, n_wave_mol

                    factor = zero
                    factor_re = zero
                    factor_im = zero

!** AG ** This assumes that the "molecule" atoms appear first in the
!         geometry list of the whole system

         ! Note the reversal of the indexes with respect to Pablo's. 
         ! j refers to the "whole system" and i to the molecule
         ! Since the absolute value of the scalar product is taken,
         ! this does not matter. The phase should also be correctly
         ! calculated

                    do i = 1, n_basis_mol
                       do jj = 1, numh(i)
                          ind = listhptr(i) + jj
                          j = listh(ind)

                          if (gamma) then
                             factor = factor + c_mol(i,n,is) * wf(1,j) * Ssparse(ind) 
                          else
                             j = indxuo(j)
                             phase=dot_product(pk(1:3,ik),xij(1:3,ind))
                             factor_re = factor_re + c_mol(i,n,is) *    &
                                      (wf(1,j) * cos(phase) -        &
                                       wf(2,j) * sin(phase)) *       &
                                      Ssparse(ind)
                             factor_im = factor_im + c_mol(i,n,is) *    &
                                      (wf(1,j) * sin(phase) +        &
                                       wf(2,j) * cos(phase)) *       &
                                      Ssparse(ind)
                          endif
                          
                       end do
                    end do

                    if (.not. gamma) factor = sqrt(factor_re*factor_re + factor_im*factor_im)
                    do ibin = 1, e_int
                       ene = e_min + ( ibin - 1 ) * delta_e
                       pdos( ibin,n,is ) = pdos( ibin,n,is ) + factor * factor * wk(ik) *              &
                            exp( - ( ene - eigval ) ** 2 / sigma ** 2 )
                    end do

                 end do    ! n (mol wfn)

              end do     ! iw 
           end do      ! is
        enddo       ! ik

   pdos = pdos / ( dsqrt( pi ) * sigma )

!===========================================================================
   if ( n_spin_sys == 1 ) then

      open( iu_out, file=fname_out )

      do i = 1, e_int 
         write(iu_out,'(10f12.5)') e_min + ( i - 1) * delta_e, pdos(i,:,1)
      end do

      close( iu_out )

   else

      fname_out_up   = 'orbmol.out.UP'
      fname_out_down = 'orbmol.out.DOWN'

      open( iu_out_up,   file=fname_out_up )

      do i = 1, e_int
         write(iu_out_up,'(20f12.5)') e_min + ( i - 1) * delta_e, pdos(i,:,1)
      end do

      open( iu_out_down, file=fname_out_down )

      do i = 1, e_int
         write(iu_out_down,'(20f12.5)') e_min + ( i - 1) * delta_e, pdos(i,:,2)
      end do


   end if

end program orbmol_proj

