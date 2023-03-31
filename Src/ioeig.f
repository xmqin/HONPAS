! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine ioeig(eo, ef, no, nspin, nk, maxo, nspinor, maxk,
     .                 kpoints, kweights)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999

      use precision, only : dp
      use siesta_cml
      use units, only : eV
      use files, only : slabel, label_length

      implicit          none

      integer,  intent(in) :: no     ! no_u: number of orbitals in unit cell
      integer,  intent(in) :: nspin  ! 'nspin_grid': 1, 2, 4 or 8
      integer,  intent(in) :: nk
      integer,  intent(in) :: maxo   ! no_u again
      integer,  intent(in) :: nspinor
      integer,  intent(in) :: maxk
      real(dp), intent(in) :: ef
      real(dp), intent(in), target :: eo(maxo, nspinor, maxk)
      real(dp), intent(in) :: kpoints(3,nk)
      real(dp), intent(in) :: kweights(nk)
      
      external          io_assign, io_close

c Internal 
      integer           ik, iu, io, is
      real(dp), pointer :: eok(:)

      character(len=label_length+4) :: fname
c -------------------------------------------------------------------

      fname = trim(slabel) // '.EIG'
      
      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(e17.9)") ef/eV
      ! The output corresponds to the number of bands.
      ! Thus it shouldn't be confused with the number of orbitals,
      ! although they are related!
      if ( nspin > nspinor ) then
        write(iu,"(tr1,i10,tr1,i0,tr1,i10)") 2*no, nspin, nk
      else
        ! This will always be nspin == nspinor with max(nspinor) == 2
        write(iu,"(tr1,i10,tr1,i0,tr1,i10)") no, nspin, nk
      end if
      do ik = 1,nk
        if ( nspin > nspinor ) then
          ! ensure we catch users doing neigwanted calculations
          ! In the NC/SOC case, the eigenvalues are written
          ! by the diag{2,3} routines to an eo(no_u*2,nk) array.
          ! So if neigwanted is, say, 0.8*no_u, there are 1.6*no_u
          ! bands, and independent loops over io and is, as in the
          ! second form below, would be wrong.
          call ravel(maxo * nspinor, eo(1,1,ik), eok)
          write(iu,"(i10,10(tr1,e17.9),/,(tr10,10(tr1,e17.9)))")
     .        ik, (eok(io)/eV,io=1,no*nspinor)
        else
          write(iu,"(i10,10(tr1,e17.9),/,(tr10,10(tr1,e17.9)))")
     .        ik, ((eo(io,is,ik)/eV,io=1,no),is=1,nspinor)
        end if
      enddo

      call io_close( iu )

      if (cml_p) then
         call cmlStartPropertyList(xf=mainXML, title="Eigenvalues")
         call cmlAddProperty(xf=mainXML, value=ef/eV, 
     .        title='Fermi Energy', dictref='siesta:E_Fermi', 
     .        fmt='r5', units='siestaUnits:ev')
         call cmlAddProperty(xf=mainXML, value=nk, 
     .        title='Number of k-points', dictRef='siesta:nkpoints',
     .        units='cmlUnits:countable')
         if ( nspin > nspinor ) then
           call cmlStartPropertyList(mainXML, dictRef='siesta:kpt_band')
           do ik = 1, nk
             call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .           weight=kweights(ik))
             ! ensure we catch users doing neigwanted calculations
             call ravel(maxo * nspinor, eo(1,1,ik), eok)
             call cmlAddProperty(xf=mainXML,
     .           value=eok(no*nspinor)/eV,
     .           dictRef='siesta:eigenenergies',
     .           units='siestaUnits:ev')
           end do
           call cmlEndPropertyList(mainXML)
         else
          do is = 1 , nspinor
            call cmlStartPropertyList(mainXML,
     .           dictRef='siesta:kpt_band')
            if ( is == 1 .and. nspinor > 1 ) then
               call cmlAddProperty(xf=mainXML, value="up", 
     .              dictRef="siesta:spin")
            else if ( nspinor > 1) then
               call cmlAddProperty(xf=mainXML, value="down", 
     .              dictRef="siesta:spin")
            end if
            do ik = 1, nk
               call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .              weight=kweights(ik))
               call cmlAddProperty(xf=mainXML, value=eo(1:no,is,ik)/eV, 
     .              dictRef='siesta:eigenenergies',
     .              units='siestaUnits:ev')
            enddo
            call cmlEndPropertyList(mainXML)
          end do
         end if
         call cmlEndPropertyList(mainXML)
      endif

      contains

      subroutine ravel(n, eo, eop)
      integer, intent(in) :: n
      real(dp), intent(in), target :: eo(n)
      real(dp), pointer :: eop(:)
      eop => eo(:)
      end subroutine
      
      end subroutine ioeig
