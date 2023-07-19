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
      subroutine ioeig(eo, ef, no, ns, nk, maxo, maxs, maxk,
     .                 kpoints, kweights)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999

      use fdf
      use precision, only : dp
      use siesta_cml
      use units, only : eV
      use files, only : slabel, label_length

      implicit          none

      integer,  intent(in) :: maxo
      integer,  intent(in) :: maxs
      integer,  intent(in) :: maxk
      real(dp), intent(in) :: eo(maxo, maxs, maxk)
      real(dp), intent(in) :: ef
      integer,  intent(in) :: no
      integer,  intent(in) :: ns
      integer,  intent(in) :: nk
      real(dp), intent(in) :: kpoints(3,nk)
      real(dp), intent(in) :: kweights(nk)
      
      external          io_assign, io_close, paste

c Internal 
      integer           ik, iu, io, is, nspin

      character(len=label_length+4), save :: fname
      logical, save                       :: frstme = .true.
c -------------------------------------------------------------------

      if (frstme) then
        fname = slabel
        fname = trim(fname) // '.EIG'
        frstme = .false.
      endif
      
      nspin = min(ns,2)

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(f14.4)") ef/eV
      write(iu,"(3i6)")   no, min(ns,2), nk
      do ik = 1,nk
        write(iu,"(i5,10f12.5,/,(5x,10f12.5))")
     .          ik, ((eo(io,is,ik)/eV,io=1,no),is=1,nspin)
      enddo

      call io_close( iu )

      if (cml_p) then
        call cmlStartPropertyList(xf=mainXML, title="Eigenvalues")
        call cmlAddProperty(xf=mainXML, value=ef/eV, 
     .       title='Fermi Energy', dictref='siesta:E_Fermi', 
     .       fmt='r5', units='siestaUnits:ev')
        call cmlAddProperty(xf=mainXML, value=nk, 
     .       title='Number of k-points', dictRef='siesta:nkpoints',
     .       units='cmlUnits:countable')
        do is = 1,nspin
          call cmlStartPropertyList(mainXML, dictRef='siesta:kpt_band')
          if (nspin.eq.2) then
            if (is.eq.1) then
              call cmlAddProperty(xf=mainXML, value="up", 
     .                            dictRef="siesta:spin")
            else
              call cmlAddProperty(xf=mainXML, value="down", 
     .                            dictRef="siesta:spin")
            endif
          endif
          do ik = 1, nk
            call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .                        weight=kweights(ik))
            call cmlAddProperty(xf=mainXML, value=eo(1:no,is,ik)/eV, 
     .                          dictRef='siesta:eigenenergies',
     .                          units='siestaUnits:ev')
!            call cmlAddBand(xf=mainXML, 
!     .           kpoint=kpoints(:, ik), kweight=kweights(ik), 
!     .           bands=eo(1:no,is,ik))
          enddo
          call cmlEndPropertyList(mainXML)
        enddo
        call cmlEndPropertyList(mainXML)
      endif

      end subroutine ioeig
