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
      subroutine pdos( NO, NSPIN, MAXSPN, NO_L, MAXNH, 
     .                 MAXO, NUMH, LISTHPTR, LISTH, H, S, 
     .                 E1, E2, SIGMA, NHIST,
     .                 GAMMA, XIJ, INDXUO, NK, KPOINT, WK, EO, NO_U )
C **********************************************************************
C Subroutine to calculate the proyected density of states on the
C atomic orbitals for a given eigenvalue spectra
C Written by J. Junquera and E. Artacho, November 1999.
C ***********  INPUT  **************************************************
C INTEGER NO                  : Number of basis orbitals in the supercell
C INTEGER NSPIN               : Spin polarizations (1 or 2)
C INTEGER MAXSPN              : Second dimension of eo and qo
C                               (Max number of different spin polarizations)
C INTEGER NO_L               : Maximum number of atomic orbitals in the unit 
C                               cell. First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
!                               IN THIS PROCESSOR
C INTEGER MAXNH               : Maximum number of orbitals interacting 
C                               with any orbital
C INTEGER MAXO                : First dimension of eo 
C INTEGER NUMH(NUO)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C INTEGER LISTH(MAXNH)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTHPTR(NUO)       : Pointer to each row (-1) of the
C                               density matrix
C REAL*8  H(MAXNH,NSPIN)      : Hamiltonian in sparse format
C REAL*8  S(MAXNH)            : Overlap in sparse format
C REAL*8  E1, E2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C REAL*8  SIGMA               : Width of the gaussian to expand the eigenvalues
C INTEGER NHIST               : Number of subdivisions of the histogram
C LOGICAL GAMMA               : Only gamma point?
C REAL*8  XIJ(3,MAXNH)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C INTEGER INDXUO(NO)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the nuber of orbitals in the unit cell
C INTEGER NK                  : Number of k points
C REAL*8  KPOINT(3,NK)        : k point vectors
C REAL*8  WK(NK)              : k point weights (must sum one)
C REAL*8  EO(NO_L,MAXSPN,NK) : Eigenvalues
C INTEGER NO_U              : Total number of orbitals in unit cell
C **********************************************************************

      use precision,    only : dp
      use parallel,     only : Node, Nodes, IOnode
      use fdf
      use siesta_geom,  only : xa, isa
      use atomlist,     only : iphorb, iaorb
      use atmfuncs,     only : zetafio, mofio, lofio, cnfigfio, labelfis
      use FoX_wxml,     only : xmlf_t, xml_OpenFile, xml_NewElement, 
     .                         xml_AddCharacters, xml_AddAttribute,
     .                         xml_EndElement, xml_Close
      use xml,          only : str, xml_dump_attribute
      use densematrix,  only : Haux, Saux, psi
      use alloc,        only : re_alloc, de_alloc
      use units,        only : eV
      use files,        only : slabel, label_length
#ifdef MPI
      use parallelsubs, only : GetNodeOrbs
      use mpi_siesta,   only : MPI_logical, MPI_Comm_World
#endif

      implicit none

      integer
     .  NO, NSPIN, MAXSPN, NO_L, MAXNH, NK, NHIST, MAXO, NO_U

      integer 
     .  NUMH(*), LISTH(MAXNH), LISTHPTR(*), INDXUO(NO)

      real(dp)
     .  H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA,  
     .  XIJ(3,MAXNH), KPOINT(3,NK), WK(NK), EO(MAXO,MAXSPN,NK)

      logical
     .  gamma

C Dynamic arrays -------------------------------------------------------
      real(dp), dimension(:,:)  , pointer :: DTOT
      real(dp), dimension(:,:,:), pointer :: DPR

C Internal variables ---------------------------------------------------
      integer
     .  nuo, nhs, npsi, iuo, ihist, ispin, 
     .  iunit1, iunit2, i

      integer iat, spec, ii, iorb

      logical, save :: ParallelOverK = .false.

      real(dp), dimension(:), pointer  :: tmp

      character*40 pos_string

      character(len=label_length+5) :: fnamepro
      character(len=label_length+4) :: fnametot

      real(dp) delta, ener

      external
     .  io_assign, io_close,
     .  pdosk, timer

#ifdef MPI
      integer :: mpierror
#endif

      type(xmlf_t) :: xf            ! For new XML output

      call timer( 'pdos', 1)

C Find the intervals between the subdivisions in the energy scale ------
      delta = (E2 - E1) / (NHIST-1)
      if (Node.eq.0) then
        ParallelOverK = fdf_boolean('Diag.ParallelOverK',.false.)
        if (Nodes.eq.1) ParallelOverK = .false.
      endif
#ifdef MPI
      call MPI_Bcast(ParallelOverK, 1, MPI_logical, 0,
     .               MPI_Comm_World, MPIerror)
#endif

C Find number of orbitals per unit cell 
#ifdef MPI
      call GetNodeOrbs(no_u,Node,Nodes,nuo)
#else
      nuo = no_u
#endif

C Check internal dimensions --------------------------------------------
       if ( nspin.le.2 .and. gamma) then
         nhs  = no_u * nuo
         npsi = no_u * no_l * nspin
       elseif ( nspin.le.2 .and. .not.gamma) then
         if (ParallelOverK) then
           nhs  = 2 * no_u * no_u
           npsi = 2 * no_u * no_u
         else
           nhs  = 2 * no_u * nuo
           npsi = 2 * no_u * nuo
         endif
       endif

C Allocate local arrays ------------------------------------------------
      call re_alloc(Haux,1,nhs,name='Haux',routine='pdos')
      call re_alloc(Saux,1,nhs,name='Saux',routine='pdos')
      call re_alloc(psi,1,npsi,name='psi',routine='pdos')
      nullify( dtot )
      call re_alloc( dtot, 1, nhist, 1, 2, name='dtot',
     &               routine='pdos' )
      nullify( dpr )
      call re_alloc( dpr, 1, nhist, 1, no_u, 1, 2, name='dpr',
     &               routine='pdos' )

C Initialize the projected density of states ---------------------------
      do ispin = 1, 2
        do ihist = 1,nhist
          dtot(ihist,ispin) = 0.d0
          do iuo = 1,no_u
            dpr(ihist,iuo,ispin) = 0.d0
          enddo
        enddo
      enddo

C  Call appropiate routine ----------------------------------------------
      if (nspin.le.2 .and. gamma) then
        call pdosg( NSPIN, NUO, NO, MAXSPN, NO_L, MAXNH,
     .              MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .              E1, E2, NHIST, SIGMA, INDXUO, EO,
     .              HAUX, SAUX, PSI, DTOT, DPR, NO_U )
      elseif ( nspin.le.2 .and. .not.gamma) then
        if (ParallelOverK) then
          call pdoskp( NSPIN, NUO, NO, MAXSPN, NO_L, MAXNH,
     .                MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .                E1, E2, NHIST, SIGMA, 
     .                XIJ, INDXUO, NK, KPOINT, WK, EO,
     .                HAUX, SAUX, PSI, DTOT, DPR, NO_U )
        else
          call pdosk( NSPIN, NUO, NO, MAXSPN, NO_L, MAXNH,
     .                MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .                E1, E2, NHIST, SIGMA, 
     .                XIJ, INDXUO, NK, KPOINT, WK, EO,
     .                HAUX, SAUX, PSI, DTOT, DPR, NO_U )
        endif
      endif

      if (IOnode) then
C Open file for write on I/O node
        fnametot = trim(slabel)//'.DOS'
        call io_assign(iunit1)
        open(unit=iunit1, file=fnametot, form='formatted', 
     .       status='unknown') 

C Output histogram
         do ihist = 1,nhist
           ENER = E1 + (ihist-1) * delta
           write(iunit1,'(3f20.5)') ener/ev,dtot(ihist,1)*ev,
     .         dtot(ihist,2)*eV
         enddo

C Close file for write
         call io_close(iunit1)
      endif

C New writing
      if (IOnode) then
         call xml_OpenFile("pdos.xml",xf)

        fnamepro = trim(slabel)//'.PDOS'
        call io_assign(iunit2)
        open(iunit2,file=fnamepro,form='formatted',status='unknown')
        call xml_NewElement(xf,"pdos")
        call xml_NewElement(xf,"nspin")
        call xml_AddCharacters(xf,nspin)
        call xml_EndElement(xf,"nspin")

        write(iunit2,'(a)') '<pdos>'
        write(iunit2,'(a,i1,a)') '<nspin>', nspin, '</nspin>'
        write(iunit2,'(a,i4,a)') '<norbitals>', no_u, '</norbitals>'
        call xml_NewElement(xf,"energy_values")
        call xml_AddAttribute(xf,"units","eV")
        write(iunit2,'(a)') '<energy_values units="eV">'

        nullify( tmp )
        call re_alloc( tmp, 1, 2*nhist, name='tmp',
     &                 routine='pdos' )
        do ihist = 1,nhist
           ENER = E1 + (ihist-1) * delta
           tmp(ihist) = ener/eV
           write(iunit2,'(f20.5)') ener/eV
        enddo
        write(iunit2,'(a)') '</energy_values>'

        call xml_AddCharacters(xf,tmp(1:nhist))
        call xml_EndElement(xf,"energy_values")

        do i = 1, no_u
          iat = iaorb(i)
          iorb = iphorb(i)
          spec = isa(iat)

          call xml_NewElement(xf,"orbital")
          write(iunit2,'(a)') '<orbital '
          call xml_AddAttribute(xf,"index",i)
          call xml_dump_attribute(iunit2,"index",str(i))
          call xml_AddAttribute(xf,"atom_index",iat)
          call xml_dump_attribute(iunit2,"atom_index",str(iat))
          call xml_AddAttribute(xf,"species",trim(labelfis(spec)))
          call xml_dump_attribute(iunit2,"species",
     .                           trim(labelfis(spec)))
          write(pos_string,'(3f11.6)') (xa(ii,iat),ii=1,3)
          call xml_AddAttribute(xf,"position",pos_string)
          call xml_dump_attribute(iunit2,"position",pos_string)
          call xml_AddAttribute(xf,"n",cnfigfio(spec,iorb))
          call xml_dump_attribute(iunit2,"n",str(cnfigfio(spec,iorb)))
          call xml_AddAttribute(xf,"l",lofio(spec,iorb))
          call xml_dump_attribute(iunit2,"l",str(lofio(spec,iorb)))
          call xml_AddAttribute(xf,"m",mofio(spec,iorb))
          call xml_dump_attribute(iunit2,"m",str(mofio(spec,iorb)))
          call xml_AddAttribute(xf,"z",zetafio(spec,iorb))
          call xml_dump_attribute(iunit2,"z",str(zetafio(spec,iorb)))
          write(iunit2,'(a)') '> '

!------------------------------------------------------------
          call xml_NewElement(xf,"data")
          if (nspin.eq.1) then
             call xml_AddCharacters(xf,dpr(1:nhist,i,1)*eV)
          elseif (nspin .eq. 2) then
             tmp(1:2*nhist-1:2) = dpr(1:nhist,i,1)*eV
             tmp(2:2*nhist:2) = dpr(1:nhist,i,2)*eV
             call xml_AddCharacters(xf,tmp(1:2*nhist),fmt="r6")
          endif
          call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"orbital")
!------------------------------------------------------------

          write(iunit2,'(a)') '<data>'
          do ihist = 1,nhist
            if (nspin.eq.1) then
              write(iunit2,'(f10.5)') dpr(ihist,i,1)*eV
            elseif (nspin .eq. 2) then
              write(iunit2,'(2f10.5)') dpr(ihist,i,1)*eV,
     .                                 dpr(ihist,i,2)*eV
            endif
          enddo
          write(iunit2,'(a)') '</data>'
          write(iunit2,'(a)') '</orbital>'
        enddo

C Close file
        call xml_EndElement(xf,"pdos")
        call xml_Close(xf)
        write(iunit2,'(a)') '</pdos>'
        call io_close(iunit2)

C Free local workspace array
        call de_alloc( tmp, name='tmp' )

      endif

C Free local arrays ----------------------------------------------------
      call de_alloc( dpr, name='dpr' )
      call de_alloc( dtot, name='dtot' )

      call timer( 'pdos', 2)

      end

