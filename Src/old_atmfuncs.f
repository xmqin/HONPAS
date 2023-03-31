! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module old_atmfuncs

C     This module contains a set of routines which provide all the information
C     about the basis set, pseudopotential, atomic mass, etc... of all the
C     chemical species present in the calculation.

!     AG: It is now a "legacy" module to interface the old and new data
!     structures.

      use precision, only: dp
      use sys
      use atmparams, only: nzetmx, lmaxd, nsemx
      use atmparams, only: maxos, nkbmx, ntbmax
      use alloc,     only: re_alloc, de_alloc
      implicit none 

      integer,  save, public     ::  nsmax            

      integer,  public, pointer  ::  izsave(:)          
      integer,  public, pointer  ::  nomax(:)           
      integer,  public, pointer  ::  nkbmax(:)          

      integer,  public, pointer  ::  lmxosave(:)        
      integer,  public, pointer  ::  nkblsave(:,:)      
      integer,  public, pointer  ::  npolorbsave(:,:,:) 
      integer,  public, pointer  ::  nsemicsave(:,:)    
      integer,  public, pointer  ::  nzetasave(:,:,:)   
      integer,  public, pointer  ::  cnfigtb(:,:,:)     

      logical,  public, pointer  ::  semicsave(:)       
     
      real(dp), public, pointer  :: zvaltb(:)
      integer,  public, pointer  :: lmxkbsave(:)
      real(dp), public, pointer  :: smasstb(:)
      real(dp), public, pointer  :: chargesave(:)
      real(dp), public, pointer  :: slfe(:)
      real(dp), public, pointer  :: lambdatb(:,:,:,:)
      real(dp), public, pointer  :: filtercuttb(:,:,:)

      real(dp), public, pointer  :: qtb(:,:)

      real(dp), public, pointer  :: table(:,:,:)
      real(dp), public, pointer  :: tabpol(:,:,:)
      real(dp), public, pointer  :: tab2(:,:,:)
      real(dp), public, pointer  :: tab2pol(:,:,:)


      real(dp), public, pointer  ::  coretab(:,:,:)
      real(dp), public, pointer  ::  chloctab(:,:,:)
      real(dp), public, pointer  ::  vlocaltab(:,:,:)
      real(dp), public, pointer  ::  rctb(:,:,:)
      real(dp), public, pointer  ::  rcotb(:,:,:,:)
      real(dp), public, pointer  ::  rcpoltb(:,:,:,:)

      character(len=20), save, public, pointer :: label_save(:)
      character(len=10), save, public, pointer :: basistype_save(:)  

!     Public routines
      public :: labelfis, izofis, zvalfis, nkbfis
      public :: massfis, lomaxfis, nofis, lmxkbfis
      public :: cnfigfio, lofio, mofio
      public :: atmpopfio , epskb, rcut
      public :: clear_tables, allocate_old_arrays
      public :: deallocate_old_arrays


      PRIVATE

      CONTAINS !================================================
!
      subroutine allocate_old_arrays()

      !allocate(rcotb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( rcotb )
      call re_alloc( rcotb, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'rcotb', 'old_atmfuncs' )

      !allocate(rcpoltb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( rcpoltb )
      call re_alloc( rcpoltb, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'rcpoltb', 'old_atmfuncs' )
      !allocate(lambdatb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( lambdatb )
      call re_alloc( lambdatb, 1, nzetmx, 0, lmaxd, 1, nsemx,
     &               1, nsmax, 'lambdatb', 'old_atmfuncs' )
      !allocate(filtercuttb(0:lmaxd,nsemx,nsmax))
      nullify( filtercuttb )
      call re_alloc( filtercuttb, 0, lmaxd, 1, nsemx,
     &               1, nsmax, 'filtercuttb', 'old_atmfuncs' )
      !allocate(qtb(maxos,nsmax))
      nullify( qtb )
      call re_alloc( qtb, 1, maxos, 1, nsmax,
     &               'qtb', 'old_atmfuncs' )
      !allocate(slfe(nsmax))
      nullify( slfe )
      call re_alloc( slfe, 1, nsmax, 'slfe', 'old_atmfuncs' )
      !allocate(rctb(nkbmx,0:lmaxd,nsmax))
      nullify( rctb )
      call re_alloc( rctb, 1, nkbmx, 0, lmaxd, 1, nsmax,
     &               'rctb', 'old_atmfuncs' )

      !allocate(smasstb(nsmax))
      nullify( smasstb )
      call re_alloc( smasstb, 1, nsmax, 'smasstb', 'old_atmfuncs' )
      !allocate(chargesave(nsmax))
      nullify( chargesave )
      call re_alloc( chargesave, 1, nsmax,
     &               'chargesave', 'old_atmfuncs' )
!
!     Table: This is a hybrid
!          negative values of the second index correspond to KB projectors
!          positive values of the second index correspond to orbitals
!          A second index of zero (0) corresponds to the local potential
! 
!          The first index has ntbmax "real points" and two extra
!          entries for bookeeping
!          The total number of angular momentum entries is lmaxd+1 (since
!          s=0, p=1, etc)
!
!
!      allocate(table((ntbmax+2),-nkbmx*(lmaxd+1):nzetmx*nsemx*
!               (lmaxd+1),nsmax))
      nullify( table )
      call re_alloc( table, 1, ntbmax+2,
     &               -nkbmx*(lmaxd+1), nzetmx*nsemx*(lmaxd+1),
     &               1, nsmax, 'table', 'old_atmfuncs' )
!      allocate(tab2(ntbmax,-nkbmx*(lmaxd+1):nzetmx*nsemx*
!        (lmaxd+1),nsmax))
      nullify( tab2 )
      call re_alloc( tab2, 1, ntbmax,
     &               -nkbmx*(lmaxd+1), nzetmx*nsemx*(lmaxd+1),
     &               1, nsmax, 'tab2', 'old_atmfuncs' )
!      allocate(tabpol((ntbmax+2),nzetmx*nsemx*(lmaxd+1),nsmax))
      nullify( tabpol )
      call re_alloc( tabpol, 1, ntbmax+2,
     &               1, nzetmx*nsemx*(lmaxd+1),
     &               1, nsmax, 'tabpol', 'old_atmfuncs' )
!      allocate(tab2pol(ntbmax,nzetmx*nsemx*(lmaxd+1),nsmax))
      nullify( tab2pol )
      call re_alloc( tab2pol, 1, ntbmax,
     &               1, nzetmx*nsemx*(lmaxd+1),
     &               1, nsmax, 'tab2pol', 'old_atmfuncs' )

!      allocate(coretab(ntbmax+1,2,nsmax))
      nullify( coretab )
      call re_alloc( coretab, 1, ntbmax+1, 1, 2, 1, nsmax,
     &               'coretab', 'old_atmfuncs' )

!      allocate(chloctab((ntbmax+1),2,nsmax))
      nullify( chloctab )
      call re_alloc( chloctab, 1, ntbmax+1, 1, 2, 1, nsmax,
     &               'chloctab', 'old_atmfuncs' )
!      allocate(vlocaltab((ntbmax+1),2,nsmax))
      nullify( vlocaltab )
      call re_alloc( vlocaltab, 1, ntbmax+1, 1, 2, 1, nsmax,
     &               'vlocaltab', 'old_atmfuncs' )

      !allocate(izsave(nsmax))
      nullify( izsave )
      call re_alloc( izsave, 1, nsmax, 'izsave', 'old_atmfuncs' )
      !allocate(lmxkbsave(nsmax))
      nullify( lmxkbsave )
      call re_alloc( lmxkbsave, 1, nsmax,
     &               'lmxkbsave', 'old_atmfuncs' )
      !allocate(lmxosave(nsmax))
      nullify( lmxosave )
      call re_alloc( lmxosave, 1, nsmax, 'lmxosave', 'old_atmfuncs' )
      !allocate(npolorbsave(0:lmaxd,nsemx,nsmax))
      nullify( npolorbsave )
      call re_alloc( npolorbsave, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'npolorbsave', 'old_atmfuncs' )
      !allocate(nsemicsave(0:lmaxd,nsmax))
      nullify( nsemicsave )
      call re_alloc( nsemicsave, 0, lmaxd, 1, nsmax,
     &               'nsemicsave', 'old_atmfuncs' )
      !allocate(nzetasave(0:lmaxd,nsemx,nsmax))
      nullify( nzetasave )
      call re_alloc( nzetasave, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'nzetasave', 'old_atmfuncs' )
      !allocate(nomax(nsmax))
      nullify( nomax )
      call re_alloc( nomax, 1, nsmax, 'nomax', 'old_atmfuncs' )
      !allocate(nkbmax(nsmax))
      nullify( nkbmax )
      call re_alloc( nkbmax, 1, nsmax, 'nkbmax', 'old_atmfuncs' )
      !allocate(zvaltb(nsmax))
      nullify( zvaltb )
      call re_alloc( zvaltb, 1, nsmax, 'zvaltb', 'old_atmfuncs' )
      !allocate(cnfigtb(0:lmaxd,nsemx,nsmax))
      nullify( cnfigtb )
      call re_alloc( cnfigtb, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'cnfigtb', 'old_atmfuncs' )
      !allocate(nkblsave(0:lmaxd,nsmax))
      nullify( nkblsave )
      call re_alloc( nkblsave, 0, lmaxd, 1, nsmax,
     &               'nkblsave', 'old_atmfuncs' )
!
      nullify (label_save)
      allocate(label_save(nsmax))
!      call re_alloc(label_save,1,nsmax,"label_save",
!     $                routine= "allocate_old_arrays")
      nullify (basistype_save)
      allocate(basistype_save(nsmax))
!      call re_alloc(basistype_save,1,nsmax,"basistype_save",
!     $                routine= "allocate_old_arrays")
      nullify (semicsave)
      call re_alloc( semicsave, 1, nsmax,
     &               'semicsave', 'old_atmfuncs' )

      end subroutine allocate_old_arrays

      subroutine deallocate_old_arrays()

      call de_alloc( rcotb,       'rcotb',       'old_atmfuncs' )
      call de_alloc( rcpoltb,     'rcpoltb',     'old_atmfuncs' )
      call de_alloc( lambdatb,    'lambdatb',    'old_atmfuncs' )
      call de_alloc( filtercuttb, 'filtercuttb', 'old_atmfuncs' )
      call de_alloc( qtb,         'qtb',         'old_atmfuncs' )
      call de_alloc( slfe,        'slfe',        'old_atmfuncs' )
      call de_alloc( rctb,        'rctb',        'old_atmfuncs' )
      call de_alloc( smasstb,     'smasstb',     'old_atmfuncs' )
      call de_alloc( chargesave,  'chargesave',  'old_atmfuncs' )
      call de_alloc( table,       'table',       'old_atmfuncs' )
      call de_alloc( tab2,        'tab2',        'old_atmfuncs' )
      call de_alloc( tabpol,      'tabpol',      'old_atmfuncs' )
      call de_alloc( tab2pol,     'tab2pol',     'old_atmfuncs' )
      call de_alloc( coretab,     'coretab',     'old_atmfuncs' )
      call de_alloc( chloctab,    'chloctab',    'old_atmfuncs' )
      call de_alloc( vlocaltab,   'vlocaltab',   'old_atmfuncs' )
      call de_alloc( izsave,      'izsave',      'old_atmfuncs' )
      call de_alloc( lmxkbsave,   'lmxkbsave',   'old_atmfuncs' )
      call de_alloc( lmxosave,    'lmxosave',    'old_atmfuncs' )
      call de_alloc( npolorbsave, 'npolorbsave', 'old_atmfuncs' )
      call de_alloc( nsemicsave,  'nsemicsave',  'old_atmfuncs' )
      call de_alloc( nzetasave,   'nzetasave',   'old_atmfuncs' )
      call de_alloc( nomax,       'nomax',       'old_atmfuncs' )
      call de_alloc( nkbmax,      'nkbmax',      'old_atmfuncs' )
      call de_alloc( zvaltb,      'zvaltb',      'old_atmfuncs' )
      call de_alloc( cnfigtb,     'cnfigtb',     'old_atmfuncs' )
      call de_alloc( nkblsave,    'nkblsave',    'old_atmfuncs' )
      call de_alloc( semicsave,   'semicsave',   'old_atmfuncs' )
      deallocate( label_save )
!      call de_alloc( label_save, 'label_save', 'old_atmfuncs' )
      deallocate( basistype_save )
!      call de_alloc( basistype_save, 'basistype_save', 'old_atmfuncs' )

      end subroutine deallocate_old_arrays

      subroutine clear_tables()

      integer is

      do is=1,nsmax
        izsave(is)=0
        lmxosave(is)=0
        lmxkbsave(is)=0
        label_save(is)='  '
        nkbmax(is)=0
        nomax(is)=0  
        semicsave(is)=.false.
              
        nsemicsave(:,is) = 0
        nzetasave(:,:,is) = 0
        rcotb(:,:,:,is) = 0.0_dp
        lambdatb(:,:,:,is) = 0.0_dp
        filtercuttb(:,:,is) = 0.0_dp
        rcpoltb(:,:,:,is) = 0.0_dp

        table(:,:,is) = 0.0_dp
        tab2(:,:,is) = 0.0_dp
        tabpol(:,:,is) = 0.0_dp
        tab2pol(:,:,is) = 0.0_dp

        qtb(1:maxos,is)=0.00_dp

      enddo 
      end subroutine clear_tables

      subroutine check_is(name,is)
      character(len=*), intent(in) :: name
      integer, intent(in) :: is

      if ((is.lt.1).or.(is.gt.nsmax)) then 
            write(6,*) trim(name),': THERE ARE NO DATA FOR IS=',IS
            write(6,*) trim(name),': ISMIN= 1, NSMAX= ',nsmax
         call die()
      endif
      end subroutine check_is
!
!
      FUNCTION IZOFIS( IS )
      integer :: izofis ! Atomic number
      integer, intent(in) :: is ! Species index

      call check_is('izofis',is)
      izofis=izsave(is)

      end function izofis
!
      FUNCTION ZVALFIS( IS )
      real(dp) :: zvalfis          ! Valence charge
      integer, intent(in) :: is            ! Species index

      call check_is('zvalfis',is)
 
      zvalfis=zvaltb(is)
      end function zvalfis
!
      FUNCTION LABELFIS (IS)
      character(len=20) ::  labelfis  ! Atomic label
      integer, intent(in) :: is            ! Species index

      call check_is('labelfis',is)
      labelfis=label_save(is)
      end function labelfis
!
      FUNCTION LMXKBFIS (IS)
      integer :: lmxkbfis    ! Maximum ang mom of the KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('lmxkbfis',is)
      lmxkbfis=lmxkbsave(is)
      end function lmxkbfis
!
      FUNCTION LOMAXFIS (IS)
      integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
      integer, intent(in) :: is            ! Species index

      integer lmx, nsm

      call check_is('lomaxfis',is)

      lomaxfis=0           
      lmx=lmxosave(is)
      do nsm=1,nsemicsave(lmx,is)+1
         if(npolorbsave(lmx,nsm,is).gt.0)   lomaxfis=lmx+1
      enddo     
      
      lomaxfis=max(lomaxfis,lmx)
      end function lomaxfis
!

      FUNCTION MASSFIS(IS)
      real(dp) :: massfis            ! Mass
      integer, intent(in) :: is            ! Species index

      call check_is('massfis',is)
      massfis=smasstb(is)
      end function massfis
!
      FUNCTION NKBFIS(IS)
      integer :: nkbfis    ! Total number of KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('nkbfis',is)
      nkbfis=nkbmax(is)
      end function nkbfis
!

      FUNCTION NOFIS(IS)
      integer :: nofis    ! Total number of Basis functions
      integer, intent(in) :: is            ! Species index

      call check_is('nofis',is)
      nofis=nomax(is)
      end function nofis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AMENoFIS

      FUNCTION ATMPOPFIO (IS,IO)
      real(dp) atmpopfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the population of the atomic basis orbitals in the atomic 
C ground state configuration.

      call check_is('atmpopfio',is)
      if((io.gt.nomax(is)).or.(io.lt.1)) then
            write(6,*) 'ATMPOPFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'ATMPOPFIO: IOMIN= 1', ' IOMAX= ',nomax(is)
         call die
      endif
 
      atmpopfio=qtb(io,is) 
      end function atmpopfio

!
!
      FUNCTION CNFIGFIO(IS,IO)
      integer cnfigfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C   INTEGER CNFIGFIO: Principal quantum number of the shell to which
C                     the orbital belongs

C               (Formerly, for polarization orbitals
C               the quantum number corresponded to the shell which
C               is polarized by the orbital io.
C               This behavior for polarization orbitals is meaningless,
C               and has been replaced.)
      
      integer l, norb, izeta, ipol, nsm

      call check_is('cnfigfio',is)
      if ((io.gt.nomax(is)).or.(io.lt.1)) then
            write(6,*) 'CNFIGFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'CNFIGFIO: IOMIN= 1',
     .           ' IOMAX= ',nomax(is)
         call die()
      endif

      ! This routine assumes that polarization orbitals are the last ones
      ! in the list indexed by 'io'
      
        norb=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,nsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) then
               cnfigfio=cnfigtb(l,nsm,is)
               return
            endif
 5        continue
 8       continue
10      continue

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) then
                 ! Deal properly with polarization orbitals
                 ! Return the highest n for l+1
                 cnfigfio=maxval(cnfigtb(l+1,:,is))
                 return
              endif
15          continue
18        continue
20      continue
           write(6,*) 'CNFIGFIO: ERROR: ORBITAL INDEX IO=',IO
           write(6,*) 'CNFIGFIO: ERROR: NOT FOUND'
        call die()

      end function cnfigfio
!
!
      FUNCTION LOFIO (IS,IO)
      integer lofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns total angular momentum quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, nsm

      call check_is('lofio',is)
      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then 
            write(6,*) 'LOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'LOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif
 
       if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
          write(6,*) 'LOFIO: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'LOFIO: ERROR: NOT FOUND'
        call die

30      lofio=l
        return

40      lofio=l+1
        return

       elseif(io.lt.0) then

        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue 

60      lofio=l       

c       elseif (io.eq.0) then
        else

        lofio=0

        endif
      end  function lofio
!
!
      FUNCTION MOFIO (IS,IO)
      integer mofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns magnetic quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number M of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, lorb, lkb, nsm

      call check_is('mofio',is)
      if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
            write(6,*) 'MOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'MOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif

      if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue 

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        write(6,*) 'MOFIO: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'MOFIO: ERROR: NOT FOUND'
        call die()

30      lorb=l 
        mofio=io-norb+lorb
        return

40      lorb=l+1 
        mofio=io-norb+lorb
        return


       elseif(io.lt.0) then


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue

60      lkb=l
        mofio=-io+nkb+lkb
c       elseif (io.eq.0) then
        else

        mofio=0

        endif
        
      end function mofio
!

!  End of FIOs 
!

      FUNCTION EPSKB (IS,IO)
      real(dp) epskb
      integer, intent(in)   ::  is   ! Species index
      integer, intent(in)   ::  io   ! KB proyector index (within atom)
                                     ! May be positive or negative 
                                     ! (only ABS(IO) is used).

C  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
C       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C  where Phi_lm is returned by subroutine PHIATM.


C  REAL(DP) EPSKB  : Kleynman-Bylander projector energy
C  Energy in Rydbergs.

      integer  ik, nkb, indx, l, ikb
C
C******************************************************************
C*****************Variables in the common blocks*******************
C
      call check_is('epskb',is)
 
         ik=-abs(io)
         if (ik.eq.0) then 
             write(6,*) 'EPSKB: FUNCTION CANNOT BE CALLED WITH'
     .         ,' ARGUMENT EQUAL TO ZERO' 
           CALL DIE 
         endif

         if(ik.lt.-nkbmax(is)) then
             write(6,*) 'EPSKB: THERE ARE NO DATA FOR IO=',IK
             write(6,*) 'EPSKB: IOMIN= ',-nkbmax(is)
           CALL DIE()
         endif

         nkb=0
         indx=0
         do 10  l=0,lmxkbsave(is)
             do 5 ikb=1,nkblsave(l,is)
                indx=indx+1
                nkb=nkb-(2*l+1)
                if(nkb.le.ik) goto 20
5            continue
10       continue 

20        indx=-indx
 
        epskb=table(2,indx,is)

      end function epskb
!
!
      function rcut(is,io)
      real(dp) rcut
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io> => basis orbitals
                                   ! io<0  => KB projectors
                                   ! io=0 : Local screened pseudopotential

C  Returns cutoff radius of Kleynman-Bylander projectors and
C  atomic basis orbitals.
C  Distances in Bohr

      integer l, norb, lorb, nzetorb, izeta, ipol, nkb,lkb,nsm
      integer  indx, nsmorb        
C
      call check_is('rcut',is)

      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'RCUT: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'RCUT: IOMIN= ',-nkbmax(is),
     .      ' IOMAX= ',nomax(is)
        CALL DIE
      endif


       if (io.gt.0) then

        norb=0 
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,nsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1 
            if(norb.ge.io) goto 30 
 5        continue      
 8       continue 
10      continue   

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)  
              indx=indx+1  
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue     
          write(6,*) 'RCUT: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'RCUT: ERROR: NOT FOUND'
        call die()

30      lorb=l
        nzetorb=izeta
        nsmorb=nsm
        rcut=rcotb(nzetorb,lorb,nsmorb,is)  
c       rcut=table(1,indx,is)*(ntbmax-1)
        return  

40      lorb=l
        nzetorb=ipol
        nsmorb=nsm
        rcut=rcpoltb(nzetorb,lorb,nsmorb,is) 
c       rcut=tabpol(1,indx,is)*(ntbmax-1)
        return 


       elseif(io.lt.0) then 


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 60
45        continue
50      continue 

60      lkb=l 
        indx=-(lkb+1)
        rcut=rctb(izeta,lkb,is) 
c       rcut=table(1,indx,is)*(ntbmax-1)
        return 

c       elseif (io.eq.0) then
        else

        rcut=table(2,0,is)

        endif

      end function rcut
!
!
      end module old_atmfuncs
