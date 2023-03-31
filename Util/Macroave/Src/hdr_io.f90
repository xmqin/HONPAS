! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_io
!! NAME
!! hdr_io
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! The parallel aspects are NOT treated here.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_clean when hdr is not used anymore.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2003 ABINIT group (XG)
!! This file is distributed under the terms of the 
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the file, 
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!  unitfi=unit number of the file (unformatted if rdwr=1 or 2, formatted if rdwr=3)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1 : will be output
!!   if rdwr=2,3,4 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!!
!! PARENTS
!!      conducti,cut3d,initaim,inwffil3,ioarr,newsp,outstat,outwf,rdem1,rdkss
!!      testepsm1,testlda,uderiv,vtorho3,wrem1
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! OUTPUT
!!
!! SOURCE

 subroutine hdr_io(fform,hdr,rdwr,unitfi)

 use defs_basis
 use defs_common
 
 implicit none
 
!Arguments ------------------------------------
 integer,intent(in) :: rdwr,unitfi
 integer :: fform
!This type is defined in defs_common
 type(hdr_type) :: hdr
 
!Local variables-------------------------------
 integer :: bantot,headform,iatom,ikpt,ipsp,isym,lloc,lmax,mmax,natom,nkpt,npsp
 integer :: nsppol,nsym,ntype
 real(dp) :: acell(3)
 character*6 :: codvsn
 character*500 :: message
 
! *************************************************************************

! -------------------------------------------------------------------------
! Reading the header of an unformatted file
! -------------------------------------------------------------------------

 
 if(rdwr==1)then
 
  rewind(unitfi)

! Reading the first record of the file ------------------------------------

  read(unitfi,err=1000)codvsn,fform

  if(fform==1   .or. fform==2   .or. &
&    fform==51  .or. fform==52  .or. &
&    fform==101 .or. fform==102       )then
!  This is the old format
   headform=22

  else

!  Format beyond 22 have a different first line, so need reading again the first line
   rewind (unitfi)
   read (unitfi)   codvsn,headform,fform
   if(headform/=23 .and. headform/=34 .and. headform/=40 )then
    write(message, '(4a,i3,3a,i8,3a)' ) ch10,&
&    ' hdr_io : ERROR -',ch10,&
&    '  The first line of the (WF, DEN or POT) file read in unit ',unitfi,' is erroneous.',ch10,&
&    '  headform is ',headform,', while it should be 23, 34 or 40.',ch10,&
&    '  Action : check the correctness of your file.'
!    call wrtout(6,message,'COLL')
!    call leave_new('COLL')
   endif

  endif

  hdr%codvsn=codvsn
  hdr%headform=headform
! fform is not a record of hdr_type

!DEBUG
! write(6,*)' hdr_io : debug '
! write(6,*)' hdr_io : codvsn,headform,fform',codvsn,headform,fform
!ENDDEBUG

! Reading the second record of the file ------------------------------------

  if(headform==22)then

    read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&    nkpt, nsppol, nsym, ntype,&
&    acell, hdr%ecut_eff, hdr%rprimd
    hdr%nspden=1
    hdr%nspinor=1
    hdr%occopt=1
    npsp=ntype

  else if(headform==23)then

!   Compared to v2.2, add nspden, nspinor, occopt
    read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, ntype, hdr%occopt,&
&    acell, hdr%ecut_eff, hdr%rprimd
    npsp=ntype

  else if(headform==34)then

!   Compared to v2.3, subtract acell, and add npsp
    read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntype, hdr%occopt,&
&    hdr%ecut_eff, hdr%rprimd

  else if(headform==40)then

!   Compared to v3.4, add ecut, ecutsm, tphysel, tsmear
    read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntype, hdr%occopt,&
&    hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

  endif  

  hdr%bantot=bantot
  hdr%natom =natom
  hdr%nkpt  =nkpt
  hdr%npsp  =npsp
  hdr%nsppol=nsppol
  hdr%nsym  =nsym
  hdr%ntype =ntype

!DEBUG
!  write(6,*)' hdr_io : before allocate '
!  write(6,*)' hdr_io : bantot,natom,nkpt,npsp,nsppol,nsym,ntype',&
!&  bantot,natom,nkpt,npsp,nsppol,nsym,ntype
!ENDDEBUG

! Allocate all parts of hdr that need to be --------------------------------

  allocate(hdr%istwfk(nkpt))
  allocate(hdr%nband(nkpt*nsppol))
  allocate(hdr%npwarr(nkpt)) ! Warning : npwarr here has only one dim
  allocate(hdr%pspcod(npsp))
  allocate(hdr%pspdat(npsp))
  allocate(hdr%pspso(npsp))
  allocate(hdr%pspxc(npsp))
  allocate(hdr%so_typat(ntype))
  allocate(hdr%symafm(nsym))
  allocate(hdr%symrel(3,3,nsym))
  allocate(hdr%type(natom))
  allocate(hdr%kptns(3,nkpt))
  allocate(hdr%occ(bantot))
  allocate(hdr%tnons(3,nsym))
  allocate(hdr%xred(3,natom))
  allocate(hdr%znuclpsp(npsp))
  allocate(hdr%znucltypat(ntype))
  allocate(hdr%zionpsp(npsp))
  allocate(hdr%title(npsp))

!DEBUG
! write(6,*)' hdr_io : after allocate '
!ENDDEBUG
  
! Reading the third record of the file ------------------------------------

  if(headform==22 .and. (fform==1 .or. fform==51 .or. fform==101))then

!   This is very old (pre-2.0) format !
    read(unitfi) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
&    hdr%type(:), hdr%kptns(:,:), hdr%occ(:), &
&    hdr%tnons(:,:), hdr%znucltypat(:)

    hdr%istwfk(:)=1

  else if(headform==22 .or. headform==23 .or. headform==34)then

!   Compared to pre v2.0, add istwfk
    read(unitfi) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
&    hdr%type(:), hdr%istwfk(:), hdr%kptns(:,:), hdr%occ(:), &
&    hdr%tnons(:,:), hdr%znucltypat(:)

  else if(headform==40)then

!   Compared to 2.3, add so_typat and symafm, and switch istwfk
    read(unitfi)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
&    hdr%so_typat(:), hdr%symafm(:), hdr%symrel(:,:,:), &
&    hdr%type(:), hdr%kptns(:,:), hdr%occ(:), &
&    hdr%tnons(:,:), hdr%znucltypat(:)

  endif

! Reading the records with psp information ---------------------------------

  if(headform==22)then

   do ipsp=1,npsp
    read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&    hdr%zionpsp(ipsp), hdr%pspdat(ipsp), hdr%pspcod(ipsp), &
&    hdr%pspxc(ipsp), lmax, lloc, mmax
    hdr%pspso(ipsp)=1
   enddo

  else if(headform==23)then

!  Compared to 2.2, add pspso
   do ipsp=1,npsp
    read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&    hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&    hdr%pspcod(ipsp), hdr%pspxc(ipsp), lmax, lloc, mmax
   enddo

  else if(headform==34 .or. headform==40)then

!  Compared to 2.3, suppress lmax, lloc, mmax
   do ipsp=1,npsp
    read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&    hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&    hdr%pspcod(ipsp), hdr%pspxc(ipsp)
   enddo

  endif

! Reading the final record of the header  ---------------------------------
 
  if(headform==22)then
   read(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot
   hdr%fermie=zero
  else if(headform==23 .or. headform==34 .or. headform==40)then
   read(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie
  endif

!DEBUG
!  write(6,*)' hdr_io : read mode, hdr%so_typat(:), hdr%symafm(:)=',&
!&                                 hdr%so_typat(:), hdr%symafm(:)
!ENDDEBUG


! -------------------------------------------------------------------------
! Writing the header of an unformatted file
! -------------------------------------------------------------------------

 else if(rdwr==2)then

! natom,nkpt,npsp,ntype... are not defined in this section : 
! always address then from hdr

  rewind(unitfi)

! Writing always use format 40
  headform=40
  write(unitfi) hdr%codvsn, headform, fform

  write(unitfi) hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, &
&  hdr%natom, hdr%ngfft(1:3), hdr%nkpt, &
&  hdr%nspden, hdr%nspinor, &
&  hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntype, hdr%occopt,&
&  hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

  write(unitfi) hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:),&
&   hdr%so_typat(:), hdr%symafm(:), hdr%symrel(:,:,:), &
&   hdr%type(:), hdr%kptns(:,:), hdr%occ(:), &
&   hdr%tnons(:,:), hdr%znucltypat(:)

!DEBUG
!   write(6,*)' hdr_io : write psp record, headform= ',hdr%headform
!ENDDEBUG

  do ipsp=1,hdr%npsp

   write(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&   hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&   hdr%pspcod(ipsp), hdr%pspxc(ipsp)

  enddo

!DEBUG
!  write(6,*)' hdr_io : write mode, hdr%so_typat(:), hdr%symafm(:)=',&
!&                                  hdr%so_typat(:), hdr%symafm(:)
!ENDDEBUG

  write(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie

! -------------------------------------------------------------------------
! Writing the header of an unformatted file
! -------------------------------------------------------------------------

 else if(rdwr==3 .or. rdwr==4)then

  write(unitfi, '(a)' )&
&  ' ==============================================================================='
  if(rdwr==3)write(unitfi, '(a)' ) ' ECHO of part of the ABINIT file header '
  if(rdwr==4)write(unitfi, '(a)' ) ' ECHO of the ABINIT file header '
  write(unitfi, '(a)' ) ' '
  write(unitfi, '(a)' ) ' First record :'
  write(unitfi, '(a,a6,2i5)' )  ' codvsn,headform,fform = ',&
&  hdr%codvsn, hdr%headform, fform     ! Do not worry about 22 format

  write(unitfi, '(a)' ) ' '
  write(unitfi, '(a)' ) ' Second record :'
  write(unitfi, '(a,4i6)') ' bantot,intxc,ixc,natom  =',&
&                            hdr%bantot, hdr%intxc, hdr%ixc, hdr%natom
  write(unitfi, '(a,4i6)') ' ngfft(1:3),nkpt         =',&
&                            hdr%ngfft(1:3), hdr%nkpt

  if(hdr%headform>=23)then
   write(unitfi, '(a,4i6)') ' nspden,nspinor          =',&
&                            hdr%nspden, hdr%nspinor
  endif

  if(hdr%headform<=23)then
   write(unitfi, '(a,4i6)' ) ' nsppol,nsym,ntype,occopt=',&
&                            hdr%nsppol,hdr%nsym,hdr%ntype,hdr%occopt
  else
   write(unitfi, '(a,5i6)' ) ' nsppol,nsym,npsp,ntype,occopt=',&
&                            hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntype,hdr%occopt
  endif

  if(hdr%headform>=40)then
   write(unitfi, '(a,2es18.10)') ' ecut,ecutsm             =',hdr%ecut, hdr%ecutsm
  endif

  write(unitfi, '(a, es18.10)' ) ' ecut_eff                =',hdr%ecut_eff
  write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,1)           =',hdr%rprimd(1:3,1)
  write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,2)           =',hdr%rprimd(1:3,2)
  write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,3)           =',hdr%rprimd(1:3,3)

  if(hdr%headform>=40)then
   write(unitfi, '(a,2es18.10)') ' tphysel,tsmear          =',hdr%tphysel, hdr%tsmear
  endif

  write(unitfi, '(a)' )
  if(rdwr==3)then
   write(unitfi, '(a,i3,a)' ) ' The header contain ',hdr%npsp+2,' additional records.'
  else

   write(unitfi, '(a)' ) ' Third record :'
   write(unitfi, '(a,(12i4,8x))') ' istwfk=',hdr%istwfk(:)
   write(unitfi, '(a,(12i4,8x))') ' nband =',hdr%nband(:)
   write(unitfi, '(a,(10i5,8x))') ' npwarr=',hdr%npwarr(:)

   if(hdr%headform>=40)then
    write(unitfi, '(a,(12i4,8x))') ' so_typat=',hdr%so_typat(:)
   endif

   if(hdr%headform>=40)then
    write(unitfi, '(a)') ' symafm='
    write(unitfi, '(8x,24i3,8x)') hdr%symafm(:)
   endif

   write(unitfi, '(a)' ) ' symrel='
   do isym=1,hdr%nsym/2
    write(unitfi, '(a,9i4,a,9i4)' ) '        ',hdr%symrel(:,:,2*isym-1),&
&                                         '  ',hdr%symrel(:,:,2*isym)
   enddo
   if(2*(hdr%nsym/2)/=hdr%nsym)write(unitfi, '(a,9i4)' ) '        ',hdr%symrel(:,:,hdr%nsym)

   write(unitfi, '(a,(12i4,8x))') ' type  =',hdr%type(:)
   write(unitfi, '(a)' ) ' kptns =                 (max 50 k-points will be written)'
   do ikpt=1,min(hdr%nkpt,50)
    write(unitfi, '(a,3es16.6)' ) '        ',hdr%kptns(:,ikpt)
   enddo
   write(unitfi, '(a,(10f6.2,8x))') ' occ   =',hdr%occ(:)
   write(unitfi, '(a)' ) ' tnons ='
   do isym=1,hdr%nsym/2
    write(unitfi, '(a,3f10.6,a,3f10.6)' ) '        ',hdr%tnons(:,2*isym-1),&
&                                               '  ',hdr%tnons(:,2*isym)
   enddo
   if(2*(hdr%nsym/2)/=hdr%nsym)write(unitfi, '(a,3f10.6)' ) '        ',hdr%tnons(:,hdr%nsym)
   write(unitfi, '(a,(10f6.2,8x))') '  znucl=',hdr%znucltypat(:)
   write(unitfi,'(a)')

   write(unitfi, '(a)' ) ' Pseudopotential info :'
   do ipsp=1,hdr%npsp
    write(unitfi,'(a,a)' ) ' title=',trim(hdr%title(ipsp))  
    write(unitfi,'(a,f6.2,a,f6.2,a,i3,a,i6,a,i3,a,i3)' ) &
&    '  znuclpsp=',hdr%znuclpsp(ipsp),    ', zionpsp=',  hdr%zionpsp(ipsp),    &
&    ', pspso=' , hdr%pspso(ipsp),  ', pspdat=',hdr%pspdat(ipsp),  &
&    ', pspcod=', hdr%pspcod(ipsp), ', pspxc=', hdr%pspxc(ipsp)
   enddo

   write(unitfi, '(a)' ) ' '
   write(unitfi, '(a)' ) ' Last record :'
   write(unitfi, '(a,es16.6,es22.12,es16.6)' ) &
&   ' residm,etot,fermie=',hdr%residm, hdr%etot, hdr%fermie
   write(unitfi, '(a)' ) ' xred ='
   do iatom=1,hdr%natom
    write(unitfi, '(a,3es16.6)' ) '        ',hdr%xred(:,iatom)
   enddo
  if(rdwr==3)write(unitfi, '(a)' ) ' End the ECHO of part of the ABINIT file header '
  if(rdwr==4)write(unitfi, '(a)' ) ' End the ECHO of the ABINIT file header '
  write(unitfi, '(a)' )&
&  ' ==============================================================================='

  endif ! rdwr is 3 or 4

 endif ! choice read/write/echo

 return
 1000 fform=0 ; return   ! This is to allow treatment of old epsm1 format

 end subroutine
!!***
