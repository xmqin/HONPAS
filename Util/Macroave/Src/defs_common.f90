! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!{\src2tex{textfont=tt}}
!!****f* ABINIT/defs_common
!! NAME
!! defs_common
!!
!! FUNCTION
!! This module contains definitions of structured datatypes.
!! In the future, it might contain all these definitions.
!! At present, it contains those for the routines in the common directory
!! for the iowfdenpot directory, and all higher directories :
!! - dataset_type : the "dataset" for the main abinit code
!! - dtfil_type : the data related to files
!! - hdr_type   : the header of wf, den and pot files
!! - bandstructure_type : different information about the band structure
!! - anaddb_dtset_type : the "dataset" for anaddb
!! - aim_dtset_type : the "dataset" for aim
!!
!! COPYRIGHT
!! Copyright (C) 2001-2003 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! (1) The dataset_type datatype
!! The dataset_type structured datatype gather all the input variables,
!! except those that are labelled NOT INTERNAL.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset.
!! The "evolving" input variables are also stored, with their
!! name appended with _orig, to make clear that this is the original
!! value, decided by the user, and not a possibly modified, intermediate value.
!! The following input variables are NOT INTERNAL, that is, they
!! are input variables used to determine other input variables,
!! after suitable processing, and do not appear anymore afterwards 
!! (so, they do not appear as components of a dataset_type variable) :
!! cpuh,cpum(but cpus is present),fband,kptbounds,ndivk,nobj,
!! objaat,objbat,objaax,objbax,objan,objbn,objarf,objbrf,objaro,objbro
!! objatr,objbtr,vaclst,vacuum 
!!
!! (2) The datafil_type datatype
!! The datafiles_type structures datatype gather all the variables related
!! to files, such as filename, and file units.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset. 
!!
!! (3) The hdr_type datatype
!! It contains all the information needed to write a header for a 
!! wf, den or pot file.
!! The structure of the header is explained in the abinis_help.htm file.
!! The datatype is considered as an object, to which are attached a whole
!! set of "methods", actually, different subroutines.
!! A few of these subroutines are : hdr_init, hdr_update, hdr_clean,
!! hdr_check, hdr_io, hdr_skip.
!!
!! (4) The bz_type datatype
!! It contains different information about the Brillouin zone to be
!! considered, according to the context : k points, occupation numbers,
!! storage mode of wavefunctions, weights ...
!! For example, the initial Brillouin zone, set up in the dataset, will be treated
!! in the response function part of the code, to give a reduced
!! Brillouin zone different from the original one, due to the
!! breaking of the symmetries related to the existence of a wavevector,
!! or the lack of time-reversal invariance
!!
!! (5) The anaddb_dataset_type datatype
!! The anaddb_dataset_type structured datatype (will in the future)
!! gather all the input variables for the anaddb code.
!!
!! (6) The aim_dataset_type datatype
!! The aim_dataset_type structured datatype 
!! gathers all the input variables for the anaddb code.
!!
!! TODO
!!
!! SOURCE

 module defs_common

 use defs_basis
 
 implicit none

!Structures

!----------------------------------------------------------------------

 type dataset_type
! Since all these input variables are described in the abinis_help.htm
! file, they are not described in length here ...
! Integer
  integer :: berryopt,brvltt,ceksph,chkexit,chkprim,delayperm,enunit,&
&  getcell,getddk,geteps,getden,getkss,getocc,getvel,getwfk,&
&  getwfkden,getwfq,&
&  getxcart,getxred,get1den,get1wf,get1wfden,&
&  ikhxc,intxc,ionmov,iprcch,iprcel,iprcfc,irdwfk,irdwfq,ird1wf,irdddk,&
&  iscf,isecur,istatr,istatshft,ixc,& 
!  jdtset contains the actual number of the dataset 
&  jdtset,kpara,&   
&  kptopt,kssform,localrdwf,&
&  mband,mffmem,mgfft,mkmem,mkqmem,mk1mem,mpw,mqgrid,&
&  natom,natrd,nbdblock,nbdbuf,&
&  nberry,nbndsto,ncomsto,nconeq,ndtset,nfft,nfreqsus,ngwpt,nkpt,nline,&
&  nnsclo,npack,npara,npsp,npspalch,npweps,npwmat,npwwfn,nqpt,&
&  nsheps,nshiftk,nshmat,nshwfn,nspden,&
&  nspinor,nsppol,nstep,nsym,ntime,&
&  ntypalch,ntype,ntyppure,occopt,optcell,optdriver,ortalg,parareel,paw,pawgratio,&
&  prtbbb,prtcml,prtden,prtdos,prtgeo,prtkpt,&
&  prtpot,prtvha,prtvhxc,prtvol,prtvxc,&
&  prtwf,prt1dm,ptgroupma,&
&  restartxf,rfasr,rfelfd,rfmeth,rfphon,rfstrs,rfthrd,&
&  rfuser,rf1elfd,rf1phon,rf2elfd,rf2phon,rf3elfd,rf3phon,&
&  signperm,spgaxor,spgorig,spgroup,td_mexcit,&
&  timopt,useylm,useria,&
&  userib,useric,userid,userie,vacnum,wfoptalg 
! Integer arrays
  integer :: bdberry(4),dsifkpt(3),kptrlatt(3,3),ngfft(8),nloalg(5),&
&  rfatpol(2),rfdir(3),rf1atpol(2),rf1dir(3),&
&  rf2atpol(2),rf2dir(3),rf3atpol(2),rf3dir(3)
! Integer pointers
  integer, pointer ::  algalch(:)    ! algalch(ntypalch)
  integer, pointer ::  bdgw(:,:)     ! bdgw(2,ngwpt)
  integer, pointer ::  iatfix(:,:)   ! iatfix(3,natom)
  integer, pointer ::  istwfk(:)     ! istwfk(nkpt)
  integer, pointer ::  kberry(:,:)   ! kberry(3,nberry)
  integer, pointer ::  nband(:)      ! nband(nkpt*nsppol)
  integer, pointer ::  so_typat(:)   ! so_typat(ntype)
  integer, pointer ::  symafm(:)     ! symafm(nsym)
  integer, pointer ::  symrel(:,:,:) ! symrel(3,3,nsym)
  integer, pointer ::  type(:)       ! type(natom)         
! Real 
  real(dp) :: charge,cpus,dedlnn,diecut,diegap,dielam,&
&  dielng,diemac,diemix,dilatmx,dtion,&
&  ecut,ecuteps,ecutgros,ecutmat,ecutsm,ecutwfn,&
&  eshift,fband,fixmom,freqsusin,freqsuslo,friction,kptnrm,kptrlen,mdftemp,&
&  mditemp,mdwall,nelect,noseinert,plasfrq,qptnrm,sciss,strfact,strprecon,&
&  td_maxene,toldfe,toldff,&
&  tolmxf,tolvrs,tolwfr,tphysel,tsmear,userra,userrb,userrc,userrd,&
&  userre,vacwidth,vis,zcut
! Real arrays
  real(dp) :: acell_orig(3),angdeg_orig(3),boxcenter(3),&
&  genafm(3),qprtrb(3),qpt(3),qptn(3),rprim_orig(3,3),&
&  rprimd_orig(3,3),strtarget(6),vprtrb(2) 
! Real pointers
  real(dp), pointer :: amu(:)         ! amu(ntype)
  real(dp), pointer :: densty(:,:)    ! densty(ntype,4)
  real(dp), pointer :: kpt(:,:)       ! kpt(3,nkpt)
  real(dp), pointer :: kptgw(:,:)     ! kptgw(3,ngwpt)
  real(dp), pointer :: kptns(:,:)     ! kptns(3,nkpt)
  real(dp), pointer :: mixalch(:,:)   ! mixalch(npspalch,ntypalch)
  real(dp), pointer :: occ_orig(:)    ! occ_orig(mband*nkpt*nsppol) 
  real(dp), pointer :: shiftk(:,:)    ! shifk(3,nshiftk)
  real(dp), pointer :: spinat(:,:)    ! spinat(3,natom)
  real(dp), pointer :: tnons(:,:)     ! tnons(3,nsym)
  real(dp), pointer :: vel_orig(:,:)  ! vel_orig(3,natom)
  real(dp), pointer :: wtatcon(:,:,:) ! wtatcon(3,natom,nconeq)
  real(dp), pointer :: wtk(:)         ! wtk(nkpt)
  real(dp), pointer :: xred_orig(:,:) ! xred_orig(3,natom) 
  real(dp), pointer :: ziontypat(:)   ! ziontypat(ntype)
  real(dp), pointer :: znucl(:)       ! znucl(npsp)
 end type dataset_type

!----------------------------------------------------------------------

 type datafiles_type

  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffknm file must be read 
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1nm file must be read 
   !   (if ird1wf non-zero or get1wf non-zero)
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: undot   ! unit number for ddk 1WF file
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: ungsc   ! unit number for <g|S|c> data (Paw only)

  character*(fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input :
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character*(fnlen) :: fildensin  
   ! if no dataset mode             : abi//'DEN' 
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN' 
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character*(fnlen) :: filkss  
   ! if no dataset mode             : abi//'KSS' 
   ! if dataset mode, and getkss==0 : abi//'_DS'//trim(jdtset)//'KSS' 
   ! if dataset mode, and getkss/=0 : abo//'_DS'//trim(jgetkss)//'KSS'

  character*(fnlen) :: filem1  
   ! if no dataset mode             : abi//'EM1' 
   ! if dataset mode, and geteps==0 : abi//'_DS'//trim(jdtset)//'EM1' 
   ! if dataset mode, and geteps/=0 : abo//'_DS'//trim(jgeteps)//'EM1'

! character*(fnlen) :: filpsp(ntype)
   ! the filenames of the pseudopotential files, from the standard input.

  character*(fnlen) :: filstat
   ! tmp//'_STATUS'

  character*(fnlen) :: wffknm
   ! the name of the ground-state wavefunction file to be read (see driver.f)

  character*(fnlen) :: wffqnm
   ! the name of the k+q ground-state wavefunction file to be read (see driver.f)
   ! only useful in the response-function case

  character*(fnlen) :: wffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.f)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character*(fnlen) :: wff1nm
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.f)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character*(fnlen) :: fildens1in   ! to be described by MVeithen
  character*(fnlen) :: wffknmden   ! to be described by MVeithen
  character*(fnlen) :: wff1nmden   ! to be described by MVeithen

 end type datafiles_type

!----------------------------------------------------------------------

 type hdr_type
  integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
  integer :: date          ! starting date
  integer :: headform      ! format of the header
  integer :: intxc,ixc,natom,nkpt,npsp,nspden        ! input variables
  integer :: nspinor,nsppol,nsym,ntype,occopt        ! input variables
  integer :: ngfft(3)      ! input variable

! This record is not a part of the hdr_type, although it is present in the
! header of the files. This is because it depends on the kind of file
! that is written, while all other information does not depend on it.
! It was preferred to let it be initialized or defined outside of hdr_type.
! integer :: fform         ! file descriptor (or file format)

  integer, pointer :: istwfk(:)    ! input variable istwfk(nkpt)
  integer, pointer :: nband(:)     ! input variable nband(nkpt*nsppol)
  integer, pointer :: npwarr(:)    ! npwarr(nkpt) array holding npw for each k point
  integer, pointer :: pspcod(:)    ! pscod(npsp) from psps
  integer, pointer :: pspdat(:)    ! psdat(npsp) from psps
  integer, pointer :: pspso(:)     ! pspso(npsp) from psps
  integer, pointer :: pspxc(:)     ! pspxc(npsp) from psps
  integer, pointer :: so_typat(:)  ! input variable so_typat(ntype)
  integer, pointer :: symafm(:)    ! input variable symafm(nsym)
  integer, pointer :: symrel(:,:,:)! input variable symrel(3,3,nsym)
  integer, pointer :: type(:)      ! input variable type(natom)

  real(dp) :: ecut                  ! input variable
  real(dp) :: ecutsm                ! input variable
  real(dp) :: ecut_eff              ! ecut*dilatmx**2 (dilatmx is an input variable)
  real(dp) :: etot,fermie,residm    ! EVOLVING variables
  real(dp) :: rprimd(3,3)           ! EVOLVING variables
  real(dp) :: tphysel               ! input variable
  real(dp) :: tsmear                ! input variable
  real(dp), pointer :: kptns(:,:)   ! input variable kptns(3,nkpt)
  real(dp), pointer :: occ(:)       ! EVOLVING variable occ(bantot)
  real(dp), pointer :: tnons(:,:)   ! input variable tnons(3,nsym)
  real(dp), pointer :: xred(:,:)    ! EVOLVING variable xred(3,natom)
  real(dp), pointer :: zionpsp(:)   ! zionpsp(npsp) from psps
  real(dp), pointer :: znuclpsp(:)  ! znuclpsp(npsp) from psps     
                                    ! Note the difference between znucl and znuclpsp !! 
  real(dp), pointer :: znucltypat(:)! znucltypat(ntype) from alchemy
  
  character*6 :: codvsn              ! version of the code
  character*132, pointer :: title(:) ! title(npsp) from psps
!Should make a list of supplementary infos

 end type hdr_type

!----------------------------------------------------------------------

 type bandstructure_type

  integer :: bantot                  ! total number of bands (sum(nband(:))
  integer :: nkpt                    ! number of k points
  integer :: nsppol                  ! number of spin-polarizations
  integer, pointer :: istwfk(:)      ! istwfk(nkpt) storage mode at each k point
  integer, pointer :: nband(:)       ! nband(nkpt*nsppol) number of bands
                                     !    at each k point and spin-polarisation
  integer, pointer :: npwarr(:)      ! npwarr(nkpt) number of plane waves at each k point
  real(dp), pointer :: kptns(:,:)    ! kptns(3,nkpt)  k-point vectors
  real(dp), pointer :: eig(:)        ! eig(bantot)  eigenvalues of each band
  real(dp), pointer :: occ(:)        ! occ(bantot)  occupation of each band
  real(dp), pointer :: doccde(:)     ! doccde(bantot)  derivative of the
                                     !    occupation of each band wrt energy (needed for RF)
  real(dp), pointer :: wtk(:)        ! wtk(nkpt)  weight of each k point

 end type bandstructure_type

!----------------------------------------------------------------------

 type anaddb_dataset_type

! Since all these input variables are described in the anaddb_help.htm
! file, they are not described in length here ...

! Real pointers
  real(dp), pointer :: qnrml1(:)  ! qnrml1(nph1l)
  real(dp), pointer :: qnrml2(:)  ! qnrml2(nph2l)
  real(dp), pointer :: qph1l(:,:) ! qph1l(3,nph1l)
  real(dp), pointer :: qph2l(:,:) ! qph2l(3,nph2l)

 end type anaddb_dataset_type

!----------------------------------------------------------------------

 type aim_dataset_type

! Since all these input variables are described in the aim_help.htm
! file, they are not described in length here ...

! Integer
  integer :: crit,denout,dltyp,gpsurf,irho,ivol,lapout,nsa,nsb,nsc
  integer :: ngrid(3)
  integer :: batom  !! Warning : corresponds to the input variable atom
  integer :: foll   !! Warning : corresponds to the input variable follow
  integer :: isurf  !! Warning : corresponds to the input variable surf
  integer :: irsur  !! Warning : corresponds to the input variable rsurf
  integer :: nph    !! Warning : corresponds to the input variable nphi
  integer :: npt    !! Warning : corresponds to the input variable inpt
  integer :: nth    !! Warning : corresponds to the input variable ntheta
  integer :: plden  !! Warning : not documented in help file ?!

! Real
  real(dp) :: atrad,coff1,coff2,dpclim,folstp,lgrad,lgrad2,lstep,lstep2,&
&  maxatd,maxcpd,phimax,phimin 
  real(dp) :: foldep(3),scal(3),vpts(3,4)
  real(dp) :: dr0    !! Warning : correspond to the input variable radstp
  real(dp) :: phi0   !! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   !! Warning : correspond to the input variable ratmin
  real(dp) :: th0    !! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax !! Warning : correspond to the input variable thetamax
  real(dp) :: themin !! Warning : correspond to the input variable thetamin

 end type aim_dataset_type

!----------------------------------------------------------------------

 end module defs_common
!!***
