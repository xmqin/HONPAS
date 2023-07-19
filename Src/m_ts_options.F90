!==========================================================================*
!                                                                          *
!  TRANSIESTA MODULE m_ts_options : Declaration of the variables           *
!  involved in a TRANSIESTA calculation.                                   *
!                                                                          *
!  Written by F.D.Novaes, May'07                                           *
!  onlyS option added by M.Paulsson May'09                                 *
!==========================================================================*
!  Contains the Subroutines:                                               *
!                                                                          *
!  1) read_ts_options : Reads the optional parameters from the fdf file    *
!                                                                          *
!==========================================================================*                                         


module m_ts_options

! SIESTA Modules used
USE precision, only : dp
USE siesta_options, only : fixspin
USE sys, only : die

implicit none
PUBLIC

!=========================================================================*
!  Arguments read from input file using the fdf package                    *
!--------------------------------------------------------------------------*

logical  :: savetshs     ! Saves the Hamiltonian and Overlap matrices if the 
                         ! the option TS.SaveHS is specified in the input file
logical  :: onlyS	 ! Option to only save overlap matrix
logical  :: mixH         ! Mixing of the Hamiltoninan instead of DM
logical  :: USEBULK      ! Use Bulk Hamiltonian in Electrodes
logical  :: TriDiag      ! true if tridiagonalization
logical  :: updatedmcr   ! Update DM values of ONLY Central Region
logical  :: FixQ         ! Fix Contact region charge
logical  :: UseVFix      ! Call the routine TSVHFix 
real(dp) :: voltfdf      ! Bias applied, Internally Volt=voltfdf/eV. 
                         ! EFermiL = voltfdf/2.0
                         ! EFermiR = -voltfdf/2.0
real(dp) :: CCEmin       ! EMin for the Complex Contour (LB)
real(dp) :: GFEta        ! Imaginary part of the Bias Contour  
real(dp) :: kT           ! Electronic Temperature
integer  :: nline        ! Number of points on the "line" segment of the Contour
integer  :: ncircle      ! Number of points on the circle part of the contour 
integer  :: npol         ! Number of poles included in the contour
integer  :: nvolt        ! Number of points for the Bias integartion part
integer  :: NBUFATL      ! Number of Left Buffer Atoms
integer  :: NBUFATR      ! Number of Right Buffer Atoms
character(20) :: smethod ! GF Numerical Integration Methods 
character(33) :: GFFileL ! Electrode Left GF File
character(33) :: GFFileR ! Electrode Right GF File
logical :: calcGF        ! Calculate the electrodes GF

!==========================================================================*
!==========================================================================*
!  Default Values for arguments read from input file                       *
!--------------------------------------------------------------------------*

logical, parameter :: savetshs_def = .true.
logical, parameter :: onlyS_def = .false.
logical, parameter :: tsdme_def = .true.
logical, parameter :: mixH_def = .false.
logical, parameter :: USEBULK_def = .true.
logical, parameter :: TriDiag_def = .false.
logical, parameter :: updatedmcr_def = .true.
logical, parameter :: FixQ_def = .false.
logical, parameter :: UseVFix_def = .true.
real(dp), parameter :: voltfdf_def = 0._dp   ! in Ry
real(dp), parameter :: CCEmin_def = -3.0_dp  ! in Ry
real(dp), parameter :: GFEta_def = 0.000001_dp  ! in Ry
real(dp), parameter :: kT_def = 0.0019_dp  ! in Ry
integer, parameter :: nline_def = 6
integer, parameter :: ncircle_def = 24
integer, parameter :: npol_def = 6
integer, parameter :: nvolt_def = 5
integer, parameter :: NBUFATL_def = 0
integer, parameter :: NBUFATR_def = 0
character(20), parameter :: smethod_def = 'gaussfermi'
character(33), parameter :: GFFileL_def = 'Left.GF'
character(33), parameter :: GFFileR_def = 'Right.GF'
logical, parameter :: calcGF_def = .true.



      CONTAINS

! *********************************************************************
! Subroutine to read the data for the TRANSIESTA program
!
!     It uses the FDF (Flexible Data Format) package
!     of J.M.Soler and A.Garcia
!
! Writen by F.D.Novaes May'07
!
! ***************************** INPUT *********************************
! integer na               : Number of atoms
! integer ns               : Number of species
! integer nspin            : Spin polarization
! **************************** OUTPUT *********************************

subroutine read_ts_options()

! SIESTA Modules Used
use parallel, only: IOnode, Nodes
use m_fdf_global, only: fdf_global_get
use units, only: eV
use m_ts_global_vars, only : ts_istep

#ifdef MPI
use mpi_siesta, only: MPI_Bcast, MPI_character, MPI_Comm_World
#endif


! Internal Variables
#ifdef MPI
integer :: MPIerror
#endif


if (IOnode) then
 write(*,*)
 write(*,'(2a)') 'ts_read_options: ', repeat('*', 62)
end if

!Set ts_istep default
ts_istep=0

! Reading TS Options from fdf ...
call fdf_global_get(savetshs,'TS.SaveHS',savetshs_def)
call fdf_global_get(onlyS,'TS.onlyS',onlyS_def)
call fdf_global_get(mixH,'TS.MixH',mixH_def)
call fdf_global_get(voltfdf,'TS.Voltage',voltfdf_def,'Ry') 
call fdf_global_get(USEBULK,'TS.UseBulkInElectrodes',USEBULK_def)
call fdf_global_get(TriDiag,'TS.TriDiag',TriDiag_def)
call fdf_global_get(updatedmcr,'TS.UpdateDMCROnly',updatedmcr_def)
call fdf_global_get(FixQ,'TS.FixContactCharge',FixQ_def)
call fdf_global_get(NBUFATL,'TS.BufferAtomsLeft',NBUFATL_def)
call fdf_global_get(NBUFATR,'TS.BufferAtomsRight',NBUFATR_def)
call fdf_global_get(CCEMin,'TS.ComplexContourEmin',CCEMin_def,'Ry')
call fdf_global_get(GFEta,'TS.biasContour.Eta',GFEta_def,'Ry')
call fdf_global_get(kT,'ElectronicTemperature',kT_def,'Ry')
call fdf_global_get(smethod,'TS.biasContour.method',smethod_def)
call fdf_global_get(npol,'TS.ComplexContour.NPoles',npol_def)
call fdf_global_get(ncircle,'TS.ComplexContour.NCircle',ncircle_def)
call fdf_global_get(nline,'TS.ComplexContour.NLine',nline_def)
call fdf_global_get(nvolt,'TS.biasContour.NumPoints',nvolt_def)
call fdf_global_get(GFFIleL,'TS.GFFileLeft',GFFileL_def)
call fdf_global_get(GFFileR,'TS.GFFileRight',GFFileR_def)
call fdf_global_get(calcGF,'TS.calcGF',calcGF_def)
call fdf_global_get(UseVFix,'TS.UseVFix',UseVFix_def)

! Output Used Options in OUT file ....
if (ionode) then
 write(*,1) 'ts_read_options: Save H and S matrices        =', savetshs
 write(*,1) 'ts_read_options: Mixing Hamiltonian           =', mixH
 write(*,6) 'ts_read_options: TranSIESTA Voltage           =', voltfdf/eV,' Volts'
! write(*,1) 'ts_read_options: Bulk Values in Elecs    =', USEBULK
 write(*,1) 'ts_read_options: TriDiag                      =', TriDiag 
 write(*,1) 'ts_read_options: Update DM Contact Reg. only  =', updatedmcr
! write(*,1) 'ts_read_options: Use VFix                     =', UseVFix
! write(*,1) 'ts_read_options: Fix Contact Charge      =', FixQ
 write(*,5) 'ts_read_options: N. Buffer At. Left           =', NBUFATL
 write(*,5) 'ts_read_options: N. Buffer At. Right          =', NBUFATR
 write(*,5) 'ts_read_options: N. Pts. Circle               =', ncircle
 write(*,5) 'ts_read_options: N. Pts. Line                 =', nline
 write(*,5) 'ts_read_options: N. Poles in Contour          =', npol
 write(*,5) 'ts_read_options: N. Pts. Bias Contour         =', nvolt
 write(*,6) 'ts_read_options: Contour E Min.               =', CCEmin,' Ry'
 write(*,7) 'ts_read_options: GFEta                        =', GFEta,' Ry'
 write(*,6) 'ts_read_options: Electronic Temperature       =', kT, ' Ry'
 write(*,10) 'ts_read_options: Bias Contour Method         =', smethod
 write(*,10) 'ts_read_options: Left GF File                =', GFFileL
 write(*,10) 'ts_read_options: Right GF File               =', GFFileR
 write(*,1) 'ts_read_options: Calculate GF                 =', calcGF
 write(*,1) 'ts_read_options: Save S and quit (onlyS)      =', onlyS
end if

if (IOnode) then
 write(*,'(2a)') 'ts_read_options: ', repeat('*', 62)
 write(*,*)
end if

if (IOnode) then
  write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24) 

! USEBULK and TriDiag
  if((.not. USEBULK) .and. TriDiag) then
    write(*,*) & 
              "WARNING: TriDiag only for UseBulkInElectrodes"
    write(*,*) "         Using normal inversion"
    TriDiag = .false. 
  end if

! Integration Method
  if( .not. (smethod .eq. 'gaussfermi' .or.   &
    smethod .eq. 'sommerfeld') ) then 
    write(*,*) &
       'WARNING: TS.biasContour.method=',smethod
    write(*,*) &
            'not defined: Using gaussfermi instead'
    smethod='gaussfermi'
  endif

  if (fixspin ) then
   write(*,*) 'Fixed Spin not possible in TS Calculations !'
   call die('Stopping code')
  end if

  write(*,'(3a)') repeat('*',24),' End: TS CHECKS AND WARNINGS ',repeat('*',26) 
  write(*,*)
end if

#ifdef MPI
call MPI_BCast(smethod,20,MPI_character,0,MPI_Comm_World,MPIerror)
#endif


1   format(a,4x,l1)
5   format(a,i5,a)
6   format(a,f10.4,a)
7   format(a,f12.6,a)
10  format(a,4x,a)
end subroutine read_ts_options

end module m_ts_options
