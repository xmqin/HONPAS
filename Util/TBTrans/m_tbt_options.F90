module m_tbt_options

use precision, only : dp
use fdf
use m_fdf_global, only : fdf_global_get
use files, only : slabel

implicit none
public

real(dp) Emin               !Calc. transmission on [Emin;Emax]
real(dp) Emax
real(dp) Volt
real(dp) GFeta              ! state broadening

integer ncontour
integer neigch            ! No. eigenchannels to calculate
integer NBUFATL,NBUFATR   ! No. buffer atoms, L/R

logical USEBULK           ! If true use bulk params. in L/R regions
logical sppol

character*33 hsfile       !name of HS-input file
integer isoat1,isoat2

integer NA1L,NA2L,nqL 
integer NA1R,NA2R,nqR

character*33 Lhsfile,Rhsfile

integer Lnucuse, Rnucuse

logical CalcIeig


! PARAMETERS

integer RNode   ! Root Node

contains

subroutine tbt_read_options()

use parallel, only : IOnode

character*33 paste
external paste

real(dp) eV
parameter ( eV = 1_dp / 13.60580_dp )

real(dp), parameter :: Emin_def = -2.0_dp*eV
real(dp), parameter :: Emax_def = 2.0d0*eV
real(dp), parameter :: Volt_def = 0._dp
integer, parameter :: ncontour_def = 100
real(dp), parameter :: GFeta_def = 0.00001_dp
logical, parameter :: UseBulk_def = .true.
integer, parameter :: neigch_def = 0
logical, parameter :: sppol_def = .false.
integer, parameter :: NBUFATL_def=0
integer, parameter :: NBUFATR_def=0
logical, parameter :: CalcIeig_def = .false.
character*33 :: hsfile_def

hsfile_def=paste(slabel,'.TSHS')



! Read from fdf

call fdf_global_get(Emin,'TS.TBT.Emin',Emin_def,'Ry')
call fdf_global_get(Emax,'TS.TBT.Emax',Emax_def,'Ry')
call fdf_global_get(Volt,'TS.Voltage',Volt_def)

call fdf_global_get(ncontour,'TS.TBT.NPoints',ncontour_def)
call fdf_global_get(GFeta,'TS.TBT.Eta',GFeta_def,'Ry')

call fdf_global_get(USEBULK,'TS.UseBulkInElectrodes',UseBulk_def)
call fdf_global_get(neigch,'TS.TBT.NEigen',neigch_def)

call fdf_global_get(sppol,'SpinPolarized',sppol_def)

call fdf_global_get(NBUFATL,'TS.BufferAtomsLeft',NBUFATL_def)
call fdf_global_get(NBUFATR,'TS.BufferAtomsRight',NBUFATR_def) 

call fdf_global_get(hsfile,'TS.TBT.HSFile',hsfile_def)

call fdf_global_get(Lhsfile,'TS.HSFileLeft','LeftELEC.TSHS')
call fdf_global_get(Rhsfile,'TS.HSFileRight','RightELEC.TSHS')

call fdf_global_get(Lnucuse,'TS.NumUsedAtomsLeft',0)
call fdf_global_get(Rnucuse,'TS.NumUsedAtomsRight',0)

call fdf_global_get(isoat1,'TS.TBT.PDOSFrom',0)
call fdf_global_get(isoat2,'TS.TBT.PDOSTo',0)

nqL=1
nqR=1
NA1L=1
NA2L=1
NA1R=1
NA2R=1

call fdf_global_get(CalcIeig,'TS.TBT.CalcIeig',CalcIeig_def)




! Print Values in output
if (IOnode) then
  write(*,'(a,f10.4,a)') &
  'tbtrans: Voltage                         = ',Volt
  write(*,'(a,2F14.6)') 'Emin,Emax = ',Emin/eV,Emax/eV
  write(*,'(a,I10)') 'Energy points in tbtrans = ',ncontour
  write(*,'(a,7L)') 'tbtrans: Spin-polarized = ', sppol
end if ! IOnode


end subroutine tbt_read_options

end module m_tbt_options
