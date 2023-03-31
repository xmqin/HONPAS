!============================================================
!= Sample program using the f90 FDF module : September 2007 =
!============================================================
!
!     Shows FDF capabilities..
!
PROGRAM SAMPLE
  USE fdf
  USE prec
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug
  character(20)              :: fname, axis, status
  character(2)               :: symbol(maxa)
  integer(sp)                :: i, j, ia, na, external_entry
  integer(sp)                :: isa(maxa)
  real(sp)                   :: wmix
  real(dp)                   :: cutoff, phonon_energy, factor
  real(dp)                   :: xa(3, maxa)
  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline


!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('sample.fdf', 'sample.out')

! Handle/Use fdf structure
  if (fdf_defined('new-style')) write(6,*) 'New-style stuff'

  na = fdf_integer('NumberOfAtoms', 0)
  write(6,*) 'Examples: na =', na

  fname = fdf_string('NameOfFile', 'whatever')
  write(6,*) 'Name of file:', fname
 
  cutoff = fdf_physical('MeshCutoff', 8.d0, 'Ry')
  write(6,*) 'MeshCutOff:', cutoff

  phonon_energy = fdf_physical('phonon-energy', 0.01d0, 'eV')
  write(6,*) 'Phonon Energy:', phonon_energy

  i = fdf_integer('SomeInt', 34)
  write(6,*) '#elems:', i

  wmix = fdf_single('WmixValue', 0.55)
  write(6,*) 'WmixValue:', wmix

  factor = fdf_double('FactorValue', 1.d-10)
  write(6,*) 'Factor:', factor

  debug = fdf_boolean('Debug', .TRUE.)
  write(6,*) 'Debug:', debug

  doit = fdf_boolean('DoIt', .FALSE.)
  write(6,*) 'Doit:', doit

  doit = fdf_defined('AtomicCoordinatesAndAtomicSpecies')
  write(6,*) 'AtomCoordsBlockDefined:', doit

  if (fdf_block('AtomicCoordinatesAndAtomicSpecies', bfdf)) then
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      do i= 1, 3
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      isa(ia) = fdf_bintegers(pline, 1)
      ia = ia + 1
    enddo
  endif

  write(6,*) 'AtomicCoordinatesAndAtomicSpecies:'
  do ia= 1, na
    write(6,'(3F10.6,I5)') (xa(i,ia),i=1,3), isa(ia)
  enddo

  if (fdf_block('AtomicInfo', bfdf)) then
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      do i= 1, 3
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo
  endif

  write(6,*) 'AtomicInfo:'
  do ia= 1, na
    write(6,'(3F10.6)') (xa(i,ia),i=1,3)
  enddo

  if (fdf_block('Other-Block', bfdf)) then

!   Forward reading
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Backward reading
    ia = 1
    do while((fdf_bbackspace(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Backward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Forward reading
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Forward reading with rewind
    call fdf_brewind(bfdf)
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward-with-rewind):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo
  endif

  if ( fdf_block('ListBlock',bfdf) ) then
     i = 0
     do while ( fdf_bline(bfdf,pline) )
        i = i + 1
        na = fdf_bnlists(pline)
        write(*,'(2(a,i0),a)') 'Listblock line: ',i,' has ',na,' lists'
        do ia = 1 , na
           j = -1
           call fdf_blists(pline,ia,j,isa)
           write(*,'(tr5,2(a,i0),a)') 'list ',ia,' has ',j,' entries'
           call fdf_blists(pline,ia,j,isa)
           write(*,'(tr5,a,1000(tr1,i0))') 'list: ',isa(1:j)
        end do
     end do
  end if

  if ( fdf_islist('MyList') ) then
     na = -1
     call fdf_list('MyList',na,isa)
     write(*,'(tr5,a,i0,a)') 'MyList has ',na,' entries'
     call fdf_list('MyList',na,isa)
     write(*,'(tr5,a,1000(tr1,i0))') 'MyList: ',isa(1:na)
  else
     write(*,*)'MyList was not recognized'
  end if

  if ( fdf_islist('externalentry') ) then
     write(*,*) 'externalentry is a list'
  else
     write(*,*) 'externalentry is not a list'
  end if

  external_entry = fdf_integer('externalentry', 60)
  write(6,*) 'ExternalEntry:', external_entry

  axis   = fdf_string('AxisXY', 'Cartesian')
  status = fdf_string('StatusXY', 'Enabled')
  write(6,*) 'Axis: ', TRIM(axis), ' | ', TRIM(status)

! Shutdown and deallocates fdf structure
  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM SAMPLE
