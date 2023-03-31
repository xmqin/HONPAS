! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program horizontal

! Reads an array of data and reshapes it into a horizontal shape
! J.M.Soler. Nov.2012

  implicit none
  integer,parameter:: dp = kind(1.d0)
  integer,parameter:: maxCols = 100
  integer,parameter:: maxRows = 100

  character(len=1000):: line
  integer :: col0, iRow, jRow, nCols, totCols, totRows
  real(dp):: cols(maxCols), data(maxRows,maxCols)

  jRow = 0
  col0 = 0
  totRows = 0
  totCols = 0
  do iRow = 1,maxRows
    read(5,'(a)',end=999) line
    call parser( line, maxCols, nCols, cols )
    if (nCols==0) then    ! blank line
      jRow = 0
      col0 = totCols
    else
      if (jRow==0) totCols = col0+nCols
      jRow = jRow+1
      totRows = max(totRows,jRow)
      data(jRow,col0+1:col0+nCols) = cols(1:nCols)
    endif
  enddo
  stop 'ERROR: parameter maxRows too small'

999 continue   ! Come here upon end of file
  do iRow = 1,totRows
    write(6,'(100e18.9)') data(iRow,1:totCols)
  end do

end program horizontal


subroutine parser( line, maxCols, nCols, cols )

  implicit none
  integer,parameter:: dp = kind(1.d0)
  character(len=*),intent(inout):: line           ! returns without comments
  integer,         intent(in)   :: maxCols        ! size(cols)
  integer,         intent(out)  :: nCols          ! number of columns
  real(dp),        intent(out)  :: cols(maxCols)  ! columns in line

  character(len=1000):: myLine
  integer:: iCol, nc

  ! Remove comments
  nc = scan(line,'#')-1          ! last character before comment
  if (nc<0) nc=len(trim(line))   ! if no comment, last nonblank character
  line = line(1:nc)
  myLine = line

  ! Loop on columns
  do iCol = 1,maxCols
    myLine = adjustl(myLine)     ! remove leading blanks
    nc = scan(myLine,' ')-1      ! last character before next blank
    if (nc<=0) then              ! no more cols
      nCols = iCol-1
      return                     ! normal return point
    else
      read(myLine(1:nc),*) cols(iCol)
      myLine = myLine(nc+1:)
    endif
  end do ! iCol
  stop 'parser ERROR: size(cols) too small'

end subroutine parser


