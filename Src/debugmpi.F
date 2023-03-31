! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module debugMpi
      logical, save     :: debugON = .false.
      integer, save     :: mpiUnit  = 23
      integer, parameter :: local_dp = selected_real_kind(10,300)
      end module

      subroutine debugMpiOn( )
      use parallel
      use debugMpi
      implicit none
      character*128         :: fname, cnode
      integer               :: istat
!------------------------------------------------------------------------- BEGIN
      if (.not.debugON) then
        debugON = .true.
        if (Node.lt.10) then
          write(cnode,'(I1)') Node
        else if (Node.lt.100) then
          write(cnode,'(I2)') Node
        else if (Node.lt.1000) then
          write(cnode,'(I3)') Node
        else if (Node.lt.10000) then
          write(cnode,'(I4)') Node
        else
          write(cnode,'(I5)') Node
        endif
        fname = TRIM('debug.'//cnode)
        open( mpiUnit, file=fname, status='unknown', iostat=istat )
        IF (ISTAT.NE.0) THEN
          write(*,*) 'cannot open file'
          STOP
        ENDIF
        call write_debug( "STARTING DEBUG FILE" )
      endif
!--------------------------------------------------------------------------- END
      end subroutine debugMpiOn

      subroutine write_debug( mesg )
      use debugMpi
      implicit none
      character*(*)        :: mesg
!------------------------------------------------------------------------- BEGIN
      if (debugON) then
        write(mpiUnit,*) mesg
        call pxfflush(mpiUnit)
      endif
!--------------------------------------------------------------------------- END
      end subroutine write_debug

      subroutine write_doublev( n, x )
      use debugMpi
      implicit none
      integer        :: n
      real(local_dp)       :: x(n)
!------------------------------------------------------------------------- BEGIN
      if (debugON) then
        write(mpiUnit,*) x
        call pxfflush(mpiUnit)
      endif
!--------------------------------------------------------------------------- END
      end subroutine write_doublev

      subroutine write_double( x )
      use debugMpi
      implicit none
      real(local_dp)       :: x
!------------------------------------------------------------------------- BEGIN
      if (debugON) then
        write(mpiUnit,*) x
        call pxfflush(mpiUnit)
      endif
!--------------------------------------------------------------------------- END
      end subroutine write_double

      subroutine write_int( a )
      use debugMpi
      implicit none
      integer         :: a
!------------------------------------------------------------------------- BEGIN
      if (debugON) then
        write(mpiUnit,*) a
        call pxfflush(mpiUnit)
      endif
!--------------------------------------------------------------------------- END
      end subroutine write_int
