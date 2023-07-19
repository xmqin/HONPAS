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
      subroutine iofa( na, fa )
c *******************************************************************
c Writes forces in eV/Ang
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer na           : Number atoms
c real*8  fa(3,na)     : Forces on the atoms
c *******************************************************************

      use fdf
      use files,     only : slabel, label_length
      use precision, only : dp

      implicit          none

      character(len=label_length+3) :: paste
      integer                       :: na
      real(dp)                      :: fa(3,*)
      external          io_assign, io_close, paste

c Internal 
      character(len=label_length+3), save :: fname
      integer                             :: ia, iu, ix
      logical,                       save :: frstme = .true.
      real(dp),                      save :: Ang, eV
c -------------------------------------------------------------------

      if (frstme) then
        Ang    = 1.d0 / 0.529177d0
        eV     = 1.d0 / 13.60580d0
        fname  = paste( slabel, '.FA' )
        frstme = .false.
      endif

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown',
     $      position='rewind')      

      write(iu,'(i6)') na
      write(iu,'(i6,3f12.6)') (ia, (fa(ix,ia)*Ang/eV,ix=1,3), ia=1,na)

      call io_close( iu )

      return
      end
