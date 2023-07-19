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
      subroutine pixmol(iza, xa, na, last)
c *******************************************************************
c Writes and accumulates coordinates to be animated by Xmol
c Written by E. Artacho. February 1999. Modified to open/close 2003
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c logical   last      : true if last time step
c *******************************************************************

      use precision,      only: dp
      use periodic_table, only: symbol
      use files,          only: slabel, label_length

      implicit          none

      character(len=label_length+4)       :: paste
      integer                             :: na
      integer                             :: iza(na)
      real(dp)                            :: xa(3,na)
      logical                             :: last
      external          io_assign, io_close, paste

c Internal variables and arrays
 
      character(len=label_length+4), save :: fname
      integer                             :: i, ia
      integer,                       save :: unit
      logical,                       save :: frstme = .true.
      real(dp),                      save :: Ang = 1.0_dp/0.529177_dp
c -------------------------------------------------------------------

      character(len=2) :: sym

      if ( frstme ) then
        fname = paste(slabel,'.ANI')
        frstme = .false.
      endif

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', position='append',
     .      status='unknown')

      write(unit,'(i5)') na
      write(unit,*)
      do ia = 1, na
         sym =  symbol(iza(ia))
         write(unit,'(a2,2x,3f12.6)') sym , (xa(i,ia)/Ang,i=1,3)
      enddo

      call io_close(unit)
      
      return
      end

