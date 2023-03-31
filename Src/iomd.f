! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iomd( na, isa, iza, xa, va, cell, vcell, varcel, istep,
     .                 istep0, temp, eks, getot, volume, Psol)
c *******************************************************************
c Saves positions, cell, and energies in a MD run (accumulative)
c J.Kohanoff August 1998, slightly modified by E. Artacho, Feb. 1999
c Modified to open and close each time by E. Artacho, Aug 2002
c *******************************************************************
c Two possibilities: |*everything into one single unformatted file
c                    | (for postprocessing programs, space saving)
c (formtt parameter) |*separated ascii files
c ********** INPUT **************************************************
c real*8  cell(3,3)  : Unit cell vectors
c real*8  vcell(3,3) : Velocities thereof
c integer na         : Number of atoms
c integer isa(na)    : Atomic species index
c integer iza(na)    : Atomic numbers
c real*8  xa(3,na)   : Atomic positions
c real*8  va(3,na)   : Velocities thereof
c logical varcel     : .true. when variable cell
c real*8  temp       : Temperature of ions
c real*8  eks        : Kohn-Sham energy
c real*8  getot      : Total energy
c integer istep      : Present time step
c integer istep0     : First time step
c real*8  volume     : cell volume in Ang**3
c real*8  Psol       : total pressure (static plus kinetik) in kBar
c *******************************************************************

      use precision,  only : dp
      use files,      only : slabel, label_length
      use units, only: eV

      implicit          none

      integer           na, isa(na), iza(na)
      integer           istep, istep0
      logical           varcel
      real(dp)          cell(3,3), xa(3,na), va(3,na), vcell(3,3),
     .                  temp, eks, getot, volume, Psol

c Internal variables and arrays
      logical, parameter :: IS_FORMATTED = .false.
      character(len=label_length+4) :: fncel
      character(len=label_length+4) :: fnene
      character(len=label_length+4) :: fnpos

      integer :: ia, iv, ix
      integer :: iucel, iuene, iupos
      
      external          io_assign, io_close

      ! Find name of file
      fnene = trim(slabel) // '.MDE'
      if ( IS_FORMATTED ) then
         fnpos = trim(slabel) // '.MDX'
         if (varcel) fncel = trim(slabel) // '.MDC'
      else
         fnpos = trim(slabel) // '.MD'
      end if

c Open file 

      call io_assign( iuene )
      open(iuene, file=fnene, form='formatted', position='append', 
     .     status='unknown')
      if ( IS_FORMATTED ) then
         call io_assign( iupos )
         open(iupos, file=fnpos, form='formatted', position='append',
     .        status='unknown')
         if ( varcel ) then 
            call io_assign( iucel )
            open(iucel,file=fncel,form='formatted',position='append',
     .           status='unknown' )
         end if
      else
         call io_assign( iupos )
         open(iupos,file=fnpos,form='unformatted',status='unknown',
     $        position="append")
      end if

      if(istep .eq. istep0) then
        write(iuene,"(6a)") '# Step','     T (K)','     E_KS (eV)',
     .     '    E_tot (eV)','   Vol (A^3)','    P (kBar)' 
      endif

C Write data on files

      write(iuene,'(i6,1x,f9.2,2(1x,f13.5),1x,f11.3,1x,f11.3)') 
     .            istep, temp, eks/eV, getot/eV, volume, Psol
      if ( IS_FORMATTED ) then
        write(iupos,*) istep
        do ia = 1,na
          write(iupos,'(i3,i6,3f13.6,3x,3f13.6)')
     .      isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
        enddo
        if ( varcel ) then
          write(iucel,*) istep
          write(iucel,'(2(3x,3f13.6))')
     .      ((cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3),iv=1,3)
        endif
      else
        write(iupos) istep, xa, va
        if ( varcel ) write(iupos) cell, vcell
      endif

C Close files

      call io_close( iuene )
      call io_close( iupos )
      if ( IS_FORMATTED .and. varcel ) call io_close( iucel )

      end
