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
      module m_iostruct

c     Reads and Saves structural information in "crystallography" format
!     
!     Cell vectors in Angstroms
!     Atomic positions in fractional coordinates

c     Alberto Garcia, Sep. 2005. Based on ioxv by J.M.Soler. July 1997.

      use precision,   only : dp
      use parallel,    only : IONode
      use fdf,         only : fdf_string
      use units,       only : Ang
      use m_mpi_utils, only : broadcast
      use siesta_geom,    only : xa, isa, cisa
      use alloc,       only : re_alloc
      use sys,         only : die
      use files,       only : slabel, label_length
      
      implicit none

      public :: write_struct
      public :: read_struct  

      private

      CONTAINS

      subroutine read_struct( na, cell)
!     
      integer, intent(out)  ::          na
      real(dp), intent(out) ::          cell(3,3)

      real(dp) :: xfrac(3)
      integer  :: dummy
      character(len=label_length+10), save :: fname
      integer                              :: ia, iu, iv
      integer                              :: ix, iostat
      logical,                        save :: frstme = .true.
      character(len=label_length+10)       :: paste
      external          io_assign, io_close, paste


      if (frstme) then
         if (IOnode) then
            fname = paste( slabel, '.STRUCT_IN' )
         endif
         frstme = .false.
      endif

      if (IOnode) then
         call io_assign( iu )
         open(iu,file=fname,form='formatted',
     $                      status='old', iostat=iostat)      
         if (iostat /= 0) call die(trim(fname) // " not found")

         read(iu,*) ((cell(ix,iv),ix=1,3),iv=1,3)
         cell = cell * Ang
         read(iu,*) na
      endif

      call broadcast(na)
      call broadcast(cell(1:3,1:3))

      nullify(isa,xa)
      call re_alloc(isa,1,na,name='isa',routine='read_struct')
      call re_alloc(xa,1,3,1,na,name='xa',routine='read_struct')
      if (IOnode) then
         do ia = 1,na
            read(iu,*) isa(ia), dummy, xfrac(1:3)
            xa(:,ia) = matmul(cell,xfrac(1:3))
         enddo
         call io_close( iu )
         write(*,*)
     $         " -- Read structural information from " // trim(fname)
      endif

      call broadcast(isa(1:na))
      call broadcast(xa(1:3,1:na))

! Construct references
      nullify(cisa)
      call re_alloc(cisa,1,na,name="cisa",routine="read_struct")
      do ia = 1, na
         write(cisa(ia), '("siesta:e",i3.3)') isa(ia)
      enddo

      end subroutine read_struct
!---------------------------------------------------------------------
      subroutine write_struct(cell, na, isa, iza, xa, moved)
!     
c     real*8  cell(3,3)  : Unit cell vectors
c     integer na         : Number of atoms
c     integer isa(na)    : Atomic species index
c     integer iza(na)    : Atomic numbers
c     real*8  xa(3,na)   : Atomic positions
c     logical moved      :  True if structure is the "predicted"
c                          one after application of forces/stress.

      integer, intent(in)  ::          na, isa(na), iza(na)
      real(dp), intent(in) ::          cell(3,3), xa(3,na)
      logical, intent(in), optional :: moved

      character(len=label_length+11)       :: paste
      external          io_assign, io_close, paste, reclat

c     Internal variables and arrays
      real(dp)                             :: celli(3,3)
      real(dp)                             :: xfrac(3)
      character(len=200)                   :: fname
      integer                              :: ia, iu, iv, ix

      logical                              :: atoms_moved_after_forces 

C     Only do reading and writing for IOnode

      if (.not. IOnode) RETURN

      atoms_moved_after_forces = .false.
      if (present(moved)) then
         atoms_moved_after_forces = moved
      endif

      if (atoms_moved_after_forces) then
         fname = paste( slabel, '.STRUCT_NEXT_ITER' )
      else
         fname = paste( slabel, '.STRUCT_OUT' )
      endif

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,'(3x,3f18.9)')
     .     ((cell(ix,iv)/Ang,ix=1,3),iv=1,3)
      write(iu,*) na
      call reclat(cell, celli, 0)
      do ia = 1,na
         xfrac(:) = matmul(transpose(celli),xa(:,ia))
         write(iu,'(i3,i6,3f18.9)') isa(ia),iza(ia),xfrac(1:3)
      enddo

      call io_close( iu )

      end subroutine write_struct

      end module m_iostruct

