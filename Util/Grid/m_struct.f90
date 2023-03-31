! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      module m_struct
!     
!     Cell vectors in Angstroms
!     Atomic positions in fractional coordinates

!     Alberto Garcia, Sep. 2005. Based on ioxv by J.M.Soler. July 1997.

      implicit none

      integer, parameter, private   :: dp = selected_real_kind(14,100)
      real(dp), parameter, private  :: Ang = 1._dp / 0.529177_dp

      type, public :: struct_t
        integer ::          na = 0
        real(dp)::          cell(3,3)
        real(dp), allocatable ::          xa(:,:)
        integer, allocatable  ::          isa(:)
        integer, allocatable  ::          iza(:)
      end type 
  
      public :: write_struct
      public :: read_struct  
      public :: clean_struct
      public :: replicate_struct

      private

      interface
         subroutine die(str)
         character(len=*), intent(in), optional  :: str
         end subroutine die
      end interface

      CONTAINS

      subroutine clean_struct(str)
      type(struct_t), intent(inout) :: str

      if (allocated(str%xa)) then
         deallocate(str%xa)
      endif
      if (allocated(str%isa)) then
         deallocate(str%isa)
      endif
      if (allocated(str%iza)) then
         deallocate(str%iza)
      endif
      str%cell = 0.0_dp
      str%na = 0
      end subroutine clean_struct


      subroutine read_struct( fname, str)
!     
      character(len=*), intent(in)       :: fname
      type(struct_t), intent(inout)      :: str

      real(dp) :: xfrac(3)
      integer                              :: ia, iv
      integer                              :: ix, iostat, na

      integer :: iu = 1

      call clean_struct(str)

         open(iu,file=fname,form='formatted', &
              status='old', iostat=iostat)      
         if (iostat /= 0) call die(trim(fname) // " not found")

         read(iu,*) ((str%cell(ix,iv),ix=1,3),iv=1,3)
         str%cell = str%cell * Ang
         read(iu,*) na
         str%na = na
         allocate(str%xa(3,na),str%isa(na),str%iza(na))

         do ia = 1,na
            read(iu,*) str%isa(ia), str%iza(ia), xfrac(1:3)
            str%xa(:,ia) = matmul(str%cell,xfrac(1:3))
         enddo
         close(iu)

      end subroutine read_struct

!---------------------------------------------------------------------
      subroutine write_struct(fname, str)
!     
      character(len=*), intent(in)     :: fname
      type(struct_t), intent(in)       :: str

!     Internal variables and arrays
      real(dp)                             :: celli(3,3)
      real(dp)                             :: xfrac(3)
      integer                              :: ia, iv, ix

      integer :: iu = 1

      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,'(3x,3f18.9)') ((str%cell(ix,iv)/Ang,ix=1,3),iv=1,3)
      write(iu,*) str%na
      call reclat(str%cell, celli, 0)
      do ia = 1,str%na
         xfrac(:) = matmul(transpose(celli),str%xa(:,ia))
         write(iu,'(i3,i6,3f18.9)') str%isa(ia),str%iza(ia),xfrac(1:3)
      enddo

      close(iu)

      end subroutine write_struct

      subroutine replicate_struct(str0,nsc,str)
!     
      type(struct_t), intent(inout)      :: str0
      integer, intent(in)                :: nsc(3)  
      type(struct_t), intent(inout)      :: str

      real(dp) :: xc(3)
      integer                              :: ia, iv, ja
      integer                              :: ix, na, na0
      integer                              :: i1, i2, i3
  
      call clean_struct(str)
      do iv = 1, 3
         str%cell(:,iv) = nsc(iv) * str0%cell(:,iv)
      enddo
      na0 = str0%na
      str%na = na0 * product(nsc(1:3))

      na = str%na
      allocate(str%xa(3,na),str%isa(na),str%iza(na))

      do ia = 1,na
         ja = mod(ia-1,na0) + 1
         str%isa(ia)  = str0%isa(ja)
         str%iza(ia)  = str0%iza(ja)
      enddo
      
      ia = 0
      
      DO I3 = 0,NSC(3)-1
         DO I2 = 0,NSC(2)-1
            DO I1 = 0,NSC(1)-1
              
               do  IX = 1,3
                  XC(IX) = str0%cell(IX,1)*I1 + &
                           str0%cell(IX,2)*I2 + &
                           str0%cell(IX,3)*I3
               enddo
               do JA = 1,NA0
                  IA = IA + 1
                  do  IX = 1,3
                     str%XA(IX,IA) = str0%XA(IX,JA) + XC(IX)
                  enddo
               enddo
              
            enddo
         enddo
      enddo

    end subroutine replicate_struct
  end module m_struct

