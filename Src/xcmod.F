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
      module xcmod

        use precision, only : dp
        use parsing, only: parse

        implicit none

        integer, parameter      :: MaxFunc = 10
        integer,           save :: nXCfunc
        character(len=10), save :: XCauth(MaxFunc)
        character(len=10), save :: XCfunc(MaxFunc)
        real(dp),          save :: XCweightX(MaxFunc)
        real(dp),          save :: XCweightC(MaxFunc)

        public :: nXCfunc, XCauth, XCfunc
        public :: XCweightX, XCweightC
        public :: setXC
        private

      contains

        subroutine setXC
C
C  This subroutine reads the exchange-correlation functional information
C
        use fdf
        use parsing
        use parallel,  only : Node
        use sys,       only : die
#ifdef MPI
        use mpi_siesta
#endif

        implicit none

C Local variables
        character(len=80)  :: names

        character(len=130) :: line
        integer            :: integs(4)
        integer            :: iu
        integer            :: lastc
        integer            :: lc(0:3)
        integer            :: n
        integer            :: ni
        integer            :: nn
        integer            :: nr
        integer            :: nv
        real(dp)           :: reals(4)
        real(dp)           :: values(4)

#ifdef MPI
        integer            :: MPIerror
#endif

        if (Node.eq.0) then
C Read XC functionals
          if (fdf_block('xc.hybrid',iu)) then
            read(iu,'(a)') line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (ni.eq.0) then
              call die('Number of functionals missing in XC.hybrid')
            endif
            nXCfunc = abs(integs(1))
            if (nXCfunc.gt.MaxFunc) then
              call die('Too many functionals in XC.hybrid')
            endif
            do n = 1,nXCfunc
              read(iu,'(a)') line
              lastc = index(line,'#') - 1
              if (lastc .le. 0) lastc = len(line)
              call parse( line(1:lastc), nn, lc, names, nv, values,
     .                    ni, integs, nr, reals )
              if (nn.gt.0) then
                XCfunc(n) = names(lc(0)+1:lc(1))
              else
                XCfunc(n) = 'LDA'
              endif
              if (nn.gt.1) then
                XCauth(n) = names(lc(1)+1:lc(2))
              else
                XCauth(n) = 'PZ'
              endif
              if (nr.gt.1) then
                XCweightX(n) = reals(1)
                XCweightC(n) = reals(2)
              elseif (nr.eq.1) then
                XCweightX(n) = reals(1)
                XCweightC(n) = reals(1)
              else
                XCweightX(n) = 1.0_dp
                XCweightC(n) = 1.0_dp
              endif
            enddo
          else
            nXCfunc = 1 
            XCfunc(1) = fdf_string('xc.functional','LDA')
            XCauth(1) = fdf_string('xc.authors','PZ')
            XCweightX(1) = 1.0_dp
            XCweightC(1) = 1.0_dp
          endif

C Output data for hybrid functionals
          if (nXCfunc.gt.1) then
            write(6,'(/,''xc:'')')
            write(6,'(''xc:  Hybrid exchange-correlation functional:'')
     .        ')
            write(6,'(''xc:'')')
            write(6,'(''xc: Number Functional Authors  '',
     .        ''   Weight(Ex)   Weight(Ec)'')')
            do n = 1,nXCfunc
            write(6,'(''xc: '',i4,7x,a3,5x,a10,3x,f5.3,8x,f5.3)')
     .          n,XCfunc(n),XCauth(n),XCweightX(n),XCweightC(n)
            enddo
            write(6,'(''xc:'')')
          endif

        endif

C Globalise data over all nodes
#ifdef MPI
        call MPI_Bcast(nXCfunc,1,MPI_integer,0,MPI_Comm_World,
     .                 MPIerror)
        call MPI_Bcast(XCfunc,10*nXCfunc,MPI_character,0,MPI_Comm_World,
     .                 MPIerror)
        call MPI_Bcast(XCauth,10*nXCfunc,MPI_character,0,MPI_Comm_World,
     .                 MPIerror)
        call MPI_Bcast(XCweightX,nXCfunc,MPI_double_precision,0,
     .                 MPI_Comm_World,MPIerror)
        call MPI_Bcast(XCweightC,nXCfunc,MPI_double_precision,0,
     .                 MPI_Comm_World,MPIerror)
#endif

        end subroutine setXC

      end module xcmod
