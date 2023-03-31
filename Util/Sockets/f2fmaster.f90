! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program f2fmaster

! Test fortran-to-fortran socket communication
! JMS, Mar.2015

  use f90sockets, only: create_socket, open_socket, writebuffer, readbuffer

  integer,parameter    :: MSGLEN=32
  ! this should be the IP address or fully-qualified host name on which the master is running  
  ! in order to use UNIX sockets, this should be the name of the socket on which master will listen
  character(len=1024)  :: host='localhost'  
  integer              :: inet = 0   ! this has to be one to use a TCP/IP socket, zero to have a UNIX socket
  integer              :: port=21212 ! this is the port number (only used for TCP/IP sockets) 
  integer              :: socket, intIn, vectorLen
  character(len=MSGLEN):: msgIn, msgOut, strIn, strOut
  real(kind=8)         :: realIn, vectorIn(3)

  integer :: i
  
  host = TRIM(host)//achar(0)
  call create_socket( socket, inet, port, host )

  msgOut = 'wait'
  do i = 1,5
    call writebuffer( socket, msgOut, MSGLEN )
  enddo

  msgOut = 'receive string'
  strOut = 'Hello server'
  call writebuffer( socket, msgOut, MSGLEN )
  call writebuffer( socket, strOut, MSGLEN )

  msgOut = 'receive integer'
  call writebuffer( socket, msgOut, MSGLEN )
  call writebuffer( socket, 2468 )

  msgOut = 'receive real'
  call writebuffer( socket, msgOut, MSGLEN )
  call writebuffer( socket, 4.d0/3 )

  msgOut = 'receive vector'
  call writebuffer( socket, msgOut, MSGLEN )
  call writebuffer( socket, (/3.5d0,4.5d0,5.5d0/), 3 )

  msgOut = 'send string'
  call writebuffer( socket, msgOut, MSGLEN )
  call readbuffer( socket, strIn, MSGLEN )
  print*,'client: strIn = ',trim(strIn)

  msgOut = 'send integer'
  call writebuffer( socket, msgOut, MSGLEN )
  call readbuffer( socket, intIn )
  print*,'client: intIn = ',intIn

  msgOut = 'send real'
  call writebuffer( socket, msgOut, MSGLEN )
  call readbuffer( socket, realIn )
  print*,'client: realIn = ',realIn

  msgOut = 'send vector'
  call writebuffer( socket, msgOut, MSGLEN )
  call readbuffer( socket, vectorIn, 3 )
  print*,'client: vectorIn = ',vectorIn

  msgOut = 'stop'
  call writebuffer( socket, msgOut, MSGLEN )

end program f2fmaster

