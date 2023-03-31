! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program f2fslave

! Test fortran-to-fortran socket communication
! JMS, Mar.2015

  use f90sockets, only: open_socket, writebuffer, readbuffer

  integer,parameter    :: MSGLEN=32
  ! this should be the IP address or fully-qualified host name on which the slave should connect
  ! in order to use UNIX sockets, this should be the name of the socket.
  ! these (and port in case of TCP sockets) should match the address and port on which the master is running
  character(len=1032)  :: host='localhost'
  integer              :: inet = 0   ! this has to be one to use a TCP/IP socket, zero to have a UNIX socket
  integer              :: port=21212 ! this is the port number (only used for TCP/IP sockets) 
  integer              :: socket, intIn, vectorLen
  character(len=MSGLEN):: msgIn, strIn, strOut
  real(kind=8)         :: doubleIn, vectorIn(3)
    
  host = TRIM(host)//achar(0)
!  call create_socket( socket, inet, port, host )
  call open_socket( socket, inet, port, host )
  
  do    ! receive-send iteration
    call readbuffer( socket, msgIn, MSGLEN )
    print*,'server: msgIn= ',trim(msgIn)

    if (trim(msgIn)=='receive string') then
      call readbuffer( socket, strIn, MSGLEN )
      print*,'server: strIn = ',trim(strIn)
    elseif (trim(msgIn)=='receive integer') then
      call readbuffer( socket, intIn )
      print*,'server: intIn =',intIn
    elseif (trim(msgIn)=='receive real') then
      call readbuffer( socket, doubleIn )
      print*,'server: doubleIn =',doubleIn
    elseif (trim(msgIn)=='receive vector') then
      call readbuffer( socket, vectorIn, 3 )
      print*,'server: vectorIn =',vectorIn
    elseif (trim(msgIn)=='send string') then
      strOut = 'Hello client'
      call writebuffer( socket, strOut, MSGLEN )
    elseif (trim(msgIn)=='send integer') then
      call writebuffer( socket, 123456 )
    elseif (trim(msgIn)=='send real') then
      call writebuffer( socket, 1.d0/3 )
    elseif (trim(msgIn)=='send vector') then
      call writebuffer( socket, (/1.1d0,2.2d0,3.3d0/), 3 )
    elseif (trim(msgIn)=='stop') then
      exit
    elseif (trim(msgIn)=='wait') then
      cycle
    endif

  enddo

end program f2fslave

