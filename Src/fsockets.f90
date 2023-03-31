!F90 ISO_C_BINGING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

! Public functions provided by this module (with overloaded names):
!
!  SUBROUTINE create_socket(psockfd, inet, port, host)  
!    Creates and opens a server socket and begins listening to it
!    INTEGER,         INTENT(OUT):: psockfd ! socket file descriptor, to be
!                                             used to call read/write_buffer
!    INTEGER,         INTENT(IN) :: inet    ! socket type (0=>unix, 1=>inet)
!    INTEGER,         INTENT(IN) :: port    ! socket port number 
!    CHARACTER(LEN=*),INTENT(IN) :: host    ! host name ('localhost' or IP addr)
!
!  SUBROUTINE open_socket(psockfd, inet, port, host)      
!    Opens a client socket and connects it to the server
!    INTEGER,         INTENT(OUT):: psockfd ! socket file descriptor, to be
!                                             used to call read/write_buffer
!    INTEGER,         INTENT(IN) :: inet    ! socket type (0=>unix, 1=>inet)
!    INTEGER,         INTENT(IN) :: port    ! server's socket port number 
!    CHARACTER(LEN=*),INTENT(IN) :: host    ! server's host name
!
!  SUBROUTINE close_socket(psockfd) 
!    Closes a socket (server or client)     
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!
!  SUBROUTINE writebuffer(psockfd, fstring, plen)
!    Sends a character string through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    CHARACTER(LEN=*),INTENT(IN) :: fstring ! string to be sent
!    INTEGER,OPTIONAL,INTENT(IN) :: plen    ! string length. If present,
!                                             it must match that at readbuffer
!
!  SUBROUTINE writebuffer(psockfd, fdata)
!    Sends an integer number through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    INTEGER,         INTENT(IN) :: fdata   ! number to be sent
!
!  SUBROUTINE writebuffer(psockfd, fdata)
!    Sends a double-precision number through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    REAL(KIND=8),    INTENT(IN) :: fdata   ! number to be sent
!
!  SUBROUTINE writebuffer(psockfd, fdata, plen)
!    Sends a double-precision vector through the socket
!    INTEGER,            INTENT(IN) :: psockfd      ! socket file descriptor
!    REAL(KIND=8),TARGET,INTENT(IN) :: fdata(plen)  ! vector to be sent
!    INTEGER,            INTENT(IN) :: plen         ! vector length
!
!  SUBROUTINE readbuffer(psockfd, fstring, plen)
!    Receives a character string through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    CHARACTER(LEN=*),INTENT(OUT):: fstring ! string to be received
!    INTEGER,OPTIONAL,INTENT(IN) :: plen    ! string length. If present,
!                                             it must match that at writebuffer
!
!  SUBROUTINE readbuffer(psockfd, fdata)
!    Receives an integer number through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    INTEGER,         INTENT(OUT):: fdata   ! number to be received
!
!  SUBROUTINE readbuffer(psockfd, fdata)
!    Receives a double-precision number through the socket
!    INTEGER,         INTENT(IN) :: psockfd ! socket file descriptor
!    REAL(KIND=8),    INTENT(OUT):: fdata   ! number to be received
!
!  SUBROUTINE readbuffer(psockfd, fdata, plen)
!    Receives a double-precision vector through the socket
!    INTEGER,            INTENT(IN) :: psockfd      ! socket file descriptor
!    REAL(KIND=8),TARGET,INTENT(OUT):: fdata(plen)  ! vector to be received
!    INTEGER,            INTENT(IN) :: plen         ! vector length

MODULE F90SOCKETS
  USE ISO_C_BINDING
   
  IMPLICIT NONE

  INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
                       writebuffer_d, writebuffer_dv, &
                       writebuffer_i
  END INTERFACE 

  INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
                       readbuffer_dv, readbuffer_d, &
                       readbuffer_i
  END INTERFACE 

  INTERFACE
    SUBROUTINE open_csocket(psockfd, inet, port, host) BIND(C, name="open_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                  :: psockfd, inet, port
      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: host
    END SUBROUTINE open_csocket

    SUBROUTINE create_csocket(psockfd, inet, port, host) BIND(C, name="create_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                  :: psockfd, inet, port
      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: host
    END SUBROUTINE create_csocket
    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                  :: psockfd
      TYPE(C_PTR), VALUE                   :: pdata
      INTEGER(KIND=C_INT)                  :: plen
    END SUBROUTINE writebuffer_csocket       

    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                  :: psockfd
      TYPE(C_PTR), VALUE                   :: pdata
      INTEGER(KIND=C_INT)                  :: plen
    END SUBROUTINE readbuffer_csocket   

    SUBROUTINE close_csocket(psockfd) BIND(C, name="close_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                  :: psockfd
    END SUBROUTINE close_csocket
  END INTERFACE

CONTAINS
   
  SUBROUTINE open_socket(psockfd, inet, port, host)      
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(IN)                       :: inet, port
    INTEGER, INTENT(OUT)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)              :: host
    CHARACTER(LEN=1,KIND=C_CHAR) :: chost(1024)
    CALL fstr2cstr(host, chost)
    CALL open_csocket(psockfd, inet, port, chost)
  END SUBROUTINE

  SUBROUTINE create_socket(psockfd, inet, port, host)      
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(IN)                       :: inet, port
    INTEGER, INTENT(OUT)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)              :: host
    CHARACTER(LEN=1,KIND=C_CHAR) :: chost(1024)
    CALL fstr2cstr(host, chost)
    CALL create_csocket(psockfd, inet, port, chost)
  END SUBROUTINE
   
  SUBROUTINE fstr2cstr(fstr, cstr, plen)
    USE ISO_C_BINDING
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)              :: fstr
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: cstr(:)
    INTEGER, INTENT(IN), OPTIONAL             :: plen
    INTEGER i,n
    IF (PRESENT(plen)) THEN
       n = plen
       DO i=1,n
          cstr(i) = fstr(i:i)
       ENDDO
    ELSE
       n = LEN_TRIM(fstr)
       DO i=1,n
          cstr(i) = fstr(i:i)
       ENDDO
    END IF
    cstr(n+1) = C_NULL_CHAR
  END SUBROUTINE

  SUBROUTINE writebuffer_d (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(IN)                 :: fdata
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE

  SUBROUTINE writebuffer_i (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata
    INTEGER(KIND=C_INT), TARGET              :: cdata
    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE

  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER,OPTIONAL, INTENT(IN)             :: plen
    INTEGER                                  :: i
    INTEGER(KIND=C_INT), TARGET              :: flen
    CHARACTER(LEN=1, KIND=C_CHAR), POINTER   :: cstring(:)
    IF (PRESENT(plen)) THEN
      flen = plen
    ELSE
      flen = LEN(fstring)
      CALL writebuffer_csocket(psockfd, c_loc(flen), 4)
    ENDIF
    allocate( cstring(flen) )
    DO i = 1,flen
      cstring(i) = fstring(i:i)
    ENDDO
    CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), flen)
    deallocate( cstring )
  END SUBROUTINE

  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(IN), TARGET         :: fdata(plen)
    CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE

  SUBROUTINE readbuffer_d (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(OUT)                :: fdata
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
    CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
    fdata=cdata
  END SUBROUTINE

  SUBROUTINE readbuffer_i (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata
    INTEGER(KIND=C_INT), TARGET              :: cdata
    CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
    fdata = cdata
  END SUBROUTINE

  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER,OPTIONAL, INTENT(IN)             :: plen
    INTEGER                                  :: i
    INTEGER(KIND=C_INT), TARGET              :: clen
    CHARACTER(LEN=1, KIND=C_CHAR), POINTER   :: cstring(:)
    IF (PRESENT(plen)) THEN
      clen = plen
    ELSE
      CALL readbuffer_csocket(psockfd, c_loc(clen), 4)
    ENDIF
    allocate( cstring(clen) )
    CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), clen)
    fstring=""   
    DO i = 1,min(clen,len(fstring))
      fstring(i:i) = cstring(i)
    ENDDO
    deallocate( cstring )
  END SUBROUTINE

  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(OUT), TARGET        :: fdata(plen)
    CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE

  SUBROUTINE close_socket(psockfd)      
    USE ISO_C_BINDING  
    IMPLICIT NONE
    INTEGER, INTENT(IN)                      :: psockfd
    CALL close_csocket(psockfd)
  END SUBROUTINE

END MODULE
