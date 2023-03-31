MODULE m_rusage

  USE ISO_C_BINDING

  IMPLICIT NONE

  integer, parameter :: int_8 = selected_int_kind(15)

  PRIVATE

  PUBLIC :: m_memory, rss_max

TYPE, BIND(C) :: timeval 
  INTEGER(C_LONG) :: tv_sec, tv_usec
END TYPE

TYPE, BIND(C)    :: rusage
!
! structure to call getrusage
!
TYPE(timeval)   :: ru_utime, ru_stime
INTEGER(C_LONG) :: ru_maxrss, ru_ixrss, ru_idrss
INTEGER(C_LONG) :: ru_isrss, ru_minflt, ru_majflt, &
    ru_nswap, ru_inblock, ru_oublock, ru_msgsnd, ru_msgrcv, &
    ru_nsignals, ru_nvcsw, ru_nivcsw 
END TYPE

INTERFACE
INTEGER(C_INT) FUNCTION getrusage (who, rusage) BIND(C, name='getrusage')
 USE ISO_C_BINDING
 INTEGER(C_INT), VALUE :: who
 TYPE(C_PTR),    VALUE :: rusage
END FUNCTION getrusage
END INTERFACE

CONTAINS

  ! Returns the maximum resident size (a rough proxy for memory) in Kbytes

FUNCTION m_memory()
INTEGER(KIND=int_8)                      :: m_memory

INTEGER(C_INT)                           :: ret
TYPE(rusage), TARGET                     :: usage

ret = getrusage(0, C_LOC(usage))
m_memory = usage%ru_maxrss 
END FUNCTION m_memory

FUNCTION rss_max()
real :: rss_max
INTEGER(KIND=int_8)    :: memk
 
memk = m_memory()
rss_max = real(memk) / 1024
end FUNCTION rss_max

end module m_rusage

