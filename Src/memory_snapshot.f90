subroutine memory_snapshot(str)
 use m_rusage, only :  rss_max

 character(len=*), intent(in) :: str

 write(6,"(a,f12.2)") " &m -- Max memory " // trim(str), rss_max()

end subroutine memory_snapshot
