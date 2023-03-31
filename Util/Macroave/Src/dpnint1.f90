! 
! Modified by Alberto Garcia, March 2015
! This routine is included in this module with permission from D.R. Hamann.
!
 subroutine dpnint1(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 from routine
! dpnint by D.R. Hamann. 
! Changes:
!   -- A single value is returned
!   -- It can extrapolate, instead of stopping,
!      when called with an abscissa outside the
!      data range.
!   -- If the number of data points is less than
!      npoly+1, npoly is implicitly reduced, without
!      error, and without warning.
!   -- Debug interface 
!
! local polynomial interpolation of data yy on nn points xx
! giving value val on point r
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output interpolated value val on point r

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp), intent(in) :: xx(*),yy(*)
 real(dp), intent(in) :: r
 real(dp), intent(out) :: val
 integer, intent(in)   ::  nn,npoly
 logical, intent(in)   ::  debug

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,kk,iend

! interval halving search for xx(ii) points bracketing r

   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(r>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=r

!   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)
   iend = min(istart+npoly,nn)

 !  if (debug) print *, "istart, iend: ", istart, iend
   sum=0.0d0
   do iy=istart,iend
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, iend
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   val=sum

 end subroutine dpnint1
