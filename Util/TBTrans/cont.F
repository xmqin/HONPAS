      module cont
      interface 
      subroutine setupcontour(joutfile,NEn,
     .     EFermi0,EFermiL,EFermiR,
     .     WGF,contour,ncontour,ncontour0,nvolt)

      integer joutfile,nen,ncontour,ncontour0,nvolt
      real*8 EFermi0,EFermiL,EFermiR  
      complex*16, dimension (:), pointer:: contour,wgf
      end subroutine setupcontour

      end interface
      end module cont

