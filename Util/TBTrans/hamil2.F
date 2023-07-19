      

      module hamil2
      interface 
      subroutine setupkham2(joutfile,tinit,kpoint,ispin,hsfile,
     &  lasto,Hk,Sk,nua,nuo,nspin) 
      integer joutfile        
      logical tinit
      real*8 kpoint(3) 
      integer ispin            
      character*33 hsfile    
      integer, dimension (:), pointer:: lasto
      complex*16, dimension (:), pointer:: Hk, Sk
      integer nua,nuo,nspin
      end subroutine setupkham2

      end interface
      end module hamil2

