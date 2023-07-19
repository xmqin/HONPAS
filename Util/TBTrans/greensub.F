      module greensub
      interface 

C FDN cell,kscell and kdispl added as dummy
      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,
     &     hsfile,nua,lasto,nuo,nspin,cell,kscell,kdispl,
     &     Hk,Sk,Hk2,Sk2)    ! ->
C FDN

      integer joutfile
      logical tkham,tinit,tlast
      real*8 kpoint(3)
      integer ispin
      character*33 hsfile
      integer, dimension (:), pointer:: lasto
      complex*16, dimension (:), pointer:: Hk, Sk, Hk2,Sk2
      integer nua,nuo,nspin
C FDN
      real*8 cell(3,3)
      integer kscell(3,3)
      real*8 kdispl(3)
C FDN
      end subroutine sethhm2

      end interface

      interface 
C FDN cell, kscell and displ added as dummy
      subroutine mkqgrid(joutfile,TRS,NA1,NA2,nq,q,wq,
     . cell, kscell, displ)

      integer joutfile          !output file unit
      logical TRS               ! Use time-reversal symmetry or not
      integer NA1,NA2           !no. repetitions of simple unitcell

      integer nq               
      real*8, pointer:: q(:,:),wq(:)
C FDN
      real*8 cell(3,3),displ(3)
      integer kscell(3,3)
C FDN

      end subroutine mkqgrid

      end interface

      end module greensub
