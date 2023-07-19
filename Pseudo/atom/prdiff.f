C
      subroutine prdiff(nconf,econf,jobold)
      implicit none
c
c   Prints out the energy differences between
c   different atomic configurations.
c
c   njtj  ***  modifications  ***
c     econf is able to handle larger numbers
c     of configurations.
c   njtj  ***  modifications  ***
c
c
C     .. Scalar Arguments ..
      integer nconf
      integer jobold  
C     ..
C     .. Array Arguments ..
      double precision econf(100)
C     ..
      integer nconfmax
      parameter (nconfmax = 20)

C     .. Local Scalars ..
      integer i, j, iu, nconf_ae, n_excitations
      double precision absval, maxabs, norm1, norm2
      logical ae_found
C     ..
C     .. Local Arrays ..
      double precision econf_ae(nconfmax,nconfmax)

C     ..
      write(6,9000) (i,i=1,nconf)
      do 10 i = 1, nconf
         write(6,9010) i, (econf(i)-econf(j),j=1,i)
   10 continue
 9000 format(/' &d total energy differences in series',
     $      //,' &d',2x,9i9)
 9010 format(' &d',1x,i2,1x,20f9.4)
      write(6,'(/,a,/)') '*----- End of series ----* spdfg &d&v'
c
c     Write external files for easier interfacing
c
      if (jobold .eq. 0) then      ! All-electron series just done

         call get_unit(iu)
         open(iu,file='AE_ECONF',form='formatted',status='unknown')
         rewind(iu)
         write(iu,"(i4)") nconf
         write(iu,*) ((econf(i)-econf(j),j=1,i-1),i=2,nconf)
         close(iu)

      else if (jobold .eq. 4) then   ! Ps test series just done

         call get_unit(iu)
         open(iu,file='PT_ECONF',form='formatted',status='unknown')
         rewind(iu)
         write(iu,"(i4)") nconf
         write(iu,*) ((econf(i)-econf(j),j=1,i-1),i=2,nconf)
         close(iu)
c
c        Check whether there is an AE_ECONF file and compute the differences
c
         inquire(file="AE_ECONF",exist=ae_found)
         if (ae_found) then
            call get_unit(iu)
            open(iu,file='AE_ECONF',form='formatted',status='old')
            rewind(iu)
            read(iu,*) nconf_ae
            if (nconf .ne. nconf_ae) STOP "nconf_ae_pt"
            if (nconf .ge. nconfmax) STOP "nconf_ae_overflow"
            read(iu,*) ((econf_ae(i,j),j=1,i-1),i=2,nconf)
            close(iu)
c
c           Write a new file
c
            call get_unit(iu)
            open(iu,file='ECONF_DIFFS',
     $             form='formatted',status='unknown')
            rewind(iu)
c
c           Compute AE-PT differences in excitation energies,
c           and print them, together with the maximum, mean_abs,
c           and root-mean-square values.
c
            n_excitations = 0
            maxabs = -1.0d0
            norm1 = 0.0d0
            norm2 = 0.0d0
            do i=2,nconf
               do j=1,i-1
                  econf_ae(i,j) = econf_ae(i,j) - (econf(i)-econf(j))
                  absval = abs(econf_ae(i,j))
                  if (absval .gt. maxabs) maxabs = absval
                  norm1 = norm1 + absval
                  norm2 = norm2 + absval*absval
                  n_excitations =  n_excitations + 1
               enddo
            enddo
            norm1 = norm1 / n_excitations
            norm2 = sqrt( norm2 / n_excitations)
            write(iu,"(i4)") n_excitations
            write(iu,"(40f12.5)") ((econf_ae(i,j),j=1,i-1),i=2,nconf)
            write(iu,"(3f12.5)") maxabs, norm1, norm2
            close(iu)
            
         endif ! found AE_ECONF
      endif
c
      end
