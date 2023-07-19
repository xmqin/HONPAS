c
      program vps2cdf
c
      use netcdf

      implicit none
c
c     This program converts pseudopotential "VPS" files (created by 
c     the ATOM program) from Binary to netCDF.
c
      integer ncid, iret
      integer nrp_id, n_up_id, n_down_id
      integer r_id, v_up_id, v_down_id, l_up_id, l_down_id
      integer core_id, val_id
      integer ptr
c
      double precision a, b, zion
      character*2 nameat, corr
      character*3 rel
      character*4 core
      character*10 ray(6), title(7)
c
      character*70 title_str
      character*60 ray_str
c
      integer i, j, lo, nrp, npotd, npotu
      double precision, allocatable  :: v(:)
c
      integer iargc, nargs
      character*70 cdf_file, binary_file
c
c     Let's get the files from the command line:
c
      nargs = iargc()
      if (nargs .ne. 2) then
         write(0,*) 'Usage: vps2cdf vps_file netCDF_file'
         stop 
      endif
c
      call getarg(1,binary_file)
      call getarg(2,cdf_file)
c
c      open files
c
      open(unit=2,file=binary_file,form='unformatted',status='old')
      iret = nf90_create(cdf_file,NF90_CLOBBER,ncid)

c
      read(2) nameat, corr, rel, core, (ray(j),j=1,6), 
     &         (title(j),j=1,7), npotd, npotu, nrp, a, b, zion

      do i=1,7
         ptr = (i-1)*10+1
         title_str(ptr:ptr+9) = title(i)
      enddo
      do i=1,6
         ptr = (i-1)*10+1
         ray_str(ptr:ptr+9) = ray(i)
      enddo

      allocate(v(nrp))
c
c     Attributes or variables??
c
      
      iret = nf90_put_att(ncid,nf90_global,'Element',nameat)
      iret = nf90_put_att(ncid,nf90_global,'Title',ray_str)
      iret = nf90_put_att(ncid,nf90_global,'Config',title_str)
      iret = nf90_put_att(ncid,nf90_global,'Correlation',corr)
      iret = nf90_put_att(ncid,nf90_global,'Relativistic',rel)
      iret = nf90_put_att(ncid,nf90_global,'Core',core)
      iret = nf90_put_att(ncid,nf90_global,'Valence_charge',zion)
      iret = nf90_put_att(ncid,nf90_global,'a_parameter',a)
      iret = nf90_put_att(ncid,nf90_global,'b_parameter',b)
c
      iret = nf90_def_dim(ncid,'nrp',nrp,nrp_id)

      if (npotd.ne.0) then
         iret =  nf90_def_dim(ncid,'n_down',npotd,n_down_id)
         iret = nf90_def_var(ncid,'l_down',nf90_int,n_down_id,l_down_id)
      endif
      if (npotu.ne.0) then
         iret = nf90_def_dim(ncid,'n_up',npotu,n_up_id)
         iret = nf90_def_var(ncid,'l_up',nf90_int,n_up_id,l_up_id)
      endif

      if (npotd.ne.0) then
         iret = nf90_def_var(ncid,'v_down',nf90_double,
     $                         (/n_down_id,nrp_id/),v_down_id)
      endif
      if (npotu.ne.0) then
         iret = nf90_def_var(ncid,'v_up',nf90_double,
     $                            (/n_up_id,nrp_id/),v_up_id)
      endif
c
      iret = nf90_def_var(ncid,'r',nf90_double,nrp_id,r_id)
      iret = nf90_def_var(ncid,'core_charge',nf90_double,
     $                          nrp_id,core_id)
      iret = nf90_def_var(ncid,'val_charge',nf90_double,
     $                          nrp_id,val_id)
c
      iret = nf90_enddef(ncid)
c
c     Read the meat
c
      read(2) (v(j),j=1,nrp)
      iret = nf90_put_var(ncid,r_id,v(1:nrp))
c 
c     "Down" potentials
c
      do 30 i = 1, npotd
         read(2) lo, (v(j),j=1,nrp)
         iret = nf90_put_var(ncid,l_down_id,lo,start=(/i/))
         iret = nf90_put_var(ncid,v_down_id,v(1:nrp),start=(/i,1/),
     $                       count=(/1,nrp/))
         call check(iret)
   30 continue
c
c     "Up" potentials
c
      do 35 i = 1, npotu
         read(2) lo, (v(j),j=1,nrp)
         iret = nf90_put_var(ncid,l_up_id,lo,start=(/i/))
         iret = nf90_put_var(ncid,v_up_id,v(1:nrp),start=(/i,1/),
     $                       count=(/1,nrp/))
         call check(iret)
   35 continue
c
c     Core and valence charge
c
      read(2) (v(j),j=1,nrp)
      iret = nf90_put_var(ncid,core_id,v(1:nrp))
      read(2) (v(j),j=1,nrp)
      iret = nf90_put_var(ncid,val_id,v(1:nrp))
c
      iret= nf90_close(ncid)

      contains

      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'Stopped'
      endif
      end subroutine check

      end program vps2cdf






