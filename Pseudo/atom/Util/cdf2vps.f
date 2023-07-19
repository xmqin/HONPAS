c
      program cdf2vps
c
      use netcdf

      implicit none
c
c     This program converts pseudopotential netCDF files to "VPS" 
c     format.
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
         write(0,*) 'Usage: cdf2vps netCDF_file vps_file'
         stop 
      endif
c
      call getarg(1,cdf_file)
      call getarg(2,binary_file)
c
c      open files
c
      iret = nf90_open(cdf_file,NF90_NOWRITE,ncid)
      open(unit=2,file=binary_file,form='unformatted',status='replace')
c
      title_str = ' '
      ray_str = ' '
      iret = nf90_get_att(ncid,nf90_global,'Element',nameat)
      iret = nf90_get_att(ncid,nf90_global,'Title',ray_str)
      iret = nf90_get_att(ncid,nf90_global,'Config',title_str)
      iret = nf90_get_att(ncid,nf90_global,'Correlation',corr)
      iret = nf90_get_att(ncid,nf90_global,'Relativistic',rel)
      iret = nf90_get_att(ncid,nf90_global,'Core',core)
      iret = nf90_get_att(ncid,nf90_global,'Valence_charge',zion)
      iret = nf90_get_att(ncid,nf90_global,'a_parameter',a)
      iret = nf90_get_att(ncid,nf90_global,'b_parameter',b)
c
      do i=1,7
         ptr = (i-1)*10 + 1
         title(i) = title_str(ptr:ptr+9)
      enddo
      do i=1,6
         ptr = (i-1)*10 + 1
         ray(i) = ray_str(ptr:ptr+9)
      enddo

      iret = nf90_inq_dimid(ncid,'nrp',nrp_id)
      iret = nf90_inquire_dimension(ncid,nrp_id,len=nrp)

      iret =  nf90_inq_dimid(ncid,'n_down',n_down_id)
      iret = nf90_inquire_dimension(ncid,n_down_id,len=npotd)

      iret =  nf90_inq_dimid(ncid,'n_up',n_up_id)
      if (iret .ne. nf90_noerr) then
         write(0,'(a)') 'No up potentials'
         npotu = 0
      else
         iret = nf90_inquire_dimension(ncid,n_up_id,len=npotu)
      endif

      write(2) nameat, corr, rel, core, (ray(j),j=1,6), 
     &         (title(j),j=1,7), npotd, npotu, nrp, a, b, zion

      allocate(v(nrp))

      iret = nf90_inq_varid(ncid,'l_down',l_down_id)
      iret = nf90_inq_varid(ncid,'v_down',v_down_id)
      if (npotu .ne. 0) then
         iret = nf90_inq_varid(ncid,'l_up',l_up_id)
         iret = nf90_inq_varid(ncid,'v_up',v_up_id)
      endif
c
      iret = nf90_inq_varid(ncid,'r',r_id)
      iret = nf90_inq_varid(ncid,'core_charge',core_id)
      iret = nf90_inq_varid(ncid,'val_charge',val_id)
c
      iret = nf90_get_var(ncid,r_id,v(1:nrp))
      write(2) (v(j),j=1,nrp)
c 
c     "Down" potentials
c
      do 30 i = 1, npotd
         iret = nf90_get_var(ncid,l_down_id,lo,start=(/i/))
         iret = nf90_get_var(ncid,v_down_id,v(1:nrp),start=(/i,1/),
     $                       count=(/1,nrp/))
         write(2) lo, (v(j),j=1,nrp)
   30 continue
c
c     "Up" potentials
c
      do i = 1, npotu
         iret = nf90_get_var(ncid,l_up_id,lo,start=(/i/))
         iret = nf90_get_var(ncid,v_up_id,v(1:nrp),start=(/i,1/),
     $                       count=(/1,nrp/))
         write(2) lo, (v(j),j=1,nrp)
      enddo
c
c     Core and valence charge
c
      iret = nf90_get_var(ncid,core_id,v(1:nrp))
      write(2) (v(j),j=1,nrp)
      iret = nf90_get_var(ncid,val_id,v(1:nrp))
      write(2) (v(j),j=1,nrp)
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

      end program cdf2vps






