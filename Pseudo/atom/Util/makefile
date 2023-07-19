#
include arch.make

#
vps2cdf: vps2cdf.f modules
	$(FC) $(FFLAGS) -o $@  $@.f  $(NETCDF_LIB)
cdf2vps: cdf2vps.f modules
	$(FC) $(FFLAGS) -o $@  $@.f  $(NETCDF_LIB)
#
clean:
	rm -f *.o vps2cdf cdf2vps









