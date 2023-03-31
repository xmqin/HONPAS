# ncdf #

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&business=NGNU2AA3JXX94&lc=DK&item_name=Papior%2dCodes&item_number=codes&currency_code=EUR&bn=PP%2dDonationsBF%3abtn_donate_SM%2egif%3aNonHosted)


A library to ease the creation of NetCDF files using a standard interface
without the need to keep track of file, variable and attribute ID's.

## Usage ##

It has proven quite useful as the necessity to switch between MPI and non-MPI [NetCDF][netcdf] is extremely easy.

## Downloading and installation ##

Installing ncdf requires a download of the library 
hosted at [github](https://github.com/) at [ncdf@git].

Note that ncdf depends on the dictionary library [fdict@git].

If `fdict` is not already installed you may need to add the submodule

    git submodule init
    git submodule update

Extract and create an `setup.make` file for compilation, a minimal
`setup.make` file can look like this

	FC=gfortran
	FFLAGS = -g
	LDFLAGS = -L<path-to-netcdf> -lnetcdff -lnetcdf

Type `make` and a library called `libncdf.a` is created.  
Subsequently the installation may be performed by:

    make PREFIX=/papth/to/ncdf install

which installs the required files (modules and libraries) to the folder.

To use the library you need to add include statements for the
modules as well as linking to the program.

To link ncdf to your program the following can be used in a `Makefile`

    FDICT_PATH  = /path/to/fdict/parent
    FDICT_LIBS  = -L$(FDICT_PATH) -lfdict
    FDICT_INC   = -I$(FDICT_PATH)

    NCDF_PATH  = /path/to/ncdf/parent
    NCDF_LIBS  = -L$(NCDF_PATH) -lncdf
    NCDF_INC   = -I$(NCDF_PATH)



### Simple API discussion ###

The compilation of the ncdf library is performed through a couple of 
preprocessor flags.
They are all prefixed with `NCDF_` to not interfere with any custom flags.

- MPI, allows for communicators and should be used for any code which
  uses MPI.
  It does not requires to have either `CDF4`, nor the `PCDF` flag.
  Flag: `MPI`
- `NCDF_4` allow the NetCDF 4 API for compression and the MPI layer with 
  the MPI flag.
  Flag: `NCDF_4`

#### Coding usage ####

To ease the conversion back and forth between the NetCDF API and the ncdf API the names are very
similar.

An easy remembering rule is that "nf90" has changed to "ncdf".

Here we provide a list which shows the details of the name-conversions
 - `ncdf_open`           (`nf90_open`)
 - `ncdf_close`          (`nf90_close`)
 - `ncdf_enddef`         (`nf90_enddef`)
 - `ncdf_redef`          (`nf90_redef`)
 - `ncdf_sync`           (`nf90_sync`)
 - `ncdf_inq`            (`nf90_inquire`)
 - `ncdf_def_dim`        (`nf90_def_dim`)
 - `ncdf_inq_dim`        (`nf90_inq_dim`, `nf90_inq_dimid`, `nf90_inquire_dimension`)
 - `ncdf_rename_dim`     (`nf90_rename_dim`)
 - `ncdf_def_var`        (`nf90_def_var`)
 - `ncdf_inq_var`        (`nf90_inq_var`, `nf90_inq_varid`, `nf90_inquire_variable`)
 - `ncdf_put_att`        (`nf90_put_att`)
 - `ncdf_put_gatt`       (`nf90_put_att(..,nf90_global,...)`)
 - `ncdf_get_att`        (`nf90_get_att`)
 - `ncdf_get_gatt`       (`nf90_get_att(..,nf90_global,...)`)
 - `ncdf_del_att`        (`nf90_del_att`)
 - `ncdf_del_gatt`       (`nf90_del_att(..,nf90_global,...)`)

Basically there was a couple of issues that was troubling the author.

1. Using inq_dim, inq_dimid, inquire_dimension didn't seem obvious to me.
   Hence everything is fetched from one inquire routine. (inq_dim)
2. The same thing applied for the inq_var, inq_varid and inquire_variable
3. Global variables I thought should be able to be denoted specifically.
   Simply because of code clarity. However, the old method of `NF90_GLOBAL`
   is still allowed (basically the gatt routines are wrappers with the 
   `NF90_GLOBAL` argument).
4. For consistency I have renamed everything which uses the inquire suffix to inq

Here is the simple example:

    use netcdf_ncdf

    type(hNCDF) :: nf

    call ncdf_create(nf,'Test.nc')
    call ncdf_def_dim(nf,'MD',100)
    call ncdf_def_var(nf,'E',NF90_DOUBLE,(/'MD'/))
    call ncdf_put_var(nf,'E',E)
    call ncdf_close(nf)



#### TODO ####

- Consider changing the parallel logical to an integer to have several layers of interactions
- Add to the original NetCDF the `VAR_FILL` routines (they are only defined in the INTEGER part)
- Add an ncdf_assert (which can take several options)
  1. check that certain variables exists
  2. check dimension sizes...
- Add a ncdf_dim function which returns a dimension array
  with correct formatting... `ncdf_dim('a', 'oesntuh') > (/'a     ','oesntuh'/)`
- Scour the NetCDF manual and add missing stuff
- Fully determine whether the dictionary is the right way to go. It appears to be 
  a great solution for the character attributes, however, it could be extended.



## Contributions, issues and bugs ##

I would advice any users to contribute as much feedback and/or PRs to further
maintain and expand this library.

Please do not hesitate to contribute!

If you find any bugs please form a [bug report/issue][issue].

If you have a fix please consider adding a [pull request][pr].


## License ##

The ncdf license is [LGPL][lgpl], see the LICENSE file.

## Thanks ##

A big thanks goes to Alberto Garcia for contributing ideas and giving
me bug reports. Without him the interface would have been much more
complex!

<!---
Links to external and internal sites.
-->
[netcdf]: http://www.unidata.ucar.edu/software/netcdf/][NetCDF
[ncdf@git]: https://github.com/zerothi/ncdf
[fdict@git]: https://github.com/zerothi/fdict
<!-- [ncdf-doc]: https://github.com/zerothi/ncdf/wiki -->
[issue]: https://github.com/zerothi/ncdf/issues
[pr]: https://github.com/zerothi/ncdf/pulls
[lgpl]: http://www.gnu.org/licenses/lgpl.html
