# Guide to STM/STS simulation in Siesta

This is work in progress. Ideas for further features are welcome.

The basic premise is that the simulation of STM images is done in the
Tersoff-Hamann approximation, using the Local Density of States
(LDOS), which can be seen as the charge density computed with
wavefunctions in a certain energy window: ldos=ldos(x,y,z;{Emin,Emax})

For STS, the energy itself is a variable: ldos=ldos(x,y,z,E), and the
spectra can be obtained as dI/dE at each point, where I (the
"current") is proportional to the ldos. A broadening function (currently
only a gaussian, but others are in the works) is applied to the data.

LDOS information is stored in the same kind of file used by Siesta to
represent grid magnitudes (charge density, potential, etc). These
files can be handled by the utilities in Util/Grid, in particular they
can be read and written by the routines in the module m_gridfunc.

An LDOS file can be produced in several ways:

  * By Siesta itself, using the LDOS options detailed in the manual.
    The LDOS file will contain LDOS(x,y,z) in the same grid used by Siesta
    in that calculation (i.e., it will depend on the MeshCutoff used).
    This LDOS function is *periodic*.  

  * By the program wfs2ldos (also known as 'ol-stm'), which processes
    a wavefunction file from Siesta (in WFSX format), plus some other
    files (see the manual for ol-stm), to generate LDOS information in
    an energy window (or, in STS mode, as a function of energy) in a
    set of planes.  The set can contain a single plane, or
    several. The LDOS file produced (with extension '.STM.LDOS') will
    be in the same 'gridfunc' format, but in this case, in general,
    the LDOS function *will not* be periodic.

    The wfs2ldos/ol-stm program can optionally perform an intermediate step
    of projection of the wave-functions from a reference plane and
    into the vacuum, before computing the LDOS. (See the manual for
    ol-stm/wfs2ldos). This is completely transparent to the user.
    
The information in an LDOS file can be processed in several ways:

   * If the file contains information for STM-image mode, the program
     'plstm' can extract a 2D section in "constant-height" or
     "constant-current" mode, optionally projected on spin components
     (see the header/manual for plstm, and note that non-collinear and
     spin-orbit modes are supported).  The 2D section is ready to be
     plotted by gnuplot. Implementations of other post-processing
     options are welcome.

   * If the file contains STS information, the energy variable is
     represented along the 'spin' dimension of the 'gridfunc' file. It
     can be argued that this is an abuse of the format, but it is
     extremely convenient, because it simplifies the coding and the
     analysis tools. Note that the "LDOS" is then restricted to the
     "charge" mode (what is implied by the 'q' spin option in plstm).
     (Again, if spin-resolved STS is needed, more coding is required
     -- suggestions and merge requests welcome).

     STS information, as produced by wfs2ldos in STS mode, is in a
     file with extension '.STS', and can be processed by ad-hoc tools
     based on functionality in Util/Grid. For example, one might be
     interested in plotting the STS spectra at a given point, or a
     given line. Then, the "value extractors" in Util/Grid can be
     leveraged to produce the given values (keeping the energy as an
     extra variable).  These kinds of tools are not written yet. The
     only available demonstrator is a simple program ('plsts') that
     gets the I(E) information at a point, and also produces
     information for a contour plot of the LDOS at a particular energy
     in the file. (Energy data is kept in an auxiliary file with
     extension .STS_AUX).

     It is also possible to convert the gridfunc files to netCDF
     format and use existing netCDF visualizers, or scripts in Python or
     other languages, to get the job done, but this is unexplored.

Alberto Garcia, June 2019

    
    