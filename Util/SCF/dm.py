#!/usr/bin/env python

#
# Example of python script for processing of DMHS netcdf files.
#
# -- plot "time evolution" of DM_in and DM_out for a particular
#    orbital pair, entered as an "element number".
# 

from Numeric import *
from Scientific.IO.NetCDF import *
import biggles
import sys

#
fname=sys.argv[1]
try:
    element=int(sys.argv[2])
except:
    element=0

print "File, element:", fname, element

of = NetCDFFile(fname)
print "Dimensions:", of.dimensions

#
dmin=of.variables['dm_in'][1:,0,element]
dmout=of.variables['dm_out'][1:,0,element]

nsteps = len(dmin)
print nsteps

step = arrayrange(nsteps)

p=biggles.FramedPlot(title=fname)
p.add(biggles.Curve(step,dmin,color="red"))
p.add(biggles.Curve(step,dmout,color="blue"))
p.show()
#p.save_as_img( "gif", 400, 400, fname+".temp.gif" )
