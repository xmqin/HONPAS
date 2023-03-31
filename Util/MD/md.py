#!/usr/bin/env python

#
# Plot temperature from MD.nc file
#
from Numeric import *
from Scientific.IO.NetCDF import *
import biggles
import sys

fname=sys.argv[1]

step=1         ## Stride over array

of = NetCDFFile(fname)
print "Dimensions:", of.dimensions

temp=of.variables['temp'][:]
index=arange(0,len(temp),step)    # To hold step number

p=biggles.FramedPlot(title=fname+" Temperature")
p.add(biggles.Curve(index[:],temp[0::step]))
p.show()
