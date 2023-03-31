#!/usr/bin/env python

#
# Example of python script for processing of DM netcdf files.
#
#  -- plot histogram of last DM_out elements

from Histo import *
from Numeric import *
from Scientific.IO.NetCDF import *
import biggles
import sys
import MLab

#
fname=sys.argv[1]
try:
    nbins=int(sys.argv[2])
except:
    nbins=0

print "File, nbins:", fname, nbins

of = NetCDFFile(fname)
print "Dimensions:", of.dimensions

last_dm=of.variables['dm_out'][-1,0,:]
print len(last_dm)
h=Histo(last_dm[:],nbins)

p=biggles.FramedPlot(title=fname)
p.add(biggles.Curve(h.bins,h.h))

p.show()
#p.save_as_img( "gif", 400, 400, fname+".temp.gif" )
