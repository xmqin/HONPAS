#!/usr/bin/env python

#
# Example of python script for processing of DM netcdf files.
#
#  -- plot histogram of last DM_out elements

from Histo import *
from Numeric import *
from Scientific.IO.NetCDF import *
##import biggles
import Gnuplot, Gnuplot.funcutils
import sys
import MLab

#
fname=sys.argv[1]
try:
    nbins=int(sys.argv[3])
except:
    nbins=100
try:
    iter=int(sys.argv[2])
except:
    iter=-1

print "File, nbins:", fname, nbins

of = NetCDFFile(fname)
print "Dimensions:", of.dimensions

dmin=of.variables['dm_in'][iter,0,:]
print len(dmin)
hin=Histo(dmin[:],nbins)
print hin.min, hin.max

dmout=of.variables['dm_out'][iter,0,:]
print len(dmout)
hout=Histo(dmout[:],nbins)
print hout.min, hout.max

g = Gnuplot.Gnuplot(debug=1)
g('set data style linespoints') # give gnuplot an arbitrary command
d1 = Gnuplot.Data(hin.bins,hin.h)
d2 = Gnuplot.Data(hout.bins,hout.h)
g.plot(d1,d2)

raw_input('Please press return to continue...\n')

#p=biggles.FramedPlot(title=fname)
#p.add(biggles.Curve(h.bins,h.h))

#p.show()
#p.save_as_img( "gif", 400, 400, fname+".temp.gif" )
