#!/usr/bin/env python

#
# Example of python script for processing of DM netcdf files.
#
#  -- plot histogram of two consecutive Hs

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
    iter=-2

print "File, nbins:", fname, nbins

of = NetCDFFile(fname)
print "Dimensions:", of.dimensions

h1=of.variables['h'][iter,0,:]
h2=of.variables['h'][iter+1,0,:]
print len(h1)
hh1=Histo(h1[:],nbins)
print hh1.min, hh1.max
hh2=Histo(h2[:],nbins)
print hh2.min, hh2.max

hdiff=h2-h1
hhdiff=Histo(hdiff[:],nbins)
print hhdiff.min, hhdiff.max

g = Gnuplot.Gnuplot(debug=1)
g('set data style lines') # give gnuplot an arbitrary command
d1 = Gnuplot.Data(hh1.bins,hh1.h)
d2 = Gnuplot.Data(hh2.bins,hh2.h)
g.plot(d1,d2)

raw_input('Please press return to continue...\n')

g = Gnuplot.Gnuplot(debug=1)
g('set data style lines') # give gnuplot an arbitrary command
d1 = Gnuplot.Data(hhdiff.bins,hhdiff.h)
g.plot(d1)

raw_input('Please press return to continue...\n')

#p=biggles.FramedPlot(title=fname)
#p.add(biggles.Curve(h.bins,h.h))

#p.show()
#p.save_as_img( "gif", 400, 400, fname+".temp.gif" )
