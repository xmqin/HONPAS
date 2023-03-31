#!/usr/bin/env python

#
import Numeric 
from Scientific.IO.NetCDF import *
from biggles import *
import sys


fname=sys.argv[1]
of = NetCDFFile(fname)

nkbs=of.dimensions['nkbs']            # Number of KB projectors
ntb=of.dimensions['ntb']

nrows, rest = nkbs/2, nkbs%2
if rest==1: nrows = nrows+1  # If odd number, we need one more row...

t=Table(nrows,2,cellspacing=0.02)  
t.title="KB PROJS from " + fname

for i in range(nkbs):
	row, column  = i/2 , i%2

        L = of.variables['pjnl_l'][i]
        N = of.variables['pjnl_n'][i]

	delta=of.variables['kbdelta'][i]

	p=FramedPlot()
	p.title="KB proj l=%2d, n=%2d" % (L,N)
	r = Numeric.arange(1,ntb+1)*float(delta)
	p.add( Curve(r,of.variables['proj'][i,:]) )
	t.set(row,column,p)

t.show()

#
t.save_as_eps( fname+".kbs.eps" )
t.save_as_img( "gif", 500, 500, fname+".kbs.gif" )
t.save_as_img( "svg", 400, 400, fname+".kbs.svg" )



# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
