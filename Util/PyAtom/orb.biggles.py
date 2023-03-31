#!/usr/bin/env python

#
from Numeric import *
from Scientific.IO.NetCDF import *
from biggles import *
import sys


fname=sys.argv[1]
of = NetCDFFile(fname)

norbs=of.dimensions['norbs']    # Number of orbitals
ntb=of.dimensions['ntb']

nrows, rest = norbs/2, norbs%2
if rest==1: nrows = nrows+1  # If odd number, we need one more row...

t=Table(nrows,2,cellspacing=0.02)  
t.title="ORBITALS from " + fname

for i in range(norbs):
	row, column  = i/2 , i%2

        L = of.variables['orbnl_l'][i]
        Z = of.variables['orbnl_z'][i]

        delta=of.variables['delta'][i]

	p=FramedPlot()
	p.title="Orbital l=%2d, z=%2d" % (L,Z)
        r = arange(1,ntb+1)*float(delta)
	p.add(Curve(r,of.variables['orb'][i,:]))
	t.set(row,column,p)

t.show()
#
t.save_as_eps( fname+".orbs.eps" )
t.save_as_img( "gif", 400, 400, fname+".orbs.gif" )
t.save_as_img( "svg", 400, 400, fname+".orbs.svg" )
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
