#!/usr/bin/env python

#
# cutoff.py : An example of Siesta scripting
#
#           Compute and plot the Energy vs. cutoff curve
#
import biggles
from Siesta.Interface import Atom      
from ASE import ListOfAtoms
from Siesta.efs import Siesta_efs as Siesta
import math, os, shutil, urllib
import Numeric as Num

p=biggles.FramedPlot(title="Energy/cutoff")

#
# Range of cutoffs
#
cutoffs = Num.arange(100.0,801.0,100.0)
energy = 0.0*cutoffs

a = Siesta(executable="$HOME/bin/siesta-2.6.9")          # Initialize object

# create a work subdirectory
orig_dir = os.getcwd()
dir = "cutoff_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

open("basis.fdf", "w").write("""
%block kgrid_Monkhorst_Pack
   7   0    0    0.5
   0   7    0    0.5
   0    0   7    0.5
%endblock kgrid_Monkhorst_Pack

DM.NumberPulay         3 
DM.MixingWeight        0.05

WriteDM F 

%Block PAO.Basis
Pb   3     -0.30296
 n=6   0   2   E     3.76509     2.94865
     5.41509     4.89944
     1.00000     1.00000
 n=6   1   2   E     2.50435     0.86601
     6.12615     5.62330
     1.00000     1.00000
 n=6   2   1   E   135.64896     4.82387
     5.14075
     1.00000
%EndBlock PAO.Basis
""")

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"
urllib.urlretrieve(URLbase + "Pb.psf", "Pb.psf") 

cell = Num.array([ [0.5,0.5,0.0], [0.5,0.0,0.5], [0.0,0.5,0.5]])
cell = 4.89*cell
atoms = ListOfAtoms([Atom('Pb', (0.0,0.0,0.0))], cell=cell, periodic=1)

for i in range(len(cutoffs)):
  cutoff = cutoffs[i]
  print i, cutoff
  a = Siesta(executable="$HOME/bin/siesta-2.6.9")          # Initialize object
###  a.SetOption("FilterCutoff"," 100.0 Ry")               # Optional Filtering
  a.SetOption("Meshcutoff", str(cutoff)+" Ry")
  a.SetOption("%include"," basis.fdf")
  energy[i], dum, dum2 = a.run(atoms)   # Run Siesta and get the (free)energy

#
# Plot
#
p.add(biggles.Curve(cutoffs[:],energy[:]))
p.show()
os.chdir(orig_dir) 
