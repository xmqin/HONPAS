#!/usr/bin/env python

#
# bond.py : An example of Siesta scripting
#
#           Compute and plot the Energy vs. bond length curve
#           for the water molecule
#
import biggles
from Siesta.Interface import Atom      
from ASE import ListOfAtoms
from Siesta.siesta import Siesta
import math, os, shutil, urllib
import Numeric as num

p=biggles.FramedPlot(title="Energy/bond")

#
# Standard bond length (units: Angstrom)
#
d0 = math.sqrt(0.757**2 +  0.586**2)

# Range of bond lengths areoun the presumed minimum
bonds = num.arange(0.8*d0,1.3*d0,0.1)
energy = 0.0*bonds
theta=0.5*(180.-105.0)*math.pi/180.  # Angle with the x axis, in radians
print bonds

a = Siesta(executable="$HOME/bin/siesta-2.4-optim")          # Initialize object

# create a work subdirectory
orig_dir = os.getcwd()
dir = "bond_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"
urllib.urlretrieve(URLbase + "H.psf", "H.psf") 
urllib.urlretrieve(URLbase + "O.psf", "O.psf") 

for i in range(len(bonds)):
  d = bonds[i]
  print i, d
  atoms = ListOfAtoms(
          [Atom('H', (d*math.cos(theta),d*math.sin(theta),0.0),magmom=1),
           Atom('H', (-d*math.cos(theta),d*math.sin(theta),0.0),magmom=1),
           Atom('O', (0, 0, 0),magmom=1)
          ])
             #                     cell=(4, 4, 4), periodic=1)

  energy[i] = a.run(atoms)   # Run Siesta and get the (free)energy

#
# Plot
#
p.add(biggles.Curve(bonds[:],energy[:]))
p.show()
os.chdir(orig_dir) 
