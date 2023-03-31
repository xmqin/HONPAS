#!/usr/bin/env python

from Siesta.Interface import  Atom, Crystal
from Siesta.calculator import SiestaCalculator

import urllib 
import os, shutil 

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"

atoms = Crystal([Atom('H', (0.757,0.586,0.)),
                     Atom('H', (-0.757,0.586,0.),magmom=1),
                     Atom('O', (0, 0, 0),magmom=1)],
                     cell=(6.0, 6.0, 6.0), periodic=1)

# create a work subdirectory
orig_dir = os.getcwd()
dir = "ase_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

urllib.urlretrieve(URLbase + "H.psf", "H.psf") 
urllib.urlretrieve(URLbase + "O.psf", "O.psf") 

b = SiestaCalculator(executable="$HOME/bin/siesta-2.4-optim")

b.SetOption("DM.Tolerance","0.005")
atoms.SetCalculator(b)

energy = atoms.GetPotentialEnergy()
forces = atoms.GetCartesianForces()
stress = atoms.GetStress()

print "The energy is: ", energy
print "Forces ", forces
print "Stress ", stress
stress = atoms.GetStress()
print "Stress ", stress

os.chdir(orig_dir) 
