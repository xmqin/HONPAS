#!/usr/bin/env python

from Siesta.Interface import  Atom, Crystal
from Siesta.calculator import SiestaCalculator
from ASE.Dynamics.QuasiNewton import QuasiNewton

import urllib 
import os, shutil 

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"

atoms = Crystal([Atom('H', (0.757,0.586,0.),label="H_test"),
                     Atom('H', (-0.757,0.586,0.),magmom=1),
                     Atom('O', (0, 0, 0),magmom=1)],
                     cell=(6.0, 6.0, 6.0), periodic=1)

# create a work subdirectory
orig_dir = os.getcwd()
dir = "ase_newton_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

urllib.urlretrieve(URLbase + "H.psf", "H_test.psf") 
urllib.urlretrieve(URLbase + "H.psf", "H.psf") 
urllib.urlretrieve(URLbase + "O.psf", "O.psf") 

b = SiestaCalculator(executable="/Users/ag/bin/siesta-xlf")
b.SetOption("DM.Tolerance","0.001")
atoms.SetCalculator(b)

energy = atoms.GetPotentialEnergy()
forces = atoms.GetCartesianForces()

print "initialEnergy ", energy
maxForce = max(abs(forces.flat))
print "maxForce ",maxForce
        
qn = QuasiNewton(atoms,fmax=0.01)
qn.Converge()

finalEnergy = atoms.GetPotentialEnergy()
print "finalEnergy",finalEnergy
    
finalForces = atoms.GetCartesianForces()
maxForce = max(abs(finalForces.flat))
print "maxForce",maxForce

b.stop()

os.chdir(orig_dir) 
