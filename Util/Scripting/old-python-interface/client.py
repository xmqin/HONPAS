#!/usr/bin/env python

from Siesta.Interface import  Atom, Crystal
from Siesta.server import Siesta

import urllib 
import os, shutil 

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"


# create a work subdirectory
orig_dir = os.getcwd()
dir = "client_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

urllib.urlretrieve(URLbase + "H.psf", "H_test.psf") 
urllib.urlretrieve(URLbase + "H.psf", "H.psf") 
urllib.urlretrieve(URLbase + "O.psf", "O.psf") 


a = Siesta(executable="/Users/ag/bin/siesta-xlf")
a.run(atoms,out="OUT")

blist = [ 0.4, 0.50,  0.60,  0.70,  0.80,  0.90]
energies = []
forces = []
energies = []
for d in blist:
   atoms = Crystal([Atom('H', (d,0.586,0.),label="H_test"),
                     Atom('H', (-d,0.586,0.),magmom=1),
                     Atom('O', (0, 0, 0),magmom=1)],
                     cell=(6.0, 6.0, 6.0), periodic=1)

   e, stress, fs = a.get_forces(atoms)
   energies.append(e)
   forces.append(fs)

print "Energies ", energies
print "Forces ", forces

a.stop()

