#!/usr/bin/env python

from Siesta.Interface import  Atom, Crystal
from Siesta.siesta import Siesta

import urllib 
import os, shutil 

URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"

atoms = Crystal([Atom('H', (0.757,0.586,0.),label="H_test"),
                     Atom('H', (-0.757,0.586,0.),magmom=1),
                     Atom('O', (0, 0, 0),magmom=1)])
             #                     cell=(4, 4, 4), periodic=1)

# create a work subdirectory
orig_dir = os.getcwd()
dir = "h2o_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

urllib.urlretrieve(URLbase + "H.psf", "H.psf") 
urllib.urlretrieve(URLbase + "H.psf", "H_test.psf") 
urllib.urlretrieve(URLbase + "O.psf", "O.psf") 

a = Siesta(executable="$HOME/bin/siesta-2.4-optim")
energy = a.run(atoms,out="OUT")

print "The energy is: ", energy

os.chdir(orig_dir) 
