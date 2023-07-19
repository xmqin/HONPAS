import amoeba
import math

from Siesta.Interface import Atom
from Siesta.siesta import Siesta
from ASE import ListOfAtoms

import urllib 
import os, shutil 

"""

Example of use of the amoeba module to drive a Simplex basis-set
optimization

"""
def afunc(var,data=None): 
  rcs, rcp  = var[0], var[1]
  f=open("base.fdf","w")
  f.write("%block PAO.Basis\n  Ga  2\n")
  f.write(" 0   1\n")
  f.write("%10.5f\n1.0\n" % (rcs,) )
  f.write("1   1\n")
  f.write("%10.5f\n1.0\n" % (rcp,) )
  f.write("%endblock PAO.Basis\n")
  f.close()
             
  atoms = ListOfAtoms([Atom('Ga', (0.0,0.0,0.0),magmom=0),
                       Atom('Ga', (2.0,0.0,0.0),magmom=0)])

  a = Siesta(executable="/Users/ag/bin/siesta-xlf")

  a.SetOption("%include","base.fdf")
  energy = a.run(atoms)

  print rcs, rcp, energy
  return  -energy


URLbase = "http://fisica.ehu.es/ag/siesta-psffiles/"

# create a work subdirectory
orig_dir = os.getcwd()
dir = "simplex_work"
if os.path.isdir(dir): # does dir exist? 
  shutil.rmtree(dir) # yes, remove old directory 
os.mkdir(dir) # make dir directory 
os.chdir(dir) # move to dir 

urllib.urlretrieve(URLbase + "Ga.lda.psf", "Ga.psf") 

print amoeba.amoeba([4.25,4.25],[2.5,2.5],afunc,xtolerance=5.0e-2)
