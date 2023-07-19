#!/usr/bin/env python

from ASE.IO.PDB import WritePDB

from Siesta.read_struct import ReadStruct
import sys, os

try:
  file = sys.argv[1]
except:
  print "Usage:",sys.argv[0], "<struct file>     [ n1 n2 n3]"; sys.exit(1) 

repeat = (1,1,1)
try:
  repeat = (int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
except:
  pass

species_map = { "1" : "Mg",
                "2" : "C" ,
                "3" : "O" }

atoms = ReadStruct(file,species_map)

pdbfile=file + ".pdb"
atoms = atoms.Repeat(repeat)
WritePDB(pdbfile,atoms=atoms)
print "Wrote " + pdbfile


