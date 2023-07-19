#!/usr/bin/env python

import os, string, shutil
from tempfile import mkdtemp

class Siesta:

  def __init__(self,executable=None):

        self.executable="/Users/ag/bin/siesta-xlf"
	if executable is not None:
           self.executable=executable
        self.redirect_output=" > OUT"
	self.options = {}

  def SetOption(self,name,value):
        self.options[name]= value

  def run(self,atoms,out=None):
       if out is not None:
          self.redirect_output=" > " + out

       hybrid_number = 201
       self.there_are_hybrids = 0
       self.species = {}
       self.hybrids = {}
       ispec = 0
       for atom in atoms:
         name = atom.GetLabel()
         z = atom.GetAtomicNumber()
	 if not self.species.has_key(name):
           ispec = ispec + 1
           if z == 0:
             z = hybrid_number
             hybrid_number = hybrid_number + 1
             self.there_are_hybrids = 1
             self.hybrids[name] = [z,atom.valence_gs]
           self.species[name] = [ispec, z]
	
  #       filename = tempfile.mktemp() + ".fdf"
       filename = "test.fdf"
       f = open(filename,"w")
            
       f.write("NumberOfSpecies")
       f.write("%3i\n" % len(self.species))
       f.write("%block ChemicalSpeciesLabel\n")
       for i in self.species.keys():
          ispec, z = self.species[i]
          f.write("%3i %3i %4s\n" % (ispec, z, i))
       f.write("%endblock ChemicalSpeciesLabel\n")

       if self.there_are_hybrids == 1:
         f.write("%block SyntheticAtoms\n")
         for i in self.species.keys():
            ispec, z = self.species[i]
	    if z > 200:
               zdum, valgs  = self.hybrids[i]
               f.write("%3i\n" % (ispec,))
               for j in valgs[0]:
                 f.write("%3i" % j )
               f.write("\n")
               for j in valgs[1]:
                 f.write("%12.8f" % j )
               f.write("\n")
         f.write("%endblock SyntheticAtoms\n")


       # see if we have periodic  boundary conditions
       bc = atoms.GetBoundaryConditions()
       if (bc[0] or bc[1] or bc[2]):
         ucell = atoms.GetUnitCell()
         f.write("LatticeConstant 1.0 Ang\n")
         f.write("%block LatticeVectors\n")
	 for i in range(3):
           for j in range(3):
             f.write("%15.8f" % ucell[i,j])
           f.write("\n")
         f.write("%endblock LatticeVectors\n")

       f.write("NumberOfAtoms")
       f.write("%5i\n" % len(atoms))
       f.write("AtomicCoordinatesFormat Ang\n")
       f.write("%block AtomicCoordinatesAndAtomicSpecies\n")

       for atom in atoms:
         name = atom.GetLabel()
         spec = self.species[name][0]
         xyz = atom.GetCartesianPosition()
         for j in range(3):
            f.write("%15.8f" % xyz[j])
         f.write("%3i\n" % spec)
       f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
#
#      Now the one-liner options
#
       f.write("# ----- Options follow ---\n")
       for key in self.options.keys():
          f.write("%s %s\n" % (key,self.options[key]))

       f.close()

       os.system(self.executable + " <  " + filename + self.redirect_output)
       os.system(" grep FreeEng  OUT | tail -1 | awk '{print $4}' > tmp.energy")
       g = open("tmp.energy","r")
       e = g.readline()
       g.close()
       os.system("rm -f tmp.energy")

       return string.atof(e)


