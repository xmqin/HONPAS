#!/usr/bin/env python

import os, string, sys
import Numeric

class Siesta:

  """ A modified version (it should be a subclass) which
    uses the SiestaAsServer features.  """

  def __init__(self,executable=None):

        self.executable="/Users/ag/bin/siesta-xlf"
	if executable is not None:
           self.executable=executable
        self.redirect_output=" > OUT"
#	self.options = { }
	self.options = { "MD.TypeOfRun" : "forces" ,
                         "SystemLabel"  : "siesta" }

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

       os.mkfifo("siesta.coords")
       os.mkfifo("siesta.forces")
       os.system(self.executable + " <  " + filename + self.redirect_output + "&")
       print "Launched Siesta server"
       self.f = open("siesta.coords","w")
       self.toserver = self.f.write
       self.toserver("wait\n")

  def get_forces(self,atoms):
       self.toserver("begin_coords\n")
       self.toserver("Ang\n")
       self.toserver("eV\n")

       ucell = atoms.GetUnitCell()
       for i in range(3):
          for j in range(3):
             self.toserver("%15.8f" % ucell[i,j])
          self.toserver("\n")

       self.toserver("%d\n" % len(atoms))

       for atom in atoms:
         xyz = atom.GetCartesianPosition()
         for i in range(3):
           self.toserver("%15.8f " % xyz[i])
         self.toserver("\n")
       self.toserver("end_coords\n")
       self.f.flush()
#      print "Sent coordinates"

#      Open pipe for reading energy, stress, and forces
#
       self.fromserver = open("siesta.forces","r").readline
       msg = self.fromserver().split()
       if  msg[0] != "begin_forces" : 
             print "Not begin_forces"
             sys.exit(1)
       energy = float(self.fromserver().split()[0])

       stress = Numeric.zeros((3,3),Numeric.Float)
       for i in range(3):
         st = self.fromserver().split()
         for j in range(3):
           stress[i,j] = float(st[j])

       na = int(self.fromserver().split()[0])
       if na != len(atoms): 
           print "number of atoms", na
           sys.exit(1)
       forces = Numeric.zeros((na,3),Numeric.Float)  # Note index order
       for i in range(na):
         force = self.fromserver().split()
         for j in range(3):
           forces[i,j] = float(force[j])

       msg = self.fromserver().split()
       if  msg[0] != "end_forces" :
             print "Not end_forces"
             sys.exit(1)

       return (energy, stress, forces)

  def stop(self):
       self.toserver("quit\n")


