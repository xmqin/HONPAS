#!/usr/bin/env python

import os, string, shutil
from tempfile import mkdtemp
import Numeric

def create_fdf_struct_file(atoms,filename="test.fdf"):
       """ Create the species and structural part in fdf
        form, starting with an Atoms object
       """
    
       hybrid_number = 201
       there_are_hybrids = 0
       species = {}
       hybrids = {}
       ispec = 0
       for atom in atoms:
         name = atom.GetLabel()
         z = atom.GetAtomicNumber()
	 if not species.has_key(name):
           ispec = ispec + 1
           if z == 0:
             z = hybrid_number
             hybrid_number = hybrid_number + 1
             there_are_hybrids = 1
             hybrids[name] = [z,atom.valence_gs]
           species[name] = [ispec, z]
	
       f = open(filename,"w")
            
       f.write("NumberOfSpecies")
       f.write("%3i\n" % len(species))
       f.write("%block ChemicalSpeciesLabel\n")
       for i in species.keys():
          ispec, z = species[i]
          f.write("%3i %3i %4s\n" % (ispec, z, i))
       f.write("%endblock ChemicalSpeciesLabel\n")

       if there_are_hybrids == 1:
         f.write("%block SyntheticAtoms\n")
         for i in species.keys():
            ispec, z = species[i]
	    if z > 200:
               zdum, valgs  = hybrids[i]
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
         spec = species[name][0]
         xyz = atom.GetCartesianPosition()
         for j in range(3):
            f.write("%15.8f" % xyz[j])
         f.write("%3i\n" % spec)
       f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")


def dump_fdf_options(options, filename):
       #
       #      Now the one-liner options
       #   The file must exist
       #
       f = open(filename,"a")
       f.write("# ----- Options follow ---\n")
       for key in options.keys():
          f.write("%s %s\n" % (key,options[key]))

       f.close()
       

def get_efs(filename="FORCE_STRESS"):
    """
    Return energy, forces, and stress from file

    """

    f = open(filename)
    energy = float(f.readline())
    #
    #  Read stress
    #
    stress = Numeric.zeros((3,3),Numeric.Float)
    for i in range(3):
      st = f.readline().split()
      for j in range(3):
             stress[i,j] = float(st[j])
#
#   Now the atoms and forces
#
    na = int(f.readline())

    forces = Numeric.zeros((na,3),Numeric.Float)  # Note index order
    for i in range(na):
      isa, iza, x, y, z, label = f.readline().split()
      force = x, y, z
      for j in range(3):
          forces[i,j] = float(force[j])

    f.close()
    
    return energy, forces, stress
#
       
class Siesta_efs:

  def __init__(self,executable=None):

        self.executable="$HOME/bin/siesta-2.4-optim"
	if executable is not None:
           self.executable=executable
        self.redirect_output=" > OUT"
	self.options = {}
        self.filename = "test.fdf"

  def SetOption(self,name,value):
        self.options[name]= value

    
  def run(self,atoms,out=None):
       Ang = 1.0 / 0.529177
       eV  = 1.0 / 13.60580

       if out is not None:
          self.redirect_output=" > " + out

       create_fdf_struct_file (atoms,filename=self.filename)
       dump_fdf_options (options=self.options, filename=self.filename)

       #
       #      Execute
       #
       os.system(self.executable + " <  " + self.filename + self.redirect_output)

       #
       # Direct reading of the FORCE_STRESS file
       #

       energy, forces, stress = get_efs()
       return energy/eV , forces*Ang/eV, stress*Ang**3/eV

       



