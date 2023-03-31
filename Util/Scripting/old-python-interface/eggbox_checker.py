#----------------------------------------------------
# Eggbox checker
# Alberto Garcia, August 2009
#
# Computes the energy and force as a function of displacement for a single
# atom in a cubic box.
#
# Usage:
#       python eggbox_checker PseudoFile.psf [Cutoff] [basis file]
#
#   The PseudoFile.psf file and the optional basis file can be
#   anywhere in the filesystem. The script will create a temporary
#   directory and do its work in it.
#
# If you move this script, you can try something along
# the following lines to be able to find the relevant modules
# (or maybe just make a symbolic link in the current directory
#  to the Siesta script directory)
#
import sys
sys.path.append("$HOME/lib/python/Siesta")
#

from Siesta.Interface import Atom, Crystal
from Siesta.efs import Siesta_efs as Siesta
import Numeric as Num
import os, urllib, shutil, math, sys

try:
  pseudo_file = sys.argv[1]
except:
  print "Usage:",sys.argv[0], "PseudoFile.psf [Cutoff] [basis file]"
  sys.exit(1) 

try:
  cutoff = float(sys.argv[2])
  print "Using a cutoff of ", cutoff
except:
  cutoff = 200.0
  print "Using default cutoff of ", cutoff
#
dx = 0.529*math.pi/math.sqrt(cutoff)
print "Predicted eggbox period, based on the nominal Meshcutoff: ", dx
#
#
try:
  basis_file = sys.argv[3]
  print "Using basis file: ", basis_file
except:
  open("tmp_basis.fdf", "w").write("#Dummy basis info")
  basis_file = "tmp_basis.fdf"

basis_file = os.path.abspath(basis_file)

#
#  Get symbol
#
f=open(pseudo_file,"r")
symbol=f.readline().split()[0]
f.close()
label = os.path.basename(pseudo_file[:-len(".psf")])

pseudo_file = os.path.abspath(pseudo_file)
print "Processing ", pseudo_file

print "Element label, symbol: ", label, symbol

a=10.0              # Default box length in Ang
cell=Num.array([(a, 0, 0),
      (0, a, 0),
      (0, 0, a)])

a = Siesta(executable="$HOME/bin/siesta-test")          # Initialize object

a.SetOption("DM.NumberPulay","3")
a.SetOption("DM.UseSaveDM","T")
a.SetOption("UseSaveDM","T")
a.SetOption("Meshcutoff", str(cutoff)+" Ry")
#
# Make sure that your fdf has been patched to allow long %include filenames
#
a.SetOption("%include",basis_file)

m=20             # Number of displacements to try
disps=Num.zeros(m,Num.Float)
energy=Num.zeros(m,Num.Float)

output_file=os.path.abspath(label+".eggbox."+str(cutoff))
print "The output will be in file: ", output_file

# create a work subdirectory
orig_dir = os.getcwd()
import tempfile
workdir = tempfile.mkdtemp(prefix="eggbox.work.",dir=".")
print "Working directory: ", workdir
os.chdir(workdir) # move to dir
#
os.system("cp -p " + pseudo_file + " .")
os.system("cp -p " + basis_file + " .")
g=open(output_file,"w")

for i in range(len(disps)):
  d = i*dx/(m-1)
  # could generalize the displacement here (1,1,1), (1,0,1)...
  atoms=Crystal([Atom(symbol=symbol,label=label,position=(0,0,d))])
  atoms.SetUnitCell(cell,fix=True)
  atoms.SetBoundaryConditions(periodic=True)
  energy[i], force, stress = a.run(atoms)   # Run Siesta - get the (free)energy
  print  i, d, energy[i]
  fforce = [ float(j) for j in force[0]]
  g.write(" %10.6f %20.6f  %15.6f%15.6f%15.6f\n" % (float(d), float(energy[i]),  fforce[0], fforce[1], fforce[2] ))

#
g.close()
os.chdir(orig_dir) 

