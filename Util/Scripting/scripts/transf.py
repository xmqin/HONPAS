import Numeric
import LinearAlgebra
import math
import sys

def get_info(filename):

    f = open(filename)
#
#   First read unit cell
#
    cell = []
    for i in range(3):
      vector = f.readline()
      x, y, z = vector.split()
      cell.append([ float(x), float(y), float(z) ])
    acell = Numeric.array(cell)
#
#   Now the atoms
#
#   Robust code to allow for broken lines. 
#   Accumulate all the data in a big list...

    crystal = []
    natoms = int(f.readline())
    strings = []
    while 1:
       line = f.readline()
       if not line: break
       strings = strings + line.split()

#   ... and extract the information for each atom in turn
#
    pos = 0
    for a in range(natoms):
        sublist = strings[pos:pos+5]
        spindex, z , x, y, z = sublist
        x = float(x) % 1.0
        y = float(y) % 1.0
        z = float(z) % 1.0
        print spindex, "-", float(x), float(y), float(z)
        crystal.append([float(x), float(y), float(z)])
        pos = pos + 5
    atoms = Numeric.array(crystal)

    return cell, atoms
#

def toparams(cell):
    va = cell[0]
    vb = cell[1]
    vc = cell[2]
    rtodeg = 180./3.1416
    ma = math.sqrt(Numeric.dot(va,va))
    mb = math.sqrt(Numeric.dot(vb,vb))
    mc = math.sqrt(Numeric.dot(vc,vc))
    alpha = rtodeg*math.acos(Numeric.dot(vb,vc)/(mb*mc))
    beta = rtodeg*math.acos(Numeric.dot(va,vc)/(ma*mc))
    gamma = rtodeg*math.acos(Numeric.dot(va,vb)/(ma*mb))
    return  ma, mb, mc, alpha, beta, gamma

"""
Basis transformation, including possible translations of coordinates,
as given by Bplot page in the BCS
"""

try:
  filename = sys.argv[1]
except:
  print "Usage:",sys.argv[0], "<struct file> "; sys.exit(1) 

cell, atoms = get_info(filename)

print "Old cell:", cell
print toparams(cell)

#
# The transformation matrix is hardwired in this version.
#
matrix = Numeric.array([[ 0., -1., 1.], [ 0., 0., 1.], [ -2., 2. , -1.]])
emat = Numeric.array([[ 0., -1., 1., 0.28488], 
                      [ 0., 0., 1.,0.12343], 
                      [ -2., 2. , -1., 0.0],
                      [ 0.0, 0.0, 0.0, 1.0] ])

invemat = LinearAlgebra.inverse(emat)

# Extend matrix to deal with translations, in the crystallographic way
#
print atoms
natoms = atoms.shape[0]
ones= 1 + Numeric.zeros((natoms,1))
print ones.shape
extatom = Numeric.concatenate((atoms,ones),axis=1)
print extatom

print "---------------"

newcell = Numeric.dot(matrix,cell)
print newcell
print toparams(newcell)

for i in range(natoms):
  xold = extatom[i]
  print Numeric.dot(invemat,xold)



