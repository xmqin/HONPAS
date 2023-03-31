#!/usr/bin/python

"""
  Struct_to_bplot: Convert the information in a Struct file
  to a format compatible with the BPLOT program in the 
  Bilbao Crystallographic Server.

  Limitations: The species_map dictionary is currently
               hardwired. This is a shortcoming of the
               Struct format.

               The space group is assumed to be 1
"""

import Numeric, math
import sys, os

try:
  file = sys.argv[1]
except:
  print "Usage:",sys.argv[0], "<struct file>"; sys.exit(1) 

species_map = { "1" : "Mg",
                "2" : "C" ,
                "3" : "O" }
print "# Space group"
print 1
#
f = open(file)
#
#   First read unit cell
#
print "# a, b, c, alpha, beta, gamma"
cell = []
for i in range(3):
  vector = f.readline()
  x, y, z = vector.split()
  cell.append([ float(x), float(y), float(z) ])

acell = Numeric.array(cell)
va = acell[0,:]
vb = acell[1,:]
vc = acell[2,:]
rtodeg = 180./3.1416
ma = math.sqrt(Numeric.dot(va,va))
mb = math.sqrt(Numeric.dot(vb,vb))
mc = math.sqrt(Numeric.dot(vc,vc))
alpha = rtodeg*math.acos(Numeric.dot(vb,vc)/(mb*mc))
beta = rtodeg*math.acos(Numeric.dot(va,vc)/(ma*mc))
gamma = rtodeg*math.acos(Numeric.dot(va,vb)/(ma*mb))
print ma, mb, mc, alpha, beta, gamma

#
#   Now the atoms
#
#   Robust code to allow for broken lines. 
#   Accumulate all the data in a big list...

natoms = int(f.readline())
print "# Number of atoms"
print natoms
strings = []
while 1:
   line = f.readline()
   if not line: break
   strings = strings + line.split()

#... and extract the information for each atom in turn

pos = 0
for a in range(natoms):
    sublist = strings[pos:pos+5]
    spindex, z , x, y, z = sublist
    x = float(x) % 1.0
    y = float(y) % 1.0
    z = float(z) % 1.0
    symbol = species_map[spindex]
    print symbol, spindex, "-", float(x), float(y), float(z)
    pos = pos + 5


