#!/usr/bin/python

"""
In Metadynamics we perform Nose-thermostat MD simulations at fixed cell.
If then one wants to relax an intermediate structure:

1. One can do a standard variable-cell CG or Broyden run.
2. It might be interesting to perform a quenched Parrinello-Rahman run
   starting from the configuration at which the Metadynamics step ended.
   For that, we need the XV file, but also a PR_RESTART file.
   We are missing: 
    a) the initial kinetic energy of the cell. We will let
       that develop with time.
    b) The "old" coordinates. Those are in the .NOSE_RESTART file, but
       this is not always kept... We can use the .STRUCT_OUT file instead.
       (but converting from fractional to cartesian in Bohrs)
    c) The "old" lattice vectors. For this we will use just the current ones

This script then takes XV and STRUCT_OUT files, and generates a PR_RESTART file.
"""

import sys, Numeric, math

try:
  xvfile = sys.argv[1] ; strfile = sys.argv[2]
except:
  print "Usage:",sys.argv[0], "<xv file> <struct file> "; sys.exit(1) 

#
xv = open(xvfile)
#
cell = []
for i in range(3):
  vector = xv.readline()
  x, y, z, vx, vy, vz  = vector.split()
  cell.append([ float(x), float(y), float(z) ])
acell = Numeric.array(cell)
va = acell[0,:]
vb = acell[1,:]
vc = acell[2,:]
for i in range(3):
  sys.stdout.write("%14.9f " % va[i])
sys.stdout.write("\n")
for i in range(3):
  sys.stdout.write("%14.9f " % vb[i])
sys.stdout.write("\n")
for i in range(3):
  sys.stdout.write("%14.9f " % vc[i])
sys.stdout.write("\n")

xv.close()
#
f = open(strfile)
# skip three lines
for i in range(3):
  f.readline()
#
#   Now the atoms
#
#   Robust code to allow for broken lines. 
#   Accumulate all the data in a big list...

natoms = int(f.readline())
print natoms
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
        x = float(x)
        y = float(y)
        z = float(z)
        fract = Numeric.array([x,y,z])
        cartesian = Numeric.dot(Numeric.transpose(cell),fract)
        sys.stdout.write("%5i " % spindex)
        for i in range(3):
         sys.stdout.write("%14.9f " % cartesian[i])
        sys.stdout.write("\n")
        pos = pos + 5





