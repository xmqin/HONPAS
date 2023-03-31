#!/usr/bin/python

"""
Converts a Siesta .STRUCT_OUT file into CIF format
"""


import sys, Numeric, math, string

try:
  file = sys.argv[1]
  species_names = sys.argv[2:]
  if  len(species_names) == 0:
    print "Usage:",sys.argv[0], "<struct file>  Symbol Symbol ... " ; sys.exit(1)
except:
  print "Usage:",sys.argv[0], "<struct file>  Symbol Symbol ... "; sys.exit(1) 
  
#
f = open(file)
#
#   First read unit cell
#
print """
data_C3D_block
loop_ """

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
print "_cell_length_a ", ma
print "_cell_length_b ", mb
print "_cell_length_c ", mc
print "_cell_angle_alpha ", alpha
print "_cell_angle_beta ", beta
print "_cell_angle_gamma ", gamma

print """
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number    '1'

loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
  1     'x, y, z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""

#
#   Now the atoms
#
#   Robust code to allow for broken lines. 
#   Accumulate all the data in a big list...

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
        spindex = string.atoi(spindex)
        x = float(x) % 1.0
        y = float(y) % 1.0
        z = float(z) % 1.0
        symbol = species_names[spindex-1]
        print symbol, symbol, float(x), float(y), float(z)
        pos = pos + 5





