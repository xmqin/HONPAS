from Siesta.Interface import Atom, Crystal

#
def ReadStruct(filename,species_map):
    """Return "Crystal" object read from Siesta STRUCT-file.
       species_map is a dictionary mapping the species index to 
       atomic symbols. For example:

          species_map = { "1" : "Mg",
                          "2" : "C" ,
                          "3" : "O" }

	In a forthcoming version this information could be obtained
        directly from the atomic number in the .STRUCT_{IN,OUT} file.

     """
    f = open(filename)
#
#   First read unit cell
#
    cell = []
    for i in range(3):
      vector = f.readline()
      x, y, z = vector.split()
      cell.append([ float(x), float(y), float(z) ])
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
    crystal = Crystal([])
    pos = 0
    for a in range(natoms):
        sublist = strings[pos:pos+5]
        spindex, z , x, y, z = sublist
        symbol = species_map[spindex]
        crystal.append(Atom(symbol, [float(x), float(y), float(z)]))
        pos = pos + 5

    crystal.SetUnitCell(cell)

    return crystal
#

