import sisl

# Read initial geometry
init = sisl.Geometry.read('flos_h2o_neb.fdf')

# Create images (also initial and final [0, 180])
for i, ang in enumerate([0, 30, 60, 90, 120, 150, 180]):
    # Rotate around atom 3 (Oxygen), and only rotate atoms
    #  [4, 5] (rotating 3 is a no-op)
    new = init.rotatec(ang, origo=3, atom=[4, 5], only='xyz')
    new.write('image_{}.xyz'.format(i))
