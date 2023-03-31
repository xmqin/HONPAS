title: Overall structure

The main program [[siesta(program)]] contains an outer loop for
geometries, that can be used in Molecular Dynamics (MD) simulations,
geometry relaxations, and vibrational-properties.  Inside the loop,
the total energy and forces are computed for each particular geometry.

Loop termination is controlled by explicit numbers of iterations
(in the case of MD or phonon calculations), or by a flag to signal convergence
(in the case of geometry relaxations).

