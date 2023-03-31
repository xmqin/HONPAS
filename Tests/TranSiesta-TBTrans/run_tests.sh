#!/bin/bash

# Run TS tests
MPI=${MPI:-mpirun -np 2}

# To run in serial mode, replace the 'mpirun' line
# by the appropriate incantation.

# Here we run the tests 
for d in ts_au \
	     ts_au_100 \
	     ts_au_100_0.25V \
	     ts_au_100_repetition \
	     ts_au_100_repetition_0.25V \
	     ts_au_repetition \
	     ts_graphene \
	     ts_term3 \
	     ts_term4
do
    cd $d
    #make clean
    make MPI="$MPI"
    cd ..
done
