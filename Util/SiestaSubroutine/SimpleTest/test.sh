#!/bin/sh

#
# Select the appropriate run below (at the end).
#
# Take care to use the appropriate (e.g., parallel or serial) copy
# of siesta. To run parallel jobs, you have to figure out how to
# set the MPI environment. For a simple single-node interactive calculation
# (where possible, such as in a Rocks cluster), you can use the parallel.sh
# script.
#
# Make sure that you are using the right version of Siesta. 
# The SIESTA setting below is quite naive and might not work
# in all cases. You can call this script as:
#
#            SIESTA=/path/to/siesta test.sh
#
#

ROOT="../../../.."
PSEUDOS=${ROOT}/Tests/Pseudos
#
if [ -z "$SIESTA" ] ; then
      SIESTA=${ROOT}/Obj/siesta
fi
echo "Using Siesta executable: $SIESTA"
#

#rm -r work
if [ -d work ] ; then
   echo "Work directory 'work' exists. Please delete it"
   exit
fi

#
mkdir work
cd work
cp -p ../*.fdf .
cp ${PSEUDOS}/H.psf  .
cp ${PSEUDOS}/O.psf  .
ln -sf ${SIESTA} ./siesta
#

echo ""; echo "simple_pipes_serial"
../Src/simple_pipes_serial    | tee simple_pipes_serial.out
mv h2o.out siesta_pipes_serial.out

echo ""; echo "simple_pipes_parallel"
../Src/simple_pipes_parallel  | tee simple_pipes_parallel.out
mv h2o.out siesta_pipes_parallel.out

echo ""; echo "simple_mpi_serial"
../Src/simple_mpi_serial      | tee simple_mpi_serial.out
mv h2o.out siesta_mpi_serial.out

echo ""; echo "simple_mpi_parallel"
mpirun -np 2 -output-filename simple_mpi_parallel.out ../Src/simple_mpi_parallel
mv h2o.out siesta_mpi_parallel.out

cat socket.fdf >> h2o.fdf

echo ""; echo "simple_sockets_serial"
../Src/simple_sockets_serial    | tee simple_sockets_serial.out
mv h2o.out siesta_sockets_serial.out

echo ""; echo "simple_sockets_parallel"
../Src/simple_sockets_parallel  | tee simple_sockets_parallel.out
mv h2o.out siesta_sockets_parallel.out

