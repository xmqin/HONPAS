#!/bin/sh

#
# Simple launcher of mpd daemons and handler for single-node 
# interactive parallel calculation
#

echo $(hostname) > mpd.hosts
NSLOTS=2
NODES=`wc mpd.hosts | awk '{print $1}'`
#
# Boot the MPI2 engine.
#
mpdboot --rsh=ssh --totalnum=$NODES --file=mpd.hosts
sleep 5
#
# Inspect if all MPI nodes have been activated.
#
mpdtrace -l
#
# Check the connectivity.
#
mpdringtest 100
#
# Check if you can run trivial non-MPI jobs.
#
mpdrun -l -n $NSLOTS hostname
#
# Execute your MPI program.
#
prog_default="../../../siesta"
prog=${SIESTA:-${prog_default}}
#
make SIESTA="mpiexec -n $NSLOTS $prog"
#
# Shut down the MPI2 engine and exit the script.
#
echo "shutting down mpd..."
mpdallexit
exit 0
