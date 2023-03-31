#!/bin/sh

#
# Simple handler for interactive single-node parallel jobs
# It just launches the mpd daemons
# To kill them: mdpallexit

echo $(hostname) > mpd.hosts
NSLOTS=2
NODES=`wc mpd.hosts | awk '{print $1}'`
#
# Boot the MPI2 engine.
#
echo "Booting mpd daemons..."
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
echo "Ready to run parallel jobs. Use 'mpiexec -n NODES program'"
