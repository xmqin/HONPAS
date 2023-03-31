#!/bin/sh
#
# Simple launcher of mpd daemons and handler for single-node 
# interactive parallel calculation
# (thanks to Javier Lopez Solano)

echo $(hostname) > mpd.hosts
NSLOTS=4
NODES=`wc mpd.hosts | awk '{print $1}'`

prog_default="../../../siesta"
prog=${SIESTA:-${prog_default}}
#
PATH=/opt/intel/impi/3.1/bin64:$PATH
make SIESTA="mpirun -r ssh -np $NSLOTS --mpd=/opt/intel/impi/3.1/bin64/mpd  $prog"
#
exit 0
