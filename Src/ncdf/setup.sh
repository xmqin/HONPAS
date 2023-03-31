#!/bin/bash

_vpath=.
[ -n "$VPATH" ] && _vpath=$VPATH

# Run the setup in the lib/fdict directory
VPATH="$_vpath/fdict" $_vpath/fdict/setup.sh $@
retval=$?

exit $retval
