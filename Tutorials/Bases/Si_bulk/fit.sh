#!/bin/sh

rm -f fit.dat
cat > fit.dat <<EOF
E vs V fit for Si
8
EOF
#
awk '{print $1/(0.529**3), $2/13.6}' ev.dat >> fit.dat
murna fit
cat fit.mout
