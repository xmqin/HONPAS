#!/bin/sh
#
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
user_specified_dir=$(dirname $0)
testdir=$(dirname $srcdir)/Tests
#
destdir=$(pwd)
#
# Replicate the hierarchy of makefiles
#
(cd $srcdir;
  for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile) ; do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp $relpath/*akefile ${destdir}/$relpath
  done
)
# Replicate any .h files
# This is needed in some systems with broken include file import heuristics
# (e.g., CSCS blanc)
#
(cd $srcdir;
  for i in $(find . -name '*.h' ); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp -f $relpath/*.h ${destdir}/$relpath
  done
)
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/Makefile > ${destdir}/Makefile

#
# Tests directory
# Create a list of files and use tar to process the list and copy the files
# to the destination directory
#
( cd ${testdir} ; cd .. ; find Tests  \
              -path *Reference -prune -o  \
              -path *Reference-xml -prune -o  \
              -path *work -prune      -o  \
              -path *.arch-ids  -prune -o -print \
              | tar -cf - --no-recursion -T- )   | ( cd ${destdir} ; tar xf -)
#
# (deactivated for now)
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file into the directory."
echo " *** These files are template arch.make files:"
echo " ***    gfortran.make (for gfortran compiler)"
echo " ***    intel.make (for intel compiler)"
echo " ***    DOCUMENTED-TEMPLATE.make (requires customization)"

