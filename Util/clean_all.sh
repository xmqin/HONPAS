#!/bin/sh
#
# clean all Utils
#
# Usage: [ OBJDIR=Objdir] sh clean_all
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
topdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
#---------------------------------------------
#
if [ -z "$OBJDIR" ] ; then
    OBJDIR=Obj
fi
echo ${OBJDIR}

for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile ) ; do
      relpath=${i%/*}
      echo "====> Processing $relpath ..." ; echo
      cd $relpath
      make OBJDIR=${OBJDIR} clean 
      cd $topdir
done
