#!/bin/bash

# This should be altered in case your
# lapack source is somewhere else
if [ -z "$LAPACK_DIR" ]; then
    LAPACK_DIR=$HOME/LA/lapack
fi
if [ ! -e $LAPACK_DIR/lapack.pc.in ]; then
    echo "Could not find lapack.pc.in in directory:"
    echo "  $LAPACK_DIR"
    echo
    echo "Files in $LAPACK_DIR:"
    ls -l $LAPACK_DIR
    echo ""
    echo "Please run:"
    echo "  LAPACK_DIR=<path-to-git-clone-lapack> $0"
    exit 1
fi

cur_dir=$(pwd)

# First gather the lapack sources
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -l $LAPACK_DIR/SRC $LAPACK_DIR/SRC/DEPRECATED $LAPACK_DIR/INSTALL \
       --list-add-file lapack_add.files \
       --list-add-routine lapack_add.routines \
       --list-remove-file lapack_remove.files \
       --list-remove-routine lapack_remove.routines \
       -o lapack_tmp.f

# Append the license to the sources
cat $LAPACK_DIR/LICENSE | sed -e 's/^/!/' > lapack_license_tmp
cat lapack_license_tmp lapack_tmp.f > lapack.F
rm lapack_license_tmp lapack_tmp.f


# Now we need to locate the lapack sources as well
python linalg2file.py -d $cur_dir/.. $cur_dir/../../Util/TS/TBtrans \
       -f lapack.F \
       -l $LAPACK_DIR/BLAS/SRC \
       --list-add-file blas_add.files -o blas_tmp.f

# Append the license to the sources
cat $LAPACK_DIR/LICENSE | sed -e 's/^/!/' > blas_license_tmp
cat blas_license_tmp blas_tmp.f > blas.F
rm blas_license_tmp blas_tmp.f
