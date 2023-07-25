    cd Obj
    sh ../Src/obj_setup.sh

    cp ../libint-1.1.5.tar.gz .
    tar -xvzf libint-1.1.5.tar.gz

    mkdir libint_install
    cd libint-1.1.5

## For GNU compuler
    ./configure --prefix="$PWD/../libint_install"  CC=gcc CXX=g++ F77=gfortran
## For INtel compuler
#    ./configure --prefix="$PWD/../libint_install"  CC=icc CXX=icpc F77=ifort
#
    make -j8
    make install
    cd ..

## For GNU compuler and openmpi
    cp ../ARCH-HONPAS/honpas-gnu-openmpi.make arch.make
## For intel compiler and intelmpi
#    cp ../ARCH-HONPAS/honpas-intelmpi.make arch.make

    make -j8


