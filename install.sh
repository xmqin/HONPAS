    cd Obj
    sh ../Src/obj_setup.sh

    cp ../libint-1.1.5.tar.gz .
    tar -xvzf libint-1.1.5.tar.gz

    mkdir libint_install
    cd libint-1.1.5

    ./configure --prefix="$PWD/../libint_install"  CC=icc CXX=icpc F77=ifort
    make -j8
    make install
    cd ..

    ln -s ../ARCH-HONPAS/honpas-intel-mpi.make arch.make
    make -j8


