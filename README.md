# HONPAS

HONPAS comes with a source tarball of the libint library in the Obj/
subdirectory of the SIESTA source tree, e.g. libint-1.1.5.tar.gz. In the
following, we will suppose that you want to install it within the SIESTA build
tree with MPI support and compiled with INTEL::

    cd ~honpas/Obj
    sh ../Src/obj_setup.sh

    tar xvzf libint-1.1.5.tar.gz
    cd libint-1.1.5
    ../configure --prefix="$PWD/../libint_install" \
       CC=icc CXX=icpc F77=ifort
    make -j8 install
    cd ..


INSTALL HONPAS with libint call

    cd ~honpas/Obj

    ln -s ../ARCH-HONPAS/honpas-intel-mpi.make arch.make
    vi arch.make ## according your system

    make -j8

TESTS ::  Sibulk Diamond  Graphene
   cd ~honpas/HONPAS_Examples
You can run following honpas.pbs for torque

