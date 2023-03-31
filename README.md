HONPAS (Hefei Order-N Packages for Ab initio Simulations) is an ab initio electronic structure calculation software package 
for linear scaling first-principles density functional theory (DFT) calculations of large-scale systems with standard norm-
conserving pseudopotentials and numerical atomic orbitals (NAOs) under the periodic boundary conditions. 
HONPAS is developed in the framework of the SIESTA methodology and focuses on the development and implementation 
of efficient linear scaling algorithms for ab initio electronic structure calculations. 

The current version is based on SIESTA-4.1.5

 Here is a brief description of the compilation and use of HONPAS.
  You can refer INSTALL（or the script install.sh） for installation, and refer to the HONPAS_Examples for HSE06/PBE0 calculations

1. Dependencies
  HONPAS requires an external library LIBINT for two-electron repulsion integrals over Gaussian-type orbitals.
 (1) LIBINT-1.1.5 is located in the home directory of HONPAS
 (2) Build LIBINT
   # cd libint-1.1.5
   # ./configure --prefix=/your-PATH/libint-1.1.5 CC=icc CXX=icpc F77=ifort
   # make -j && make install

2. Build honpas
 (1) Go to honpas/Obj
   $ cd honpas/Obj
   $ sh ../Src/obj_setup.sh
 (2) Modify the file "arch.make" according to your system configuration as SIESTA.
   The only difference is that HONPAS requires an additionally calling path of LIBINT.
    LIBINT_LIBS=/your-PATH/libint-1.1.5/lib/libderiv.a /your-PATH/libint-1.1.5/lib/libint.a

 (3) Compile it by "make" and obtain an executable file "honpas"
   $ make -j

 
How to use :
 (1) All parameters and calculation modules of SIESTA 4.1.5 are supported, thus you can also use it for conventional DFT calculations.

 (2) For HSE06 calculations, you only need to change one parameter of the SIESTA-PBE calculation：

        xc.authors   HSE06

       PS: Other parameters have default values, no additional settings are required for most systems.

 (3) The directory "HONPAS_Examples" contains HSE06 bandstructure examples (PBE pseudopotentials, 
     input and output files) for different semiconductors, which you can refer to for the desired calculation settings.


