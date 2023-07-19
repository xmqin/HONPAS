dnl
dnl Copyright (c) Toby White <tow21@cam.ac.uk>
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
dnl "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
dnl LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
dnl A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
dnl OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
dnl SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
dnl LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
dnl DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
dnl THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
dnl (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
dnl OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl                                                                           
AC_INIT([siesta],[2.0],[tow21@cam.ac.uk])
AC_CONFIG_SRCDIR([siesta.F])
dnl
dnl All macros are kept in Src/m4:
dnl AC_CONFIG_MACRO_DIR([m4])
dnl FIXME that doesn't work. Why not?
dnl
dnl Find system type.
dnl Uses rewritten AC_CONFIG_AUX_DIR to look for config.sub
dnl instead of install script
AC_CONFIG_AUX_DIR(.)
AC_CANONICAL_HOST
dnl
dnl FIXME: just for the moment to work with NG's fpp detection:
dnl m4_pattern_allow([AM_CPPFLAGS])
dnl this is no longer necessary in my current fortran.m4. Upstream needs it, though.
dnl
dnl Get Fortran compiler and check characteristics.
dnl Must check first if we are using MPI, since compiler
dnl choice depends on this.
dnl
AC_ARG_ENABLE([mpi],
  AC_HELP_STRING([--enable-mpi], [Compile the parallel version of SIESTA]),
  [use_mpi=$enableval], [use_mpi=no])
if test $use_mpi = yes; then
  AC_DEFINE(MPI, 1, [Define if compiling parallel version.])
fi
dnl
AC_LANG(Fortran)
AC_PROG_FC
dnl
AC_FC_FIXEDFORM([f])
AC_SUBST(FCFLAGS_fixed_f)
dnl
AC_FC_FREEFORM([f90])
AC_SUBST(FCFLAGS_free_f90)
dnl
AC_LANG_PUSH(Preprocessed Fortran)
AC_PROG_FPP
dnl
AC_FPP_FIXEDFORM([F])
AC_SUBST(FPPFLAGS_fixed_F)
dnl
AC_FPP_FREEFORM([F90])
AC_SUBST(FPPFLAGS_free_F90)
AC_LANG_POP(Preprocessed Fortran)
dnl
if test $ac_cv_fpp_build_rule = direct; then
     F_RULE=['$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< ']
     F90_RULE=['$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< ']
elif test $ac_cv_fpp_build_rule = indirect; then
     F_RULE=['$(FPP) $(FPPFLAGS) $< $(FPP_OUTPUT) aux_$*.f && $(FC) -c -o $*.o $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f) aux_$*.f']
     F90_RULE=['$(FPP) $(FPPFLAGS) $< $(FPP_OUTPUT) aux_$*.f90 && $(FC) -c -o $*.o $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90) aux_$*.f90']
fi
dnl Generate dependency rules according to whether we need
dnl separate preprocessing stage.
AC_SUBST(F_RULE)
AC_SUBST(F90_RULE)
dnl
AC_PROG_RANLIB
dnl
dnl Compile-type options.
dnl
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [Compile with debugging support]),
  [use_debug=$enableval], [use_debug=no])
if test $use_debug = yes; then
  FCFLAGS="$FCFLAGS $FFLAGS_DEBUG"
fi
dnl
AC_ARG_ENABLE([fast],
  AC_HELP_STRING([--enable-fast], [Compile with best known optimization flags]),
  [use_debug=$enableval], [use_debug=no])
if test $use_debug = yes; then
  FCFLAGS="$FCFLAGS $FFLAGS_FAST"
fi
dnl
dnl Generate host/compiler ID string to encode in program
TW_FC_ID
SIESTA_ARCH="$host--$FC_ID"
AC_SUBST(SIESTA_ARCH)
dnl
dnl Check for any extra compiler-specific flags necessary
dnl Also set FFLAGS_DEBUG and FFLAGS_FAST
dnl SYS will be "nag" (ie standard f95) for almost all cases
SYS=nag
TW_FC_ID_FLAGS
dnl
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [Compile with debugging support]),
  [use_debug=$enableval], [use_debug=no])
if test $use_debug = yes; then
  FCFLAGS="$FCFLAGS $FFLAGS_DEBUG"
fi

AC_ARG_ENABLE([fast],
  AC_HELP_STRING([--enable-fast], [Compile with best known optimization flags]),
  [use_debug=$enableval], [use_debug=no])
if test $use_debug = yes; then
  FCFLAGS="$FCFLAGS $FFLAGS_FAST"
fi
dnl
dnl Find how to compile MPI
if test $use_mpi = yes; then
  ACX_MPI([:],AC_MSG_ERROR([Don't know how to compile MPI program]))
  FC=$MPIFC
  MPI_INTERFACE=libmpi_f90.a
  MPI_INCLUDE=.
  DEFS="$DEFS MPI"
fi
AC_SUBST(MPI_INTERFACE)
AC_SUBST(MPI_INCLUDE)
dnl Check compiler characteristics; it must support F95
dnl to use the cpu_time intrinsic. Otherwise, SYS must
dnl be set appropriately. Only pgf90 should fail this.
dnl Actually, we should simply check to see if cpu_time is supported.
dnl TW_CHECK_FC_95([SYS=nag],dnl
dnl   TW_CHECK_FC_90([SYS=cputim],[])
dnl                   )
AC_SUBST(SYS)
dnl
dnl Check for compiler kind flags
dnl
AC_FC_GET_REAL_KINDS
SP_KIND=$ac_fc_real_kind_sp
DP_KIND=$ac_fc_real_kind_dp
AC_SUBST(SP_KIND)
AC_SUBST(DP_KIND)
dnl
dnl Check whether we need to pass a flag to get at the flush intrinsic.
TW_FC_CHECK_FLUSH
TW_FC_CHECK_ABORT
dnl this next flag only actually needed if compiling siesta blas & lapack.
TW_FC_CHECK_DCFUNS
dnl
dnl Check for each relevant library in turn
dnl in each case, assume that a system-supplied version is better;
dnl otherwise fallback to SIESTA-supplied version.
dnl We need NetCDF and BLAS for both serial and parallel builds.
dnl
dnl 1: netcdf
AC_ARG_WITH(netcdf,
        [AC_HELP_STRING([--with-netcdf=<lib>], [use NetCDF library <lib>])])
if test x$with_netcdf != xno; then
  TW_PATH_NETCDF([DEFS="$DEFS CDF";NETCDF_INTERFACE="libnetcdf_f90.a"],
               [if test ! -z $with_netcdf; then AC_MSG_ERROR([Could not find NetCDF library.]); fi])
fi
dnl 
dnl Linear Algebra libraries
dnl
# Check for BLAS
AC_ARG_WITH(siesta-blas,
           [AC_HELP_STRING([--with-siesta-blas], [use BLAS library packaged with SIESTA])])
if test "x$with_siesta_blas" = xyes; then
   AC_MSG_NOTICE([Using BLAS library packaged with SIESTA])
   COMP_LIBS="libblas.a $COMP_LIBS"
   # if we're using SIESTA's blas, we have to use SIESTA's lapack as well, since we can't
   # test for lapack without a blas library already built. Would not make much sense anyway.
   with_siesta_lapack=yes
else
   TW_FIND_FC_BLAS([],
                   [AC_MSG_NOTICE(Using BLAS packaged with SIESTA); BLAS_LIBS="libblas.a"; COMP_LIBS="libblas.a $COMP_LIBS"])
fi
dnl
# Check for LAPACK
AC_ARG_WITH(siesta-lapack,
      [AC_HELP_STRING([--with-siesta-lapack], [use LAPACK library packaged with SIESTA])])
if test "x$with_siesta_lapack" = xyes; then
     AC_MSG_NOTICE([Using LAPACK packaged with SIESTA])
     COMP_LIBS="dc_lapack.a liblapack.a $COMP_LIBS"
else
     siesta_lapack=no
     TW_FIND_LAPACK([siesta_lapack=yes],
                    [AC_MSG_NOTICE([Using LAPACK packaged with SIESTA])
                     COMP_LIBS="dc_lapack.a liblapack.a $COMP_LIBS";
                     LAPACK_LIBS="dc_lapack.a liblapack.a"])
#FIXME: also need to test if we have all necessary calls. eg dsygst
     if test $siesta_lapack = yes; then
        TW_CHECK_DC_LAPACK(
         [AC_MSG_NOTICE([using DC_LAPACK routines packaged with SIESTA due to bug in library. Linker flag might be needed to avoid duplicate symbols])
          COMP_LIBS="dc_lapack.a $COMP_LIBS"],
         [AC_MSG_NOTICE([using additional LAPACK routines packaged with SIESTA])
          COMP_LIBS="dc_lapack.a $COMP_LIBS"]
        )
     fi
fi
dnl If we're configuring an MPI version, we need BLACS and SCALAPACK
dnl Don't error out if we can't find them - they're usually in some
dnl obscure place.
dnl We check if BLACS_LIBS or SCALAPACK_LIBS has been set; and if so
dnl we check if it works. We also check for a couple of the most common
dnl names for libraries. Otherwise, the user must sort it out themselves.
if test $use_mpi = yes; then
   TW_CHECK_BLACS([],
                  [AC_MSG_NOTICE([Don't know how to compile with BLACS - edit arch.make by hand])])
   TW_CHECK_SCALAPACK([],
                      [AC_MSG_NOTICE([Don't know how to compile with ScaLAPACK - edit arch.make by hand])])
fi
dnl
dnl 4: F90 MPI interface
dnl FIXME: no check here - always assume SIESTA version.
dnl ACX_MPI should do something for us.
dnl
dnl
AC_SUBST(NETCDF_LIBS)
AC_SUBST(NETCDF_INTERFACE)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(BLACS_LIBS)
AC_SUBST(SCALAPACK_LIBS)
AC_SUBST(COMP_LIBS)
dnl
dnl Fix up DEFS -> FPPFLAGS conversion.
for i in $DEFS; do
   FPPFLAGS="$FPPFLAGS $FPPFLAGS_DEF$i"
done
dnl
AC_SUBST(FPPFLAGS)
AC_SUBST(FFLAGS)
dnl
AC_SUBST(FCFLAGS)
AC_SUBST(LDFLAGS)
dnl
AC_CONFIG_FILES([arch.make])
AC_OUTPUT
dnl
