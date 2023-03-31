title: Contributing to Siesta

The Siesta code is heavily dependent on external contributions and we
welcome the help.

Any size of contribution is welcome, from small bug-fixes to fully new implementations.

We do have some guidelines for contributing code to Siesta which will be
discussed in the following.

## Development platform

Any contribution (bug or addition) *must* be transmitted through the Launchpad
development site [here](http://launchpad.net/siesta).
To begin development please follow the following FAQs:

1. [Setting up a shared repository](https://answers.launchpad.net/siesta/+faq/2751)
2. [Creating a bugfix merge-request](https://answers.launchpad.net/siesta/+faq/2747)
3. [Creating merge-request for new feature](https://answers.launchpad.net/siesta/+faq/2748)


## Code development

All code *must* adhere to the following rules.

### Fortran version

The fortran code version is the [2003 standard](ftp://ftp.nag.co.uk/sc22wg5/N1601-N1650/N1601.pdf.gz).
You can check compliance of your code by using the `-std=f2003` flag in `gfortran` (other compilers might have similar diagnostic capabilities).

One feature that should not be used is `reallocation on assignment`:

        real, allocatable :: r(:)
        r = (/2., 3./)
        r = (/2., 3./, 4./)

This "whole array" assignments can trigger checks by the compiler for the need of re-allocation, impacting
performance.

When using allocatables *always* use array-sections on the left-hand-side of assignments.  Formally, an array section is not "allocatable", and this has the effect of avoiding the check for reallocation:
		
        real, allocatable :: r(:)
        allocate(r(2))
        r(:) = (/2., 3./)

To learn more about this, see [this post](https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran). To check that your code will not trigger checks for reallocation, you can use the `-Wrealloc-lhs-all` flag in `gfortran`.

### Parallel environments

* OpenMP, 4.0 is allowed.

* MPI, use MPI-2. If you need MPI-3 please contact the developers.

### Modules

All code should be put in modules. The name of the modules should reflect the code contained.

* Filenames must *not* contain `m_*` or `*_m`.
* Modules *must* be named `*_m` such that the module name may be reused for routine names.
* Derived types *must* be named `*_t` such that the type name may be reused for variable names.

Note that a large portion of the current code base does not obey the above admonitions.
The conversion is a long-term work-in-progress.


### Declarations

* All fortran declarations *must* use `implicit none` at the top level. Note that `implicit none` is
	global when declared at the module level, hence it is not necessary to declare it in subroutines or
	functions.

* Always declare precisions using `real(dp)`, the older `real*8` is not an ANSI standard. The module [[precision]]
  contains the necessary complex, real and integer precisions required.

* Assigning constants to variables should specify precision via the suffix `double = 0._dp`

* Defining pointers should declare them as `=> null()` upon declaration

* The module [[alloc]] contains routines to allocate pointers (must be nullified before call)

### File format

* Free form files (`.f90` or `.F90`)
* Maximum line-length of 132 characters
* Use `&` for continuation lines
* Use `!` for comments
* Follow indentation in steps of 2
* Use Doxygen like documentation `!>` for starting documentation

### Casting data-types to different precision/data-types

Casting should *always* have a precision specifier (if the interface allows)

	! Single precisison
	c_var = cmplx(0._sp, 0._sp, sp)
	! Taking real/imaginary parts
	s_var1 = real(c_var, sp)
	s_var2 = aimag(c_var)
	
	z_var = cmplx(0._dp, 0._dp, dp)
	! Taking real/imaginary parts
	d_var1 = real(z_var, dp)
	d_var2 = aimag(z_var)
	
	! Casting from same type
	s_var1 = real(d_var1, sp)
	s_var2 = real(d_var2, sp)
	d_var1 = real(s_var1, dp)
	d_var2 = real(s_var2, dp)

Do *not* use any of these:

* `dcmplx`, use `cmplx(.., .., dp)`
* `dimag`, use `aimag`
* `dreal`, use `real(.., dp)`


### Logical operators

Use logical operators

	== ! instead of .eq.
	/= ! instead of .ne.
	<  ! instead of .lt.
	<= ! instead of .le.
	>  ! instead of .gt.
	>= ! instead of .ge.


### Controlling extensions via preprocessors

Implementation of optional dependencies are done via preprocessor statements.

Preprocessor names *must* be of the form:

    SIESTA__MRRR ! NOTE: two underscores
    SIESTA__ELPA

where the suffix is unique and has meaning for the intent of the optional dependency.
Note the two underscores to ensure unique names.
