! We define a concatenating function which allows to
! "abstratize" names of routines.
! This is required as
! #define TYPE_NAME foo
! #define PREC bar
!  subroutine TYPE_NAME_PREC
! will fail due to unambiguous naming scheme.

!
! TODO
! 
! This directive does not behave as expected.
! The expansion is performed too late. Hence, the 
! concatenation is errorneous!?!
#define CCAT(x,y) x##y