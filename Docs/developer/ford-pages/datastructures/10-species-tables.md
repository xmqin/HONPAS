title: Tables for species information

Lots of operations in SIESTA depend on the knowledge of the specifics
of the atomic orbitals, KB projectors, ionic potentials, etc. By
design from the earliest versions of the program, this is done through
an independent interface. Routines in the [[atmfuncs(module)]] such as
[[atmfuncs:rcut]] are called everytime the information is needed (in
this case, the cutoff radius for a specific orbital, projector, or
neutral-atom potential).

The routines in [[atmfuncs(module)]] take their information from
the [[atm_types:species]] array, whose elements (one per species) are instances of the
[[atm_types:species_info]] derived type.

`species_info` records are filled by the routines that generate
PAOs,KBs, etc in module [[atom(module)]], or the information can be read
from `.ion` files produced by earlier runs of the program.


The typical interface for a routine in [[atmfuncs]] makes use of two indices:

```fortran
      function rcut(is,io)
      real(dp) rcut
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io> => basis orbitals
                                   ! io<0  => KB projectors
                                   ! io=0 : Local screened pseudopotential
```

The way in which the second index (`io`) is overloaded means that
there is no possibility to extend this interface to new types of
radial functions associated to the species in the system. When new
functions were needed (when LDA+U was implemented) an extended
interface was introduced, but only for a limited number of operations.

The new interface uses the concept of a `registry` of radial functions
(see [[m_matel_registry(module)]]).
