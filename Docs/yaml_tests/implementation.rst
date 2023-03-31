==========================================
Implementation of the YAML variable dumper
==========================================

Changes in SIESTA
-----------------

The YAML variable dumper has been implemented in as an orthogonal way as
possible from the rest of SIESTA. It consists in:

- a Fortran module producing valid YAML 1.2, codename *m_io_yaml*;
- the inclusion of the corresponding object file into the MAkefile of SIESTA;
- a call to the main routine of the module, *siesta_write_yaml* at the end of
  *siesta.F*.

The call to the I/O routine is both compatible with serial and parallel runs
of SIESTA:

.. code::

         if (ionode) call siesta_write_yaml()

This strategy is expected to be harmless to any other ongoing developments
within SIESTA and allows to add the YAML output capability into already
published versions of SIESTA. Patching such versions will indeed provide
useful information on possible issues that were previously undetected.

.. note::

   Contrary to XML which has tags and attributes, YAML is just plain text
   using indenting (think of Python) to define data structures. As a
   consequence, emitting YAML does **not** introduce any additional dependency
   of SIESTA on another external library and is much easier to read by humans.


Source code
-----------

The complete source code of the YAML dumper can be found in
*~siesta/Src/m_io_yaml.F90*. The file is self-documented using the Doxygen
format.

.. note::

   The Doxygen infrastructure is not merged yet...


Documentation
-------------

To generate the Doxygen documentation, just run *doxygen* from the top
source directory of SIESTA. The resulting files can then be found within the
*~siesta/Docs/developers/* directory.

To access the HTML website, look for *index.html* inside the *source/*
subdirectory. To build the PDF manual, look at the *refman.tex* main document
inside the *latex/* subdirectory.

.. tip::

   The provided *Doxyfile* has been tested with Doxygen 1.8.11. If you have an
   earlier version and experience crashes, just add the problematic files to
   the *EXCLUDE* variable within the *Doxyfile* and re-run doxygen. Repeat
   this procedure as many times as necessary.

