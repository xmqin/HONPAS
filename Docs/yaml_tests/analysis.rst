===============================================
Analyzing tests using the YAML output of SIESTA
===============================================

Minimum requirements
--------------------

Performing the analysis of the YAML data produced by SIESTA requires the
following packages:

- Python 3.4 or later;
- The *ruamel.yaml* Python package.

Optionally, the *termcolor* Python package can be installed to provide
colorized on-screen output (recommended).

To install Python, please refer to the documentation of your system
distribution. Once Python is installed, you can install in turn the additional
Python packages using PIP:

.. code::

   pip install ruamel.yaml termcolor

Depending on how you have installed Python, the PIP executable may have a
different name, e.g. *pip3* or *pip3.6*.


Data workflow
-------------

The data of every test case goes through two steps: generation by a SIESTA
executable and post-processing by a Python script. The workflow is open to
further post-processing thanks to the generation of a YAML report of the
analysis.

The following diagram illustrates how the data is processed when running a
YAML-enabled SIESTA executable:

.. graphviz::

   digraph yaml_workflow {

     graph [fontname="Helvetica"];
     node [fontname="Helvetica"];
     edge [fontname="Helvetica"];

     node [shape="box", style="filled", fillcolor="skyblue"];
     siesta [label="SIESTA"];
     yaml_cmp [label="yaml_compare.py"];

     node [shape="oval", style="filled", fillcolor="gold"];
     config [label="siesta-testsuite.yml"];
     input [label="Input file", fillcolor="mistyrose"];
     outvars [label="OUTVARS.yml"];
     refs [label="YAML_Refs/"];
     report [label="tests-report.yml"];

     node [shape="octagon", style="filled", fillcolor="white"];
     summary [label="Screen summary"];

     input -> siesta;
     siesta -> outvars;
     config -> yaml_cmp;
     outvars -> yaml_cmp;
     refs -> yaml_cmp;
     yaml_cmp -> report;
     yaml_cmp -> summary;
   }

At the end of the calculation, SIESTA writes down a file named *OUTVARS.yml*,
which contains information about both the build parameters of SIESTA and the
values of selected variables (energies for the moment). The *yaml_compare.py*
script then analyzes the content of *OUTVARS.yml* for each test case and
compares it to the corresponding reference file located in the
*Tests/YAML_Refs/* directory. The behaviour of the script can be tuned by
editing the *siesta-testuite.yml* configuration file, which contains the full
specifications of the SIESTA Test Suite, including fine-grained tolerances for
the comparisons. Once the output of all tests has been processed, the
*yaml_compare.py* script generates a detailed report in YAML format,
*tests-report.yml*, as well as a condensed summary on the terminal.

The YAML report can be read by humans but is optimized for further automated
analysis and problem diagnosis, while the on-screen summary is intended
primarily for the user and can be colorized if the terminal allows it (see the
*--fancy* option of *yaml_compare.py*).


How to use yaml_compare.py
--------------------------

The *yaml_compare.py* script is loacted in *Tests/Scripts/* and expected to
run from the top directory of the SIESTA Test Suite (*Tests/*), through the
following command: *./Scripts/yaml_compare.py [options]*. It goes through all
test cases described in *siesta-testsuite.yml* and follows the specifications
available therein.

Several parameters of the scripts can be tuned through command-line options.
For detailed information about them, just run:

.. code::

   ./Scripts/yaml_compare.py --help

The *termcolor* Python package is necessary to get fancy colorized output
through the *--fancy* option, which will simply be ignored otherwise. If this
option is not used, the on-screen summary displayed by *yaml_compare.py* is
valid RestructuredText and can be further processed to produce HTML and/or PDF
summaries.


Configuring the tests
---------------------

All configuration parameters related to the processing of the tests are stored
in the self-documented *Tests/siesta-testsuite.yml* file. All information on
how to tune the various parameters can thus be found in the file itself.


Generating new references
-------------------------

The procedure to generate new references is the following:

#. Build SIESTA using the reference configuration specifications.
#. Run the test case and carefully check the output, both on screen and in
   *OUTVARS.yml*.
#. Copy *OUTVARS.yml* to *~siesta/Tests/YAML_Refs/test_name.yml*, where you
   replace *test_name* by the actual name of the test.
#. Add the test case specifications into *~siesta/Tests/siesta-testsuite.yml*.
#. Run *yaml_compare.py* and check that the reference passes the test against
   itself.

.. note:: 

   TODO: Define the reference specifications for the build of SIESTA.

