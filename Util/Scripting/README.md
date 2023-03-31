## Support for scripting in Siesta

* Compiling support for the embedded Lua interpreter (-DSIESTA__FLOOK
  pre-processor option, plus linking to the flook library) enables
  'internal' scripting of some operations in Siesta using Lua
  scripts. See the documentation for more information.

* Siesta is interfaced to the [AiiDA simulation
  framework](https://www.aiida.net), through the use of the
  [aiida-siesta](https://aiida_siesta_plugin.readthedocs.io) package
  of plugins and workflows.

* The [ASE (Atomic Simulation
  Environment)](https://wiki.fysik.dtu.dk/ase/) has Siesta as one of
  its supported calculators.

* Some old experiments for a python interface to Siesta built on top
  of ASE can be found in the directory `old-python-interface`.

