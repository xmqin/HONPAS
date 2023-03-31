#!/usr/bin/env python

from distutils.core import setup

setup(name='Siesta',
      version='0.1',
      description='Siesta Scripting Utilities',
      author='Alberto Garcia',
      author_email='albertog@icmab.es',
      url='http://fisica.ehu.es/ag/',
      packages=['Siesta'],
      scripts=['scripts/plot_struct_rasmol.py',
               'scripts/plot_struct_jmol.py',
               'scripts/struct_to_pdb.py',
               'scripts/struct_to_cif.py',
               'scripts/struct_to_bplot.py',]
     )

