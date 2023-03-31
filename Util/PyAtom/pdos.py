#
# PDOS processor (alpha stage)
#
# Uses Python XML package (Follow XML SIG link in www.python.org)
# Python 2.1 and 2.2 might have this built-in.
#
from xml.sax import saxutils
from Numeric import *
import biggles

#def normalize_whitespace(text):
#  "Remove redundant whitespace from a string"
#   return string.join(string.split(text), ' ')


#
# User-configurable function to decide which orbitals to process
#
def orbitalp(seq_no,n,l,m,z):
    if string.atoi(seq_no) == 1:
       return 1
    else:
       return 0

class PDOS(saxutils.DefaultHandler):
   def __init__(self):
	self.inData = 0
	self.inOrbital = 0
	self.inEnergyGrid = 0
	self.data = ""
	self.energies = ""
       
   def startElement(self, name, attrs):

      if name == 'energy_values':
         self.inEnergyGrid = 1

      if name == 'orbital':

	 n = attrs.get('n', None)
         l = attrs.get('l', None)
         m = attrs.get('m', None)
         z = attrs.get('z', None)
         seq_no = attrs.get('index',None)
         self.processOrbital = orbitalp(seq_no,n,l,m,z)
         print "Orbital", n, l, m, z, seq_no

      if name == 'data':
         self.inData = 1

   def characters(self, ch):
         if self.inData:
            if self.processOrbital:       # if not, disregard
              self.data = self.data + ch
         if self.inEnergyGrid:
            self.energies = self.energies + ch

   def endElement(self, name):

       if name == 'energy_values':
         self.inEnergyGrid = 0
         ll = string.split(self.energies)
         llfloat = map(string.atof,ll)
         self.energy = array(llfloat)
         self.npts = len(self.energy)
         self.dos = 0*self.energy

       if name == 'data':
         self.inData = 0
         if self.processOrbital:
           ll = string.split(self.data)
	   llfloat = map(string.atof,ll)
           dos = array(llfloat)

           p = biggles.FramedPlot(title="PDOS")
           p.add(biggles.Curve(self.energy,dos))
           p.show()

         self.data=""

#------------------------------
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces

if __name__ == '__main__':
        # Create a parser
        parser = make_parser()
        # Tell the parser we are not interested in XML namespaces
        parser.setFeature(feature_namespaces, 0)

        # Create the handler
        dh = PDOS()

        # Tell the parser to use our handler
        parser.setContentHandler(dh)

        # Parse the input
        parser.parse('pdos.xml')



# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
