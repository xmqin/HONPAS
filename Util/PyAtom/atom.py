#!/usr/bin/python

#
from Numeric import *
from Scientific.IO.NetCDF import *

class rfunc:
	def __init__(self,f,delta,cutoff):
		self.delta = delta
		self.cutoff = cutoff
		self.f = f
        def __str__(self):
           return ("(rfunc) n: %s delta: %f cutoff: %f" % 
                    (len(self.f), self.delta, self.cutoff))
        def dump(self):
           for i in range(len(self.f)):
             print i*self.delta, self.f[i] 


class orbital(rfunc):
  def __init__(self,f,delta,cutoff,n,l,z):
     rfunc.__init__(self,f,delta,cutoff)
     self.l = l
     self.z = z
     self.n = n

  def __str__(self):
     return ("orb n:%s l:%s z:%s " % (self.n,self.l,self.z) + rfunc.__str__(self))
  def dump(self):
     print "# ", str(self)
     rfunc.dump(self)
	
class atom:
  def __init__(self,filename):
    f=NetCDFFile(filename)
    self.element=f.Element
    vars = f.variables
    p = vars['vna']
    self.vna=rfunc(p[:],p.Vna_delta[0], p.Vna_cutoff[0])
    norbs = f.dimensions['norbs']
    self.base= []
    for i in range(norbs):
      p = vars['orb'][i]
      l = vars['orbnl_l'][i]
      n = vars['orbnl_n'][i]
      z = vars['orbnl_z'][i]
      delta = vars['delta'][i]
      cutoff = vars['cutoff'][i]
      orb=orbital(p[:],delta,cutoff,n,l,z)
      self.base.append(orb)
    del p, l, n, z, delta, cutoff, norbs, vars
    del f
  def dump(self):
    print "Orbitals:"
    for i in range(len(self.base)):
       print self.base[i]
    print "Vna: ", self.vna
 

if __name__ == "__main__":
  print """

  If you had used the "-i" flag to python
  you could now type things such as:
 
    a = atom("H.ion.nc")
    print a.base[0]
    a.base[0].dump()

  """
a = atom("O.ion.nc")
print a.base[0]
#    a.base[0].dump()

# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
