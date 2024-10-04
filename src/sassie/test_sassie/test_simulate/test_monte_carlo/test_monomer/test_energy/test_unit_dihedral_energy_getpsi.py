'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
	 Core-Testing: Copyright (C) 2011 Hailiang Zhang, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from sassie.core_testing.util import env

'''
unit test:

contract

test for all possible aa pairs

'''

import os
import locale

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy




DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../../../data/simulate/monte_carlo/monomer/energy/'


class Test_dihedral_energy_getpsi(MockerTestCase): 


   def setUp(self):
      self.lp = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HSD','HIS','HSE','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','NULLN','MAN']
      self.lt = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HSD','HIS','HSE','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','NULLC','MAN']
      self.data = {}
      for p in self.lp:
         self.data[p] = {}
      f = os.path.join(DataPath,'standard_psi_types.txt')
      for line in open(f,'r').readlines():
         dum = line.split()
         p = dum[0]
         t = dum[1]
         v = locale.atoi(dum[2])
         self.data[p][t]=v


   def test_all(self):
      '''
      test for all possible aa pairs
      '''
      #
      for p in self.lp:
         for t in self.lt:
            value = energy.getpsi(p,t)
            self.assertEqual(value,self.data[p][t])


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

