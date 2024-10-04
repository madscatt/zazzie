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
contract

test for all possible types 

'''

import os
import locale

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy





class Test_dihedral_energy_getparampsi(MockerTestCase): 


   def setUp(self):
      self.data = {}
      self.data[1] = (0.6,1.0,0.0,0.0,0.0,0.0)
      self.data[2] = (1E90,0.0,0.0,0.0,0.0,0.0)
      self.data[3] = (0.3,1.0,0.0,-0.3,4.0,0.0)
      self.data[4] = (0.3,1.0,0.0,-0.3,4.0,0.0)
      self.data[5] = (1E90,0.0,0.0,0.0,0.0,0.0)
      self.data[6] = (0.6,1.0,0.0,0.0,0.0,0.0)
      self.data[7] = (0.4,1.0,0.0,0.0,0.0,0.0)
      self.data[8] = (1E90,0.0,0.0,0.0,0.0,0.0)
      self.data[9] = (0.4,1.0,0.0,0.0,0.0,0.0)
      self.data[97] = (0.0,0.0,0.0,0.0,0.0,0.0)
      self.data[98] = (0.0,0.0,0.0,0.0,0.0,0.0)
      self.data[99] = (0.0,0.0,0.0,0.0,0.0,0.0)


   def test_all(self):
      '''
      test for all possible types
      '''
      #
      for typ in [1,2,3,4,5,6,7,8,9,97,98,99]:
         result = energy.getparampsi(typ)
         expected = self.data[typ]
         self.assertEqual(result,expected)


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

