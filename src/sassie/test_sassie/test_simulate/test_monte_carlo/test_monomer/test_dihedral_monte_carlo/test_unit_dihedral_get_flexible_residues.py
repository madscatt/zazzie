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

test if numranges==1
test if numranges==2


'''

import os

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import dihedral_monte_carlo as dihedral



class Test_dihedral_get_flexible_residues(MockerTestCase): 

   def setUp(self):
      pass


   def test_numranges1(self):
      '''
      test if numranges==1
      '''
      #
      numranges=1
      reslow = [10]
      numcount =[3]
      result = dihedral.get_flexible_residues(numranges,reslow,numcount)
      self.assertEqual(result,[10,11,12])

   def test_numranges2(self):
      '''
      test if numranges==2
      '''
      #
      numranges=2
      reslow = [10,20]
      numcount =[3,2]
      result = dihedral.get_flexible_residues(numranges,reslow,numcount)
      self.assertEqual(result,[10,11,12,20,21])

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

