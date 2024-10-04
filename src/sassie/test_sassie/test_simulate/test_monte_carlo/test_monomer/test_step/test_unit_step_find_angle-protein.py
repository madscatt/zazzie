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


test for a 3-aa system; find psi angle of the first residue
test for a 3-aa system; find phi angle of the 2nd residue
test for a 3-aa system; find psi angle of the 2nd residue
test for a 3-aa system; find phi angle of the 3rd residue


'''

import os
import locale
import numpy
from numpy.random import RandomState

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import step


class Test_unit_step_find_angle(MockerTestCase): 

   def assert_list_almost_equal(self,a,b,places=3):
      if (len(a)!=len(b)):
         raise TypeError
      else:
         for i in range(len(a)):
            if isinstance(a[i],(int,float,numpy.generic)):
               if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
               self.assertAlmostEqual(a[i],b[i],places)
            else:
               self.assert_list_almost_equal(a[i],b[i],places)

   def setUp(self):
      self.coor_3aa=numpy.array([[[13.725, 11.174, 16.425],[13.257, 10.745, 15.081],[14.275, 9.687, 14.612],[14.342, 8.64, 15.422],[15.445, 7.667, 15.246],[15.171, 6.533, 14.28 ], [13.966, 6.502, 13.739],[13.512, 5.395, 12.878],[13.311, 5.853, 11.455]]],numpy.float) # N-CA-C...N-CA-C...N-CA-C
      self.first_last_resid_3aa = [1,3]
      self.indices_3aa = None #This is an unecessary parameter
      self.seed_object = RandomState(1)
      self.molecule_type = 'protein'
      self.o = step.Setup()
      self.anglelist=['phi','psi']
      self.parm = [[0.2,1.0,180.0,0.0,0.0,0.0],[0.6,1.0,0.0,0.0,0.0,0.0]]
      self.nonbondflag = 0

   def test_3aa_firstres_psi(self):
      '''
      test for a 3-aa system; find psi angle of the first residue
      '''
      #
      this_dihedral = 1 #psi angle
      this_mask = [1,1,1,1,0,0,0,0,0]
      q0 = 1
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 0.900, 0.907,-0.830)
      results = self.o.find_angle(self.coor_3aa, self.anglelist[this_dihedral], self.indices_3aa, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3aa_secondres_phi(self):
      '''
      test for a 3-aa system; find phi angle of the 2nd residue
      '''
      #
      this_dihedral = 0 #psi angle
      this_mask = [0,0,1,1,1,1,1,0,0]
      q0 = 2
      this_dtheta = 15.0
      beta = 1.0/(1.380658E-23*100)
      #
      expected = ( 0.196, 0.204,-2.489)
      results = self.o.find_angle(self.coor_3aa, self.anglelist[this_dihedral], self.indices_3aa, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3aa_secondres_psi(self):
      '''
      test for a 3-aa system; find psi angle of the 2nd residue
      '''
      #
      this_dihedral = 1 #psi angle
      this_mask = [0,0,1,1,1,1,1,0,0]
      q0 = 2
      this_dtheta = 25.0
      beta = 1.0/(1.380658E-23*1000)
      #
      expected = ( 1.199, 1.196,-4.149)
      results = self.o.find_angle(self.coor_3aa, self.anglelist[this_dihedral], self.indices_3aa, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)


   def test_3aa_thirdres_psi(self):
      '''
      test for a 3-aa system; find phi angle of the 3rd residue
      '''
      #
      this_dihedral = 0 #psi angle
      this_mask = [0,0,0,0,0,1,1,1,1]
      q0 = 3
      this_dtheta = 35.0
      beta = 1.0/(1.380658E-23*10)
      #
      expected = ( 0.278, 0.206,21.052)
      results = self.o.find_angle(self.coor_3aa, self.anglelist[this_dihedral], self.indices_3aa, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)


   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

