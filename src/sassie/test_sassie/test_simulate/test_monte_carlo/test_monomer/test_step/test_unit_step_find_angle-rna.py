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

test for a 3-na system; find beta angle of na-1
test for a 3-na system; find gamma angle of na-1
test for a 3-na system; find epsilon angle of na-1
test for a 3-na system; find eta angle of na-1
test for a 3-na system; find alpha angle of na-2
test for a 3-na system; find beta angle of na-2
test for a 3-na system; find gamma angle of na-2
test for a 3-na system; find epsilon angle of na-2
test for a 3-na system; find eta angle of na-2
test for a 3-na system; find alpha angle of na-3
test for a 3-na system; find beta angle of na-3
test for a 3-na system; find gamma angle of na-3



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
      self.coor_3na=numpy.array([[[-12.978, -3.894, 12.958], [-14.435, -3.372, 12.66], [-15.582, -4.232, 12.858], [-16.887, -3.567, 12.384], [-17.019, -3.394, 10.864], [-17.415, -4.642, 10.261], [-17.389, -4.908, 8.719], [-18.267, -3.756, 8.092], [-19.632, -4.041, 7.692], [-20.403, -2.789, 7.233], [-20.027, -2.262, 5.836], [-20.688, -3.044, 4.81], [-20.323, -2.862, 3.284], [-20.517, -1.333, 2.992], [-21.78, -0.842, 2.489], [-21.785, 0.687, 2.353], [-20.923, 1.265, 1.226], [-21.576, 1.109, -0.049]]],numpy.float) # P-O5'-C5'-C4'-C3'-O3'...P-O5'-C5'-C4'-C3'-O3'...P-O5'-C5'-C4'-C3'-O3'
      self.first_last_resid_3aa = [1,3]
      self.indices_3na = None #This is an unecessary parameter
      self.seed_object = RandomState(1)
      self.molecule_type = 'rna'
      self.o = step.Setup()
      self.nonbondflag = 0
      self.anglelist=['alpha','beta','gamma','delta','epsilon','eta']
      talpha = [1.20,1,180.0,0.1,2,180.0,0.1,3,180.0,0.0,6,0.0]
      tbeta = [0.2, 1, 120.0]
      tgamma = [0.2,4,180.0,0.8,3,180.0,0.4,2,0.0,2.5,1,180.0]
      tdelta = [0.2,4,0.0,0.8,3,180.0]
      tepsilon = [0.6,5,0.0,0.2,4,0.0,0.0,3,180.0,0.4,2,0.0,1.9,1,180.0]
      teta = [1.20,1,180.0,0.10,2,180.0,0.1,3,180.0,0.0,6,0.0]
      self.parm = [talpha, tbeta, tgamma, tdelta, tepsilon, teta]

   def test_3na_firstres_beta(self):
      '''
      test for a 3-na system; find beta angle of na-1
      '''
      #
      this_dihedral = 1 #beta angle
      this_mask = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
      q0 = 1
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 0.284, 0.286,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_firstres_gamma(self):
      '''
      test for a 3-na system; find gamma angle of na-1
      '''
      #
      this_dihedral = 2 #gamma angle
      this_mask = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
      q0 = 1
      this_dtheta = 15.0
      beta = 1.0/(1.380658E-23*100)
      #
      expected = ( 3.398, 3.399,-2.489)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_firstres_epsilon(self):
      '''
      test for a 3-na system; find epsilon angle of na-1
      '''
      #
      this_dihedral = 4 #epsilon angle
      this_mask = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
      q0 = 1
      this_dtheta = 25.0
      beta = 1.0/(1.380658E-23*1000)
      #
      expected =( 5.152, 5.069,-4.149)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_firstres_eta(self):
      '''
      test for a 3-na system; find eta angle of na-1
      '''
      #
      this_dihedral = 5 #gamma angle
      this_mask = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
      q0 = 1
      this_dtheta = 10.0
      beta = 1.0/(1.380658E-23*100)
      #
      expected = ( 0.852, 0.888,-1.660)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res2_alpha(self):
      '''
      test for a 3-na system; find alpha angle of na-2
      '''
      #
      this_dihedral = 0 #alpha angle
      this_mask = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0]
      q0 = 2
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 1.670, 1.682,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res2_beta(self):
      '''
      test for a 3-na system; find beta angle of na-2
      '''
      #
      this_dihedral = 1 #beta angle
      this_mask = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0]
      q0 = 2
      this_dtheta = 15.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 0.317, 0.324,-2.489)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res2_gamma(self):
      '''
      test for a 3-na system; find gamma angle of na-2
      '''
      #
      this_dihedral = 2 #gamma angle
      this_mask = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0]
      q0 = 2
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 3.372, 3.377,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res2_epsilon(self):
      '''
      test for a 3-na system; find epsilon angle of na-2
      '''
      #
      this_dihedral = 4 #epsilon angle
      this_mask = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0]
      q0 = 2
      this_dtheta = 15.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 5.094, 5.051,-2.489)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res2_eta(self):
      '''
      test for a 3-na system; find eta angle of na-2
      '''
      #
      this_dihedral = 5 #epsilon angle
      this_mask = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0]
      q0 = 2
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 0.843, 0.861,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res3_alpha(self):
      '''
      test for a 3-na system; find alpha angle of na-3
      '''
      #
      this_dihedral = 0 #alpha angle
      this_mask = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
      q0 = 3
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 1.543, 1.556,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res3_beta(self):
      '''
      test for a 3-na system; find beta angle of na-3
      '''
      #
      this_dihedral = 1 #beta angle
      this_mask = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
      q0 = 3
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 0.310, 0.313,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def test_3na_res3_gamma(self):
      '''
      test for a 3-na system; find gamma angle of na-3
      '''
      #
      this_dihedral = 2 #gamma angle
      this_mask = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
      q0 = 3
      this_dtheta = 5.0
      beta = 1.0/(1.380658E-23*300)
      #
      expected = ( 3.393, 3.396,-0.830)
      results = self.o.find_angle(self.coor_3na, self.anglelist[this_dihedral], self.indices_3na, this_dihedral, q0, this_dtheta, beta, self.parm, this_mask, self.nonbondflag, self.first_last_resid_3aa, self.molecule_type, self.seed_object)
      print '(%6.3f,%6.3f,%6.3f)'%results
      self.assert_list_almost_equal(results,expected)

   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

