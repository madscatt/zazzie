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

test for a 3-na system; measure beta angle of the first residue
test for a 3-na system; measure gamma angle of the first residue
test for a 3-na system; measure delta angle of the first residue
test for a 3-na system; measure epsilon angle of the first residue
test for a 3-na system; measure eta angle of the first residue
test for a 3-na system; measure alpha angle of the second residue
test for a 3-na system; measure beta angle of the second residue
test for a 3-na system; measure gamma angle of the second residue
test for a 3-na system; measure delta angle of the second residue
test for a 3-na system; measure epsilon angle of the second residue
test for a 3-na system; measure eta angle of the second residue
test for a 3-na system; measure alpha angle of the third residue
test for a 3-na system; measure beta angle of the third residue
test for a 3-na system; measure gamma angle of the third residue



'''

import os
import locale
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import dihedral_rotate


class Test_dihedral_rotate_measure(MockerTestCase): 


   def setUp(self):
      self.coor_3na=numpy.array([[[-12.978, -3.894, 12.958], [-14.435, -3.372, 12.66], [-15.582, -4.232, 12.858], [-16.887, -3.567, 12.384], [-17.019, -3.394, 10.864], [-17.415, -4.642, 10.261], [-17.389, -4.908, 8.719], [-18.267, -3.756, 8.092], [-19.632, -4.041, 7.692], [-20.403, -2.789, 7.233], [-20.027, -2.262, 5.836], [-20.688, -3.044, 4.81], [-20.323, -2.862, 3.284], [-20.517, -1.333, 2.992], [-21.78, -0.842, 2.489], [-21.785, 0.687, 2.353], [-20.923, 1.265, 1.226], [-21.576, 1.109, -0.049]]],numpy.float)  # P-O5'-C5'-C4'-C3'-O3'...P-O5'-C5'-C4'-C3'-O3'...P-O5'-C5'-C4'-C3'-O3'
      self.first_last_resid_3na = [1,3]
      self.indices_3na = None



   def test_3na_firstres_beta(self):
      '''
      test for a 3-na system; measure beta angle of the first residue
      '''
      #
      an = 'beta'
      this_mask = [1]*6 + [1,1]+[0]*4 + [0]*6
      q0 = 1
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -174.70 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_firstres_gama(self):
      '''
      test for a 3-na system; measure gamma angle of the first residue
      '''
      #
      an = 'gamma'
      this_mask = [1]*6 + [1,1]+[0]*4 + [0]*6
      q0 = 1
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 69.55 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_firstres_delta(self):
      '''
      test for a 3-na system; measure delta angle of the first residue
      '''
      #
      an = 'delta'
      this_mask = [1]*6 + [1,1]+[0]*4 + [0]*6
      q0 = 1
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 79.34 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_firstres_epsilon(self):
      '''
      test for a 3-na system; measure epsilon angle of the first residue
      '''
      #
      an = 'epsilon'
      this_mask = [1]*6 + [1,1]+[0]*4 + [0]*6
      q0 = 1
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -168.28 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_firstres_eta(self):
      '''
      test for a 3-na system; measure eta angle of the first residue
      '''
      #
      an = 'eta'
      this_mask = [1]*6 + [1,1]+[0]*4 + [0]*6
      q0 = 1
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -55.43 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_secondres_alpha(self):
      '''
      test for a 3-na system; measure alpha angle of the second residue
      '''
      #
      an = 'alpha'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -101.16 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_secondres_beta(self):
      '''
      test for a 3-na system; measure beta angle of the second residue
      '''
      #
      an = 'beta'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 174.15 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_secondres_gamma(self):
      '''
      test for a 3-na system; measure gamma angle of the second residue
      '''
      #
      an = 'gamma'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 75.05 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_secondres_delta(self):
      '''
      test for a 3-na system; measure delta angle of the second residue
      '''
      #
      an = 'delta'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 80.03 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)


   def test_3na_secondres_epsilon(self):
      '''
      test for a 3-na system; measure epsilon angle of the second residue
      '''
      #
      an = 'epsilon'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -171.05 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_secondres_eta(self):
      '''
      test for a 3-na system; measure eta angle of the second residue
      '''
      #
      an = 'eta'
      this_mask = [0]*5+[1] + [1]*6 + [1,1]+[0]*4
      q0 = 2
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -55.01 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_thirdres_alpha(self):
      '''
      test for a 3-na system; measure alpha angle of the third residue
      '''
      #
      an = 'alpha'
      this_mask = [0]*6 + [0]*5+[1] + [1]*6
      q0 = 3
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = -92.79 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_thirdres_beta(self):
      '''
      test for a 3-na system; measure beta angle of the third residue
      '''
      #
      an = 'beta'
      this_mask = [0]*6 + [0]*5+[1] + [1]*6
      q0 = 3
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 176.55 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3na_thirdres_gamma(self):
      '''
      test for a 3-na system; measure gamma angle of the third residue
      '''
      #
      an = 'gamma'
      this_mask = [0]*6 + [0]*5+[1] + [1]*6
      q0 = 3
      molecule_type = 'rna'
      result = dihedral_rotate.measure(self.coor_3na, self.indices_3na, an, this_mask, q0, self.first_last_resid_3na, molecule_type)
      expected = 71.05 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

