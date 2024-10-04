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

test for a 3-aa system; measure phi angle of the first residue (it will return psi angle)
test for a 3-aa system; measure psi angle of the first residue
test for a 3-aa system; measure the phi angle of the second residue
test for a 3-aa system; measure the psi angle of the second residue
test for a 3-aa system; measure the phi angle of the last residue
test for a 3-aa system; measure the psi angle of the last residue (it will measure the psi angle)


'''

import os
import locale
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import dihedral_rotate


class Test_dihedral_rotate_measure(MockerTestCase): 


   def setUp(self):
      self.coor_3aa=numpy.array([[[13.725, 11.174, 16.425],[13.257, 10.745, 15.081],[14.275, 9.687, 14.612],[14.342, 8.64, 15.422],[15.445, 7.667, 15.246],[15.171, 6.533, 14.28 ], [13.966, 6.502, 13.739],[13.512, 5.395, 12.878],[13.311, 5.853, 11.455]]],numpy.float) # N-CA-C...N-CA-C...N-CA-C
      self.first_last_resid_3aa = [1,3]
      self.indices_3aa = None




   """
   #Will not test this case, since phi angle will never be picked up for the first residue
   def test_3aa_firstres_phi(self):
      '''
      test for a 3-aa system; measure phi angle of the first residue (it will return psi angle)
      '''
      #
      an = 'phi'
      this_mask = [1,1,1,1,0,0,0,0,0]
      q0 = 1
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = 60.008 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=3)
   """

   def test_3aa_firstres_psi(self):
      '''
      test for a 3-aa system; measure psi angle of the first residue
      '''
      #
      an = 'psi'
      this_mask = [1,1,1,1,0,0,0,0,0]
      q0 = 1
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = 60.008 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=3)


   def test_3aa_secondres_phi(self):
      '''
      test for a 3-aa system; measure the phi angle of the second residue
      '''
      #
      an = 'phi'
      this_mask = [0,0,1,1,1,1,1,0,0]
      q0 = 2
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = -88.77 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3aa_secondres_psi(self):
      '''
      test for a 3-aa system; measure the psi angle of the second residue
      '''
      #
      an = 'psi'
      this_mask = [0,0,1,1,1,1,1,0,0]
      q0 = 2
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = -2.53 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   def test_3aa_lastres_phi(self):
      '''
      test for a 3-aa system; measure the phi angle of the last residue
      '''
      #
      an = 'phi'
      this_mask = [0,0,0,0,0,1,1,1,1]
      q0 = 3
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = -112.85 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)

   """
   # Will not test this case, since psi will never be picked for the last residue
   def test_3aa_lastres_psi(self):
      '''
      test for a 3-aa system; measure the psi angle of the last residue (it will measure the psi angle)
      '''
      #
      an = 'psi'
      this_mask = [0,0,0,0,0,1,1,1,1]
      q0 = 3
      molecule_type = 'protein'
      result = dihedral_rotate.measure(self.coor_3aa, self.indices_3aa, an, this_mask, q0, self.first_last_resid_3aa, molecule_type)
      expected = -112.85 #from vmd dihedral
      self.assertAlmostEqual(result,expected,places=2)
   """

   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

