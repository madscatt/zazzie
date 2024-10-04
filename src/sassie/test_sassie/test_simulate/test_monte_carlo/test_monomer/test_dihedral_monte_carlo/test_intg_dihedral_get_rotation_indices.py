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

test a 2-aa pdb file; get the rotation indices for residue-1
test a 2-aa pdb file; get the rotation indices for residue-2
test a 2-aa pdb file; get the rotation indices for residue-1/2
test a 3-aa pdb file; get the rotation indices for residue-1/2
test a 3-aa pdb file; get the rotation indices for residue-2/3
test a 3-aa pdb file; get the rotation indices for residue-1/2/3
test a regular file; get the rotation indices for residue-1-2,10-12,21-23


'''

import os
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol
from sassie.simulate.monte_carlo.monomer import dihedral_monte_carlo



class Test_dihedral_monte_carlo_get_rotation_indices(MockerTestCase): 

   def assert_list_almost_equal(self,a,b,places=5):
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
      self.PdbDataPath = os.path.dirname(os.path.realpath(__file__))+'/../../../../data/simulate/monte_carlo/monomer/dihedral_monte_carlo/'
      self.o = sasmol.SasMol(0)
      self.txtOutput = None

   def test_2aa_res1(self):
      '''
      test a 2-aa pdb file; get the rotation indices for residue-1
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      flexible_residues = [1]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {1:[0,1,6,8]})
      self.assertEqual(residue_rotation_mask, {1:[1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]})


   def test_2aa_res2(self):
      '''
      test a 2-aa pdb file; get the rotation indices for residue-2
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      flexible_residues = [2]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {2:[6,8,9,13]})
      self.assertEqual(residue_rotation_mask, {2:[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0]})

   def test_2aa_res1and2(self):
      '''
      test a 2-aa pdb file; get the rotation indices for residue-1/2
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      flexible_residues = [1,2]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {1:[0,1,6,8],2:[6,8,9,13]})
      self.assertEqual(residue_rotation_mask, {1:[1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],2:[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0]})


   def test_3aa_res1and2(self):
      '''
      test a 3-aa pdb file; get the rotation indices for residue-1/2
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      flexible_residues = [1,2]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {1:[0,1,6,8],2:[6,8,9,13,15]})
      self.assertEqual(residue_rotation_mask, {1:[1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0],2:[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0,1,0,0,0,0,0,0]})


   def test_3aa_res2and3(self):
      '''
      test a 3-aa pdb file; get the rotation indices for residue-2/3
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      flexible_residues = [2,3]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {2:[6,8,9,13,15],3:[13,15,16,20]})
      self.assertEqual(residue_rotation_mask, {2:[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0,1,0,0,0,0,0,0],3:[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,1,1,0,0,0,1,0]})

   def test_3aa_res1to3(self):
      '''
      test a 3-aa pdb file; get the rotation indices for residue-1/2/3
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      flexible_residues = [1,2,3]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {1:[0,1,6,8],2:[6,8,9,13,15],3:[13,15,16,20]})
      self.assertEqual(residue_rotation_mask, {1:[1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0],2:[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0,1,0,0,0,0,0,0],3:[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,1,1,0,0,0,1,0]})


   def test_1CRN_3flexibleregions(self):
      '''
      test a regular file; get the rotation indices for residue-1-2,10-12,21-23
      '''
      #to be modified
      molecule_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'1CRN.pdb'))
      flexible_residues = [1,2,10,11,12,21,22,23]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      self.assertEqual(residue_rotation_indices, {1: [0, 1, 2, 7], 2: [2, 7, 8, 9, 14], 10: [56, 59, 60, 61, 70], 11: [61, 70, 71, 72, 76], 12: [72, 76, 77, 78, 84], 21: [144, 146, 147, 148, 153], 22: [148, 153, 154, 155, 160], 23: [155, 160, 161, 162, 169]})
      expected_residue_rotation_mask = {}
      for (key,value) in zip(residue_rotation_indices.keys(),residue_rotation_indices.values()):
         self.assert_list_almost_equal(map(int,residue_rotation_mask[key]),[1 if i in value else 0 for i in range(327)])

   def test_rna_2flexibleregions(self):
      '''
      test a rna file; get the rotation indices for residue-5-6,10-12
      '''
      #to be modified
      molecule_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      flexible_residues = [5,6,10,11,12]
      #run
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(self.o,molecule_type,flexible_residues,self.txtOutput)
      #verify
      print residue_rotation_indices
      self.assertEqual(residue_rotation_indices, {10: [289, 290, 293, 294, 297, 317, 319, 320, 323], 11: [319, 320, 323, 324, 327, 351, 353, 354, 357], 12: [353, 354, 357, 358, 361, 384, 386, 387, 390], 5: [131, 132, 135, 136, 139, 159, 161, 162, 165], 6: [161, 162, 165, 166, 169, 190, 192, 193, 196]})
      expected_residue_rotation_mask = {}
      for (key,value) in zip(residue_rotation_indices.keys(),residue_rotation_indices.values()):
         self.assert_list_almost_equal(map(int,residue_rotation_mask[key]),[1 if i in value else 0 for i in range(518)])

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

