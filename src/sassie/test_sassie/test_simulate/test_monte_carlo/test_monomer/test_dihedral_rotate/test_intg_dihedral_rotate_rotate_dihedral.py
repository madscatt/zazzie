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


test a 2-aa pdb file; rotate the psi angle of aa-1 for 0deg
test a 2-aa pdb file; rotate the psi angle of aa-1 for 90deg
test a 2-aa pdb file; rotate the psi angle of aa-1 for 180deg
test a 2-aa pdb file; rotate the phi angle of aa-2 for 0deg
test a 2-aa pdb file; rotate the phi angle of aa-2 for 90deg
test a 2-aa pdb file; rotate the phi angle of aa-2 for 180deg
test a 3-aa pdb file; rotate the psi angle of aa-1 for 72deg
test a 3-aa pdb file; rotate the phi angle of aa-2 for 86deg
test a 3-aa pdb file; rotate the psi angle of aa-2 for 108deg
test a 3-aa pdb file; rotate the phi angle of aa-3 for 168deg
test a regular pdb file (1CRN) ; rotate the phi angle of aa-30 for 68deg
test a single strand rna; rotate the beta angle of na-1 for 68deg
test a single strand rna; rotate the gamma angle of na-1 for 168deg
test a single strand rna; rotate the epsilon angle of na-1 for 10deg
test a single strand rna; rotate the eta angle of na-1 for 170deg
test a single strand rna; rotate the alpha angle of na-10 for -10deg
test a single strand rna; rotate the eta angle of na-10 for -170deg
test a single strand rna; rotate the eta angle of na-last for 375deg
test a single strand rna; rotate the eta angle of na-last for 365deg



'''

import os
import locale
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol
from sassie.simulate.monte_carlo.monomer import dihedral_rotate



class Test_energy_dihedral_rotate_rotate_dihedral(MockerTestCase): 


   def setUp(self):
      self.PdbDataPath = os.path.dirname(os.path.realpath(__file__))+'/../../../../data/simulate/monte_carlo/monomer/dihedral_rotate/'
      self.o = sasmol.SasMol(0)

   def test_2aa_psi1_0deg(self):
      '''
      test a 2-aa pdb file; rotate the psi angle of aa-1 for 0deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 1
      itheta = 0.0
      an = 'psi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2407.676, places=3) #The new pdb has been visually checked in vmd

   def test_2aa_psi1_90deg(self):
      '''
      test a 2-aa pdb file; rotate the psi angle of aa-1 for 90deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 1
      itheta = 90.0
      an = 'psi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2388.930, places=3) #The new pdb has been visually checked in vmd

   def test_2aa_psi1_180deg(self):
      '''
      test a 2-aa pdb file; rotate the psi angle of aa-1 for 180deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 1
      itheta = 180.0
      an = 'psi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2380.309, places=3) #The new pdb has been visually checked in vmd


   def test_2aa_phi2_0deg(self):
      '''
      test a 2-aa pdb file; rotate the phi angle of aa-2 for 0deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 2
      itheta = 0.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2407.676, places=3) #The new pdb has been visually checked in vmd

   def test_2aa_phi2_90deg(self):
      '''
      test a 2-aa pdb file; rotate the phi angle of aa-2 for 90deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 2
      itheta = 90.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2401.509, places=3) #The new pdb has been visually checked in vmd

   def test_2aa_phi2_180deg(self):
      '''
      test a 2-aa pdb file; rotate the phi angle of aa-2 for 180deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'2AAD.pdb'))
      q0 = 2
      itheta = 180.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([1,2],mol_type)
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,2396.937, places=3) #The new pdb has been visually checked in vmd

   def test_3aa_psi1_72deg(self):
      '''
      test a 3-aa pdb file; rotate the psi angle of aa-1 for 72deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      q0 = 1
      itheta = 72.0
      an = 'psi'
      mask = self.o.get_dihedral_subset_mask([1,2,3],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,3488.521, places=3) #The new pdb has been visually checked in vmd

   def test_3aa_phi2_72deg(self):
      '''
      test a 3-aa pdb file; rotate the phi angle of aa-2 for 86deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      q0 = 2
      itheta = 86.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([1,2,3],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,3511.518, places=3) #The new pdb has been visually checked in vmd

   def test_3aa_psi2_108deg(self):
      '''
      test a 3-aa pdb file; rotate the psi angle of aa-2 for 108deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      q0 = 2
      itheta = 108.0
      an = 'psi'
      mask = self.o.get_dihedral_subset_mask([1,2,3],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,3557.158, places=3) #The new pdb has been visually checked in vmd

   def test_3aa_phi3_168deg(self):
      '''
      test a 3-aa pdb file; rotate the phi angle of aa-3 for 168deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'3AAD.pdb'))
      q0 = 3
      itheta = 168.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([1,2,3],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[q0-1]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,3548.450, places=3) #The new pdb has been visually checked in vmd

   def test_1CRN_phi30_68deg(self):
      '''
      test a regular pdb file (1CRN) ; rotate the phi angle of aa-30 for 68deg
      '''
      #to be modified
      mol_type = 0
      molecular_type = 'protein'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'1CRN.pdb'))
      q0 = 30
      itheta = 68.0
      an = 'phi'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,8663.450, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_beta1_68deg(self):
      '''
      test a single strand rna; rotate the beta angle of na-1 for 68deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 1
      itheta = 68.0
      an = 'beta'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,8024.782, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_gamma1_168deg(self):
      '''
      test a single strand rna; rotate the gamma angle of na-1 for 168deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 1
      itheta = 168.0
      an = 'gamma'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,10942.108, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_epsilon1_10deg(self):
      '''
      test a single strand rna; rotate the epsilon angle of na-1 for 10deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 1
      itheta = 10.0
      an = 'epsilon'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-310.530, places=3) #The new pdb has been visually checked in vmd


   def test_rna_single_eta1_170deg(self):
      '''
      test a single strand rna; rotate the eta angle of na-1 for 170deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 1
      itheta = 170.0
      an = 'eta'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-10876.761, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_alpha10_n10deg(self):
      '''
      test a single strand rna; rotate the alpha angle of na-10 for -10deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 10
      itheta = -10.0
      an = 'alpha'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-2098.086, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_eta10_n170deg(self):
      '''
      test a single strand rna; rotate the eta angle of na-10 for -170deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 10
      itheta = -170.0
      an = 'eta'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-1663.490, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_alpha16_375deg(self):
      '''
      test a single strand rna; rotate the eta angle of na-last for 375deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 16
      itheta = 375.0
      an = 'alpha'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-1772.984, places=3) #The new pdb has been visually checked in vmd

   def test_rna_single_gama16_300365deg(self):
      '''
      test a single strand rna; rotate the eta angle of na-last for 365deg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      self.o.read_pdb(os.path.join(self.PdbDataPath,'rna_single.pdb'))
      q0 = 16
      itheta = 300365.0
      an = 'gamma'
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      print 'mask \n',mask
      #may not modify
      first_last_resid = [1,self.o.number_of_resids()]
      print first_last_resid
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      dihedral_rotate.rotate_dihedral(coor,self.o,frame,q0,itheta,an,indices,this_mask,first_last_resid,molecular_type)
      #self.o.write_pdb('out.pdb',frame,'w')
      result_static = coor.sum()
      self.assertAlmostEqual(result_static,-1615.201, places=3) #The new pdb has been visually checked in vmd


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

