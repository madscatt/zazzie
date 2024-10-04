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

test a rna file with acceptable conditions
test rna with bad overlap
test rna with bad rg
test rna with bad z


'''

import os
import locale
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol
from sassie.simulate.monte_carlo.monomer import dihedral_rotate as nrotate


class Test_energy_nrotate_rotate(MockerTestCase): 


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
      self.PdbDataPath = os.path.dirname(os.path.realpath(__file__))+'/../../../../data/simulate/monte_carlo/monomer/dihedral_rotate/'
      self.o = sasmol.SasMol(0)


   def test_rna_reg(self):
      '''
      test a rna file with acceptable conditions
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      pdbid = 'rna_single'
      pdbFile = os.path.join(self.PdbDataPath,pdbid+'.pdb')
      dcdFile = os.path.join(self.PdbDataPath,pdbid+'.dcd')
      self.o.read_pdb(pdbFile)
      q0 = 16
      th = 5.0
      an = 'beta'
      cut = 3.0
      lowrg = 0.0
      highrg = 160.0
      re = [0,0,0,0,0,10.0,50.0,0,0,0]
      zflag = 0
      zval = 0.0
      cflag = 0
      basis = 'P'
      lowres1 = 1 
      highres1 = 16
      #may not modify
      if(os.path.isfile(dcdFile)):
         os.remove(dcdFile)
      dcdoutfile = self.o.open_dcd_write(dcdFile)
      taccepted = 0
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      first_last_resid = [1,self.o.number_of_resids()]
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      error,basis_mask = self.o.get_subset_mask('name[i] == "'+basis+'"')
      align_filter = 'name[i] == "'+basis+'" and (resid[i] >= '+str(lowres1)+' and resid[i] <= '+str(highres1)+')'   
      error,align_mask = self.o.get_subset_mask(align_filter)
      sub_m1=sasmol.SasMol(2)
      error = self.o.copy_molecule_using_mask(sub_m1,align_mask,0)
      com_sub_m1 = sub_m1.calccom(0)
      sub_m1.center(0)
      coor_sub_m1 = sub_m1.coor()[0]
      sub_m2 = sasmol.SasMol(4)
      error = self.o.copy_molecule_using_mask(sub_m2,align_mask,0)
      mask_a_array,mask_b_array,distance_array,type_array=[],[],[],[]
      #run
      filename = nrotate.rotate(coor,self.o,q0,th,an,cut,lowrg,highrg,re,taccepted,zflag,zval,cflag,dcdoutfile,indices,this_mask,basis_mask,sub_m2,align_mask,coor_sub_m1,com_sub_m1,mask_a_array,mask_b_array,distance_array,type_array,first_last_resid,molecular_type)
      print(re)
      self.o.close_dcd_write(dcdoutfile)
      #verify
      self.assertEqual(filename, 'winner')
      re_expected = [1, 0, 0, 14.999, 14.999, 10.0, 50.0, 0, 0, [numpy.array([-22.828, -15.981, -12.435]), numpy.array([ 12.806, 9.599, 24.177])]]
      self.assert_list_almost_equal(re, re_expected,3)

   def test_rna_overlap(self):
      '''
      test rna with bad overlap
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      pdbid = 'rna_single'
      pdbFile = os.path.join(self.PdbDataPath,pdbid+'.pdb')
      dcdFile = os.path.join(self.PdbDataPath,pdbid+'.dcd')
      self.o.read_pdb(pdbFile)
      q0 = 10
      th = 180.0
      an = 'beta'
      cut = 30.0
      lowrg = 0.0
      highrg = 160.0
      re = [0,0,0,0,0,10.0,50.0,0,0,0]
      zflag = 0
      zval = 0.0
      cflag = 0
      basis = 'P'
      lowres1 = 1 
      highres1 = 16
      #may not modify
      if(os.path.isfile(dcdFile)):
         os.remove(dcdFile)
      dcdoutfile = self.o.open_dcd_write(dcdFile)
      taccepted = 0
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      first_last_resid = [1,self.o.number_of_resids()]
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      error,basis_mask = self.o.get_subset_mask('name[i] == "'+basis+'"')
      align_filter = 'name[i] == "'+basis+'" and (resid[i] >= '+str(lowres1)+' and resid[i] <= '+str(highres1)+')'   
      error,align_mask = self.o.get_subset_mask(align_filter)
      sub_m1=sasmol.SasMol(2)
      error = self.o.copy_molecule_using_mask(sub_m1,align_mask,0)
      com_sub_m1 = sub_m1.calccom(0)
      sub_m1.center(0)
      coor_sub_m1 = sub_m1.coor()[0]
      sub_m2 = sasmol.SasMol(4)
      error = self.o.copy_molecule_using_mask(sub_m2,align_mask,0)
      mask_a_array,mask_b_array,distance_array,type_array=[],[],[],[]
      #run
      filename = nrotate.rotate(coor,self.o,q0,th,an,cut,lowrg,highrg,re,taccepted,zflag,zval,cflag,dcdoutfile,indices,this_mask,basis_mask,sub_m2,align_mask,coor_sub_m1,com_sub_m1,mask_a_array,mask_b_array,distance_array,type_array,first_last_resid,molecular_type)
      self.o.close_dcd_write(dcdoutfile)
      #verify
      self.assertEqual(filename, '')
      self.assert_list_almost_equal(re, [0,1,0,12.139,0.0,10.0,50.0,0,0,[]],3)


   def test_rna_badrg(self):
      '''
      test rna with bad rg
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      pdbid = 'rna_single'
      pdbFile = os.path.join(self.PdbDataPath,pdbid+'.pdb')
      dcdFile = os.path.join(self.PdbDataPath,pdbid+'.dcd')
      self.o.read_pdb(pdbFile)
      q0 = 10
      th = 180.0
      an = 'beta'
      cut = 3.0
      lowrg = 0.0
      highrg = 10.0
      re = [0,0,0,0,0,10.0,50.0,0,0,0]
      zflag = 0
      zval = 0.0
      cflag = 0
      basis = 'P'
      lowres1 = 1 
      highres1 = 16
      #may not modify
      if(os.path.isfile(dcdFile)):
         os.remove(dcdFile)
      dcdoutfile = self.o.open_dcd_write(dcdFile)
      taccepted = 0
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      first_last_resid = [1,self.o.number_of_resids()]
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      error,basis_mask = self.o.get_subset_mask('name[i] == "'+basis+'"')
      align_filter = 'name[i] == "'+basis+'" and (resid[i] >= '+str(lowres1)+' and resid[i] <= '+str(highres1)+')'   
      error,align_mask = self.o.get_subset_mask(align_filter)
      sub_m1=sasmol.SasMol(2)
      error = self.o.copy_molecule_using_mask(sub_m1,align_mask,0)
      com_sub_m1 = sub_m1.calccom(0)
      sub_m1.center(0)
      coor_sub_m1 = sub_m1.coor()[0]
      sub_m2 = sasmol.SasMol(4)
      error = self.o.copy_molecule_using_mask(sub_m2,align_mask,0)
      mask_a_array,mask_b_array,distance_array,type_array=[],[],[],[]
      #run
      filename = nrotate.rotate(coor,self.o,q0,th,an,cut,lowrg,highrg,re,taccepted,zflag,zval,cflag,dcdoutfile,indices,this_mask,basis_mask,sub_m2,align_mask,coor_sub_m1,com_sub_m1,mask_a_array,mask_b_array,distance_array,type_array,first_last_resid,molecular_type)
      self.o.close_dcd_write(dcdoutfile)
      #verify
      self.assertEqual(filename, '')
      self.assert_list_almost_equal(re, [0,1,0,12.139,0.0,10.0,50.0,0,0,[]],3)

   def test_rna_badz(self):
      '''
      test rna with bad z
      '''
      #to be modified
      mol_type = 1
      molecular_type = 'rna'
      pdbid = 'rna_single'
      pdbFile = os.path.join(self.PdbDataPath,pdbid+'.pdb')
      dcdFile = os.path.join(self.PdbDataPath,pdbid+'.dcd')
      self.o.read_pdb(pdbFile)
      q0 = 10
      th = 10.0
      an = 'beta'
      cut = 3.0
      lowrg = 0.0
      highrg = 50.0
      re = [0,0,0,0,0,10.0,50.0,0,0,0]
      zflag = 1
      zval = 1000.0
      cflag = 0
      basis = 'P'
      lowres1 = 1 
      highres1 = 16
      #may not modify
      if(os.path.isfile(dcdFile)):
         os.remove(dcdFile)
      dcdoutfile = self.o.open_dcd_write(dcdFile)
      taccepted = 0
      mask = self.o.get_dihedral_subset_mask([q0],mol_type)
      first_last_resid = [1,self.o.number_of_resids()]
      this_mask = mask[0]
      coor = self.o.coor()
      indices = self.o.get_indices_from_mask(this_mask)
      frame = 0
      error,basis_mask = self.o.get_subset_mask('name[i] == "'+basis+'"')
      align_filter = 'name[i] == "'+basis+'" and (resid[i] >= '+str(lowres1)+' and resid[i] <= '+str(highres1)+')'   
      error,align_mask = self.o.get_subset_mask(align_filter)
      sub_m1=sasmol.SasMol(2)
      error = self.o.copy_molecule_using_mask(sub_m1,align_mask,0)
      com_sub_m1 = sub_m1.calccom(0)
      sub_m1.center(0)
      coor_sub_m1 = sub_m1.coor()[0]
      sub_m2 = sasmol.SasMol(4)
      error = self.o.copy_molecule_using_mask(sub_m2,align_mask,0)
      mask_a_array,mask_b_array,distance_array,type_array=[],[],[],[]
      #run
      filename = nrotate.rotate(coor,self.o,q0,th,an,cut,lowrg,highrg,re,taccepted,zflag,zval,cflag,dcdoutfile,indices,this_mask,basis_mask,sub_m2,align_mask,coor_sub_m1,com_sub_m1,mask_a_array,mask_b_array,distance_array,type_array,first_last_resid,molecular_type)
      self.o.close_dcd_write(dcdoutfile)
      #verify
      self.assertEqual(filename, 'winner')
      self.assert_list_almost_equal(re, [0,0,0,14.646,0.0,10.0,50.0,1,0,[]],3)


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

