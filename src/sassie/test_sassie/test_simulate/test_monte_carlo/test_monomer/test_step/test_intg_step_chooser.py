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

choose from a 3aa system with every flexible
choose from a small protein (1CRN) with 1 flexible region
choose from a small protein (1CRN) with 2 flexible regions
choose from a small protein (1CRN) with 3 flexible regions
choose from an rna molecule with 1 flexible region
choose from an rna molecule with 2 flexible regions
choose from an rna molecule with 3 flexible regions



'''

import os
import numpy
from numpy.random import RandomState

import multiprocessing

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol
from sassie.simulate.monte_carlo.monomer import dihedral_monte_carlo,step
from sassie.simulate.energy import dihedral_energy




class Test_step_chooser(MockerTestCase): 

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
      self.ochooser = step.Setup()
      self.nonbondflag = 0
      self.seed_object = RandomState(1)
      self.pairdat=['dum',1,1.0]


   def init_dihedral_parameters(self,o,basis,numranges,reslow,numcont,first_last_resid):
      respsi=[] ; resphi=[]
      basis_filter = 'name[i] == "'+basis+'"'
      error,basis_mask = o.get_subset_mask(basis_filter)
      basis_m1=sasmol.SasMol(1)
      error = o.copy_molecule_using_mask(basis_m1,basis_mask,0)
      basis_resid = basis_m1.resid()
      basis_resname = basis_m1.resname()
      txtQueue=multiprocessing.JoinableQueue()
      dihedral_energy.protein_initialization(respsi,resphi,basis_resid,basis_resname,numranges,reslow,numcont,first_last_resid,txtQueue)
      dihedral_parameters = [respsi,resphi]
      return dihedral_parameters

   def init_rna_dihedral_parameters(self,o,basis,numranges,reslow,numcont,first_last_resid):
      resalpha = [] ; resbeta = [] ; resgamma = [] ; resdelta = [] ; resepsilon = [] ; reseta = []
      basis_filter = 'name[i] == "'+basis+'"'
      error,basis_mask = o.get_subset_mask(basis_filter)
      basis_m1=sasmol.SasMol(1)
      error = o.copy_molecule_using_mask(basis_m1,basis_mask,0)
      basis_resid = basis_m1.resid()
      basis_resname = basis_m1.resname()
      txtQueue=multiprocessing.JoinableQueue()
      dihedral_energy.rna_initialization(resalpha,resbeta,resgamma,resdelta,resepsilon,reseta,basis_resid,basis_resname,numranges,reslow,numcont,first_last_resid,txtQueue)
      dihedral_parameters = [resalpha,resbeta,resgamma,resdelta,resepsilon,reseta]
      return dihedral_parameters


   def get_rotation_indices_masks(self,m1,molecule_type,numranges,reslow,numcont,flexible_residues):
      flexible_residues = dihedral_monte_carlo.get_flexible_residues(numranges,reslow,numcont)
      txtOutput = None
      residue_rotation_indices,residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(m1,molecule_type,flexible_residues,txtOutput)
      return residue_rotation_indices,residue_rotation_mask


   def test_3aa(self):
      '''
      choose from a 3aa system with every flexible
      '''
      # to-be-modified
      pdb = '3AAD.pdb'
      numranges = 1
      reslow = [1]
      numcont = [2]
      flexible_residues = range(reslow[0],numcont[0]+2)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'CA'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'protein'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      expected = [0.429, 0.409 , [6, 8, 9, 13, 15]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)


   def test_1CRN_1flexible(self):
      '''
      choose from a small protein (1CRN) with 1 flexible region
      '''
      # to-be-modified
      pdb = '1CRN.pdb'
      numranges = 1
      reslow = [1]
      numcont = [6]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+2)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'CA'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'protein'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices
      expected = [-0.099,-0.142, [22, 26, 27, 28, 33]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)


   def test_1CRN_2flexible(self):
      '''
      choose from a small protein (1CRN) with 2 flexible regions
      '''
      # to-be-modified
      pdb = '1CRN.pdb'
      numranges = 2
      reslow = [1,16]
      numcont = [6,6]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+1)+range(reslow[1],reslow[1]+numcont[1]+1)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'CA'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'protein'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices
      expected = [-0.099,-0.142, [22, 26, 27, 28, 33]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)

   def test_1CRN_3flexible(self):
      '''
      choose from a small protein (1CRN) with 3 flexible regions
      '''
      # to-be-modified
      pdb = '1CRN.pdb'
      numranges = 3
      reslow = [3,16,30]
      numcont = [6,6,3]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+1)+range(reslow[1],reslow[1]+numcont[1]+1)+range(reslow[2],reslow[2]+numcont[2]+1)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'CA'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'protein'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices,flexible_residues
      expected = [0.432, 0.452, [137, 142, 143, 144, 146]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)


   def test_rna_1flexible(self):
      '''
      choose from an rna molecule with 1 flexible region
      '''
      # to-be-modified
      pdb = 'rna_single.pdb'
      numranges = 1
      reslow = [1]
      numcont = [6]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+2)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'P'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'rna'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_rna_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices
      expected = [0.914, 0.956, [131, 132, 135, 136, 139, 159, 161, 162, 165]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)


   def test_rna_2flexible(self):
      '''
      choose from an rna molecule with 2 flexible regions
      '''
      # to-be-modified
      pdb = 'rna_single.pdb'
      numranges = 2
      reslow = [1,16]
      numcont = [6,6]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+1)+range(reslow[1],reslow[1]+numcont[1]+1)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'P'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'rna'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_rna_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices
      expected = [0.914, 0.956, [131, 132, 135, 136, 139, 159, 161, 162, 165]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)


   def test_rna_3flexible(self):
      '''
      choose from an rna molecule with 3 flexible regions
      '''
      # to-be-modified
      pdb = 'rna_single.pdb'
      numranges = 2
      reslow = [1,16]
      numcont = [6,6]
      flexible_residues = range(reslow[0],reslow[0]+numcont[0]+1)+range(reslow[1],reslow[1]+numcont[1]+1)
      # may not need modify
      dtheta = [5.0]*numranges
      basis = 'P'
      beta = 1.0/(1.380658E-23*300)
      molecule_type = 'rna'
      # DO NOT modify
      self.o.read_pdb(os.path.join(self.PdbDataPath,pdb))
      coor=self.o.coor()
      first_last_resid = [1,self.o.number_of_resids()]
      dihedral_parameters = self.init_rna_dihedral_parameters(self.o,basis,numranges,reslow,numcont,first_last_resid)
      residue_rotation_indices,residue_rotation_mask = self.get_rotation_indices_masks(self.o, molecule_type,numranges,reslow,numcont,flexible_residues)
      vdi,vdf,indices,this_mask = self.ochooser.chooser(coor,self.o,self.pairdat,dtheta,numranges,reslow,numcont,dihedral_parameters,beta,residue_rotation_indices,residue_rotation_mask,self.nonbondflag,first_last_resid,molecule_type,self.seed_object)
      print '%6.3f,%6.3f,'%(vdi,vdf),indices
      expected = [0.914, 0.956, [131, 132, 135, 136, 139, 159, 161, 162, 165]]
      self.assert_list_almost_equal([vdi,vdf,indices],expected,3)

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

