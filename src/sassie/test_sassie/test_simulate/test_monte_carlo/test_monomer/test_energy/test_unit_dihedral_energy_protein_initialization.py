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

test for a 18-aa system with no PRO/GLY, 1 flexible range in middle
test for a 18-aa system with no PRO/GLY, 1 flexible range at nter
test for a 18-aa system with no PRO/GLY, 1 flexible range at cter
test for a 18-aa system with PRO in the middle, and 1 flexible range in middle
test for a 18-aa system with PRO at nter, 1 flexible range at nter
test for a 18-aa system with PRO at cter, 1 flexible range at cter
test for a 18-aa system with GLY in the middle, and 1 flexible range in middle
test for a 18-aa system with GLY at nter, 1 flexible range at nter
test for a 18-aa system with GLY at cter, 1 flexible range at cter
test for a 18-aa system with 3 flexible ranges at n-ter, middle and cter respectively



'''

import os
import locale

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy







class Test_dihedral_energy_protein_initialization(MockerTestCase): 


   def setUp(self):
      pass

   def test_18aa_1middle(self):
      '''
      test for a 18-aa system with no PRO/GLY, 1 flexible range in middle
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0],[0.2,1.0,180.0,0.0,0.0,0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0],[0.6,1.0,0.0,0.0,0.0,0.0]]]]
      print 'exp ',expected_resphi
      print 'res ',resphi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_1nter(self):
      '''
      test for a 18-aa system with no PRO/GLY, 1 flexible range at nter
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [1]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[1E90,0.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0],[0.2,1.0,180.0,0.0,0.0,0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0],[0.6,1.0,0.0,0.0,0.0,0.0]]]]
      print 'exp ',expected_resphi
      print 'res ',resphi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_1cter(self):
      '''
      test for a 18-aa system with no PRO/GLY, 1 flexible range at cter
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [first_last_resid[1]-3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0],[0.2,1.0,180.0,0.0,0.0,0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0],[1e+90, 0.0, 0.0, 0.0, 0.0, 0.]]]]
      print 'exp ',expected_respsi
      print 'res ',respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)

   def test_18aa_PRO_1middle(self):
      '''
      test for a 18-aa system with PRO in the middle, and 1 flexible range in middle
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','PRO','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.8, 3.0, 0.0, 0.0, 0.0, 0.0],[0.2,1.0,180.0,0.0,0.0,0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0], [0.4, 1.0, 0.0, 0.0, 0.0, 0.0], [0.3, 1.0, 0.0, -0.3, 4.0, 0.0],[0.6,1.0,0.0,0.0,0.0,0.0]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_PRO_1nter(self):
      '''
      test for a 18-aa system with PRO at nter, 1 flexible range at nter
      '''
      #
      resname = ['PRO','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [1]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[1e+90, 0.0, 0.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]]]
      expected_respsi = [[0, [[0.3, 1.0, 0.0, -0.3, 4.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_PRO_1cter(self):
      '''
      test for a 18-aa system with PRO at cter, 1 flexible range at cter
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','PRO']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [first_last_resid[1]-3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0],[0.8, 3.0, 0.0, 0.0, 0.0, 0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0], [0.6,1.0,0.0,0.0,0.0,0.0], [0.4, 1.0, 0.0, 0.0, 0.0, 0.0], [1e+90, 0.0, 0.0, 0.0, 0.0, 0.]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)

   def test_18aa_GLY_1middle(self):
      '''
      test for a 18-aa system with GLY in the middle, and 1 flexible range in middle
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','GLY','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2,1.0,180.0,0.0,0.0,0.0], [0.2,1.0,180.0,0.0,0.0,0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0],[0.2,1.0,180.0,0.0,0.0,0.0]]]]
      expected_respsi = [[0, [[0.6,1.0,0.0,0.0,0.0,0.0],[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0],[0.6,1.0,0.0,0.0,0.0,0.0]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_GLY_1nter(self):
      '''
      test for a 18-aa system with GLY at nter, 1 flexible range at nter
      '''
      #
      resname = ['GLY','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [1]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[1e+90, 0.0, 0.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]]]
      expected_respsi = [[0, [[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_GLY_1cter(self):
      '''
      test for a 18-aa system with GLY at cter, 1 flexible range at cter
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','GLY']
      resid = range(1,len(resname))
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [first_last_resid[1]-3]
      rnum = [3]
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]]]
      expected_respsi = [[0, [[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [1e+90, 0.0, 0.0, 0.0, 0.0, 0.]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def test_18aa_GLY_1nter_1middle_1cter(self):
      '''
      test for a 18-aa system with 3 flexible ranges at n-ter, middle and cter respectively
      '''
      #
      resname = ['ALA','ARG','ASP','ASN','CYS','GLU','GLN','HIS','ILE','LEU','LYS','MET','PHE','SER','THR','TRP','TYR','VAL']
      resid = range(1,len(resname))
      first_last_resid = [1,len(resname)]
      rlow = [first_last_resid[0],5,first_last_resid[1]-3]
      rnum = [3,3,3]
      numranges=len(rlow)      
      txtOutput=None
      resphi=[]
      respsi=[]
      energy.protein_initialization(respsi,resphi,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      expected_resphi = [[0, [[1e+90, 0.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]], [1, [[0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]], [2, [[0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0], [0.2, 1.0, 180.0, 0.0, 0.0, 0.0]]]]
      expected_respsi = [[0, [[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0]]], [1, [[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0]]], [2, [[0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [0.6, 1.0, 0.0, 0.0, 0.0, 0.0], [1e+90, 0.0, 0.0, 0.0, 0.0, 0.]]]]
      print 'exp ',expected_resphi,expected_respsi
      print 'res ',resphi,respsi
      self.assertEqual(resphi,expected_resphi)
      self.assertEqual(respsi,expected_respsi)


   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

