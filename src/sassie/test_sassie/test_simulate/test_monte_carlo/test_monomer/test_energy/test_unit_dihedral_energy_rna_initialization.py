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

test for a 8-base system, 1 flexible base in middle
test for a 8-base system, 1 flexible base at 5'-ter
test for a 8-base system, 1 flexible base at 3'-ter
test for a 8-base system, 4 flexible bases in middle
test for a 8-base system,  4 flexible bases range at 5'-ter
test for a 8-base system,  4 flexible bases range at 3'-ter
test for a 8-base system, 2 flexible bases at 3-ter, 2 flexible bases in middle, and 2 flexible basesat 5-ter



'''

import os
import numpy

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy





class Test_dihedral_energy_protein_initialization(MockerTestCase): 

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
      self.talpha = [1.20,1,180.0,0.1,2,180.0,0.1,3,180.0,0.0,6,0.0]
      self.talpha_ter = [1E90,0,180.0,1E90,2,180.0,1E90,3,180.0,1E90,0,0.0]
      self.tbeta = [0.2, 1, 120.0]
      self.tgama = [0.2,4,180.0,0.8,3,180.0,0.4,2,0.0,2.5,1,180.0]
      self.tdelta = [0.2,4,0.0,0.8,3,180.0]
      self.tepsilon = [0.6,5,0.0,0.2,4,0.0,0.0,3,180.0,0.4,2,0.0,1.9,1,180.0]
      self.tepsilon_ter = [1E90,5,0.0,1E90,4,0.0,1E90,3,180.0,1E90,2,0.0,1E90,1,180.0]
      self.teta = [1.20,1,180.0,0.10,2,180.0,0.1,3,180.0,0.0,6,0.0]
      self.teta_ter = [1E90,1,180.0,1E90,2,180.0,1E90,3,180.0,1E90,6,0.0]


   def test_8b_1middle(self):
      '''
      test for a 8-base system, 1 flexible base in middle
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [3]
      rnum = [0]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0]+1)]]])


   def test_8b_1at5ter(self):
      '''
      test for a 8-base system, 1 flexible base at 5'-ter
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [1]
      rnum = [0]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha_ter]+[self.talpha for i in range(rnum[0])]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0]+1)]]])


   def test_8b_1at3ter(self):
      '''
      test for a 8-base system, 1 flexible base at 3'-ter
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [8]
      rnum = [0]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0])]+[self.tepsilon_ter]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0])]+[self.teta_ter]]])

   def test_8b_4middle(self):
      '''
      test for a 8-base system, 4 flexible bases in middle
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [3]
      rnum = [3]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0]+1)]]])


   def test_8b_4at5ter(self):
      '''
      test for a 8-base system,  4 flexible bases range at 5'-ter
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [1]
      rnum = [3]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha_ter]+[self.talpha for i in range(rnum[0])]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0]+1)]]])


   def test_8b_4at3ter(self):
      '''
      test for a 8-base system,  4 flexible bases range at 3'-ter
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      numranges=1
      rlow = [5]
      rnum = [3]
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta for i in range(rnum[0]+1)]]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon for i in range(rnum[0])]+[self.tepsilon_ter]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta for i in range(rnum[0])]+[self.teta_ter]]])

   def test_8b_2at5ter_2middle_2at3ter(self):
      '''
      test for a 8-base system, 2 flexible bases at 3-ter, 2 flexible bases in middle, and 2 flexible basesat 5-ter
      '''
      #
      resname = ['RNUS','RNUA','RUUG','RNUC','RNUS','RNUA','RUUG','RNUC']
      resid = range(1,len(resname)+1)
      first_last_resid = [1,len(resname)]
      rlow = [1,3,7]
      rnum = [1,1,1]
      numranges=len(rlow)
      txtOutput=None
      resalpha, resbeta, resgama, resdelta, resepsilon, reseta =[],[],[],[],[],[]
      energy.rna_initialization(resalpha, resbeta, resgama, resdelta, resepsilon, reseta,resid,resname,numranges,rlow,rnum,first_last_resid,txtOutput)
      print 'res \n',resalpha, '\n',resbeta, '\n',resdelta,'\n', resepsilon,'\n', reseta
      self.assert_list_almost_equal(resalpha,[[0,[self.talpha_ter,self.talpha]],[1,[self.talpha]*2],[2,[self.talpha]*2]])
      self.assert_list_almost_equal(resbeta,[[0,[self.tbeta]*2],[1,[self.tbeta]*2],[2,[self.tbeta]*2]])
      self.assert_list_almost_equal(resdelta,[[0,[self.tdelta]*2],[1,[self.tdelta]*2],[2,[self.tdelta]*2]])
      self.assert_list_almost_equal(resepsilon,[[0,[self.tepsilon]*2],[1,[self.tepsilon]*2],[2,[self.tepsilon,self.tepsilon_ter]]])
      self.assert_list_almost_equal(reseta,[[0,[self.teta]*2],[1,[self.teta]*2],[2,[self.teta,self.teta_ter]]])

   def tearDown(self):
     pass

if __name__ == '__main__': 
   main() 

