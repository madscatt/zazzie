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
contract

test for regular alpha angle rotating a small angle at 300K
test for regular alpha angle rotating a large angle at 300K
test for regular alpha angle rotating a small angle at 30K
test for regular alpha angle rotating a large angle at 30K
test for regular alpha angle rotating a small angle at 3000K
test for regular alpha angle rotating a large angle at 3000K
test for beta angle rotating a small angle at 300K
test for beta angle rotating a large angle at 300K
test for beta angle rotating a small angle at 30K
test for beta angle rotating a large angle at 30K
test for beta angle rotating a small angle at 3000K
test for beta angle rotating a large angle at 3000K
test for delta angle rotating a small angle at 300K
test for delta angle rotating a large angle at 300K
test for delta angle rotating a small angle at 30K
test for delta angle rotating a large angle at 30K
test for delta angle rotating a small angle at 3000K
test for delta angle rotating a large angle at 3000K
test for epsilon angle rotating a small angle at 300K
test for epsilon angle rotating a large angle at 300K
test for epsilon angle rotating a small angle at 30K
test for epsilon angle rotating a large angle at 30K
test for epsilon angle rotating a small angle at 3000K
test for epsilon angle rotating a large angle at 3000K
test for regular teta angle rotating a small angle at 300K
test for regular teta angle rotating a large angle at 300K
test for regular teta angle rotating a small angle at 30K
test for regular teta angle rotating a large angle at 30K
test for regular teta angle rotating a small angle at 3000K
test for regular teta angle rotating a large angle at 3000K


'''

import os
import locale
import numpy
from numpy.random import RandomState

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy





class Test_dihedral_energy_calc(MockerTestCase): 

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
      talpha = [1.20,1,180.0,0.1,2,180.0,0.1,3,180.0,6.0,0,0.0]
      talpha_ter = [1E90,0,180.0,1E90,2,180.0,1E90,3,180.0,1E90,0,0.0]
      tbeta = [0.2, 1, 120.0]
      tdelta = [0.2,4,180.0,0.8,3,180.0,0.4,2,0.0,2.5,1,180.0]
      tepsilon = [0.2,4,0.0,0.8,3,180.0]
      teta = [0.6,5,0.0,0.2,4,0.0,0.0,3,180.0,0.4,2,0.0,1.9,1,180.0]
      teta_ter = [1E90,0,0.0,1E90,4,0.0,1E90,3,180.0,1E90,2,0.0,1E90,1,180.0]
      self.para = [talpha, tbeta, tdelta, tepsilon, teta]
      self.para_alphater = [talpha_ter, tbeta, tdelta, tepsilon, teta]
      self.para_tetater = [talpha, tbeta, tdelta, tepsilon, teta_ter]
      self.para_bothter = [talpha_ter, tbeta, tdelta, tepsilon, teta_ter]

   def test_regular_alpha_small_angle_300K(self):
      '''
      test for regular alpha angle rotating a small angle at 300K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0, 12.311, 12.330]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_regular_alpha_small_angle_300K(self):
      '''
      test for regular alpha angle rotating a large angle at 300K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_regular_alpha_small_angle_30K(self):
      '''
      test for regular alpha angle rotating a small angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0, 12.311, 12.330]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_regular_alpha_small_angle_30K(self):
      '''
      test for regular alpha angle rotating a large angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)


   def test_regular_alpha_small_angle_3000K(self):
      '''
      test for regular alpha angle rotating a small angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0, 12.311, 12.330]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_regular_alpha_small_angle_3000K(self):
      '''
      test for regular alpha angle rotating a large angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,12.311,13.5]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_small_angle_300K(self):
      '''
      test for beta angle rotating a small angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.2,0.203]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_large_angle_300K(self):
      '''
      test for beta angle rotating a large angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.2,0.373]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_small_angle_30K(self):
      '''
      test for beta angle rotating a small angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.2,0.203]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_large_angle_30K(self):
      '''
      test for beta angle rotating a large angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_small_angle_3000K(self):
      '''
      test for beta angle rotating a small angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.2,0.203]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_beta_large_angle_3000K(self):
      '''
      test for beta angle rotating a large angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.2,0.373]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_small_angle_300K(self):
      '''
      test for delta angle rotating a small angle at 300K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,2.035,2.099]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_large_angle_300K(self):
      '''
      test for delta angle rotating a large angle at 300K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_small_angle_30K(self):
      '''
      test for delta angle rotating a small angle at 30K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_large_angle_30K(self):
      '''
      test for delta angle rotating a large angle at 30K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_small_angle_3000K(self):
      '''
      test for delta angle rotating a small angle at 3000K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,2.035,2.099]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_delta_large_angle_3000K(self):
      '''
      test for delta angle rotating a large angle at 3000K
      '''
      #
      angle_index = 2
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,2.035,3.3]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)


   def test_epsilon_small_angle_300K(self):
      '''
      test for epsilon angle rotating a small angle at 300K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.900,0.930]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_epsilon_large_angle_300K(self):
      '''
      test for epsilon angle rotating a large angle at 300K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.900,1.200]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_epsilon_small_angle_30K(self):
      '''
      test for epsilon angle rotating a small angle at 30K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.900,0.930]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_epsilon_large_angle_30K(self):
      '''
      test for epsilon angle rotating a large angle at 30K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_epsilon_small_angle_3000K(self):
      '''
      test for epsilon angle rotating a small angle at 3000K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.900,0.930]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_epsilon_large_angle_3000K(self):
      '''
      test for epsilon angle rotating a large angle at 3000K
      '''
      #
      angle_index = 3
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,0.900,1.200]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_small_angle_300K(self):
      '''
      test for regular teta angle rotating a small angle at 300K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,1.035,1.004]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_large_angle_300K(self):
      '''
      test for regular teta angle rotating a large angle at 300K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_small_angle_30K(self):
      '''
      test for regular teta angle rotating a small angle at 30K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,1.035,1.004]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_large_angle_30K(self):
      '''
      test for regular teta angle rotating a large angle at 30K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [1,1E10,1E10]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_small_angle_3000K(self):
      '''
      test for regular teta angle rotating a small angle at 3000K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,1.035,1.004]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      #print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def test_reg_teta_large_angle_3000K(self):
      '''
      test for regular teta angle rotating a large angle at 3000K
      '''
      #
      angle_index = 4
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [0,1.035,2.900]
      search, vdi, vdf = energy.calc(angle_index,itheta,theta,self.para,beta,nonbondflag,seed_object)
      result=[search,vdi,vdf]
      print 'result ',result
      self.assert_list_almost_equal(expected,result)

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

