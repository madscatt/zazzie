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

all the following tests are performed for all possible combinations of dihedral types
test for all regular phi angle rotating a small angle at 300K
test for all regular phi angle rotating a large angle at 30K
test for all regular phi angle rotating a small angle at 30K
test for all regular phi angle rotating a large angle at 30K
test for all regular phi angle rotating a small angle at 3000K
test for all regular phi angle rotating a large angle at 3000K
test for a nter phi angle rotating a small angle at 300K
test for a nter phi angle rotating a large angle at 300K
test for a nter phi angle rotating a small angle at 30K
test for a nter phi angle rotating a large angle at 30K
test for a nter phi angle rotating a small angle at 3000K
test for a nter phi angle rotating a large angle at 3000K
test for all regular psi angle rotating a small angle at 300K
test for all regular psi angle rotating a large angle at 300K
test for all regular psi angle rotating a small angle at 30K
test for all regular psi angle rotating a large angle at 30K
test for all regular psi angle rotating a small angle at 3000K
test for all regular psi angle rotating a large angle at 3000K
test for a cter psi angle rotating a small angle at 300K
test for a cter psi angle rotating a large angle at 300K
test for a cter psi angle rotating a small angle at 30K
test for a cter psi angle rotating a large angle at 30K
test for a cter psi angle rotating a small angle at 3000K
test for a cter psi angle rotating a large angle at 3000K


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
      self.l_phi = [[0.2,1.0,180.0,0.0,0.0,0.0],[0.8,3.0,0.0,0.0,0.0,0.0],[1E90,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0]]
      self.l_psi = [[0.6,1.0,0.0,0.0,0.0,0.0],[0.3,1.0,0.0,-0.3,4.0,0.0],[0.4,1.0,0.0,0.0,0.0,0.0],[1E90,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0]]
      self.l_phipsi = []
      for phi in self.l_phi:
         for psi in self.l_psi:
            self.l_phipsi.append([phi,psi])
      self.l_phipsi2 = []
      for psi in self.l_psi:
         for phi in self.l_phi:
            self.l_phipsi2.append([phi,psi])

   def test_regular_phi_small_angle_300K(self):
      '''
      test for all regular phi angle rotating a small angle at 300K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_phi_large_angle_300K(self):
      '''
      test for all regular phi angle rotating a large angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 1.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0267, 3.046e-05], [0, 0.0267, 3.046e-05], [0, 0.0267, 3.046e-05], [0, 0.0267, 3.046e-05], [0, 0.0267, 3.046e-05], [1,1E10,1E10], [1,1E10,1E10], [0,0.8,1.599], [1,1E10,1E10], [0,0.8,1.599]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_phi_small_angle_30K(self):
      '''
      test for all regular phi angle rotating a small angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_phi_large_angle_30K(self):
      '''
      test for all regular phi angle rotating a large angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 1.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [1,1E10,1E10], [1,1E10,1E10], [1,1E10,1E10], [1,1E10,1E10], [1,1E10,1E10]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_phi_small_angle_3000K(self):
      '''
      test for all regular phi angle rotating a small angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0, 0.0268, 0.0286],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758],[0,0.8,0.758]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_phi_large_angle_3000K(self):
      '''
      test for all regular phi angle rotating a large angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 1.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0, 0.0268, 3.046e-05], [0,0.8,1.599],[0,0.8,1.599], [0,0.8,1.599], [0,0.8,1.599], [0,0.8,1.599]]
      for i in range(len(self.l_phipsi[0:10])):
         parm = self.l_phipsi[0:10][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_small_angle_300K(self):
      '''
      test for a nter phi angle rotating a small angle at 300K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_large_angle_300K(self):
      '''
      test for a nter phi angle rotating a large angle at 300K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 0.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_small_angle_30K(self):
      '''
      test for a nter phi angle rotating a small angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_large_angle_30K(self):
      '''
      test for a nter phi angle rotating a large angle at 30K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 9.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_small_angle_3000K(self):
      '''
      test for a nter phi angle rotating a small angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_nter_phi_large_angle_3000K(self):
      '''
      test for a nter phi angle rotating a large angle at 3000K
      '''
      #
      angle_index = 0
      itheta = 30.0
      theta = 9.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 2E90, 2E90] for i in range(5)]
      for i in range(len(self.l_phipsi[10:15])):
         parm = self.l_phipsi[10:15][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_psi_small_angle_300K(self):
      '''
      test for all regular psi angle rotating a small angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_psi_large_angle_300K(self):
      '''
      test for all regular psi angle rotating a large angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_psi_small_angle_30K(self):
      '''
      test for all regular psi angle rotating a small angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_psi_large_angle_30K(self):
      '''
      test for all regular psi angle rotating a large angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)


   def test_regular_psi_small_angle_3000K(self):
      '''
      test for all regular psi angle rotating a small angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0, 1.120,1.114],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.410,0.425],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743],[0,0.746,0.743]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_regular_psi_large_angle_3000K(self):
      '''
      test for all regular psi angle rotating a large angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0, 1.120,0.6],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.410,-0.3],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4],[0,0.746,0.4]]
      for i in range(len(self.l_phipsi2[0:12])):
         parm = self.l_phipsi2[0:12][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_small_angle_300K(self):
      '''
      test for a cter psi angle rotating a small angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_large_angle_300K(self):
      '''
      test for a cter psi angle rotating a large angle at 300K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*300)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_small_angle_30K(self):
      '''
      test for a cter psi angle rotating a small angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_large_angle_30K(self):
      '''
      test for a cter psi angle rotating a large angle at 30K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*30)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_small_angle_3000K(self):
      '''
      test for a cter psi angle rotating a small angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 31.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)

   def test_cter_psi_large_angle_3000K(self):
      '''
      test for a cter psi angle rotating a large angle at 3000K
      '''
      #
      angle_index = 1
      itheta = 30.0
      theta = 90.0
      beta = 1.0/(1.380658E-23*3000)
      nonbondflag = 0
      seed_object = RandomState(1)
      expected = [[0,2E90,2E90] for i in range(4)]
      for i in range(len(self.l_phipsi2[12:16])):
         parm = self.l_phipsi2[12:16][i]
         search, vdi, vdf = energy.calc(angle_index,itheta,theta,parm,beta,nonbondflag,seed_object)
         result=[search,vdi,vdf]
         #print 'result ',result
         #print 'expected ',expected[i]
         self.assert_list_almost_equal(expected[i],result)


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

