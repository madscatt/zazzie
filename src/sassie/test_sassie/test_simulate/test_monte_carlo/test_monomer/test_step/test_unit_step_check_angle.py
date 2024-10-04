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

test for theta=90.0 and dtheta=0.0
test for theta=90.0 and dtheta=3.0
test for theta=90.0 and dtheta=-3.0
test for theta=90.0 and dtheta=90.0
test for theta=90.0 and dtheta=-90.0
test for theta=90.0 and dtheta=180.0
test for theta=90.0 and dtheta=180.0
test for theta=90.0 and dtheta=360.0
test for theta=90.0 and dtheta=-360.0
test for theta=0.0 and dtheta=0.0
test for theta=0.0 and dtheta=90.0
test for theta=0.0 and dtheta=180.0
test for theta=0.0 and dtheta=360.0
test for theta=0.0 and dtheta=-90.0
test for theta=0.0 and dtheta=-180.0
test for theta=0.0 and dtheta=-360.0
test for theta=180.0 and dtheta=0.0
test for theta=180.0 and dtheta=90.0
test for theta=180.0 and dtheta=180.0
test for theta=180.0 and dtheta=360.0
test for theta=180.0 and dtheta=-90.0
test for theta=180.0 and dtheta=-180.0
test for theta=180.0 and dtheta=-360.0
test for theta=-172.6 and dtheta=0.0
test for theta=-172.6 and dtheta=92.8
test for theta=-172.8 and dtheta=178.6
test for theta=-172.8 and dtheta=358.8
test for theta=-172.8 and dtheta=-92.8
test for theta=-172.8 and dtheta=-178.8
test for theta=-172.8 and dtheta=-358.8


'''

import os

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import step



class Test_step_check_angle(MockerTestCase): 

   def setUp(self):
      self.o = step.Setup()

   def test_theta90_d(self):
      '''
      test for theta=90.0 and dtheta=0.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = 0.0
      expected_theta = 90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d3(self):
      '''
      test for theta=90.0 and dtheta=3.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = 3.0
      expected_theta = 93.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d3m(self):
      '''
      test for theta=90.0 and dtheta=-3.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = -3.0
      expected_theta = 87.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d90(self):
      '''
      test for theta=90.0 and dtheta=90.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = 90.0
      expected_theta = -180.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d90m(self):
      '''
      test for theta=90.0 and dtheta=-90.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = -90.0
      expected_theta = 0.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d180(self):
      '''
      test for theta=90.0 and dtheta=180.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = 180.0
      expected_theta = -90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d180m(self):
      '''
      test for theta=90.0 and dtheta=180.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = -180.0
      expected_theta = 270.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d360(self):
      '''
      test for theta=90.0 and dtheta=360.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = 360.0
      expected_theta = 90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta90_d360m(self):
      '''
      test for theta=90.0 and dtheta=-360.0
      '''
      # to-be-modified
      angle_value = 90.0
      trial_value = -360.0
      expected_theta = 90.0
      # may not modify
      with self.assertRaises(Exception):
         result_theta = self.o.check_angle(angle_value,trial_value)

   def test_theta0_d0(self):
      '''
      test for theta=0.0 and dtheta=0.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = 0.0
      expected_theta = 0.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta0_d90(self):
      '''
      test for theta=0.0 and dtheta=90.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = 90.0
      expected_theta = 90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta0_d180(self):
      '''
      test for theta=0.0 and dtheta=180.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = 180.0
      expected_theta = -180.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta0_d360(self):
      '''
      test for theta=0.0 and dtheta=360.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = 360.0
      expected_theta = 0.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta0_d90m(self):
      '''
      test for theta=0.0 and dtheta=-90.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = -90.0
      expected_theta = 90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta0_d180m(self):
      '''
      test for theta=0.0 and dtheta=-180.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = -180.0
      expected_theta = -180.0
      # may not modify
      with self.assertRaises(Exception):
         result_theta = self.o.check_angle(angle_value,trial_value)

   def test_theta0_d360m(self):
      '''
      test for theta=0.0 and dtheta=-360.0
      '''
      # to-be-modified
      angle_value = 0.0
      trial_value = -360.0
      expected_theta = 0.0
      # may not modify
      with self.assertRaises(Exception):
         result_theta = self.o.check_angle(angle_value,trial_value)

   def test_theta180_d0(self):
      '''
      test for theta=180.0 and dtheta=0.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = 0.0
      expected_theta = -180.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d90(self):
      '''
      test for theta=180.0 and dtheta=90.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = 90.0
      expected_theta = -90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d180(self):
      '''
      test for theta=180.0 and dtheta=180.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = 180.0
      expected_theta = 0.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d360(self):
      '''
      test for theta=180.0 and dtheta=360.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = 360.0
      expected_theta = 180.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d90m(self):
      '''
      test for theta=180.0 and dtheta=-90.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = -90.0
      expected_theta = 90.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d180m(self):
      '''
      test for theta=180.0 and dtheta=-180.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = -180.0
      expected_theta = 0.0
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta180_d360m(self):
      '''
      test for theta=180.0 and dtheta=-360.0
      '''
      # to-be-modified
      angle_value = 180.0
      trial_value = -360.0
      expected_theta = 0.0
      # may not modify
      with self.assertRaises(Exception):
         result_theta = self.o.check_angle(angle_value,trial_value)

   def test_theta172m_d0(self):
      '''
      test for theta=-172.6 and dtheta=0.0
      '''
      # to-be-modified
      angle_value = -172.6
      trial_value = 0.0
      expected_theta = -172.6
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta172m_d92(self):
      '''
      test for theta=-172.6 and dtheta=92.8
      '''
      # to-be-modified
      angle_value = -172.6
      trial_value = 92.8
      expected_theta = -79.8
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta172m_d180(self):
      '''
      test for theta=-172.8 and dtheta=178.6
      '''
      # to-be-modified
      angle_value = -172.8
      trial_value = 178.6
      expected_theta = 178.6
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta172m_d358(self):
      '''
      test for theta=-172.8 and dtheta=358.8
      '''
      # to-be-modified
      angle_value = -172.8
      trial_value = 358.8
      expected_theta = 180.0
      # may not modify
      with self.assertRaises(Exception):
         result_theta = self.o.check_angle(angle_value,trial_value)

   def test_theta172_d92m(self):
      '''
      test for theta=-172.8 and dtheta=-92.8
      '''
      # to-be-modified
      angle_value = -172.8
      trial_value = -92.8
      expected_theta = 94.4
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta172_d178m(self):
      '''
      test for theta=-172.8 and dtheta=-178.8
      '''
      # to-be-modified
      angle_value = -172.8
      trial_value = -178.8
      expected_theta = 8.4
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def test_theta172_d358m(self):
      '''
      test for theta=-172.8 and dtheta=-358.8
      '''
      # to-be-modified
      angle_value = -172.8
      trial_value = -358.8
      expected_theta = -171.6
      # may not modify
      result_theta = self.o.check_angle(angle_value,trial_value)
      self.assertAlmostEqual(expected_theta, result_theta)

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

