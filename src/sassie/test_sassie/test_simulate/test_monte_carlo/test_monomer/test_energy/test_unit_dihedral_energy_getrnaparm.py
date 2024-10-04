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


test for regular types
test for nullalpha
test for nullepsilon
test for nulleta
test for nullbeta and nullepsilon
test for nullbeta and nulleta
test for nullepsilon and nulleta
test for nullalpha, nullepsilon and nulleta



'''

import os
import locale

from unittest import main
from mocker import Mocker, MockerTestCase

import sassie.simulate.energy.dihedral_energy as energy





class Test_dihedral_energy_getparamphi(MockerTestCase): 


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


   def test_regular(self):
      '''
      test for regular types
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(0,0,0)
      self.assertEqual(talpha, self.talpha)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon)
      self.assertEqual(teta, self.teta)

   def test_nullalpha(self):
      '''
      test for nullalpha
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(1,0,0)
      self.assertEqual(talpha, self.talpha_ter)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon)
      self.assertEqual(teta, self.teta)
   
   def test_nullepsilon(self):
      '''
      test for nullepsilon
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(0,1,0)
      self.assertEqual(talpha, self.talpha)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon_ter)
      self.assertEqual(teta, self.teta)
   
   def test_nulleta(self):
      '''
      test for nulleta
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(0,0,1)
      self.assertEqual(talpha, self.talpha)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon)
      self.assertEqual(teta, self.teta_ter)

   def test_nullalpha_nullepsilon(self):
      '''
      test for nullbeta and nullepsilon
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(1,1,0)
      self.assertEqual(talpha, self.talpha_ter)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon_ter)
      self.assertEqual(teta, self.teta)

   def test_nullalpha_nulleta(self):
      '''
      test for nullbeta and nulleta
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(1,0,1)
      self.assertEqual(talpha, self.talpha_ter)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon)
      self.assertEqual(teta, self.teta_ter)

   def test_nullepsilon_nulleta(self):
      '''
      test for nullepsilon and nulleta
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(0,1,1)
      self.assertEqual(talpha, self.talpha)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon_ter)
      self.assertEqual(teta, self.teta_ter)

   def test_nullalpha_nullepsilon_nulleta(self):
      '''
      test for nullalpha, nullepsilon and nulleta
      '''
      #
      talpha, tbeta, tgama, tdelta, tepsilon, teta = energy.getrnaparm(1,1,1)
      self.assertEqual(talpha, self.talpha_ter)
      self.assertEqual(tbeta, self.tbeta)
      self.assertEqual(tgama, self.tgama)
      self.assertEqual(tdelta, self.tdelta)
      self.assertEqual(tepsilon, self.tepsilon_ter)
      self.assertEqual(teta, self.teta_ter)


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

