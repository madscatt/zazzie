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

test unpack


'''

import os

from unittest import main
from mocker import Mocker, MockerTestCase

from sassie.simulate.monte_carlo.monomer import dihedral_monte_carlo



class Test_dihedral_monte_carlo_unpack_variables(MockerTestCase): 


   def setUp(self):
      pass


   def test_all(self):
      '''
      test unpack
      '''
      #
      svariables={}

      svariables['runname']		= ['run']
      svariables['dcdfile']		= ['a.dcd']
      svariables['path']		= ['./']
      svariables['pdbfile']		= ['a.pdb']
      svariables['trials']		= [1000]
      svariables['goback']		= [5]
      svariables['temp']		= [300.0]
      svariables['moltype']		= ['protein']
      svariables['numranges']		= [2]
      svariables['dtheta']		= [[10.0,10.0]]
      svariables['reslow']		= [[10,20]]
      svariables['numcont']		= [[5,3]]
      svariables['lowres1']		= [1]
      svariables['highres1']		= [30]
      svariables['basis']		= ['CA']
      svariables['cutoff']		= [0.0]
      svariables['lowrg']		= [30.0]
      svariables['highrg']		= [100.0]
      svariables['zflag']		= [1]
      svariables['zcutoff']		= [0.0]
      svariables['cflag']		= [0]
      svariables['confile']		= ['a.con']
      svariables['nonbondflag']	= [0]
      svariables['nonbondscale']	= [1.0]
      svariables['psffilepath']	= ['./']
      svariables['psffilename']	= ['a.psf']
      svariables['parmfilepath']	= ['./']
      svariables['parmfilename']	= ['para.par']
      svariables['plotflag']		= [0]
      svariables['seed']		= [[1,123]] # set this to '1,123' if you want to set the seed

      result = dihedral_monte_carlo.unpack_variables(svariables)
      expected = 'run','a.dcd','./','a.pdb',1000,5,300.0,'protein',2,[10.0,10.0],[10,20],[5,3],1,30,'CA',0.0,30.0,100.0,1,0.0,0,'a.con',0,1.0,'./','a.psf','./','para.par',0,[1,123]
      self.assertEqual(result, expected)


   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

