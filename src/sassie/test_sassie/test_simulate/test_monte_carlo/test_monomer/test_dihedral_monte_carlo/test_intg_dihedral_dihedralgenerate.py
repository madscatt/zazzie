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

import multiprocessing

from sassie.core_testing.util import env
import sassie.interface.input_filter as input_filter

'''
intg test:

contract

test a regular file; get the rotation indices for residue-1-2,10-12,21-23; the output files will be comppared against files under a provided expected output directory



'''

import os, shutil, filecmp, glob
import numpy

from unittest import main,skipIf
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol
from sassie.core_testing.util import FileCmp
from sassie.simulate.monte_carlo.monomer import dihedral_monte_carlo



class Test_dihedral_monte_carlo_get_rotation_indices(MockerTestCase): 

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
      self.DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../../../data/simulate/monte_carlo/monomer/dihedral_monte_carlo/'
      self.o = sasmol.SasMol(0)

   def verify_me(self,expected_path,result_path):
      result_files = glob.glob(result_path+'*')
      result_files = [os.path.split(result_files[i])[1] for i in range(len(result_files))]
      expected_files = glob.glob(expected_path+'*')
      expected_files = [os.path.split(expected_files[i])[1] for i in range(len(expected_files))]
      print '\nresult_files:   ',result_path, '\n', result_files, '\n\nexepected_files:',expected_path, '\n',expected_files
      flag = True
      for ifile in expected_files:
         if ifile in result_files:
            print '\ndiff ',expected_path+ifile, result_path+ifile,
            if ifile[-4:]=='.dcd':
               flag = (flag and FileCmp.cmp_dcd_skip_date(expected_path+ifile, result_path+ifile))
            else:
               flag = (flag and filecmp.cmp(expected_path+ifile, result_path+ifile))
            print '\n...to be...',flag
            if flag==False:
               return False
         else:
            return False
      return flag



   def test_1CRN_3flexibleregions(self):
      '''
      test a regular file; get the rotation indices for residue-1-2,10-12,21-23; the output files will be comppared against files under a provided expected output directory
      '''
      #to-be-modified
      svariables={}
      svariables['pdbfile']		= ('1CRN.pdb','string')
      svariables['trials']		= ('2','int')
      svariables['goback']		= ('10','int')
      svariables['temp']		= ('300.0','float')
      svariables['moltype']		= ('protein','string')
      svariables['numranges']		= ('2','int')
      svariables['dtheta']		= ('30.0,30.0','float_array')
      svariables['reslow']		= ('10, 20','int_array')
      svariables['numcont']		= ('5, 5','int_array')
      svariables['lowres1']		= ('1','int')
      svariables['highres1']		= ('10','int')
      svariables['cutoff']		= ('3.0','float')
      svariables['lowrg']		= ('0.0','float')
      svariables['highrg']		= ('400.0','float')
      svariables['zflag']		= ('0','int')
      svariables['zcutoff']		= ('0.0','float')
      #may not modify
      svariables['runname']		= ('run_0','string')
      svariables['dcdfile']		= ('run_0.dcd','string')
      svariables['path']		= (self.DataPath,'string')
      svariables['basis']		= ('CA','string')
      svariables['cflag']		= ('0','int')
      svariables['confile']		= ('constraints.txt','string')
      svariables['nonbondflag']	= ('0','int')
      svariables['nonbondscale']	= ('1.0','float')
      svariables['psffilepath']	= (self.DataPath,'string')
      svariables['psffilename']	= ('refgag.psf','string')
      svariables['parmfilepath']	= ('/usr/local/bin/sassie/simulate/namd/','string')
      svariables['parmfilename']	= ('par_all27_prot_na.inp','string')
      svariables['plotflag']		= ('0', 'int')
      svariables['seed']		= ('1,123', 'int_array') # set this to '1,123' if you want to set the seed
      error,self.variables=input_filter.type_check_and_convert(svariables)
      #run
      if os.path.isdir(self.variables['runname'][0]):
         shutil.rmtree(self.variables['runname'][0])
      self.txtQueue=multiprocessing.JoinableQueue()
      dihedral_monte_carlo.dihedralgenerate(self.variables, self.txtQueue)
      self.verify_me(os.path.join(self.DataPath,'dihedral_monte_carlogenerate_expected','1CRN/'), os.path.join(self.variables['runname'][0],'generate/'))
      shutil.rmtree(self.variables['runname'][0])

   """
   def test_rna_3flexibleregions(self):
      '''
      test a regular file; get the rotation indices for residue-1-2,10-12,21-23; the output files will be comppared against files under a provided expected output directory
      '''
      #to-be-modified
      svariables={}
      svariables['pdbfile']		= ('rna_single.pdb','string')
      svariables['trials']		= ('2','int')
      svariables['goback']		= ('10','int')
      svariables['temp']		= ('300.0','float')
      svariables['moltype']		= ('rna','string')
      svariables['numranges']		= ('2','int')
      svariables['dtheta']		= ('0.0,0.0','float_array')
      svariables['reslow']		= ('5, 10','int_array')
      svariables['numcont']		= ('3, 3','int_array')
      svariables['lowres1']		= ('1','int')
      svariables['highres1']		= ('3','int')
      svariables['cutoff']		= ('3.0','float')
      svariables['lowrg']		= ('0.0','float')
      svariables['highrg']		= ('400.0','float')
      svariables['zflag']		= ('0','int')
      svariables['zcutoff']		= ('0.0','float')
      #may not modify
      svariables['runname']		= ('run_0','string')
      svariables['dcdfile']		= ('run_0.dcd','string')
      svariables['path']		= (self.DataPath,'string')
      svariables['basis']		= ('CA','string')
      svariables['cflag']		= ('0','int')
      svariables['confile']		= ('constraints.txt','string')
      svariables['nonbondflag']	= ('0','int')
      svariables['nonbondscale']	= ('1.0','float')
      svariables['psffilepath']	= (self.DataPath,'string')
      svariables['psffilename']	= ('refgag.psf','string')
      svariables['parmfilepath']	= ('/usr/local/bin/sassie/simulate/namd/','string')
      svariables['parmfilename']	= ('par_all27_prot_na.inp','string')
      svariables['plotflag']		= ('0', 'int')
      svariables['seed']		= ('1,123', 'int_array') # set this to '1,123' if you want to set the seed
      error,self.variables=input_filter.type_check_and_convert(svariables)
      #run
      if os.path.isdir(self.variables['runname'][0]):
         shutil.rmtree(self.variables['runname'][0])
      self.txtQueue=multiprocessing.JoinableQueue()
      dihedral_monte_carlo.dihedralgenerate(self.variables, self.txtQueue)
      self.verify_me(os.path.join(self.DataPath,'dihedral_monte_carlogenerate_expected','1CRN/'), os.path.join(self.variables['runname'][0],'generate/'))
      shutil.rmtree(self.variables['runname'][0])
   """

   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_gag(self):
      '''
      test a regular file; get the rotation indices for residue-1-2,10-12,21-23
      '''
      #to-be-modified
      svariables={}
      svariables['pdbfile']		= ('min3.pdb','string')
      svariables['trials']		= ('2','int')
      svariables['goback']		= ('10','int')
      svariables['temp']		= ('300.0','float')
      svariables['moltype']		= ('protein','string')
      svariables['numranges']		= ('5','int')
      svariables['dtheta']		= ('30.0,30.0,30.0,30.0,30.0','float_array')
      svariables['reslow']		= ('123, 278, 354, 378, 408','int_array')
      svariables['numcont']		= ('21, 5, 24, 11, 4','int_array')
      svariables['lowres1']		= ('284','int')
      svariables['highres1']		= ('350','int')
      svariables['cutoff']		= ('3.0','float')
      svariables['lowrg']		= ('0.0','float')
      svariables['highrg']		= ('400.0','float')
      svariables['zflag']		= ('0','int')
      svariables['zcutoff']		= ('0.0','float')
      #may not modify
      svariables['runname']		= ('run_0','string')
      svariables['dcdfile']		= ('run_0.dcd','string')
      svariables['path']		= (self.DataPath,'string')
      svariables['basis']		= ('CA','string')
      svariables['cflag']		= ('0','int')
      svariables['confile']		= ('constraints.txt','string')
      svariables['nonbondflag']	= ('0','int')
      svariables['nonbondscale']	= ('1.0','float')
      svariables['psffilepath']	= (self.DataPath,'string')
      svariables['psffilename']	= ('refgag.psf','string')
      svariables['parmfilepath']	= ('/usr/local/bin/sassie/simulate/namd/','string')
      svariables['parmfilename']	= ('par_all27_prot_na.inp','string')
      svariables['plotflag']		= ('0', 'int')
      svariables['seed']		= ('1,123', 'int_array') # set this to '1,123' if you want to set the seed
      error,self.variables=input_filter.type_check_and_convert(svariables)
      #run
      if os.path.isdir(self.variables['runname'][0]):
         shutil.rmtree(self.variables['runname'][0])
      self.txtQueue=multiprocessing.JoinableQueue()
      dihedral_monte_carlo.dihedralgenerate(self.variables, self.txtQueue)
      self.verify_me(os.path.join(self.DataPath,'dihedral_monte_carlogenerate_expected','gag/'), os.path.join(self.variables['runname'][0],'generate/'))
      shutil.rmtree(self.variables['runname'][0])

   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

