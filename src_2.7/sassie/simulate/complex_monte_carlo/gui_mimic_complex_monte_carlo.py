'''
Driver method to run the complex_monte_carlo module
'''

import sys
import string
import os
import shutil
import time

import sassie.util.sasconfig as sasconfig
import sassie.interface.input_filter as input_filter
import sassie.simulate.complex_monte_carlo.complex_monte_carlo as complex_monte_carlo
#import complex_monte_carlo as complex_monte_carlo
import sassie.interface.complex_monte_carlo.complex_filter as complex_filter
#import complex_filter as complex_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.dcdfile = 'run_0.dcd'
    self.path = './'
    self.pdbfile = 'pai_vn_start.pdb'
#    self.pdbfile = 'rna_protein_complex.pdb'
    self.trials = '50'
    self.goback = '10'
    self.nsegments = '2'
#    self.nsegments = '7'   #rna-prot
    self.npsegments = '1'
#    self.npsegments = '7'  #rna-prot
#    self.flpsegname = 'HFQ1,HFQ2,HFQ3,HFQ4,HFQ5,HFQ6,RNA1'   #rna-prot
#    self.flpsegname = 'RNA1' 
    self.flpsegname = 'VN1'   
#    self.segbasis = 'CA,CA'
#    self.segbasis = 'all'
    self.segbasis = 'heavy'
#    self.segbasis = 'backbone'
    self.seglow = '1'
#    self.seglow = '20'
#    self.seglow = '20,20,20,20,20,20,20'     #rna-prot
    self.seghigh = '39'
#    self.seghigh = '30'
#    self.seghigh = '30,30,30,30,30,30,30'    #rna-prot
    self.temp = '300.0'
    self.lowrg = '0.0'
#    self.lowrg = '46.86'
#    self.highrg = '46.89'
    self.highrg = '50.0'
#    self.highrg = '400.0'   #rna-prot
    self.zflag = '0'
    self.zcutoff = '-19.0'
#    self.zcutoff = '-60.0'
    self.cflag = '0'
    self.confile = 'pai_vn_constraints.txt'
#    self.confile = 'rna_protein_constraints2.txt'
    self.directedmc = '0'
    self.psffilepath='./'   #not currently used
    self.psffilename = 'refgag.psf'     #not currently used
    self.parmfilepath = sasconfig.__bin_path__ + 'toppar'  #not currently used
    self.parmfilename = 'par_all27_prot_na.inp' #not currently used
    self.plotflag = '1'
    self.seed = '0,123'
#    self.seed = '1,123' #to use a specific random sequence each time -- use to create files for tests
#    self.seed = '1,321' #use for rna-prot
    
    self.psegvariables= [['1', '30', '40', '89', 'protein']]
#    self.psegvariables = [['1', '30', '128', '2', 'rna']]
#    self.psegvariables= [['2', '30,30', '40,75', '32,54', 'protein']]
#    self.psegvariables = [['2', '30,30', '128,276', '1,1', 'rna']]
#    self.psegvariables= [['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '67', '33', 'protein'],['1', '30', '128', '1', 'rna']]    

    self.testflag = False

    ### END USER INPUT ###
    ### END USER INPUT ###
    ### END USER INPUT ###
    

def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    other_data_path = paths['other_data_path']
    module_data_path = paths['module_data_path']

    self.runname = 'run_0'
    self.dcdfile = 'run_0.dcd'
    self.path = ''
    self.pdbfile = os.path.join(pdb_data_path,'pai_vn_start.pdb')
    self.trials = '20'
    self.goback = '10'
    self.nsegments = '2'
    self.npsegments = '1'
    self.flpsegname = 'VN1'
    self.segbasis = 'heavy'
    self.seglow = '1'
    self.seghigh = '39'
    self.temp = '300.0'
    self.lowrg = '0.0'
    self.highrg = '50.0'
    self.zflag = '0'
    self.zcutoff = '0.0'
    self.cflag = '0'
    self.confile = os.path.join(other_data_path,'complex_mc_constraints.txt')
    self.directedmc = '0'
    self.psffilepath = ''   #not currently used
    self.psffilename = os.path.join(other_data_path,'gag_start.psf')    #not currently used
    self.parmfilepath = sasconfig.__bin_path__ + 'toppar'  #not currently used
    self.parmfilename = 'par_all27_prot_na.inp'         #not currently used
    self.plotflag = '0'
    self.seed = '1,123' #the '1' is needed for testing so that a specific random sequence is used each time

    self.psegvariables= [['1', '30', '40', '89', 'protein']]
    
    self.precision = 3
    self.testflag = True
    
def run_module(self, **kwargs):

    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables={}

    svariables['cflag'] = (self.cflag,'int')
    svariables['confile'] = (self.confile,'string')
    svariables['dcdfile'] = (self.dcdfile,'string')
    svariables['directedmc'] = (self.directedmc,'float')
    svariables['flpsegname'] = (self.flpsegname, 'string')
    svariables['goback'] = (self.goback,'int')
    svariables['highrg'] = (self.highrg,'float')
    svariables['lowrg'] = (self.lowrg,'float')
    svariables['npsegments'] = (self.npsegments,'int')
    svariables['nsegments'] = (self.nsegments,'int')
    svariables['parmfilename'] = (self.parmfilename,'string')
    svariables['path'] = (self.path,'string')
    svariables['pdbfile'] = (self.pdbfile,'string')
    svariables['plotflag'] = (self.plotflag,'int')
    svariables['psffilename'] = (self.psffilename,'string')
    svariables['runname'] = (self.runname,'string')
    svariables['seed'] = (self.seed,'int_array')
    svariables['segbasis'] = (self.segbasis,'string')
    svariables['seghigh'] = (self.seghigh,'int_array')
    svariables['seglow'] = (self.seglow,'int_array')
    svariables['temp'] = (self.temp,'float')
    svariables['trials'] = (self.trials,'int')
    svariables['zcutoff'] = (self.zcutoff,'float')
    svariables['zflag'] = (self.zflag, 'int')


    error,self.variables=input_filter.type_check_and_convert(svariables)

#    print 'variables: ', self.variables
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error
    
    try:
        if kwargs['file_check']:
            error = complex_filter.check_complex(self.variables, self.psegvariables)
    except:
            error = complex_filter.check_complex(self.variables, self.psegvariables, no_file_check="true")

    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass
          
    runname=self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue=multiprocessing.JoinableQueue()
    this_complex_monte_carlo = complex_monte_carlo.complex_monte_carlo()
    this_complex_monte_carlo.main(self.variables, self.psegvariables,txtQueue)



class gui_mimic_complex_monte_carlo():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'complex_monte_carlo'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 1.0.
    paths = None


# We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
# That way, the test files will always be available to the user.
    if test:
        pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
        dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
        other_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'other_common') + os.path.sep
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'simulate', 'monomer_monte_carlo') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_complex_monte_carlo(test, paths)
    print "time used: ", time.time() - start






    
