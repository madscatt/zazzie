'''
Driver method to run the monomer_monte_carlo module
'''

import sys
import string
import os
import shutil
import time

import sassie.util.sasconfig as sasconfig
import sassie.interface.input_filter as input_filter
import sassie.simulate.monomer_monte_carlo.monomer_monte_carlo as monomer_monte_carlo
import sassie.interface.monomer_monte_carlo.monomer_monte_carlo_filter as monomer_monte_carlo_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.dcdfile = 'run_0.dcd'
    self.moltype = 'protein'
#    self.moltype = 'rna'
    self.path = './'
    self.pdbfile = 'min3.pdb'
    #self.pdbfile = 'gag_start.pdb'
#    self.pdbfile = 'hiv1_gag_ma.pdb'
#    self.pdbfile = 'trunc2a_min.pdb'   #rna
    self.trials = '200'
    self.goback = '20'
    self.temp = '300.0'
    self.numranges = '5'
#    self.numranges = '1'
#    self.numranges = '2'               #rna
    self.dtheta = '30.0,30.0,30.0,30.0,30.0'
#    self.dtheta = '30.0' 
#    self.dtheta = '30.0,30.0'  
    self.reslow = '123,277,354,378,408'
    self.numcont = '22,6,21,12,5'
    self.lowres1 = '284'
    self.highres1 = '350'
#    self.reslow = '24,47'              #rna
#    self.numcont = '7,9'               #rna
#    self.lowres1 = '34'                #rna
#    self.highres1 = '45'               #rna
#    self.reslow = '115'
#    self.numcont = '24'
#    self.lowres1 = '1'
#    self.highres1 = '110'    
#    self.basis = 'backbone'
    self.basis = 'heavy'
#    self.basis = 'all'     
#    self.basis = 'CA'      #do not use for rna
#    self.cutoff = '1.0'    #for backbone
    self.cutoff = '0.8'     #for heavy and all
#    self.cutoff = '3.0'    #for 'CA'
    self.lowrg = '0.0'
    self.highrg = '400.0'
    self.zflag = '0'
    self.zcutoff = '-35.0'
    self.cflag = '0'
    self.confile = 'mmc_constraints.txt'    
    self.directedmc = '0'
    self.nonbondedflag = '0'        #nonbonded feature not currently implemented
    self.nonbondedscale = '1.0'     #used for nonbonded feature not currently implemented
    self.psffilepath = './'         #used for nonbonded feature not currently implemented
    self.psffilename = 'gag_start.psf'  #used for nonbonded feature not currently implemented
    self.parmfilepath = sasconfig.__bin_path__ + 'toppar'  #used for nonbonded feature not currently implemented
    self.parmfilename = 'par_all27_prot_na.inp' #used for nonbonded feature not currently implemented
    self.plotflag = '1'
#    self.seed = '0,123'
    self.seed = '1,123' #to use a specific random sequence each time -- use to create files for tests
    
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
    self.moltype = 'protein'
    self.path = ''
    self.pdbfile = os.path.join(pdb_data_path,'gag_start.pdb')
    self.trials = '50'
    self.goback = '10'
    self.temp = '300.0'
    self.numranges = '5'
    self.dtheta = '30.0,30.0,30.0,30.0,30.0'
    self.reslow = '123,277,354,378,408'
    self.numcont = '22,6,21,12,5'
    self.lowres1 = '284'
    self.highres1 = '350'
    self.basis = 'heavy'
    self.cutoff = '0.8'
    self.lowrg = '0.0'
    self.highrg = '400.0'
    self.zflag = '0'
    self.zcutoff = '0.0'
    self.cflag = '0'
    self.confile = os.path.join(other_data_path,'mmc_constraints.txt')
    self.directedmc = '0'
    self.nonbondedflag = '0'      
    self.nonbondedscale = '1.0'    
    self.psffilepath = ''
    self.psffilename = os.path.join(other_data_path,'gag_start.psf')
    self.parmfilepath = sasconfig.__bin_path__ + 'toppar'
    self.parmfilename = 'par_all27_prot_na.inp'
    self.plotflag = '0'
    self.seed = '1,123' #the '1' is needed for testing so that a specific random sequence is used each time

    self.precision = 3
    self.testflag = True
    
def run_module(self, **kwargs):

    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname, 'string')
    svariables['dcdfile'] = (self.dcdfile, 'string')
    svariables['moltype'] = (self.moltype, 'string')
    svariables['path'] = (self.path, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['trials'] = (self.trials, 'int')
    svariables['goback'] = (self.goback, 'int')
    svariables['temp'] = (self.temp, 'float')
    svariables['numranges'] = (self.numranges, 'int')
    svariables['dtheta'] = (self.dtheta, 'float_array')
    svariables['reslow'] = (self.reslow, 'int_array')
    svariables['numcont'] = (self.numcont, 'int_array')
    svariables['lowres1'] = (self.lowres1, 'int')
    svariables['highres1'] = (self.highres1, 'int')
    svariables['basis'] = (self.basis, 'string')
    svariables['cutoff'] = (self.cutoff, 'float')
    svariables['lowrg'] = (self.lowrg, 'float')
    svariables['highrg'] = (self.highrg, 'float')
    svariables['zflag'] = (self.zflag, 'int')
    svariables['zcutoff'] = (self.zcutoff, 'float')
    svariables['cflag'] = (self.cflag, 'int')
    svariables['confile'] = (self.confile, 'string')
    svariables['directedmc'] = (self.directedmc, 'float')
    svariables['nonbondflag'] = (self.nonbondedflag, 'int')
    svariables['nonbondscale'] = (self.nonbondedscale, 'float')
    svariables['psffilepath'] = (self.psffilepath, 'string')
    svariables['psffilename'] = (self.psffilename, 'string')
    svariables['parmfilepath'] = (self.parmfilepath, 'string')
    svariables['parmfilename'] = (self.parmfilename, 'string')
    svariables['plotflag'] = (self.plotflag, 'int')
    svariables['seed'] = (self.seed, 'int_array')


    error, self.variables = input_filter.type_check_and_convert(svariables)

#    print 'variables: ', self.variables
    if len(error) > 0:
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

    eflag = 0
    monflag = 1
    
    try:
        if kwargs['file_check']:
            error = monomer_monte_carlo_filter.check_protein(self.variables, eflag, monflag)
    except:
            error = monomer_monte_carlo_filter.check_protein(self.variables, eflag, monflag, no_file_check="true")

    if len(error) > 0:
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))


    txtQueue = multiprocessing.JoinableQueue()
    this_monomer_monte_carlo = monomer_monte_carlo.monomer_monte_carlo()
    this_monomer_monte_carlo.main(self.variables, txtQueue)

class gui_mimic_monomer_monte_carlo():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'monomer_monte_carlo'

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
    run_gui = gui_mimic_monomer_monte_carlo(test, paths)
    print("time used: ", time.time() - start)






    
