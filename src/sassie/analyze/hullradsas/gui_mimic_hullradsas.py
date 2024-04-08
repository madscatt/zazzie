'''
Driver method to run the hullradsas module
'''

import sys
import os
import shutil
import time

import sassie.analyze.hullradsas.hullradsas as hullradsas
import sassie.interface.input_filter as input_filter
import sassie.interface.hullradsas.hullradsas_filter as hullradsas_filter
import multiprocessing

def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.run_name = 'run_0'
    self.pdbfile = 'c.pdb'
    self.ofile = 'hullradsas.dat'
    self.testflag = False

    #### end user input ####
    #### end user input ####
    #### end user input ####


def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_hullradsas class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    module_data_path = paths['module_data_path']
    #other_data_path = paths['other_data_path']

    self.run_name = 'run_0'
    #self.pdbfile = os.path.join(other_data_path, 'c.pdb')
    self.pdbfile = 'c.pdb'
    self.ofile = 'hullradsas.dat'

    self.testflag = True

def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['run_name'] = (self.run_name, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['ofile'] = (self.ofile, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

    print("pdbfile = ", self.pdbfile)

    #import sasmol.system as system
    #test_mol = system.Molecule(0)
    #test_mol.read_pdb(self.pdbfile, fastread=True)
    #flag = input_filter.check_binary(self.pdbfile)
    #print('flag = ',flag)

    #resname = test_mol.resname()
    #print(test_mol.resname()[0])

    #error = hullradsas_filter.check_hullradsas(self.variables)

    #if(len(error) > 0):
    #    print('error = ', error)
    #    if not(self.testflag):
    #        sys.exit()
    #    return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    run_name = self.variables['run_name'][0]

    if os.path.exists(os.path.join(run_name, self.module)):
        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_hullradsas = hullradsas.hullradsas()
    this_hullradsas.main(self.variables, txtQueue)

    return

class gui_mimic_hullradsas():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'hullradsas'

    def __init__(self, test, paths):

        if not test:
            user_variables(self)
        else:
            test_variables(self, paths)

        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 1.0.
    test = True  # option to run with test variables not implemented in 1.0.
    paths = None

# We are thinking of defining the install path so the gui mimic can be run from anywhere as long as it is called from that particular python
# That way, the test files will always be available to the user.
    if test:
        pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
        dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_hullradsas(test, paths)
    print("time used: ", time.time() - start)
