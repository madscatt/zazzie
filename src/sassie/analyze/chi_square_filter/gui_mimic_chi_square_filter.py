'''
Driver method to run the chi_square filter module
'''

import sys
import string
import os
import shutil
import time

import sassie.analyze.chi_square_filter.chi_square_filter as chi_square_filter
#import chi_square_filter as chi_square_filter
import sassie.interface.input_filter as input_filter
import sassie.interface.chi_square_filter.chi_square_filter_filter as chi_square_filter_filter
#import chi_square_filter_filter as chi_square_filter_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_2'
    self.path = './'
    self.saspaths = os.path.join(self.path, 'diUB', 'run_0', 'sascalc', 'neutron_D2Op_100')
#    self.saspaths = os.path.join(self.path, 'pai-vn', 'run_2','sascalc', 'neutron_D2Op_85')+','+os.path.join(self.path, 'pai-vn', 'run_2', 'sascalc','neutron_D2Op_0')
    self.sasintfiles = os.path.join(self.path, 'K48_UBA2_org.dat')
#    self.sasintfiles = os.path.join(self.path,'trunc.dat')
#    self.sasintfiles = os.path.join(self.path,'85p1.dat')+','+os.path.join(self.path,'0p.dat')
    self.io = '0.1229'
#    self.io = '0.031'
#    self.io = '0.013,0.85'
#    self.number_of_weight_files = '0'
    self.number_of_weight_files = '2'
#    self.basis_string = ''
    self.basis_string = '(rg<19.4) and (x2<27)' + ',' + 'x2<25'
#    self.basis_string = 'x2<50'+','+'(rg<40) and (rg>20) and (x2<8)'
#    self.weight_file_names = ''
#    self.weight_file_names = 'x2_lt_50.txt'+','+'0_x2_8_20_rg_40.txt'
    self.weight_file_names = 'x2_19_27.txt' + ',' + 'x2_lt_25.txt'
    self.sastype = '0'
    self.reduced_x2 = '1'
    self.plotflag = '0'

    self.testflag = False
    
    ### END USER INPUT ###
    ### END USER INPUT ###
    ### END USER INPUT ###


def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_align class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    other_data_path = paths['other_data_path']
    module_data_path = paths['module_data_path']

    self.runname = 'run_0'
    self.path = ''
    self.saspaths = os.path.join(
        other_data_path, 'diUb', 'run_0', 'sascalc', 'neutron_D2Op_100')
    self.sasintfiles = os.path.join(other_data_path, 'K48_UBA2_org.dat')
    self.io = '0.1229'
    self.number_of_weight_files = '0'
    self.basis_string = ''
    self.weight_file_names = ''
    self.sastype = '0'
    self.reduced_x2 = '1'
    self.plotflag = '0'

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
    svariables['saspaths'] = (self.saspaths, 'string')
    svariables['sasintfiles'] = (self.sasintfiles, 'string')
    svariables['io'] = (self.io, 'float_array')
    svariables['number_of_weight_files'] = (self.number_of_weight_files, 'int')
    svariables['basis_string'] = (self.basis_string, 'string')
    svariables['weight_file_names'] = (self.weight_file_names, 'string')
    svariables['sastype'] = (self.sastype, 'int')
    svariables['reduced_x2'] = (self.reduced_x2, 'int')
    svariables['plotflag'] = (self.plotflag, 'int')
    svariables['path'] = (self.path, 'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

#    print 'variables: ', self.variables
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()

    try:
        if kwargs['file_check']:
            error = chi_square_filter_filter.check_chi_square_filter(
                self.variables)
    except:
        error = chi_square_filter_filter.check_chi_square_filter(
            self.variables, no_file_check="true")

    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_chi_square_filter = chi_square_filter.chi_square_filter()
    this_chi_square_filter.main(self.variables, txtQueue)


class gui_mimic_chi_square_filter():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'chi_square_filter'

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
            __file__)), '..', '..', 'data', 'interface', 'chi_square_filter') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_chi_square_filter(test, paths)
    print "time used: ", time.time() - start
