'''
Driver method to run the interpolate module
'''

import sys
import os
import shutil
import time

import sassie.tools.data_interpolation.data_interpolation as data_interpolation
#import data_interpolation as data_interpolation
import sassie.interface.input_filter_sasmol as input_filter
import sassie.interface.data_interpolation.interpolate_filter as interpolate_filter
import multiprocessing


def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.run_name = 'run_0'
    self.expdata = 'sans_data.sub'
#    self.expdata = 'trunc2a_saxs.sub'
    self.ofile = 'sans_data.dat'
#    self.ofile = 'trunc2a.dat'
    self.io = '0.04'
#    self.io = '0.031'
    self.ioe = '0.001'
    self.dq = '.02'
#    self.dq = '0.007'
    self.maxpoints = '16'
#    self.maxpoints = '72'
    self.plotflag = '0'

    self.testflag = False

    #### end user input ####
    #### end user input ####
    #### end user input ####


def test_variables(self, paths):
    '''
    users of gui_mimic as a driver script to run this module should not edit the values below as they
    are used for development tests

    this module defines variables that will be used to test the module as well as its input filter
    variables are defined outside the gui_mimic_data_interpolation class so that they can be used by these other programs

    '''

    pdb_data_path = paths['pdb_data_path']
    dcd_data_path = paths['dcd_data_path']
    module_data_path = paths['module_data_path']
    other_data_path = paths['other_data_path']

    self.run_name = 'run_0'
    self.expdata = os.path.join(other_data_path, 'sans_data.sub')
    self.ofile = 'sans_data.dat'
    self.io = '0.04'
    self.ioe = '0.001'
    self.dq = '.02'
    self.maxpoints = '16'
    self.plotflag = '0'

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
    svariables['expdata'] = (self.expdata, 'string')
    svariables['ofile'] = (self.ofile, 'string')
    svariables['io'] = (self.io, 'float')
    svariables['ioe'] = (self.ioe, 'float')
    svariables['dq'] = (self.dq, 'float')
    svariables['maxpoints'] = (self.maxpoints, 'int')
    svariables['plotflag'] = (self.plotflag, 'int')

    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = interpolate_filter.check_interpolate(self.variables)
    except:
        error = interpolate_filter.check_interpolate(
            self.variables, no_file_check="true")

    if(len(error) > 0):
        print('error = ', error)
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['test_filter']:
            return error
    except:
        pass

    run_name = self.variables['run_name'][0]

    if os.path.exists(os.path.join(run_name, self.module)):
        shutil.rmtree(os.path.join(run_name, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_interpolate = data_interpolation.data_interpolation()
    this_interpolate.main(self.variables, txtQueue)


class gui_mimic_data_interpolation():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'data_interpolation'

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
        module_data_path = os.path.join(os.path.dirname(os.path.realpath(
            __file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_data_interpolation(test, paths)
    print("time used: ", time.time() - start)
