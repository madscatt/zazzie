'''
Driver method to run the density_plot module
'''

import sys
import string
import os
import shutil
import time

import sassie.interface.input_filter as input_filter
import sassie.analyze.density_plot as density_plot
#import density_plot as density_plot
import sassie.interface.density_plot_filter as density_plot_filter
#import density_plot_filter as density_plot_filter
import multiprocessing

def user_variables(self, **kwargs):

    #### BEGIN USER EDIT
    #### BEGIN USER EDIT
    #### BEGIN USER EDIT

    self.runname = 'run_0'
    self.path = ('./')
    self.dcdfile = os.path.join(self.path,'hiv1_gag_20_frames.dcd')
    self.pdbfile = os.path.join(self.path,'hiv1_gag.pdb')
#    self.dcdfile = os.path.join(self.path,'pai_vn_20_frames.dcd')
#    self.pdbfile = os.path.join(self.path,'pai_vn_start.pdb')
#    self.dcdfile = os.path.join(self.path,'trunc2a_20_frames.dcd')
#    self.pdbfile = os.path.join(self.path,'trunc2a_min.pdb')            
    self.ofile = 'test'
    self.nsegments = '1'
#    self.nsegments = '2'
    self.xlength = '100.0'
    self.ylength = '100.0'
    self.zlength = '100.0'
    self.gridsp = '5.0'
    self.equalweights = '1'
#    self.equalweights = '0'
    self.weightsfile = ''
#    self.weightsfile = 'weights_file_density_plot.txt'
#    self.save_occupancy = 'Y'    
    self.save_occupancy = 'N'

    self.segvariables = [[u'1', '6', '123', u'CA', u'GAG']]
#    self.segvariables = [[u'1', '1', '379', u'CA', u'PAI1'],[u'2', '1,40', '39,130', u'CA', u'VN1']] 
#    self.segvariables = [[u'3', '1,31,47', '30,46,80', u'P', u'TR2A']]   
    self.testflag = False

    #### END USER EDIT
    #### END USER EDIT
    #### END USER EDIT

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
    self.dcdfile = os.path.join(dcd_data_path,'hiv1_gag_20_frames.dcd')
    self.pdbfile = os.path.join(pdb_data_path,'hiv1_gag.pdb')
    self.ofile = 'test'
    self.nsegments = '1'
    self.xlength = '100.0'
    self.ylength = '100.0'
    self.zlength = '100.0'
    self.gridsp = '5.0'
    self.equalweights = '1'
    self.weightsfile = ''    
    self.save_occupancy = 'N'

    self.segvariables = [[u'1', '6', '123', u'CA', u'GAG']]
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

    svariables['runname']      = (self.runname,'string')
    svariables['path']         = (self.path, 'string')
    svariables['dcdfile']  	   = (self.dcdfile,'string')
    svariables['pdbfile']  	   = (self.pdbfile,'string')
    svariables['ofile']        = (self.ofile,'string')
    svariables['nsegments']      = (self.nsegments,'int')    
    svariables['xlength']      = (self.xlength,'float')
    svariables['ylength']      = (self.ylength,'float')
    svariables['zlength']      = (self.zlength,'float')
    svariables['gridsp']       = (self.gridsp,'float')
    svariables['equalweights'] = (self.equalweights,'int')
    svariables['weightsfile']  = (self.weightsfile,'string')
    svariables['save_occupancy']  = (self.save_occupancy,'string')

    error,self.variables=input_filter.type_check_and_convert(svariables)

    if(len(error)>0):
        print 'error = ',error
        if not(self.testflag):
            sys.exit()
        return error

    try:
        if kwargs['file_check']:
            error = density_plot_filter.check_density_plot(self.variables,self.segvariables)
    except:
            error = density_plot_filter.check_density_plot(self.variables,self.segvariables,no_file_check="true")
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

    txtQueue = multiprocessing.JoinableQueue()
    density_plot.density(self.variables, self.segvariables,txtQueue)


class gui_mimic_density_plot():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'density_plot'

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
    run_gui = gui_mimic_density_plot(test, paths)
    print "time used: ", time.time() - start







