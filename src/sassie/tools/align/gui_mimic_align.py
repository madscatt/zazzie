'''
Driver method to run the align module
'''

import sys
import time
import os
import shutil

#import align as align

import sassie.tools.align.align as align
import sassie.interface.input_filter as input_filter
import sassie.interface.align.align_filter as align_filter
import multiprocessing

def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.inpath = ''
    self.pdbmol1 = os.path.join('./', 'hiv1_gag.pdb')
    self.pdbmol2 = os.path.join('./', 'hiv1_gag.pdb')
#    self.infile = os.path.join('./', 'hiv1_gag.pdb')
    self.infile = os.path.join('./', 'hiv1_gag_20_frames.dcd')
    self.ofile = 'aligned_hiv1_gag_20_frames.dcd'
#    self.ofile = 'aligned_hiv1_gag_zcutoff.dcd'
#    self.ofile = 'aligned_hiv1_gag_20_frames.pdb'
#    self.ofile = 'aligned_hiv1_gag_zcutoff.pdb'
    self.ofile = 'aligned_hiv1_gag.pdb'
    self.basis1 = 'CA'
    self.basis2 = 'CA'
    self.lowres1 = '145'
    self.lowres2 = '145'
    self.highres1 = '350'
    self.highres2 = '350'
    self.ebasis1 = 'None'
    self.ebasis2 = 'None'
    self.zflag = False
#    self.zflag = True
    self.zcutoff = '-66.0' #If there are ANY atoms with a z-value less than the cutoff the frame will not be written to disk.

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
    module_data_path = paths['module_data_path']
    other_data_path = paths['other_data_path']    

    self.runname = 'run_0'
    self.inpath = ''
    self.pdbmol1 = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
    self.pdbmol2 = os.path.join(pdb_data_path, 'hiv1_gag.pdb')
    self.infile = os.path.join(dcd_data_path, 'hiv1_gag_20_frames.dcd')
    self.ofile = 'aligned_hiv1_gag_20_frames.dcd'
    self.basis1 = 'CA'
    self.basis2 = 'CA'
    self.lowres1 = '145'
    self.lowres2 = '145'
    self.highres1 = '350'
    self.highres2 = '350'
    self.ebasis1 = 'None'
    self.ebasis2 = 'None'

    self.zflag = False
    self.zcutoff = '0.0'

    self.precision = 3


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the Test_Align class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname, 'string')
    svariables['path'] = (self.inpath, 'string')
    svariables['pdbmol1'] = (self.pdbmol1, 'string')
    svariables['pdbmol2'] = (self.pdbmol2, 'string')
    svariables['infile'] = (self.infile, 'string')
    svariables['ofile'] = (self.ofile, 'string')
    svariables['basis1'] = (self.basis1, 'string')
    svariables['basis2'] = (self.basis2, 'string')
    svariables['lowres1'] = (self.lowres1, 'int')
    svariables['lowres2'] = (self.lowres2, 'int')
    svariables['highres1'] = (self.highres1, 'int')
    svariables['highres2'] = (self.highres2, 'int')
    svariables['ebasis1'] = (self.ebasis1, 'string')
    svariables['ebasis2'] = (self.ebasis2, 'string')

    svariables['zflag'] = (self.zflag, 'boolean')
    svariables['zcutoff'] = (self.zcutoff, 'float')
    
    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
#        print 'error = ', error
#        sys.exit()
        return error
    try:
        if kwargs['file_check']:
            error = align_filter.check_align(self.variables)
    except:
        error = align_filter.check_align(self.variables, no_file_check="true")

    if(len(error) > 0):
#        print 'error = ', error
#        sys.exit()
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
    this_align = align.align()
    this_align.main(self.variables, txtQueue)


class gui_mimic_align():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'align'

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
            __file__)), '..', '..', 'data', 'interface', 'align') + os.path.sep

        paths = {'pdb_data_path': pdb_data_path,
                 'dcd_data_path': dcd_data_path, 'other_data_path': other_data_path, 'module_data_path': module_data_path}

    start = time.time()
    run_gui = gui_mimic_align(test, paths)
    print "time used: ", time.time() - start


###Clean up json and log files that are left when runs fail!


