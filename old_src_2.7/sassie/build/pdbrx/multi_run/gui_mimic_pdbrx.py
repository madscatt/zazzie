'''
Driver method to run the align module
'''

import sys
import os
import shutil
import time
import sassie.build.pdbrx.pdb_rx as pdb_rx
import sassie.util.sasconfig as sasconfig
import sassie.interface.input_filter as input_filter
#import sassie.interface.pdbrx_filter as pdbrx_filter
import multiprocessing

import logging

sys.path.append('./')

def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname = 'run_0'
    self.pdbfile = 'testing/data/5E3L.pdb'
    self.topfile = os.path.join(sasconfig.__bin_path__, \
                           'toppar', 'top_all27_prot_na.inp')
    self.use_defaults = False

    #### end user input ####
    #### end user input ####
    #### end user input ####

    try:
        if kwargs['runname']:
            self.runname = kwargs['runname']
        if kwargs['pdbfile']:
            self.pdbfile = kwargs['pdbfile']
        if kwargs['use_defaults']:
            self.use_defaults = kwargs['use_defaults']

    except:
        pass


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname, 'string')
    svariables['pdbfile'] = (self.pdbfile, 'string')
    svariables['topfile'] = (self.topfile, 'string')
    svariables['defaults'] = (self.use_defaults, 'boolean')

    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
        return error

    # try:
    #    if kwargs['file_check']:
    #        error = align_filter.check_align(self.variables)
    # except:
    #    error = align_filter.check_align(self.variables, no_file_check="true")

    # if(len(error) > 0):
    #    print 'error = ', error
    #    if not(self.testflag):
    #        sys.exit()
    #    return error

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()

    logging.basicConfig()
    scan = pdb_rx.PDBRx()
    scan.main(self.variables, txtQueue)

    this_text = txtQueue.get(True, timeout=0.1)


class gui_mimic_pdbrx():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'pdbrx'

    def __init__(self, test, **kwargs):

        if not test:
            user_variables(self, **kwargs)
        else:
            test_variables(self, paths)

        run_module(self)

if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 2.0.
    paths = None

    start = time.time()
    run_gui = gui_mimic_pdbrx(test)
    print "time used: ", time.time() - start
