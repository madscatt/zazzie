'''
Driver method to run the contrast calculator module
'''

import sys
import os
import shutil
import time

import sassie.tools.contrast_calculator as contrast_calculator
#import contrast_calculator as contrast_calculator
import sassie.interface.input_filter as input_filter
import sassie.interface.contrast_calculator_filter as contrast_calculator_filter
#import contrast_calculator_filter as contrast_calculator_filter
import multiprocessing

def user_variables(self, **kwargs):

    #### user input ####
    #### user input ####
    #### user input ####

    self.runname = 'run_0'
    self.inpath = './'
    self.outfile = 'test'
    self.numfiles = '2'
#    self.numfiles = '1'
#    self.numfiles = '0'    
    self.solute_conc = '1.0'
    self.d2ostep = '5'
    self.fexchp = '0.95'
    self.fexchn = '1.0'
#    self.seqfiles = ['protein_sequence.txt', 'dna_sequence.txt']
#    self.seqfiles = ['pai_seq.txt']
#    self.seqfiles = ['protein_sequence.txt']
#    self.seqfiles = ['dna_sequence.txt']
#    self.seqfiles = ['rna_sequence.txt']
#    self.seqfiles = ['trunc2a_min.pdb']
#    self.seqfiles = ['hiv1_gag.pdb']
#    self.seqfiles = ['c36_dsDNA60_min.pdb']
#    self.seqfiles = ['pai_seq.txt','vn_seq.txt']
#    self.seqfiles = ['skp_trimer.pdb','ompA.pdb']
    self.seqfiles = ['pai_seq.txt','c36_dsDNA60_min.pdb']
    self.numunits = ['1', '1']
#    self.numunits = ['1']
#    self.numunits = ['2'] 
    self.fracdeut = ['0', '0']
#    self.fracdeut = ['0']
#    self.fracdeut = ['0.0','0.6']
    self.moltype = ['protein', 'dna']
#    self.moltype = ['dna']
#    self.moltype = ['protein']
#    self.moltype = ['rna']
#    self.moltype = ['protein','protein']
#    self.isFasta = ['1', '1']
#    self.isFasta = ['1']
#    self.isFasta = ['0']
#    self.isFasta = ['1', '1']
#    self.isFasta = ['0', '0']
    self.isFasta = ['1', '0']
    self.plotflag = '1'            


#    self.numsolv = '0'
    self.numsolv = '2'
    self.solv_comp = ['NaCl','KCl']
    self.solv_conc = ['0.15','0.05']

 #   self.number_of_chemicals = '0'
#    self.number_of_chemicals = '1'
    self.number_of_chemicals = '2'
    self.formula_array = ['(C3H4O3)12', '(C3H4O3)12']
#    self.formula_array = ['(C42H82NO8P)130']    
    
    self.number_exchangeable_hydrogens = ['12', '5']
    self.fraction_exchangeable_hydrogens = ['0.95', '0.45']
    self.mass_density = ['1.1', '1.3']
#    self.number_exchangeable_hydrogens = ['0']
#    self.fraction_exchangeable_hydrogens = ['0.0']
#    self.mass_density = ['1.0']    

    self.testflag = False 

    #### end user input ####
    #### end user input ####
    #### end user input ####

def test_variables(self,paths):
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
    self.inpath = other_data_path
    self.outfile = 'test'
    self.numfiles = '1' 
    self.solute_conc = '1.0'
    self.d2ostep = '5'
    self.fexchp = '0.95'
    self.fexchn = '1.0'
    self.seqfiles = ['pai_seq.txt']
    self.numunits = ['1']
    self.fracdeut = ['0']
    self.moltype = ['protein']
    self.isFasta = ['1']
    self.plotflag = '0'            

    self.numsolv = '0'
    self.solv_comp = []
    self.solv_conc = []

    self.number_of_chemicals = '0'
    self.formula_array = []
    self.number_exchangeable_hydrogens = []
    self.fraction_exchangeable_hydrogens = []
    self.mass_density = []

    self.testflag = True
    self.precision = 3

def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables={}

    svariables['runname'] = (self.runname,'string')
    svariables['inpath'] = (self.inpath,'string')
    svariables['outfile'] = (self.outfile,'string')
    svariables['numfiles'] = (self.numfiles,'int')
    svariables['solute_conc'] = (self.solute_conc,'float')
    svariables['d2ostep'] = (self.d2ostep,'int')
    svariables['numsolv'] = (self.numsolv,'int')
    svariables['fexchp'] = (self.fexchp,'float')
    svariables['fexchn'] = (self.fexchn,'float')
    svariables['number_of_chemicals'] = (self.number_of_chemicals,'int')
    svariables['plotflag'] = (self.plotflag,'int')

    error, self.variables = input_filter.type_check_and_convert(svariables)
#    print 'error, length: ', error, len(error)
#    print 'variables in run module: ', self.variables

    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error


    error1 = contrast_calculator_filter.check_numfiles(self.variables['numfiles'][0])
    error2 = contrast_calculator_filter.check_numsolvcomp(self.variables['numsolv'][0])
    error3 = contrast_calculator_filter.check_numchemcomp(self.variables['number_of_chemicals'][0])
    if len(error1) > 0 or len(error2) > 0 or len(error3) > 0:
        error = error1 + error2 + error3
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error    

    else:
    
        error = contrast_calculator_filter.check_contrast(self.variables)
        if len(error) > 0:
            print 'error = ', error
            if not(self.testflag):
                sys.exit()
            return error


    self.ivariables = []
    if(int(self.numfiles) > 0):
        for i in xrange(int(self.numfiles)):       
            self.ivariables.append([self.seqfiles[i], self.numunits[i], self.fracdeut[i], self.moltype[i], self.isFasta[i]])
        error = contrast_calculator_filter.check_ivariables(self.inpath,self.ivariables)   
        if len(error) > 0:
            print 'error = ', error
            if not(self.testflag):
                sys.exit()
            return error


    self.solvvariables = []
    if(int(self.numsolv) > 0):         
        error, self.solv_formula = input_filter.check_and_convert_formula(self.solv_comp)
        if(len(error) > 0):
            print 'error = ', str(error)
            if not(self.testflag):
                sys.exit()
            return error

        else:

            for i in xrange(int(self.numsolv)):
#                self.solvvariables.append([self.solv_comp[i],self.solv_conc[i]])
                self.solvvariables.append([self.solv_formula[i], self.solv_conc[i]])
            error = contrast_calculator_filter.check_solvvariables(self.solvvariables)
            if len(error) > 0:
                print 'error = ', error
                if not(self.testflag):
                    sys.exit()
                return error
            

    self.chemvariables = []            
    if(int(self.number_of_chemicals) > 0):
        error, self.formulas = input_filter.check_and_convert_formula(self.formula_array)
        if(len(error) > 0):
            print 'error = ', str(error)
            if not(self.testflag):
                sys.exit()
            return error
           
        else:

            for i in xrange(int(self.number_of_chemicals)):
                this_chemical_formula = self.formulas[i]
                this_number_exchangeable_hydrogens = self.number_exchangeable_hydrogens[i]
                this_fraction_exchangeable_hydrogens = self.fraction_exchangeable_hydrogens[i]
                this_mass_density = self.mass_density[i]
                self.chemvariables.append([this_chemical_formula, this_number_exchangeable_hydrogens, this_fraction_exchangeable_hydrogens, this_mass_density])
            error = contrast_calculator_filter.check_chemvariables(self.chemvariables)
            if(len(error) > 0):
                print 'error = ', str(error)
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
    contrast_calculator.contrast(self.variables,self.ivariables,self.solvvariables,self.chemvariables,txtQueue)


class gui_mimic_contrast_calculator():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'contrast_calculator'

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
    run_gui = gui_mimic_contrast_calculator(test, paths)
    print 'time used: ', time.time() - start


