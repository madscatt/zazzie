'''
Driver method to run the align module
'''

import sys
import os
import shutil
import time
import locale
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.simulate.torsion_angle_md.torsion_angle_md as torsion_angle_md
#import torsion_angle_md as torsion_angle_md
import sassie.interface.input_filter as input_filter
import sassie.interface.torsion_angle_md.torsion_angle_md_filter as torsion_angle_md_filter
#import torsion_angle_md_filter as torsion_angle_md_filter
import multiprocessing


def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname         = 'run_0'
    self.infile          = 'hiv1_gag_ma.pdb'
    self.pdbfile         = 'hiv1_gag_ma.pdb'
    self.outfile         = 'hiv1_gag_ma.dcd'
#    self.infile          = 'pai_vn_start.pdb'
#    self.pdbfile         = 'pai_vn_start.pdb'
#    self.outfile         = 'pai_vn.dcd'
#    self.infile          = 'trunc2a_min.pdb'
#    self.pdbfile         = 'trunc2a_min.pdb'
#    self.outfile         = 'trunc2a.dcd'
#    self.infile          = 'ssDNA.pdb'
#    self.pdbfile         = 'ssDNA.pdb'
#    self.outfile         = 'ssDNA.dcd'
#    self.infile          = 'c36_dsDNA60_min.pdb'
#    self.pdbfile         = 'c36_dsDNA60_min.pdb'
#    self.outfile         = 'c36_dsDNA60_min.dcd'
    self.nsteps          = '10'
    self.topfile         = sasconfig.__bin_path__ + '/toppar/top_all27_prot_na.inp'
    self.parmfile         = sasconfig.__bin_path__ + '/toppar/par_all27_prot_na.inp'
    self.keepout         = '1'
    self.dcdfreq         = '2'
    self.charmmexe       = sasconfig.__bin_path__ + '/charmm.exe'
    self.temperature     = '300.0'
    self.rgforce         = '0.0'
    self.rgvalue         = '0.0'
    self.number_flexible_segments = '1'
    self.pretamd_min_steps = '100'
    self.poll_frequency = '10'
#    self.all_flexible_segnames=['PAI1','VN1']
#    self.all_snumranges=['1','1']
#    self.all_srlow=['3','40']       ##avoid making first and last residues flexible
#    self.all_srnum=['3','89']
#    self.all_moltype=['protein','protein']
    self.all_flexible_segnames=['MA']
#    self.all_flexible_segnames=['DNA1']    
    self.all_snumranges=['1']
#    self.all_srlow=['11']       ##avoid making first and last residues flexible
#    self.all_srnum=['9']
#    self.all_moltype=['dna']
    self.all_srlow=['114']       ##avoid making first and last residues flexible
    self.all_srnum=['20']
    self.all_moltype=['protein']
#    self.all_flexible_segnames=['DNA1','DNA2']
#    self.all_snumranges=['1','1']
#    self.all_srlow=['11','101']       ##avoid making first and last residues flexible
#    self.all_srnum=['9','9']
#    self.all_moltype=['dna','dna']
    self.dna_segnames = ''
    for i in xrange(locale.atoi(self.number_flexible_segments)):
        if self.all_moltype[i] == 'dna':
            self.dna_segnames += self.all_flexible_segnames[i] + ','
    if self.dna_segnames and self.dna_segnames[-1] ==',':
        self.dna_segnames = self.dna_segnames[:-1]
    print 'dna_segnames: ', self.dna_segnames
    self.psegvariables = []
    for i in xrange(locale.atoi(self.number_flexible_segments)):
        self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])

    self.path = ''

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

    self.runname         = 'run_0'
    self.infile          = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
    self.pdbfile         = os.path.join(pdb_data_path,'hiv1_gag_ma.pdb')
    self.outfile         = 'hiv1_gag_ma.dcd'
    self.nsteps          = '10'
    self.topfile         = sasconfig.__bin_path__ + '/toppar/top_all27_prot_na.inp'
    self.parmfile        = sasconfig.__bin_path__ + '/toppar/par_all27_prot_na.inp'
    self.keepout         = '0'
    self.dcdfreq         = '2'
    self.charmmexe       = sasconfig.__bin_path__ + '/charmm.exe'
    self.temperature     = '300.0'
    self.rgforce         = '0.0'
    self.rgvalue         = '0.0'
    self.number_flexible_segments = '1'
    self.pretamd_min_steps = '100'
    self.poll_frequency = '10'
    self.all_flexible_segnames=['MA']
    self.all_snumranges=['1']
    self.all_srlow=['114']
    self.all_srnum=['20']
    self.all_moltype=['protein']
    self.dna_segnames = ''
    for i in xrange(locale.atoi(self.number_flexible_segments)):
        if self.all_moltype[i] == 'dna':
            self.dna_segnames += self.all_flexible_segnames[i] + ','
    if self.dna_segnames and self.dna_segnames[-1] ==',':
        self.dna_segnames = self.dna_segnames[:-1]
    print 'dna_segnames: ', self.dna_segnames
    self.psegvariables = []
    for i in xrange(locale.atoi(self.number_flexible_segments)):
        self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])
    self.path = ''
    self.testflag = True


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname'] = (self.runname,'string')
    svariables['infile'] = (self.infile,'string')
    svariables['pdbfile'] = (self.pdbfile,'string')
    svariables['outfile'] = (self.outfile,'string')
    svariables['nsteps']  = (self.nsteps,'int')
    svariables['topfile'] = (self.topfile,'string')
    svariables['parmfile'] = (self.parmfile,'string')
    svariables['keepout'] = (self.keepout,'int')
    svariables['dcdfreq'] = (self.dcdfreq,'int')
    svariables['charmmexe'] = (self.charmmexe,'string')
    svariables['temperature'] = (self.temperature,'float')
    svariables['rgforce'] = (self.rgforce,'float')
    svariables['rgvalue'] = (self.rgvalue,'float')
    svariables['dna_segnames']  = (self.dna_segnames,'string')
    svariables['number_flexible_segments']  = (self.number_flexible_segments,'int')
    svariables['pretamd_min_steps'] = (self.pretamd_min_steps,'string')
    svariables['poll_frequency']    = (self.poll_frequency,'float')
    svariables['path'] = (self.path,'string')


    error, self.variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print 'error = ', error
        if not(self.testflag):
            sys.exit()
        return error

#    print 'variables: ', self.variables
#    print 'psegvariables: ',self.psegvariables

    error = torsion_angle_md_filter.check_torsion_angle_md(self.variables, self.psegvariables)

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

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    this_torsion_angle_md = torsion_angle_md.torsion_angle_md()
    this_torsion_angle_md.main(self.variables,self.psegvariables,txtQueue)

class gui_mimic_torsion_angle_md():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'torsion_angle_md'

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
    run_gui = gui_mimic_torsion_angle_md(test, paths)
    print "time used: ", time.time() - start

