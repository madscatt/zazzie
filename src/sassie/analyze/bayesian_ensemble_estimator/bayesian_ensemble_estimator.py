'''
    ENSEMBLE_FIT is a module that will read in a list of scattering
    profiles, along with an optional second set of restraints, and it will
    reweight populations in subsets of these scattering profiles to find
    the best fit to the experimental data. Here, "best fit" means a balance
    between a model with a low chi^2 value that simultaneously uses as few
    fitting parameters (number of profiles) as possible, thereby trying to
    avoid overfitting. This module will output information regarding this
    "best mode", and it has the option to also present the relative
    ability of all other models in comparison. If you use this module,
    please be sure to cite the references below.

    REFERENCE:

    Bowerman et al.
    Journal of Chemical Theory and Computation, 13, pp 2418-2429 (2017)

    Bowerman et al.
    Placeholder for Release Paper
'''
import numpy as np
import os
import io
import sys
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import zipfile
import itertools
import time
import pickle
import subprocess

if sasconfig.__level__ == 'DEBUG':
    DEBUG = True
else:
    DEBUG = False

app = 'bayesian_ensemble_estimator'

# Location of python version
pythexec = '/share/apps/local/anacondaz/bin/python'
#pythexec = '/share/apps/local/anaconda2/bin/python'
# Location of the ensemble modeling runtime
executable = '/share/apps/local/anacondaz/lib/python2.7/site-packages/sassie/analyze/bayesian_ensemble_estimator/bayesian_ensemble_estimator_parallel_routine.py'
#executable = '/home/sbowerma/sassie_module/alpha/for_jec/renamed/plotQueue_dict/bayesian_ensemble_estimator_parallel_routine.py'
# Location of the proper mpiexec commmand
mpiexec = '/share/apps/local/anacondaz/bin/mpirun'
#mpiexec = '/share/apps/local/anaconda2/bin/mpirun'


class module_variables():
    def __init__(self, parent=None):
        self.app = app


class efvars():
    ''' _E_nsemble _F_it _VAR_iable_S_ '''

    def __init__(self, parent=None):
        pass

######## MPI Environment ########
# comm=MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.Get_size()
#################################


# This is the top-level object for managing the Bayesian Ensemble Fit Routine
'''
    This object builds up all the necessary variables and stores them to a 
    "pickle" file object, which is passed to the parallel runtime, which
    actually conducts the complete MC fitting protocol for each sub-basis.
    This object also handles likelihood comparisons and final outputs.
'''


class ensemble_routine(object):
    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput, plotQueues):
        self.mvars = module_variables()
        self.efvars = efvars()
        self.run_utils = module_utilities.run_utils(app, txtOutput)
        self.run_utils.setup_logging(self)
        pgui = self.run_utils.print_gui
        self.log.debug('In main()')
        self.UnpackVariables(input_variables)
        self.run_utils.general_setup(self)
        self.Initialize()
        self.PickleVars()
        pgui('Beginning iterative Bayesian sub-basis fitting.\n\n')
        pgui('STATUS\t0.0001\n\n')
        self.EnsembleFit()
        self.epilogue(plotQueues)

    def UnpackVariables(self, variables):
        '''
        extract variables from gui/mimic into passable mvars object
        '''
        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui
        log.debug('in UnpackVariables()')
        # Standard Options
        mvars.runname = variables['runname'][0]
        mvars.sas_data = variables['sas_data'][0]
        mvars.theoretical_profiles_zip = variables['theoretical_profiles_zip'][0]
        mvars.posterior_burn = variables['posterior_burn'][0]
        mvars.max_iterations = variables['max_iterations'][0]
        mvars.number_of_MCs = variables['number_of_MCs'][0]
        mvars.d_max = variables['d_max'][0]
        mvars.nproc = variables['nproc'][0]

        # Advanced Options
        mvars.use_all = variables['use_all'][0]
        mvars.every = variables['every'][0]
        mvars.shansamp = variables['shansamp'][0]
        mvars.auxiliary_data = variables['auxiliary_data'][0]
        mvars.use_bic = variables['use_bic'][0]
        mvars.zeroing_threshold = variables['zeroing_threshold'][0]
        mvars.walk_one = variables['walk_one'][0]
        mvars.sigma = variables['sigma'][0]

        # Pass debugging along to parallel process
        mvars.debug = DEBUG

        # Need to fix the mvars logic for pickling for some reason...
        if mvars.use_all.lower() == 'true':
            mvars.use_all = True
        else:
            mvars.use_all = False
        if mvars.every.lower() == 'true':
            mvars.use_all = True
        else:
            mvars.use_all = False
        if mvars.use_bic.lower() == 'true':
            pgui("Using BIC to evaluate model quality.\n\n")
            mvars.use_bic = True
        else:
            pgui("Using AIC to evaluate model quality.\n\n")
            mvars.use_bic = False
        if mvars.walk_one.lower() == 'true':
            pgui("Will increment one population at a time.\n\n")
            mvars.walk_one = True
        else:
            mvars.walk_one = False

        if (mvars.d_max == 0) and (mvars.shansamp.lower() == 'true'):
            log.error("ERROR: Please set D_max to a non-zero value.\n\n")
            quit()
        elif (mvars.d_max > 0) and (mvars.shansamp.lower() == 'false'):
            mvars.shansamp = False
            pgui("WARNING: D_max value non-zero but Shannon Sampling turned off\n\n")
            mvars.d_max = 0.0
        else:
            mvars.shansamp = True
        return

    def Initialize(self):
        '''
        Prepare efvars object
        '''
        mvars = self.mvars
        efvars = self.efvars
        pgui = self.run_utils.print_gui

        efvars.output_folder = os.path.join(mvars.runname, app)
        if not os.path.isdir(efvars.output_folder):
            os.mkdir(efvars.output_folder)

        # Clean out old log files, backs-up most recent version
        existingfiles = os.listdir(efvars.output_folder)
        for item in existingfiles:
            if item.endswith('BayesMC.log'):
                oldname = os.path.join(efvars.output_folder, item)
                newname = os.path.join(efvars.output_folder, item+'.BAK')
                os.rename(oldname, newname)

        efvars.Q, efvars.I, efvars.ERR = np.genfromtxt(mvars.sas_data,
                                                       unpack=True,
                                                       usecols=(0, 1, 2))

        # Set up _B_ayesian _E_nsemble _F_it _LOG_ file for progress, etc.
        efvars.status_file = os.path.join(efvars.output_folder,
                                          '._status.txt')
        efvars.aic_file = os.path.join(efvars.output_folder,
                                       mvars.runname+'_IC_vs_Size.dat')
        if os.path.isfile(efvars.aic_file):
            aic_bak = os.path.join(efvars.output_folder,
                                   mvars.runname+'_IC_vs_Size.dat.BAK')
            os.rename(efvars.aic_file, aic_bak)

        if mvars.shansamp:
            pgui('Using Shanning Sampling (chi^2 free).\n\n')
            statf = open(efvars.status_file, 'w')
            statf.write('STATUS\t0.0001\n')
            statf.close()
            self.BuildShannonSamples()
            pgui('There are '+str(efvars.number_of_channels)
                 + ' Shannon channels.\n\n')
            efvars.samples_Q = efvars.Q[efvars.shannon_samples]
            efvars.samples_I = efvars.I[efvars.shannon_samples]
            efvars.samples_ERR = efvars.ERR[efvars.shannon_samples]
            efvars.num_q = efvars.number_of_channels
        else:
            pgui('Using standard chi^2 (no Shannon Sampling).\n\n')
            efvars.samples_Q = efvars.Q
            efvars.samples_I = efvars.I
            efvars.samples_ERR = efvars.ERR
            efvars.num_q = len(efvars.Q)

        if mvars.auxiliary_data != '':
            pgui('Including second data set.\n\n')
            efvars.include_second_dimension = True
            efvars.aux_data = np.array([], dtype=float)
            efvars.aux_error = np.array([], dtype=float)
            efvars.aux_data, efvars.aux_error = np.genfromtxt(
                mvars.auxiliary_data, usecols=(0, 1), unpack=True, dtype=float)
            #auxiliary_data_file = open(mvars.auxiliary_data,'r')
            # for line in auxiliary_data_file:
            #    line_string = auxiliary_data_file.readline()
            #    line_string_split = line_string.split()
            #    efvars.aux_data = np.append(efvars.aux_data,
            #                                float(line_string_split[0]))
            #    efvars.aux_error= np.append(efvars.aux_error,
            #                                float(line_string_split[1]))
            efvars.num_aux = len(efvars.aux_data)
        else:
            pgui('Not using secondary data set.\n\n')
            efvars.include_second_dimension = False
            efvars.num_aux = 0

        efvars.num_points = efvars.num_q + efvars.num_aux
        self.UnpackTheoreticalProfiles()
        self.BuildBasis()
        efvars.aic_tracker = np.array([], dtype=float)
        efvars.chi_tracker = np.array([], dtype=float)

        return

    def BuildShannonSamples(self):
        '''
        Build Shannon Channels for Calculating Chi_free
        '''
        mvars = self.mvars
        efvars = self.efvars

        shannon_width = np.pi/mvars.d_max
        n_channels = int(np.floor(efvars.Q[-1]/shannon_width))

        shannon_points = np.random.choice(np.where(
            np.logical_and(efvars.Q >= 0,
                           efvars.Q < shannon_width))[0], 3001)
        for idx in range(1, n_channels):
            this_channel = np.random.choice(np.where(
                np.logical_and(efvars.Q >= idx * shannon_width,
                               efvars.Q < (idx+1) * shannon_width))[0], 3001)
            shannon_points = np.vstack((shannon_points, this_channel))

        efvars.shannon_samples = np.copy(shannon_points)
        efvars.number_of_channels = n_channels
        return

    def UnpackTheoreticalProfiles(self):
        mvars = self.mvars
        efvars = self.efvars

        efvars.profiles_directory = os.path.join(efvars.output_folder,
                                                 'TheoreticalProfiles')
        if os.path.isdir(efvars.profiles_directory):
            os.system('rm -r '+efvars.profiles_directory)
        os.system('echo \"efvars.profiles_directory - '
                  + efvars.profiles_directory+'\" >> ' + self.logfile)
        archive = zipfile.ZipFile(mvars.theoretical_profiles_zip)
        archive.extractall(efvars.profiles_directory)
        archive.close()

        return

    def BuildBasis(self):
        '''
        Stores individual theoretical profiles into single array(s) for
        faster ensemble averaging.
        '''
        mvars = self.mvars
        efvars = self.efvars
        pgui = self.run_utils.print_gui
        log = self.log
        pgui('Building Theoretical Basis set.\n\n')

        flist_as_str = os.path.join(efvars.profiles_directory, 'filelist.txt')
        efvars.name_array = np.genfromtxt(flist_as_str, usecols=0, dtype=str)
        efvars.number_of_profiles = len(efvars.name_array)
        # Below is used to quickly build subsets, even if it seems redundant
        efvars.id_array = np.arange(efvars.number_of_profiles, dtype=int)
        efvars.subspace_dict = {}
        efvars.full_scattering_basis = np.zeros((efvars.number_of_profiles,
                                                 len(efvars.Q)), dtype=float)

        if efvars.include_second_dimension:
            try:
                efvars.aux_name_array = np.genfromtxt(flist_as_str,
                                                      usecols=1, dtype=str)
                efvars.number_of_aux_profiles = len(efvars.aux_name_array)
                if efvars.number_of_aux_profiles != efvars.number_of_profiles:
                    log.error('ERROR: A different number of extra profiles'
                              + ' ('+str(efvars.number_of_aux_profiles)+')'
                              + ' have been supplied than the number of'
                              + ' scattering profiles ('
                              + str(efvars.number_of_profiles)+').')
                size_tuple = (efvars.number_of_profiles, efvars.num_aux)
                efvars.full_extra_basis = np.zeros(size_tuple, dtype=float)
            except:
                log.error('ERROR: In order to use a second dataset,'
                          + ' users must include the list of files in'
                          + ' the second column of the'
                          + ' \"filelist.txt\" file.\"')
        else:
            efvars.full_extra_basis = None

        for idx in efvars.id_array:
            try:
                scatter_file = os.path.join(efvars.profiles_directory,
                                            efvars.name_array[idx])
                the_q, the_i = np.genfromtxt(scatter_file, dtype=float,
                                             usecols=(0, 1), unpack=True)

                if not np.all(the_q == efvars.Q):
                    log.error('ERROR: Q-values of '+scatter_file
                              + ' don\'t match the experimental values.')
                efvars.full_scattering_basis[idx] = the_i
            except:
                log.error('ERROR: Data file ('+scatter_file+')'
                          + ' not found. Please inspect the .zip archive'
                          + ' and re-upload.')

            if efvars.include_second_dimension:
                extra_file = os.path.join(efvars.profiles_directory,
                                          efvars.aux_name_array[idx])
                try:
                    extra_data = np.genfromtxt(extra_file, usecols=0,
                                               dtype=float)
                except:
                    log.error('ERROR: Data file ('
                              + ' '+extra_file+') not found.'
                              + ' Please inspect the .zip archive'
                              + ' and re-upload.')
                if len(extra_data) != efvars.num_aux:
                    log.error('ERROR: The number of points in '
                              + extra_file + ' does not match'
                              + ' the number of experimental values'
                              + ' ('+str(efvars.num_aux)+').')
                efvars.full_extra_basis[idx] = extra_data

    def PickleVars(self):
        '''
        Saves self.mvars and self.efvars to serialized files that are
        loaded by the parallel runtime
        '''
        log = self.log
        mvars = self.mvars
        efvars = self.efvars
        log.debug('In PickleVars.')

        self.mvarspickle = os.path.join(
            efvars.output_folder, mvars.runname+'_mvars.p')
        pickle.dump(mvars, open(self.mvarspickle, 'wb'))
        log.debug('mvars have been written to pickle: '+self.mvarspickle)

        self.efvarspickle = os.path.join(
            efvars.output_folder, mvars.runname+'_efvars.p')
        pickle.dump(efvars, open(self.efvarspickle, 'wb'))
        log.debug('efvars have been written to pickle: '+self.efvarspickle)

        return

    def EnsembleFit(self):
        '''
        This initializes the parallel iterative fitting routine
        '''
        mvars = self.mvars
        efvars = self.efvars
        log = self.log

        runcommand = mpiexec + ' -np ' + str(mvars.nproc)\
            + ' ' + pythexec + ' ' + executable\
            + ' -mpick ' + self.mvarspickle\
            + ' -efpick ' + self.efvarspickle

        log.debug('Executable command: '+runcommand)
        p = subprocess.Popen(runcommand, shell=True, executable='/bin/bash')
        sts = os.waitpid(p.pid, 0)[1]
        log.debug('Executable PID = '+str(p.pid))

    def epilogue(self, plotQueues):

        efvars = self.efvars
        mvars = self.mvars
        log = self.log
        pgui = self.run_utils.print_gui

        try:
            auxQueue = plotQueues['bestAUXplot']
            auxResQueue = plotQueues['bestAUXresplot']
        except:
            log.debug('Not using auxiliary plots')

        try:
            pgui('Best model found: \n\n')
            parallel_outf = os.path.join(efvars.output_folder,
                                         mvars.runname+'_final_model.dat')
            sas_outf = os.path.join(efvars.output_folder,
                                    mvars.runname+'_ensemble_sas.int')
            all_model_outf = os.path.join(efvars.output_folder,
                                          mvars.runname+'_all_models.dat')
            plots_outf = os.path.join(efvars.output_folder,
                                      mvars.runname+'_plots.html')
            model_dataf = open(parallel_outf, 'r')
            for line in model_dataf:
                pgui(line)
                pgui('\n\n')
            model_dataf.close()

        except:
            log.error('ERROR: Unable to locate parallel routine output')

        pgui('Best IC model populations saved to '
             + parallel_outf + '\n\n')
        pgui('Best IC model scattering profile saved to '
             + sas_outf + '\n\n')
        if efvars.include_second_dimension:
            aux_outf = os.path.join(efvars.output_folder,
                                    mvars.runname+'_ensemble_aux.dat')
            pgui('Best IC model auxiliary profile saved to '
                 + aux_outf + '\n\n')
        pgui('All sub-ensemble model populations saved to '
             + all_model_outf+'\n\n')
        pgui('Analysis plots saved to '+plots_outf+'\n\n')

        try:
            sas_pickle = os.path.join(efvars.output_folder,
                                      mvars.runname+'_SAS_bokeh.p')
            sas_script = pickle.load(open(sas_pickle, 'rb'))
            res_pickle = os.path.join(efvars.output_folder,
                                      mvars.runname+'_SASres_bokeh.p')
            res_script = pickle.load(open(res_pickle, 'rb'))

            plotQueues['bestSASplot'].put(sas_script)
            plotQueues['bestSASresplot'].put(res_script)

            if efvars.include_second_dimension:
                aux_pickle = os.path.join(efvars.output_folder,
                                          mvars.runname+'_AUX_bokeh.p')
                aux_script = pickle.load(open(aux_pickle, 'rb'))
                aux_res_pickle = os.path.join(efvars.output_folder,
                                              mvars.runname+'_AUXres_bokeh.p')
                aux_res_script = pickle.load(open(aux_res_pickle, 'rb'))

                plotQueues['bestAUXplot'].put(aux_script)
                plotQueues['bestAUXresplot'].put(aux_res_script)

        except:
            log.error('ERROR: Unable to locate bokeh plot pickle objects')

        os.remove(efvars.status_file)

        pgui('STATUS\t1.0\n\n')

        time.sleep(2)
