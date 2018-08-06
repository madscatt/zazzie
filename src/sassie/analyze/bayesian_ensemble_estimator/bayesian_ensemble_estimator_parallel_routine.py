'''
    ENSEMBLE_FIT is the module that contains the functions
    that are used to fit an ensemble of structres to 
    experimental scattering data with optional additional data

    REFERENCE:

    Bowerman et al.
    Journal of Chemical Theory and Computation, 13, pp 2418-2429 (2017)

    Bowerman et al.
    Placeholder for Methods Paper...
'''
import numpy as np
from math import ceil
import os
import sys
import zipfile
import itertools
import time
import math
import pickle
import argparse
from scipy.special import binom as binom
from mpi4py import MPI
from bokeh.plotting import figure
from bokeh.plotting import output_file as bokeh_output
from bokeh.plotting import save as bokeh_save
from bokeh.io import curdoc
from bokeh.embed import components
from bokeh.layouts import column, widgetbox, gridplot
from bokeh.models.widgets import Panel, Tabs, DataTable, TableColumn, Select, Div
from bokeh.models import Range1d, FixedTicker, Whisker, CustomJS, Legend
from bokeh.models.sources import ColumnDataSource
from bokeh.models.glyphs import VBar
from bokeh.palettes import d3
from bokeh.core.properties import value

app = 'bayesian_ensemble_estimator'

######## MPI Environment ########
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#################################


class module_variables():
    def __init__(self, parent=None):
        self.app = app


class efvars():
    ''' _E_nsemble _F_it _VAR_iable_S_ '''

    def __init__(self, parent=None):
        pass


# This is the larger object for managing the Ensemble Fit Routine at top-level
'''
    This object builds up all necessary variables and facilitates the forking 
    of each sub-basis proc. After the sub-basis routines are complete, this 
    object also handles likelihood comparisons and printing of the 'best' model
    information.
'''


class ensemble_routine(object):
    def __init__(self, parent=None):
        pass

    def main(self, mpickle='', efpickle=''):
        '''
        main method to handle iterative Bayesian MC routine
        '''
        self.UnpackVariables(mpickle, efpickle)
        self.EnsembleFit()
        self.Epilogue()
        os.system('echo \"Rank '+str(rank) +
                  ' completed all tasks.\" >>'+self.logfile)
        quit()

    def UnpackVariables(self, mpickle, efpickle):
        '''
        extract variables from gui/mimic into system wide class instance
        '''
        self.mvars = pickle.load(open(mpickle, 'rb'))
        self.efvars = pickle.load(open(efpickle, 'rb'))
        mvars = self.mvars
        efvars = self.efvars
        self.logfile = os.path.join(
            efvars.output_folder, mvars.runname+'_runtime.log')

        # Fix Boolean -> string bug...
        try:
            if mvars.use_all.lower() == 'true':
                mvars.use_all = True
            else:
                mvars.use_all = False
        except:
            pass

        try:
            if mvars.every.lower() == 'true':
                mvars.every = True
            else:
                mvars.every = False
        except:
            pass

        try:
            if mvars.use_bic.lower() == 'true':
                mvars.use_bic = True
            else:
                mvars.use_bic = False
        except:
            pass

        try:
            if mvars.walk_one.lower() == 'true':
                mvars.walk_one = True
            else:
                mvars.walk_one = False
        except:
            pass

        try:
            if mvars.shansamp.lower() == 'true':
                mvars.shansamp = True
            else:
                mvars.shansamp = False
        except:
            pass

        return

    def PickleVars(self):
        mvars = self.mvars
        efvars = self.efvars
        if rank == 0:
            os.system('echo \"In PickleVars.\" >> '+self.logfile)

            mvarspickle = os.path.join(
                efvars.output_folder, mvars.runname+"_mvars.p")
            pickle.dump(mvars, open(mvarspickle, 'wb'))
            os.system('echo \"mvars have been written to pickle: ' +
                      mvarspickle+'\" >> '+self.logfile)

            efvarspickle = os.path.join(
                efvars.output_folder, mvars.runname+'_efvars.p')
            pickle.dump(efvars, open(efvarspickle, 'wb'))
            os.system('echo \"efvars have been written to pickle: ' +
                      efvarspickle+'\" >> '+self.logfile)
            for thread in range(1, size):
                comm.send(1, dest=thread)
        else:
            comm.recv(source=0)
        return

    def EnsembleFit(self):
        '''
        This is the actual iterative routine
        '''
        mvars = self.mvars
        efvars = self.efvars
        self.tic = time.time()
        self.min_aic = 99999.9
        self.Individuals()
        self.FindBest()
        efvars.best_sas_chi2 = efvars.subspace_dict[str(
            efvars.best_model)].sas_chi2
        efvars.best_aux_chi2 = efvars.subspace_dict[str(
            efvars.best_model)].aux_chi2
        efvars.best_total_chi2 = efvars.subspace_dict[str(
            efvars.best_model)].total_chi2

        if rank == 0:
            trackf = open(efvars.aic_file, 'a')
            trackf.write(str(1)+'\t'
                         + str(self.min_aic)+'\t'
                         + str(efvars.best_total_chi2)+'\n')
            trackf.close()

        if mvars.use_all == True:
            efvars.best_model = str(efvars.id_array)
            efvars.subspace_dict[efvars.best_model] =\
                simulated_basis(self, efvars.id_array)

            efvars.subspace_dict[efvars.best_model].BayesMC()
            model_object = efvars.subspace_dict[efvars.best_model]
            trackf = open(efvars.aic_file, 'a')
            trackf.write(str(efvars.number_of_profiles)+'\t'
                         + str(model_object.aic)+'\t'
                         + str(model_object.total_chi2)+'\n')
            trackf.close()
        else:
            for subsize in range(2, efvars.number_of_profiles + 1):
                self.current_subsize = subsize
                if rank == 0:
                    self.Progress()
                    sets = list(itertools.combinations(
                        efvars.id_array, subsize))
                    for thread in range(1, size):
                        comm.send(sets, dest=thread)
                else:
                    sets = comm.recv(source=0)
                all_sets_as_list = []
                self.Nsets = len(np.atleast_1d(sets))

                for tracker in range(rank, self.Nsets, size):
                    as_list = np.array([], dtype=int)
                    for submember in sets[tracker]:
                        as_list = np.append(as_list, submember)
                    os.system("echo \"Rank "+str(rank)
                              + ': Working on sub-basis '+str(as_list)
                              + '\" >> '+self.logfile)
                    all_sets_as_list.append(str(as_list))
                    efvars.subspace_dict[str(as_list)] = simulated_basis(
                        self, as_list)
                    efvars.subspace_dict[str(as_list)].BayesMC()

                # The second boolean is to protect from case of sets < threads
                if (rank != 0 and rank < np.min([size, self.Nsets])):
                    os.system('echo \"Rank '+str(rank)
                              + ': Sharing sub-basis information'
                              + ' with head node\" >> '+self.logfile)
                    comm.send(efvars.subspace_dict, dest=0)
                    comm.recv(source=0)
                    os.system('echo \"Rank '+str(rank)
                              + ': Received \'continue\' signal'
                              + ' from head node\" >> '+self.logfile)
                elif (rank >= self.Nsets):
                    os.system('echo \"Rank '+str(rank)
                              + ': Waiting for \'continue\' signal'
                              + ' from head node\" >> '+self.logfile)
                    comm.recv(source=0)
                else:
                    for thread in range(1, np.min([size, self.Nsets])):
                        #dum_subspace_list = comm.recv(source=thread)
                        #Nsubsets = comm.recv(source=thread)
                        # for subset in range(int(Nsubsets)):
                        #    efvars.subspace_dict[dum_subspace_list[subset]] = comm.recv(source=thread)
                        dum_subspace_dict = comm.recv(source=thread)
                        efvars.subspace_dict.update(dum_subspace_dict)
                    for thread in range(1, size):
                        comm.send(1, dest=thread)

                self.FindBest(sublist=all_sets_as_list,
                              subsize=self.current_subsize)
                if self.subset_min_aic > self.min_aic:
                    if rank == 0:
                        os.system('echo \"Best subset of size '+str(subsize)
                                  + ' is a poorer fit than '+str(subsize-1)
                                  + ' according to AIC.\"'
                                  + ' >> '+self.logfile)
                    if not mvars.every:
                        break
                else:
                    if rank == 0:
                        os.system('echo \"Best subset of size '+str(subsize)
                                  + ' is an improvement over '+str(subsize-1)
                                  + ' according to AIC.'
                                  + ' Expanding subset size.\"'
                                  + ' >> '+self.logfile)
                        self.min_aic = self.subset_min_aic
                        efvars.best_model = self.subset_best_model
                        efvars.best_sas_chi2 = efvars.subspace_dict[str(
                            efvars.best_model)].sas_chi2
                        efvars.best_aux_chi2 = efvars.subspace_dict[str(
                            efvars.best_model)].aux_chi2
                        efvars.best_total_chi2 = efvars.subspace_dict[str(
                            efvars.best_model)].total_chi2
                        for thread in range(1, size):
                            comm.send(self.min_aic, dest=thread)
                    else:
                        self.min_aic = comm.recv(source=0)
        if rank == 0:
            if mvars.use_bic:
                ICstring = 'BIC'
            else:
                ICstring = 'AIC'
            os.system('echo \"Best model is: '
                      + str(efvars.best_model)+', '
                      + ICstring+' = '+str(self.min_aic)+', '
                      + 'sas_chi2 = '+str(efvars.best_sas_chi2)+', '
                      + 'aux_chi2 = '+str(efvars.best_aux_chi2)+'\"'
                      + ' >> '+self.logfile)
            with file(efvars.status_file, 'r') as old_status:
                stat_data = old_status.read()
            with file(efvars.status_file, 'w') as new_status:
                new_status.write('STATUS\t0.9999\n' + stat_data)
            #statf = open(efvars.status_file,'w')
            # statf.write('STATUS\t0.9999')
            # statf.close()

        bokehRoutine = Bokeh_and_Save(self)
        bokehRoutine.plotModels()
        bokehRoutine.saveModels()

        return

    def Individuals(self):
        '''
        This will find the fitting capabilities of each individual basis member.
        '''
        mvars = self.mvars
        efvars = self.efvars
        os.system('echo \"Finding Individual Fits\" >> '+self.logfile)
        if rank == 0:
            for idx in range(efvars.number_of_profiles):
                efvars.subspace_dict[str(idx)] = simulated_basis(self, [idx])
                efvars.subspace_dict[str(idx)].BayesMC()
            for thread in range(1, size):
                comm.send(efvars.subspace_dict, dest=thread)
        else:
            efvars.subspace_dict = comm.recv(source=0)

        return

    def FindBest(self, sublist=np.array([], dtype=str), subsize=None):
        '''
        Identifies the best model from the subset, and its relative likelihood to the next best fit.
        '''
        mvars = self.mvars
        efvars = self.efvars
        if rank == 0:
            tic = time.time()
            os.system('echo \"Finding best available model\" >> '
                      + self.logfile)
        if len(sublist) == 0:
            if rank == 0:
                for key in efvars.subspace_dict:
                    key_AIC = efvars.subspace_dict[key].aic
                    if key_AIC < self.min_aic:
                        self.min_aic = key_AIC
                        efvars.best_model = key
                for thread in range(1, size):
                    comm.send(self.min_aic, dest=thread)
                    comm.send(efvars.best_model, dest=thread)
            else:
                if subsize == None:
                    self.min_aic = comm.recv(source=0)
                    efvars.best_model = comm.recv(source=0)
                else:
                    self.subset_min_aic = comm.recv(source=0)
                    self.subset_best_model = comm.recv(source=0)
                    self.subset_best_chi = comm.recv(source=0)
        else:
            self.subset_min_aic = 9999999.9
            for key in sublist:
                key_AIC = efvars.subspace_dict[str(key)].aic
                os.system('echo \"'+str(key_AIC)+','
                          + str(self.subset_min_aic)+'\" >> '
                          + self.logfile)
                if key_AIC < self.subset_min_aic:
                    self.subset_min_aic = key_AIC
                    self.subset_best_model = key
                    self.subset_best_chi = efvars.subspace_dict[key].total_chi2
            if (rank != 0 and rank < np.min([size, self.Nsets])):
                comm.send(self.subset_min_aic, dest=0)
                comm.send(self.subset_best_model, dest=0)
                comm.send(self.subset_best_chi, dest=0)
            if rank == 0:
                for thread in range(1, np.min([size, self.Nsets])):
                    dum_min_aic = comm.recv(source=thread)
                    dum_best_model = comm.recv(source=thread)
                    dum_best_chi = comm.recv(source=thread)
                    if dum_min_aic < self.subset_min_aic:
                        self.subset_min_aic = dum_min_aic
                        self.subset_best_model = dum_best_model
                        self.subset_best_chi = dum_best_chi
                efvars.aic_tracker = np.append(efvars.aic_tracker,
                                               self.subset_min_aic)
                efvars.chi_tracker = np.append(efvars.chi_tracker,
                                               self.subset_best_chi)
                trackf = open(efvars.aic_file, 'a')
                trackf.write(str(self.current_subsize)+'\t'
                             + str(self.subset_min_aic)+'\t'
                             + str(self.subset_best_chi)+'\n')
                trackf.close()
                for thread in range(1, size):
                    comm.send(self.subset_min_aic, dest=thread)
                    comm.send(self.subset_best_model, dest=thread)
                    comm.send(self.subset_best_chi, dest=thread)
                os.system('echo \"Best sub-basis is: '
                          + str(self.subset_best_model)
                          + ', AIC: '+str(self.subset_min_aic)+'\"'
                          + ' >> '+self.logfile)

            else:
                self.subset_min_aic = comm.recv(source=0)
                self.subset_best_model = comm.recv(source=0)
                self.subset_best_chi = comm.recv(source=0)

        return

    def Epilogue(self):
        '''
        Saves files and prints info
        '''
        mvars = self.mvars
        efvars = self.efvars

        if rank == 0:
            print_model = efvars.subspace_dict[efvars.best_model]
            outfile = open(os.path.join(efvars.output_folder,
                                        mvars.runname+'_final_model.dat'), 'w')
            outfile.write(
                '#ScatterFile avg_weight std_weight lone_sas_chi2 lone_aux_chi2 lone_total_chi2\n')
            if print_model.subset_size > 1:
                for idx in range(print_model.subset_size):
                    ##  FIXED?  ##
                    subid = print_model.subset_members[idx]
                    outfile.write(
                        efvars.name_array[print_model.subset_members[idx]]
                        + ' '
                        + str(np.around(print_model.total_mean_weights[idx], decimals=3))
                        + ' '
                        + str(np.around(print_model.weights_std[idx], decimals=3))
                        + ' '
                        + str(np.around(efvars.subspace_dict[str(subid)].sas_chi2, decimals=3))
                        + ' '
                        + str(np.around(efvars.subspace_dict[str(subid)].aux_chi2, decimals=3))
                        + ' '
                        + str(np.around(efvars.subspace_dict[str(subid)].total_chi2, decimals=3))
                        + '\n')
            else:
                outfile.write(
                    efvars.name_array[print_model.subset_members[0]]+'\n')
            outfile.write('SAXS_chi^2  (normalized): '
                          + str(np.around(print_model.sas_chi2, decimals=3))
                          + ' ('+str(np.around(print_model.sas_chi2 /
                                               efvars.num_q, decimals=3))+')'
                          + '\n')
            if efvars.include_second_dimension:
                outfile.write('AUX_chi^2   (normalized): '
                              + str(np.around(print_model.aux_chi2, decimals=3))
                              + ' ('+str(np.around(print_model.aux_chi2 /
                                                   efvars.num_aux, decimals=3))+')'
                              + '\n')
            else:
                outfile.write('AUX_chi^2   (normalized): 0.000 (0.000)\n')
            outfile.write('Total_chi^2 (normalized): '
                          + str(np.around(print_model.total_chi2, decimals=3))
                          + ' ('+str(np.around(print_model.total_chi2 /
                                               efvars.num_points, decimals=3))+')'
                          + '\n')
            if not mvars.use_bic:
                outfile.write('AIC                     : '
                              + str(np.around(print_model.aic, decimals=3))
                              + '\n')
            else:
                outfile.write('BIC                     : '
                              + str(np.around(print_model.aic, decimals=3))
                              + '\n')
            outfile.close()

            SAS_ensemble = efvars.subspace_dict[efvars.best_model].ensemble_I
            ensemblefile = os.path.join(efvars.output_folder,
                                        mvars.runname+'_ensemble_sas.int')
            ofile = open(ensemblefile, 'w')
            for qidx in range(len(efvars.Q)):
                ofile.write(str(efvars.Q[qidx])+'\t'
                            + str(SAS_ensemble[qidx])+'\n')
            ofile.close()

            if efvars.include_second_dimension:
                extra_ensemble = efvars.subspace_dict[efvars.best_model].ensemble_aux
                extrafile = os.path.join(
                    efvars.output_folder, mvars.runname+'_ensemble_aux.dat')
                ofile2 = open(extrafile, 'w')
                for eidx in range(len(efvars.aux_data)):
                    ofile2.write(str(extra_ensemble[eidx])+'\n')
                ofile2.close()

            for thread in range(1, size):
                comm.send(1, dest=thread)

            if not mvars.debug:
                outputfiles = os.listdir(efvars.output_folder)
                for item in outputfiles:
                    if item.endswith('BayesMC.log'):
                        os.remove(os.path.join(efvars.output_folder, item))
                    if item.endswith('.p'):
                        os.remove(os.path.join(efvars.output_folder, item))
            self.toc = time.time()
            runtime = self.toc - self.tic
            os.system("echo \"Bayesian Ensemble Fit routine complete in " +
                      str(runtime/3600.)+" hours.\" >> "+self.logfile)
            self.CheckPoint()
        else:
            comm.recv(source=0)

    def CheckPoint(self):
        '''
        Pickles the current efvars object.
        '''
        mvars = self.mvars
        efvars = self.efvars
        efvarspickle = os.path.join(efvars.output_folder,
                                    mvars.runname+'_efvars_checkpoint.p')
        if os.path.isfile(efvarspickle):
            os.remove(efvarspickle)
        pickle.dump(efvars, open(efvarspickle, 'wb'))

        return

    def Progress(self):
        efvars = self.efvars
        mvars = self.mvars

        try:
            Ncomplete = self.cumul_sum[self.current_subsize-1]
            Npercent = (Ncomplete/self.total_combos)
        except:
            self.cumul_sum = np.zeros(efvars.number_of_profiles, dtype=int)
            for idx in range(1, efvars.number_of_profiles+1):
                count = binom(efvars.number_of_profiles, idx)
                if idx > 1:
                    self.cumul_sum[idx-1] = self.cumul_sum[idx-2] + count
                else:
                    self.cumul_sum[idx-1] = count

            self.total_combos = float(self.cumul_sum[-1]) + 4

            Ncomplete = self.cumul_sum[self.current_subsize-1]
            Npercent = (Ncomplete/self.total_combos)

        # Don't let it count backwards if Nperc too low
        if Npercent < 0.0001:
            Npercent = Npercent + 0.0001
        # Convert 100% to something lower to reflect epilogue processes
        if Npercent == 1.0:
            Npercent = 0.99

        with file(efvars.status_file, 'r') as old_status:
            stat_data = old_status.read()
        with file(efvars.status_file, 'w') as new_status:
            new_status.write('STATUS\t'+str(Npercent)+'\n' + stat_data)
        #statf = open(efvars.status_file, 'w')
        #statf.write('STATUS\t%f' % Npercent)
        # statf.close()

        return


#######################################################################

########################  Subset MonteCarlo object class  ##########################################
''' 
    This object is the one that actually does the Monte Carlo for each subset.
    Relevant metrics (Posterior, Convergence, AIC/BIC, etc.) are stored herein.
'''


class simulated_basis(object):

    def __init__(self, parent, subset_ids):
        self.mvars = parent.mvars
        self.efvars = parent.efvars
        mvars = self.mvars
        efvars = self.efvars

        self.subset_members = subset_ids
        self.as_string = ''
        for ID in self.subset_members:
            self.as_string = self.as_string+'_'+str(ID)
        self.subset_basis = efvars.full_scattering_basis[subset_ids]
        self.subset_size = len(subset_ids)

        if efvars.include_second_dimension:
            self.subset_extra = efvars.full_extra_basis[subset_ids]
        else:
            self.aux_chi2 = 0.0

    def BayesMC(self):
        '''
        This is the Monte Carlo routine.
        '''
        mvars = self.mvars
        efvars = self.efvars
        self.logfile = os.path.join(
            efvars.output_folder, mvars.runname+'_Rank'+str(rank)+"_BayesMC.log")
        os.system('echo \"Running on sub-basis: ' +
                  str(self.subset_members)+'\" >> '+self.logfile)
        # The posterior array has form [id_1,id_2,...,id_n,SAXSchi^2,Dimension2chi^2,Likelihood]
        if self.subset_size == 1:
            self.ensemble_I = np.copy(self.subset_basis[0])
            if efvars.include_second_dimension:
                self.ensemble_aux = np.copy(self.subset_extra[0])
            self.Chi2()
            self.max_likelihood = self.Likelihood()
            if not mvars.use_bic:
                self.AIC()
            else:
                self.BIC()
        else:
            self.posterior_array = np.zeros((
                mvars.number_of_MCs,
                mvars.max_iterations,
                self.subset_size+3
            ), dtype=float)
            for self.run in range(mvars.number_of_MCs):
                os.system('echo \"Beginning run ' +
                          str(self.run)+'\" >> '+self.logfile)
                self.weights = np.random.random(size=self.subset_size)
                self.weights = self.weights/float(np.sum(self.weights))
                self.ensemble_I = np.dot(self.weights, self.subset_basis)
                if efvars.include_second_dimension:
                    self.ensemble_aux = np.dot(self.weights, self.subset_extra)
                self.Chi2()
                for member in range(self.subset_size):
                    self.posterior_array[self.run, 0,
                                         member] = self.weights[member]
                self.posterior_array[self.run, 0,
                                     self.subset_size] = self.sas_chi2
                self.posterior_array[self.run, 0,
                                     self.subset_size + 1] = self.aux_chi2
                self.likelihood = self.Likelihood()
                self.posterior_array[self.run, 0,
                                     self.subset_size + 2] = self.likelihood

                self.iteration = 1

                while self.iteration < mvars.max_iterations:
                    self.Walk()
                    self.iteration += 1

            self.WeightsFromPosterior()
            if not mvars.use_bic:
                self.AIC()
            else:
                self.BIC()
            if mvars.number_of_MCs > 1:
                self.Convergence()
            self.Epilogue()

    def Walk(self):
        efvars = self.efvars
        mvars = self.mvars

        proceed = False

        self.prev_weights = np.copy(self.weights)
        self.prev_sas_chi2 = np.copy(self.sas_chi2)
        self.prev_likelihood = np.copy(self.likelihood)
        self.prev_ensemble_I = np.copy(self.ensemble_I)
        if efvars.include_second_dimension:
            self.prev_aux_chi2 = np.copy(self.aux_chi2)
            self.prev_ensemble_aux = np.copy(self.ensemble_aux)

        if mvars.walk_one:
            while proceed == False:
                delta = np.random.normal(scale=mvars.sigma)
                increment_member = np.random.randint(self.subset_size,
                                                     size=1)
                if ((self.weights[increment_member] + delta) >= 0.0):
                    self.weights[increment_member] += delta
                    self.weights = self.weights/np.sum(self.weights)
                    proceed = True
        else:
            deltas = np.random.normal(scale=mvars.sigma,
                                      size=self.subset_size)
            self.weights = self.weights + deltas
            self.weights = self.weights/np.sum(self.weights)
        zeroers = np.where(self.weights < mvars.zeroing_threshold)[0]
        self.weights[zeroers] = 0.0
        if np.sum(self.weights) == 0:
            os.system('echo \"Sum of weights is 0!'
                      + ' Try reducing the \'sigma\' or'
                      + ' \'zeroing_threshold\' parameters\"'
                      + ' >> ' + self.logfile)
        self.weights = self.weights/np.sum(self.weights)

        self.ensemble_I = np.dot(self.weights, self.subset_basis)
        if efvars.include_second_dimension:
            self.ensemble_aux = np.dot(self.weights, self.subset_extra)
        self.Chi2()
        self.likelihood = self.Likelihood()
        if self.prev_likelihood == 0.0:
            self.prev_likelihood = 1e-8
        acceptance_ratio = self.likelihood/self.prev_likelihood
        draw = np.random.uniform(low=0.0, high=1.0)
        if not ((acceptance_ratio >= 1.0) or (draw < acceptance_ratio)):
            self.weights = np.copy(self.prev_weights)
            self.sas_chi2 = np.copy(self.prev_sas_chi2)
            self.likelihood = np.copy(self.prev_likelihood)
            self.ensemble_I = np.copy(self.prev_ensemble_I)
            if efvars.include_second_dimension:
                self.aux_chi2 = np.copy(self.prev_aux_chi2)
                self.ensemble_aux = np.copy(self.prev_ensemble_aux)
        for member in range(self.subset_size):
            self.posterior_array[self.run, self.iteration,
                                 member] = self.weights[member]
        self.posterior_array[self.run, self.iteration,
                             self.subset_size] = self.sas_chi2
        self.posterior_array[self.run, self.iteration,
                             self.subset_size + 1] = self.aux_chi2
        self.posterior_array[self.run, self.iteration,
                             self.subset_size + 2] = self.likelihood

    def Chi2(self):
        mvars = self.mvars
        efvars = self.efvars
        if mvars.shansamp:
            self.ensemble_I_shan = self.ensemble_I[efvars.shannon_samples]
            sum1 = np.sum(np.divide(np.multiply(
                efvars.samples_I, self.ensemble_I_shan), np.power(efvars.samples_ERR, 2)), axis=0)
            sum2 = np.sum(
                np.power(np.divide(self.ensemble_I_shan, efvars.samples_ERR), 2), axis=0)
            c_array = np.reshape(np.divide(sum1, sum2), (1, 3001))
            self.ensemble_I_shan = np.multiply(c_array, self.ensemble_I_shan)
            chi_array = np.sum(np.power(np.divide(np.subtract(
                efvars.samples_I, self.ensemble_I_shan), efvars.samples_ERR), 2), axis=0)
            chi = np.median(chi_array)
            chi_idx = np.where(chi_array == chi)[0][0]
            c = c_array[0, chi_idx]
            self.ensemble_I = c*self.ensemble_I
            self.sas_chi2 = chi
        else:
            sum1 = np.sum(np.divide(np.multiply(efvars.samples_I,
                                                self.ensemble_I), np.power(efvars.samples_ERR, 2)))
            sum2 = np.sum(
                np.power(np.divide(self.ensemble_I, efvars.samples_ERR), 2))
            c = np.divide(sum1, sum2)
            self.ensemble_I = c*self.ensemble_I
            chi = np.sum(np.power(np.divide(np.subtract(
                efvars.samples_I, self.ensemble_I), efvars.samples_ERR), 2))
            self.sas_chi2 = chi

        if efvars.include_second_dimension:
            chi2 = np.sum(np.power(np.divide(np.subtract(
                efvars.aux_data, self.ensemble_aux), efvars.aux_error), 2))
            self.aux_chi2 = chi2
            self.total_chi2 = self.sas_chi2 + self.aux_chi2
        else:
            self.total_chi2 = self.sas_chi2

    def Likelihood(self):

        likelihood = math.exp(-self.total_chi2/2.0)

        return likelihood

    def AIC(self):

        if self.subset_size > 1:
            self.maxlike = np.max(self.posterior_array[:, :, -1])
            self.aic = 2*(self.subset_size-1) \
                - 2 * np.log(self.maxlike)
        else:
            self.aic = self.total_chi2
        return

    def BIC(self):
        efvars = self.efvars
        if self.subset_size > 1:
            self.maxlike = np.max(self.posterior_array[:, :, -1])
            self.aic = ((self.subset_size - 1) * np.log(efvars.num_points))\
                - 2 * np.log(self.maxlike)
        else:
            self.aic = self.total_chi2
        return

    def WeightsFromPosterior(self):
        mvars = self.mvars
        efvars = self.efvars
        log = open(self.logfile, 'a')

        log.write('Posterior Array:\n'+str(self.posterior_array)+'\n')
        self.mean_weights = np.average(
            self.posterior_array[:, mvars.posterior_burn:, :self.subset_size], axis=1)
        self.mean_weights = self.mean_weights / \
            np.sum(self.mean_weights, axis=1, keepdims=True)
        log.write('mean_weights: '+str(self.mean_weights)+'\n')
        self.total_mean_weights = np.average(self.mean_weights, axis=0)
        log.write('total_mean_weights: '+str(self.total_mean_weights)+'\n')
        self.weights_std = np.zeros(
            np.shape(self.total_mean_weights), dtype=float)
        for idx in range(len(self.weights_std)):
            self.weights_std[idx] = np.std(
                self.posterior_array[:, mvars.posterior_burn:, idx])
        self.ensemble_I = np.dot(self.total_mean_weights, self.subset_basis)
        if efvars.include_second_dimension:
            self.ensemble_aux = np.dot(
                self.total_mean_weights, self.subset_extra)
        self.Chi2()
        # if mvars.print_maxlike:
        #    burnt  = self.posterior_array[:,mvars.posterior_burn:,:]
        #    maxlike= np.max(burnt[:,:,-1])
        #    maxid  = np.where(burnt == maxlike)
        #    minchi = burnt[maxid[0][0],maxid[1][0],-3]
        #    maxlike_weights = burnt[maxid[0][0],maxid[1][0],:self.subset_size]
        #    log.write('Weights of Maximum Observed Likelihood: (Chi^2 = '\
        #              +str(minchi)+')\n')
        #    log.write(str(maxlike_weights)+'\n')
        log.close()

    def Epilogue(self):
        mvars = self.mvars
        efvars = self.efvars
        del self.subset_basis
        if mvars.d_max > 0.0:
            del self.ensemble_I_shan
        # Tracking the full posteriors of all the subsets is very
        # memory-intensive, so we have to delete them as we go.
        del self.posterior_array
        if efvars.include_second_dimension:
            del self.subset_extra
        del self.efvars
        del self.mvars
        return

    def Convergence(self):
        '''Gelman-Rubin Convergence Diagnostic'''
        mvars = self.mvars
        efvars = self.efvars
        log = open(self.logfile, 'a')
        log.write('Calculating Convergence Criterion.\n')

        N = float(mvars.max_iterations - mvars.posterior_burn)

        self.within_variance = np.average(np.var(
            self.posterior_array[:, mvars.posterior_burn:, :self.subset_size], axis=1), axis=0)

        self.between_variance = (N/(mvars.number_of_MCs - 1))*np.sum(
            np.power(np.subtract(self.mean_weights, self.total_mean_weights), 2), axis=0)

        self.pooled_var = ((N-1)/N)*self.within_variance + (
            (mvars.number_of_MCs + 1)/(N*mvars.number_of_MCs))*self.between_variance

        self.PSRF = np.divide(self.pooled_var, self.within_variance)
        self.Rc = np.sqrt(
            ((self.subset_size + 2.0)/self.subset_size)*self.PSRF)

        log.write('PSRF: '+str(self.PSRF)+'\n')
        log.write('Rc: '+str(self.Rc)+'\n')
        log.close()
        del self.PSRF
        del self.Rc
        return
########################################################################################

##### Bokeh Plotting Object  ###########################################################


class Bokeh_and_Save(object):

    def __init__(self, parent):
        self.mvars = parent.mvars
        self.efvars = parent.efvars
        self.logfile = parent.logfile
        self.DataStructures()

    def plotModels(self):
        mvars = self.mvars
        efvars = self.efvars
        plotFile = os.path.join(efvars.output_folder, mvars.runname
                                + '_plots.html')
        if rank == 0:
            bokeh_output(plotFile)
            self.ensembleTab = self.PrepBestModelTab()
            self.top10Tab = self.PrepTop10ModelsTab()
            self.compareTab = self.PrepCompareAllTab()
            AllTabsTogether = Tabs(tabs=[self.ensembleTab,
                                         self.top10Tab,
                                         self.compareTab])
            bokeh_save(AllTabsTogether)
            CurrentBokehDoc = curdoc()
            try:
                script, div = components(CurrentBokehDoc)
            except:
                os.system('echo "Unable to create Bokeh components"'
                          + ' >> '+self.logfile)

            bokeh_pickle = os.path.join(efvars.output_folder,
                                        mvars.runname+'_bokeh.p')
            pickle.dump(script, open(bokeh_pickle, 'wb'))

            for thread in range(1, size):
                comm.send(1, dest=thread)
        else:
            comm.recv(source=0)

        return

    def DataStructures(self):
        mvars = self.mvars
        efvars = self.efvars
        names = efvars.name_array

        ICarray = np.array([], dtype=float)
        ChiArray = np.array([], dtype=float)
        SASchiArray = np.array([], dtype=float)
        AUXchiArray = np.array([], dtype=float)
        SubsizeArray = np.array([], dtype=int)
        MemberArray = np.array([], dtype=str)
        PopDict = dict()

        for name in names:
            PopDict[name] = np.array([], dtype=float)
            PopDict[name+'STD'] = np.array([], dtype=float)

        if rank == 0:
            keys = [str(key) for key in efvars.subspace_dict]
            mykeys = keys[0::size]
            for thread in range(1, size):
                comm.send(efvars.subspace_dict, dest=thread)
                comm.send(keys[thread::size], dest=thread)
        else:
            efvars.subspace_dict = comm.recv(source=0)
            mykeys = comm.recv(source=0)
        for key in mykeys:
            model = efvars.subspace_dict[key]
            ICarray = np.append(ICarray, model.aic)
            ChiArray = np.append(ChiArray,
                                 model.total_chi2/efvars.num_points)
            SASchiArray = np.append(SASchiArray,
                                    model.sas_chi2/efvars.num_q)
            if efvars.include_second_dimension:
                AUXchiArray = np.append(AUXchiArray,
                                        model.aux_chi2/efvars.num_aux)
            else:
                AUXchiArray = np.append(AUXchiArray, 0.000)

            SubsizeArray = np.append(SubsizeArray, model.subset_size)
            try:
                SASProf = np.vstack((SASProf, model.ensemble_I))
                if efvars.include_second_dimension:
                    AUXProf = np.vstack((AUXProf, model.ensemble_aux))
            except:
                SASProf = np.reshape(model.ensemble_I,
                                     (1, len(model.ensemble_I)))
                if efvars.include_second_dimension:
                    AUXProf = np.reshape(model.ensemble_aux,
                                         (1, len(model.ensemble_aux)))
            nameChecklist = np.copy(names)
            if model.subset_size > 1:
                MemberString = ''
                for member in range(model.subset_size):
                    profile = names[model.subset_members[member]]
                    MemberString += profile
                    if (member+1) != model.subset_size:
                        MemberString += ','

                    try:
                        PopDict[profile] = np.append(PopDict[profile],
                                                     model.total_mean_weights[member])
                        PopDict[profile+'STD'] = np.append(PopDict[profile+'STD'],
                                                           model.weights_std[member])
                    except:
                        PopDict[profile] = model.total_mean_weights[member]
                        PopDict[profile+'STD'] = model.weights_std[member]
                    whichProfile = np.where(nameChecklist == profile)[0][0]
                    nameChecklist = np.delete(nameChecklist,
                                              np.where(
                                                  nameChecklist == profile)[0][0]
                                              )
                MemberArray = np.append(MemberArray, MemberString)
            else:
                profile = names[model.subset_members[0]]
                MemberArray = np.append(MemberArray, profile)
                try:
                    PopDict[profile] = np.append(PopDict[profile], 1.000)
                    PopDict[profile +
                            'STD'] = np.append(PopDict[profile+'STD'], 0.000)
                except:
                    PopDict[profile] = 1.0
                    PopDict[profile+'STD'] = 0.0
                nameChecklist = np.delete(nameChecklist,
                                          np.where(
                                              nameChecklist == profile)[0][0]
                                          )
            for nonmember in nameChecklist:
                try:
                    PopDict[nonmember] = np.append(PopDict[nonmember], 0.000)
                    PopDict[nonmember +
                            'STD'] = np.append(PopDict[nonmember+'STD'], 0.000)
                except:
                    PopDict[nonmember] = 0.00
                    PopDict[nonmember+'STD'] = 0.00

        if rank == 0:
            for thread in range(1, size):
                RecvDict = comm.recv(source=thread)
                ICarray = np.append(ICarray,
                                    RecvDict['IC'])
                ChiArray = np.append(ChiArray,
                                     RecvDict['Chi'])
                SASchiArray = np.append(SASchiArray,
                                        RecvDict['SASChi'])
                AUXchiArray = np.append(AUXchiArray,
                                        RecvDict['AUXChi'])
                SubsizeArray = np.append(SubsizeArray,
                                         RecvDict['Subsize'])
                MemberArray = np.append(MemberArray,
                                        RecvDict['Members'])
                SASProf = np.vstack((SASProf,
                                     RecvDict['SASProf']))
                if efvars.include_second_dimension:
                    AUXProf = np.vstack((AUXProf,
                                         RecvDict['AUXProf']))
                for name in names:
                    PopDict[name] = np.append(PopDict[name],
                                              RecvDict['Pops'][name])
                    PopDict[name+'STD'] = np.append(PopDict[name+'STD'],
                                                    RecvDict['Pops'][name+'STD'])

            sort_idx = np.argsort(ICarray)
            self.ICarray = ICarray[sort_idx]
            self.ChiArray = ChiArray[sort_idx]
            self.SASchiArray = SASchiArray[sort_idx]
            self.AUXchiArray = AUXchiArray[sort_idx]
            self.SubsizeArray = SubsizeArray[sort_idx]
            self.MemberArray = MemberArray[sort_idx]
            self.SASProf = SASProf[sort_idx]
            if efvars.include_second_dimension:
                self.AUXProf = AUXProf[sort_idx]
            self.PopDict = PopDict
            for ID in PopDict:
                self.PopDict[ID] = self.PopDict[ID][sort_idx]

            minIC = self.ICarray[0]
            self.RelPerform = np.exp((minIC-self.ICarray)/2.0)

            # Round values for printing
            self.RelPerform = np.around(self.RelPerform, decimals=2)
            self.ICarray = np.around(self.ICarray, decimals=2)
            self.ChiArray = np.around(self.ChiArray, decimals=2)
            self.SASchiArray = np.around(self.SASchiArray, decimals=2)
            self.AUXchiArray = np.around(self.AUXchiArray, decimals=2)
            for name in names:
                self.PopDict[name] = np.around(self.PopDict[name],
                                               decimals=2)
                self.PopDict[name+'STD'] = np.around(self.PopDict[name+'STD'],
                                                     decimals=2)
                self.PopDict[name+'_hiErr'] = self.PopDict[name]\
                    + self.PopDict[name+'STD']
                self.PopDict[name+'_loErr'] = self.PopDict[name]\
                    - self.PopDict[name+'STD']

            self.MakeColumnDataSources()

        else:
            SendDict = dict()
            SendDict['IC'] = ICarray
            SendDict['Chi'] = ChiArray
            SendDict['SASChi'] = SASchiArray
            SendDict['AUXChi'] = AUXchiArray
            SendDict['Subsize'] = SubsizeArray
            SendDict['Pops'] = PopDict
            SendDict['Members'] = MemberArray
            SendDict['SASProf'] = SASProf
            if efvars.include_second_dimension:
                SendDict['AUXProf'] = AUXProf
            comm.send(SendDict, dest=0)

        return

    def MakeColumnDataSources(self):
        mvars = self.mvars
        efvars = self.efvars
        names = efvars.name_array
        if efvars.include_second_dimension:
            auxnames = efvars.aux_name_array
        bestmodel = efvars.subspace_dict[efvars.best_model]

        # For "Best Model" tab
        bestdat = dict()
        bestaux = dict()

        # For "Top 10 Models" tab
        top10dat = dict()
        top10plots = dict()
        top10aux = dict()

        # For "Compare All Models" tab
        PerfHists = dict()
        fullTable = dict()

        ##### Build Best Model Column Data Sources #######################
        bestdat['q'] = efvars.Q
        bestdat['i'] = efvars.I
        bestdat['err'] = efvars.ERR
        bestdat['hiErr'] = efvars.I + efvars.ERR
        bestdat['loErr'] = efvars.I - efvars.ERR
        bestdat['model'] = bestmodel.ensemble_I
        bestdat['resid'] = bestmodel.ensemble_I - efvars.I
        if efvars.include_second_dimension:
            bestaux['x'] = np.arange(efvars.num_aux)
            bestaux['y'] = efvars.aux_data
            bestaux['hiErr'] = efvars.aux_data + efvars.aux_error
            bestaux['loErr'] = efvars.aux_data - efvars.aux_error
            bestaux['model'] = bestmodel.ensemble_aux
            bestaux['resid'] = bestmodel.ensemble_aux - efvars.aux_data

        self.bestdatCDS = ColumnDataSource(data=bestdat)
        self.bestauxCDS = ColumnDataSource(data=bestaux)
        ##################################################################

        ##### Build Top 10 Tab Column Data Sources #######################
        top10dat['Subsize'] = self.SubsizeArray[:10]
        top10dat['IC'] = self.ICarray[:10]
        top10dat['RelPerf'] = self.RelPerform[:10]
        top10dat['Chi'] = self.ChiArray[:10]
        top10dat['ChiSAS'] = self.SASchiArray[:10]
        top10dat['ChiAUX'] = self.AUXchiArray[:10]
        top10dat['Members'] = self.MemberArray[:10]

        top10plots['q'] = efvars.Q
        top10plots['i'] = efvars.I
        top10plots['hiErr'] = efvars.I + efvars.ERR
        top10plots['loErr'] = efvars.I - efvars.ERR
        if efvars.include_second_dimension:
            top10aux['x'] = np.arange(efvars.num_aux)
            top10aux['y'] = efvars.aux_data
            top10aux['hiErr'] = efvars.aux_data + efvars.aux_error
            top10aux['loErr'] = efvars.aux_data - efvars.aux_error

        for idx in range(len(top10dat['IC'])):  # Don't assume 10 models
            top10plots['model'+str(idx)] = np.copy(self.SASProf[idx])
            top10plots['model'+str(idx)+'res'] = self.SASProf[idx]\
                - efvars.I
            if efvars.include_second_dimension:
                top10aux['model'+str(idx)] = np.copy(self.AUXProf[idx])
                top10aux['model'+str(idx)+'res'] = self.AUXProf[idx]\
                    - efvars.aux_data

        for idx in range(len(names)):  # Plots include individual members
            oneSAS = efvars.full_scattering_basis[idx]
            top10plots[names[idx]] = oneSAS
            top10plots[names[idx]+'res'] = oneSAS - efvars.I
            if efvars.include_second_dimension:
                oneAUX = efvars.full_extra_basis[idx]
                top10aux[auxnames[idx]] = oneAUX
                top10aux[auxnames[idx]+'res'] = oneAUX - efvars.aux_data

        self.top10datCDS = ColumnDataSource(data=top10dat)
        self.top10plotCDS = ColumnDataSource(data=top10plots)
        if efvars.include_second_dimension:
            self.top10auxCDS = ColumnDataSource(data=top10aux)
        ##################################################################

        ##### Build Compare All Models Tab Column Data Sources ###########
        # Relative Performance Bar Plot
        maxsub = np.max(self.SubsizeArray)
        self.basis_sizes = []
        Nmodels = len(self.SubsizeArray)
        PerfBins = np.arange(0, 1.01, 0.05)
        for subsize in range(1, maxsub+1):
            subMembers = np.where(self.SubsizeArray == subsize)[0]
            subPerfs = self.RelPerform[subMembers]
            try:
                PerfHists[str(subsize)] = np.histogram(
                    subPerfs, bins=PerfBins)[0]
                self.basis_sizes.append(str(subsize))
            except:  # If "use_all" option (only size 1 and maxsub)
                continue

        self.RelPerfCDS = ColumnDataSource(data=PerfHists)
        self.RelPerfCDS.data['Bins'] = PerfBins[:-1] + 0.025

        # Populations and model info for table and pop. bar plot
        fullTable['model'] = np.arange(len(self.ICarray))
        fullTable['IC'] = self.ICarray
        fullTable['Chi'] = self.ChiArray
        fullTable['SASchi'] = self.SASchiArray
        fullTable['AUXchi'] = self.AUXchiArray
        fullTable['Subsize'] = self.SubsizeArray
        fullTable['Members'] = self.MemberArray
        fullTable['RelPerf'] = self.RelPerform
        fullTable.update(self.PopDict)

        self.FullTableCDS = ColumnDataSource(data=fullTable)
        ##################################################################
        return

    def PrepBestModelTab(self):
        efvars = self.efvars
        mvars = self.mvars

        TabTitle = 'Best Model'
        SAStitle = 'Scattering Profile, Best Model vs Experiment'\
                   + ' (Chi^2 = '+str(self.SASchiArray[0])+')'
        SASxLabel = 'q'
        SASyLabel = 'I(q)'
        REStitle = 'Best Model, Scattering Residuals'
        if efvars.include_second_dimension:
            AUXtitle = 'Auxiliary Profile, Best Model vs Experiment'\
                + ' (Chi^2 = '+str(self.AUXchiArray[0])+')'
            AUXxLabel = 'Auxiliary Measurement'
            AUXyLabel = 'Measurement Value'
            AUXresTitle = 'Best Model, Auxiliary Residuals'

        SASplot = figure(title=SAStitle, x_axis_label=SASxLabel,
                         y_axis_label=SASyLabel, y_axis_type='log')
        resPlot = figure(title=REStitle, x_axis_label=SASxLabel,
                         y_axis_label=SASyLabel, plot_height=300)
        resPlot.x_range = SASplot.x_range

        SASplot.add_layout(Whisker(source=self.bestdatCDS, base='q',
                                   lower='loErr', upper='hiErr'))
        SASplot.circle(source=self.bestdatCDS, x='q', y='i', color='black',
                       legend='Experiment')
        SASplot.line(source=self.bestdatCDS, x='q', y='model',
                     line_width=4.0, legend='Best Model')

        resPlot.line(source=self.bestdatCDS, x='q', y=0, line_dash='dashed',
                     line_color='black')
        resPlot.line(source=self.bestdatCDS, x='q', y='resid')

        if efvars.include_second_dimension:
            AUXplot = figure(title=AUXtitle, x_axis_label=AUXxLabel,
                             y_axis_label=AUXyLabel)
            AUXresPlot = figure(title=AUXresTitle, x_axis_label=AUXxLabel,
                                y_axis_label=AUXyLabel, height=300)

            AUXresPlot.x_range = AUXplot.x_range
            AUXplot.xaxis.ticker = np.arange(efvars.num_aux)
            AUXresPlot.xaxis.ticker = np.arange(efvars.num_aux)

            AUXlowY = 0.9*np.min([np.min(self.bestauxCDS.data['loErr']),
                                  np.min(self.bestauxCDS.data['model'])])
            AUXhighY = 1.1*np.max([np.max(self.bestauxCDS.data['hiErr']),
                                   np.min(self.bestauxCDS.data['model'])])
            AUXplot.y_range = Range1d(AUXlowY, AUXhighY)

            AUXplot.add_layout(Whisker(source=self.bestauxCDS, base='x',
                                       lower='loErr', upper='hiErr'))
            AUXplot.circle(source=self.bestauxCDS, x='x', y='y', size=12,
                           color='black', legend='Experiment')
            AUXplot.circle(source=self.bestauxCDS, x='x', y='model',
                           size=12, legend='Best Model')

            AUXresPlot.line(source=self.bestauxCDS, x='x', y=0,
                            line_dash='dashed', line_color='black')
            AUXresPlot.circle(source=self.bestauxCDS, x='x', y='resid',
                              size=12)

            plotSpace = gridplot([[SASplot, AUXplot],
                                  [resPlot, AUXresPlot]])

        else:
            plotSpace = column([SASplot, resPlot])
        bestModelTab = Panel(child=plotSpace, title=TabTitle)
        return bestModelTab

    def PrepTop10ModelsTab(self):
        mvars = self.mvars
        efvars = self.efvars
        names = efvars.name_array
        if efvars.include_second_dimension:
            auxnames = efvars.aux_name_array

        TabTitle = 'Top 10 Models'
        ############  The data table  ####################################
        tableColumns = [
            TableColumn(field='RelPerf', sortable=True,
                        default_sort='descending',
                        title='Relative Model Performance'),
            TableColumn(field='Subsize', sortable=True,
                        title='Ensemble Size'),
            TableColumn(field='Members', title='Ensemble Members'),
            TableColumn(field='Chi', title='Chi^2')]
        if efvars.include_second_dimension:
            tableColumns.append(TableColumn(field='ChiSAS',
                                            title='SAS Chi^2'))
            tableColumns.append(TableColumn(field='ChiAUX',
                                            title='AUX Chi^2'))
        top10Table = DataTable(source=self.top10datCDS,
                               columns=tableColumns,
                               sortable=True)
        ##################################################################

        ############ The Spectra Plots  ##################################
        palette = 'Category10'
        plotLocation = 'top_left'
        SAStitle = 'Scattering Profile, Selected Model vs Experiment'
        SASxLabel = 'q'
        SASyLabel = 'I(q)'
        resTitle = 'Selected Model, Scattering Residuals'

        SASlegendItems = []
        ResLegendItems = []

        if efvars.include_second_dimension:
            AUXtitle = 'Auxiliary Profile, Selected Model vs Experiment'
            AUXxLabel = 'Auxiliary Measurement'
            AUXyLabel = 'Measurement Value'
            AUXresTitle = 'Selected Model, Auxiliary Residuals'

            AUXlegendItems = []
            AUXresLegendItems = []

            AUXplot = figure(title=AUXtitle, x_axis_label=AUXxLabel,
                             y_axis_label=AUXyLabel)

            AUXresPlot = figure(title=AUXresTitle, x_axis_label=AUXxLabel,
                                y_axis_label=AUXyLabel)
            AUXplot.add_layout(Whisker(source=self.top10auxCDS, base='x',
                                       upper='hiErr', lower='loErr'))
            AUXplot.circle(x='x', y='y', size=12, color='black',
                           source=self.top10auxCDS)
            AUXresPlot.line(x='x', y=0, line_dash='dashed',
                            line_color='black', source=self.top10auxCDS)
            AUXresPlot.x_range = AUXplot.x_range
            AUXplot.xaxis.ticker = np.arange(efvars.num_aux)
            AUXresPlot.xaxis.ticker = np.arange(efvars.num_aux)

            NumModels = len(self.top10datCDS.data['Chi'])
            minIndivid = np.min([self.top10auxCDS.data[key] for key in
                                 auxnames])
            maxIndivid = np.max([self.top10auxCDS.data[key] for key in
                                 auxnames])

            minModel = np.min([self.top10auxCDS.data['model'+str(idx)]
                               for idx in range(NumModels)])
            maxModel = np.max([self.top10auxCDS.data['model'+str(idx)]
                               for idx in range(NumModels)])
            AUXlowY = 0.9*np.min([np.min(self.top10auxCDS.data['loErr']),
                                  minModel, minIndivid])
            AUXhighY = 1.1*np.max([np.max(self.top10auxCDS.data['hiErr']),
                                   maxModel, maxIndivid])
            AUXplot.y_range = Range1d(AUXlowY, AUXhighY)

        SASplot = figure(title=SAStitle, x_axis_label=SASxLabel,
                         y_axis_label=SASyLabel, y_axis_type='log')
        resPlot = figure(title=resTitle, x_axis_label=SASxLabel,
                         y_axis_label=SASyLabel)
        resPlot.line(x='q', y=0, line_dash='dashed', line_color='black',
                     source=self.top10plotCDS)
        resPlot.x_range = SASplot.x_range

        SASplot.add_layout(Whisker(source=self.top10plotCDS, base='q',
                                   upper='hiErr', lower='loErr'))

        ######### JS Callback for Plotted Model Selection ############
        numModels = len(self.top10datCDS.data['Chi'])
        dropOptions = ['model'+str(idx) for idx in range(numModels)]
        dropSelect = Select(title='Selected Model:',
                            options=dropOptions, value=dropOptions[0])

        if efvars.include_second_dimension:
            JScall = """
            var SASdat = SAS.data;
            var AUXdat = AUX.data;

            SASdat['Selected Model'] = SASdat[cb_obj.value];
            SASdat['Selected Resid'] = SASdat[cb_obj.value+'res'];
            AUXdat['Selected Model'] = AUXdat[cb_obj.value];
            AUXdat['Selected Resid'] = AUXdat[cb_obj.value+'res'];

            SAS.change.emit();
            AUX.change.emit();"""
            callback = CustomJS(args={'SAS': self.top10plotCDS,
                                      'AUX': self.top10auxCDS},
                                code=JScall)
        else:
            JScall = """
            var SASdat  = SAS.data;

            SASdat['Selected Model'] = SASdat[cb_obj.value];
            SASdat['Selected Resid'] = SASdat[cb_obj.value+'res'];

            SAS.change.emit();
            """
            callback = CustomJS(args={'SAS': self.top10plotCDS},
                                code=JScall)

        dropSelect.callback = callback
        ##############################################################

        # Start by plotting the best/selected model
        startModel = self.top10plotCDS.data['model0']
        startResid = self.top10plotCDS.data['model0res']
        self.top10plotCDS.data['Selected Model'] = startModel
        self.top10plotCDS.data['Selected Resid'] = startResid
        line = SASplot.line(x='q', y='Selected Model', line_width=4.0,
                            source=self.top10plotCDS)
        rline = resPlot.line(x='q', y='Selected Resid',
                               source=self.top10plotCDS)
        SASlegendItems.append(('Selected Model', [line]))
        ResLegendItems.append(('Selected Model', [rline]))

        if efvars.include_second_dimension:
            AUXstartModel = self.top10auxCDS.data['model0']
            AUXresStartModel = self.top10auxCDS.data['model0res']
            self.top10auxCDS.data['Selected Model'] = AUXstartModel
            self.top10auxCDS.data['Selected Resid'] = AUXresStartModel
            scat = AUXplot.circle(x='x', y='Selected Model', size=12,
                                  source=self.top10auxCDS)
            rscat = AUXresPlot.circle(x='x', y='Selected Resid', size=12,
                                        source=self.top10auxCDS)
            AUXlegendItems.append(('Selected Model', [scat]))
            AUXresLegendItems.append(('Selected Model', [rscat]))

        # Then add all the individuals
        colors = d3['Category10'][10]
        styles = ['dashed', 'dotted', 'dotdash', 'dashdot']
        for idx in range(len(names)):
            line_idx = idx/len(colors)
            key = names[idx]
            line = SASplot.line(x='q', y=key, line_color=colors[idx % 10],
                                line_width=4, line_dash=styles[line_idx],
                                source=self.top10plotCDS)
            rline = resPlot.line(x='q', y=key+'res',
                                   line_color=colors[idx % 10],
                                   line_dash=styles[line_idx],
                                   source=self.top10plotCDS)

            line.visible = False
            rline.visible = False

            SASlegendItems.append((key, [line]))
            ResLegendItems.append((key, [rline]))

            if efvars.include_second_dimension:
                auxkey = auxnames[idx]
                if line_idx == 0:
                    scat = AUXplot.circle(x='x', y=auxkey,
                                          fill_color=colors[idx % 10],
                                          line_color='black',
                                          size=12, line_width=3,
                                          source=self.top10auxCDS)
                    rscat = AUXresPlot.circle(x='x', y=auxkey+'res',
                                                color=colors[idx],
                                                line_color='black',
                                                line_width=3, size=12,
                                                source=self.top10auxCDS)
                elif line_idx == 1:
                    scat = AUXplot.circle(x='x', y=auxkey,
                                          fill_color='black',
                                          line_color=colors[idx % 10],
                                          size=12, line_width=3,
                                          source=self.top10auxCDS)
                    rscat = AUXresPlot.circle(x='x', y=auxkey+'res',
                                                fill_color='black',
                                                line_color=colors[idx],
                                                line_width=3, size=12,
                                                source=self.top10auxCDS)
                else:
                    idx2 = (idx + line_idx + 2) % 10
                    scat = AUXplot.circle(x='x', y=auxkey,
                                          fill_color=colors[idx],
                                          line_color=colors[idx2],
                                          size=12, line_width=3,
                                          source=self.top10auxCDS)
                    rscat = AUXresPlot.circle(x='x', y=auxkey+'res',
                                                fill_color=colors[idx],
                                                line_color=colors[idx2],
                                                line_width=3, size=12,
                                                source=self.top10auxCDS)

                scat.visible = False
                rscat.visible = False

                AUXlegendItems.append((auxkey, [scat]))
                AUXresLegendItems.append((auxkey, [rscat]))

        legend = Legend(items=SASlegendItems)
        legend.inactive_fill_alpha = 0.0
        legend.click_policy = 'hide'
        SASplot.add_layout(legend, 'right')

        rlegend = Legend(items=ResLegendItems)
        rlegend.inactive_fill_alpha = 0.0
        rlegend.click_policy = 'hide'
        resPlot.add_layout(rlegend, 'right')

        if efvars.include_second_dimension:
            legend2 = Legend(items=AUXlegendItems)
            legend2.inactive_fill_alpha = 0.0
            legend2.click_policy = 'hide'
            AUXplot.add_layout(legend2, 'right')

            rlegend2 = Legend(items=AUXresLegendItems)
            rlegend2.inactive_fill_alpha = 0.0
            rlegend2.click_policy = 'hide'
            AUXresPlot.add_layout(rlegend2, 'right')

            plotSpace = gridplot([[widgetbox(top10Table), None],
                                  [dropSelect, None],
                                  [SASplot, AUXplot],
                                  [resPlot, AUXresPlot]],
                                 toolbar_location='right',
                                 merge_tools='True')
        else:
            plotSpace = column([widgetbox(top10Table),
                                dropSelect,
                                SASplot,
                                resPlot])

        top10Tab = Panel(child=plotSpace, title='Top 10 Models')
        return top10Tab

    def PrepCompareAllTab(self):
        mvars = self.mvars
        efvars = self.efvars
        names = efvars.name_array
        if efvars.include_second_dimension:
            auxnames = efvars.aux_name_array

        Nprof = len(names)
        tableTitle = """<b>Sub-ensemble Model Details</b>"""
        Tools = 'ypan,ywheel_zoom,save,reset'
        HistTitle = """<b>Relative Model Performance Histogram,\
                         Sorted by Model Ensemble Size</b>"""
        PopTitle = """<b>Member Populations, Selected Model</b>"""

        if mvars.use_bic:
            IClabel = 'BIC'
        else:
            IClabel = 'AIC'

        ##### Rel Peform Histogram ###############################
        histPlot = figure(plot_height=800, plot_width=800, tools=Tools,
                          sizing_mode='scale_width',
                          x_range=(0, 1.0), y_range=(0, 30.0),
                          x_axis_label='Relative Model Performance',
                          y_axis_label='Number of Models')
        colors = d3['Category20c'][len(self.basis_sizes)]
        legendItems = [value(x) for x in self.basis_sizes]
        renderer = histPlot.vbar_stack(self.basis_sizes,
                                       x='Bins', width=0.05,
                                       name=self.basis_sizes,
                                       legend=legendItems,
                                       color=colors, line_width=0.0,
                                       source=self.RelPerfCDS)
        histPlot.legend.click_policy = 'hide'
        ##########################################################

        ##### Full Data Table ####################################
        tableColumns = [
            TableColumn(field='RelPerf',
                        title='Relative Model Performance'),
            TableColumn(field='IC', title=IClabel),
            TableColumn(field='Subsize', title='Ensemble Size'),
            TableColumn(field='Chi', title='Model Chi^2'),
        ]
        if efvars.include_second_dimension:
            tableColumns.append(TableColumn(field='SASchi',
                                            title='SAS Chi^2'))
            tableColumns.append(TableColumn(field='AUXchi',
                                            title='AUX Chi^2'))
        fullTable = DataTable(source=self.FullTableCDS,
                              columns=tableColumns,
                              sortable=True)
        # Initialize the Selection to the Best Model
        self.FullTableCDS.selected['1d'].indices = [0]
        ##########################################################

        ##### Population Bars ####################################
        popFig = figure(plot_height=400, plot_width=800,
                        sizing_mode='scale_width', tools=Tools,
                        y_range=(0, 1.0), x_axis_label='Ensemble Member',
                        y_axis_label='Member Populations')

        tickdict = {}
        popFig.xaxis.ticker = FixedTicker(ticks=list(np.arange(Nprof)))
        for idx in range(Nprof):
            tickdict[str(idx)] = names[idx]
        popFig.xaxis.major_label_overrides = tickdict
        popFig.xaxis.major_label_orientation = math.pi/4

        nonselect = VBar(fill_alpha=0.0, line_alpha=0.0)
        for nameIDX in range(Nprof):
            profName = names[nameIDX]
            popglyph = VBar(x=nameIDX, top=profName, width=0.5,
                            line_width=0.0, line_alpha=0.0)
            eglyph = VBar(x=nameIDX, top=profName+'_hiErr',
                          bottom=profName+'_loErr', width=0.1,
                          fill_color='black', line_color='black')
            renderer = popFig.add_glyph(self.FullTableCDS, popglyph)
            erenderer = popFig.add_glyph(self.FullTableCDS, eglyph)
            renderer.nonselection_glyph = nonselect
            erenderer.nonselection_glyph = nonselect
        ##########################################################

        plotSpace = column(Div(text=HistTitle), histPlot,
                           widgetbox(Div(text=tableTitle),
                                     fullTable,
                                     Div(text=PopTitle)),
                           popFig)
        compareAllTab = Panel(child=plotSpace,
                              title='Compare All Models')

        return compareAllTab

    def saveModels(self):
        efvars = self.efvars
        mvars = self.mvars
        names = efvars.name_array

        if rank == 0:
            allfile = open(os.path.join(
                           efvars.output_folder,
                           mvars.runname+'_all_models.dat'),
                           'w')
            allfile.write('#Relative Performance,')
            if mvars.use_bic:
                allfile.write('BIC,')
            else:
                allfile.write('AIC,')
            if efvars.include_second_dimension:
                allfile.write('model chi^2,sas chi^2,aux chi^2,')
            else:
                allfile.write('model chi^2,')
            for idx in range(len(names)):
                name = names[idx]
                allfile.write(name+' weight,'+name+' stdev,')
            allfile.write('\n')

            for model in range(len(self.FullTableCDS.data['IC'])):
                modPerf = self.FullTableCDS.data['RelPerf'][model]
                modIC = self.FullTableCDS.data['IC'][model]
                modChi = self.FullTableCDS.data['Chi'][model]
                allfile.write(str(modPerf)+','+str(modIC)+','+str(modChi)+',')
                if efvars.include_second_dimension:
                    modSAS = self.FullTableCDS.data['SASchi'][model]
                    modAUX = self.FullTableCDS.data['AUXchi'][model]
                    allfile.write(str(modSAS)+','+str(modAUX)+',')

                for idx in range(len(names)):
                    name = names[idx]
                    pop = self.FullTableCDS.data[name][model]
                    std = self.FullTableCDS.data[name+'STD'][model]
                    allfile.write(str(pop)+','+str(std)+',')
                allfile.write('\n')


##### Command Line Parsing #####
'''
    What follows below is the command-line parsing from the main
    bayesian_ensemble_estimator module.  It builds and runs the
    above classes.
'''
parser = argparse.ArgumentParser()
parser.add_argument("-mpick", dest='mpick', type=str,
                    help='Pickled module_variables() object.')
parser.add_argument("-efpick", dest='efpick', type=str,
                    help='Pickled ensemble_fitting_variables() object.')
args = parser.parse_args()

reweighting = ensemble_routine()
reweighting.main(mpickle=args.mpick, efpickle=args.efpick)
