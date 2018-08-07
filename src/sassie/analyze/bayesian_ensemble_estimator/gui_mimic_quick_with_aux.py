#!/share/apps/local/anacondaz/bin/python

import sys,os
sys.path.append(".")
import numpy as np
import multiprocessing
import sassie.interface.input_filter as input_filter
import bayesian_ensemble_estimator as ensemble_fit
from mpi4py import MPI
svariables = {}

#### MPI Environment ####
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

#### user input ####

#Standard Input
runname = 'UBQtest_quick'
sas_data = './ubq_saxs.dat'
d_max    = '83.6' #If d_max is non-zero, do Shannon sampling, otherwise use standard chi^2
theoretical_profiles_zip = './QuickTest_withAuxData.zip' #'./full_K63.zip' 
number_of_MCs = '1'
max_iterations = '10' #Maximum number of MC steps per walker
posterior_burn = '1' #Number of iterations per walker to remove from beginning (lose initial conformation information
nproc = '2' #Number of processors to run on

#Advanced Input
auxiliary_data = './dist_and_angle.dat' #Experimental data for the extra dimension
use_all = 'False' #Do the entire set, do not conduct iterative search for AIC-identified subset
every = 'False' #Build every possible sub-basis model for comparison instead of stopping at IC-best model
shansamp= 'True' #'True' #Conduct Shannon Sampling
use_bic = 'True'  #Use BIC instead of AIC
walk_one = 'False' #Adjust one population at a time (pre-normalization)
sigma = '0.10' #Gives standard deviation of Gaussian used to pull walk increments.
zeroing_threshold = '0.00' #Convert populations with w < thresh to 0.0.
## end user input ##

#### Parse inputs
svariables['runname'] = (runname, 'string')
svariables['sas_data'] = (sas_data, 'string')
svariables['theoretical_profiles_zip'] = (theoretical_profiles_zip, 'string')
svariables['posterior_burn'] = (posterior_burn, 'int')
svariables['max_iterations'] = (max_iterations, 'int')
svariables['number_of_MCs'] = (number_of_MCs, 'int')
svariables['nproc'] = (nproc, 'int')
svariables['d_max'] = (d_max, 'float')
svariables['auxiliary_data'] = (auxiliary_data, 'string')
svariables['use_all'] = (use_all, 'boolean')
svariables['use_bic'] = (use_bic, 'boolean')
svariables['sigma'] = (sigma, 'float')
svariables['zeroing_threshold'] = (zeroing_threshold, 'float')
svariables['walk_one'] = (walk_one, 'boolean')
svariables['every'] = (every, 'boolean')
svariables['shansamp'] = (shansamp, 'boolean')
####

error,variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
        print("error = ",error)
        sys.exit()

txtQueue = multiprocessing.JoinableQueue()
plotQueues = dict()

#Initialize all the plotting object Queues()
#Best plot
plotQueues['bestSASplot'] = multiprocessing.JoinableQueue()
plotQueues['bestSASresplot'] = multiprocessing.JoinableQueue()
if svariables['auxiliary_data'] != '':
    plotQueues['bestAUXplot'] = multiprocessing.JoinableQueue()
    plotQueues['bestAUXresplot'] = multiprocessing.JoinableQueue()


#run it
reweighting = ensemble_fit.ensemble_routine()
reweighting.main(variables, txtQueue, plotQueues)
