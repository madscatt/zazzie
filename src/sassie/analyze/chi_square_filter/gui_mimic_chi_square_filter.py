'''
Driver method to run the chi_square filter module
'''

import sys, string

import sassie.analyze.chi_square_filter as chi_square_filter
import sassie.interface.input_filter as input_filter
import sassie.interface.chi_square_filter_filter as chi_square_filter_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0'

runname='run_1'
saspath = './test/run_0/sascalc/neutron_D2Op_100/'
saspath = './diUb/run_0/sascalc/neutron_D2Op_100/'
saspath = './K48_UBA2_corr_1H/sascalc/neutron_D2Op_100/'
sasintfile = './sans_data.dat'
sasintfile = './K48_UBA2.dat'
sasintfile = './K48_UBA2_long.dat'
io = '0.04'
io = '0.1229'
number_of_weight_files = '0'
basis_string  = ''
weight_file_names  = ''
x2highcut = '10.0'
x2highweightfile = 'x2highweights.txt'
x2lowcut = '1.0'
x2lowweightfile = 'x2lowweights.txt'
rghighcut = '60.0'
rghighweightfile = 'rghighweights.txt'
rglowcut = '40.0'
rglowweightfile = 'rglowweights.txt'
sastype = '0'
reduced_x2 = '1'
plotflag = '0'
folder_flag = False

path='./'
data_path = path

#### end user input ####
#### end user input ####
#### end user input ####


svariables['runname'] = (str(runname),'string')

svariables['saspath']           = (str(saspath),'string')
svariables['sasintfile']        = (str(string.strip(sasintfile," ")),'string')
#svariables['sasintfile']       = (str(sasintfile),'string')
svariables['io']                = (str(io),'float')
svariables['number_of_weight_files'] = (str(number_of_weight_files),'int')
svariables['basis_string']      = (str(basis_string),'string')
svariables['weight_file_names'] = (str(weight_file_names),'string')
svariables['x2highcut']         = (str(x2highcut),'float')
svariables['x2highweightfile']  = (str(x2highweightfile),'string')
svariables['x2lowcut']          = (str(x2lowcut),'float')
svariables['x2lowweightfile']   = (str(x2lowweightfile),'string')
svariables['rghighcut']         = (str(rghighcut),'float')
svariables['rghighweightfile']  = (str(rghighweightfile),'string')
svariables['rglowcut']          = (str(rglowcut),'float')
svariables['rglowweightfile']   = (str(rglowweightfile),'string')
svariables['sastype']           = (str(sastype),'int')
svariables['reduced_x2']        = (str(reduced_x2),'int')
svariables['plotflag']          = (str(plotflag),'int')

svariables['path']    = (path,'string')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()

error=chi_square_filter_filter.check_chi_square_filter(variables,no_file_check="true")
if len(error) > 0:
    print 'error = ', error
    sys.exit()

txtQueue = multiprocessing.JoinableQueue()
process=multiprocessing.Process(target=chi_square_filter.find_best,args=(variables,txtQueue))
process.start()


#print 'in GUI and txtOutput = ', this_text, '\n'





