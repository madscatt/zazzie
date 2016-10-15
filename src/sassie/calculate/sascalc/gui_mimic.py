'''
Driver method to run the SasCalc module
'''

import sys

#import sassie.calculate.sascalc as sascalc
import sascalc
import sassie.interface.input_filter as input_filter
import sassie.interface.sascalc_filter as sascalc_filter
import multiprocessing

svariables = {}

#### user input ####
#### user input ####
#### user input ####

runname = 'run_0' #_monica'
#pdbfile = 'new_lysozyme.pdb'
#dcdfile = 'new_lysozyme.dcd'
#dcdfile = 'new_lysozyme.pdb'
pdbfile = 'min3.pdb'
#dcdfile = 'c7.dcd'
dcdfile = 'c7_3.dcd'
#dcdfile = 'c7_2000.dcd'
#pdbfile = 'new_nist_mab.pdb'
#dcdfile = '10000_new_nist_mab.dcd'
#pdbfile = 'small.pdb'
#dcdfile = '4000.dcd'

xon = 'neutron'
#xon = 'xray'
xon = 'neutron_and_xray'

number_contrast_points = "1"
D2O_percentage_array = "100.0" #,0.0"
I0_array = "22.0" #,1.0"
number_exH_regions = "2"
exH_basis_string_array = ["moltype protein","not (moltype protein)"]
fraction_exH_array = "0.0,0"
number_deuteration_regions = "2"
deuterated_basis_string_array = ["not (moltype protein)","moltype protein"]
fraction_deuterated_array = "0.0,0.0"
xray_number_contrast_points = "1"
xray_D2O_percentage_array = "1.0"
xray_I0_array = "8.0"

number_q_values='21'
q_max='0.2'
number_r_values='51'
solvent_volume='20.0'
VDW_scaling_factor='0.77'
#golden_vector_method_option='converge'
golden_vector_method_option='fixed'
number_golden_vectors='35'
golden_vector_method_converge_tolerance='0.01'

pRDF = 'off'
directed_mc = '0'

#### end user input ####
#### end user input ####
#### end user input ####


svariables={}

svariables['runname'] = (runname,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['dcdfile'] = (dcdfile,'string')
svariables['xon'] = (xon,'string')
svariables['number_contrast_points'] = (number_contrast_points,'int')
svariables['D2O_percentage_array'] = (D2O_percentage_array,'float_array')
svariables['I0_array'] = (I0_array,'float_array')
svariables['xray_number_contrast_points'] = (xray_number_contrast_points,'int')
svariables['xray_D2O_percentage_array'] = (xray_D2O_percentage_array,'float_array')
svariables['xray_I0_array'] = (xray_I0_array,'float_array')
if svariables['xon'][0] in ['neutron','neutron_and_xray']:
        svariables['number_exH_regions'] = (number_exH_regions,'int')
        svariables['exH_basis_string_array'] = (exH_basis_string_array,'string')
        svariables['fraction_exH_array'] = (fraction_exH_array,'float_array')
        svariables['number_deuteration_regions'] = (number_deuteration_regions,'int')
        svariables['deuterated_basis_string_array'] = (deuterated_basis_string_array,'string')
        svariables['fraction_deuterated_array'] = (fraction_deuterated_array,'float_array')
            
svariables['number_q_values'] = (number_q_values,'int')
svariables['q_max'] = (q_max,'float')
svariables['number_r_values'] = (number_r_values,'int')
svariables['golden_vector_method_option'] = (golden_vector_method_option,'string')
if svariables['golden_vector_method_option'][0] == 'fixed':
    svariables['number_golden_vectors'] = (number_golden_vectors,'int')
elif svariables['golden_vector_method_option'][0] == 'converge':
    svariables['golden_vector_method_converge_tolerance'] = (golden_vector_method_converge_tolerance,'float')

svariables['solvent_volume'] = (solvent_volume,'float')
svariables['VDW_scaling_factor'] = (VDW_scaling_factor,'float')

if(pRDF == "on"):
    svariables['surf_r'] = (surfr,'float')
    svariables['surf_d'] = (surfd,'float')
    svariables['surf_rou'] = (surfrou,'float')
    svariables['surf_eta'] = (surfeta,'float')
    svariables['surf_N'] = (surfn,'int')
    svariables['surf_coverage'] = (surfcoverage,'float')
    svariables['prdfpath'] = (prdfpath,'string')
    svariables['prdffile'] = (prdffile,'string')
    svariables['cube_length'] = (prdfcubelength,'float')
    svariables['cube_cutoff'] = (prdfcubecutoff,'float')

#variables['directed_mc'] = (directed_mc,'string')
svariables['complex_amplitudes'] = ('0.0','string')

svariables[''] = ('0.0','string')

error, variables = input_filter.type_check_and_convert(svariables)
if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    error=sascalc_filter.check_sascalc(variables)
    if(len(error) != 0):
        print 'error = ',error
        sys.exit()

#import pprint; pprint.pprint(variables); exit()

txtQueue = multiprocessing.JoinableQueue()
sascalc = sascalc.sascalc()
sascalc.main(variables, txtQueue)
this_text = txtQueue.get(True, timeout=0.1)


#print 'in GUI and txtOutput = ', this_text, '\n'


