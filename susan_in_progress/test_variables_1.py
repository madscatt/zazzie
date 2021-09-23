# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

def process_input_variables(stoichiometry_variables):

    print('stoichiometry_variables: %s' %(str(stoichiometry_variables))) 

    fraction_d2o = stoichiometry_variables[0]
    partial_specific_volume = stoichiometry_variables[1]
    izero = stoichiometry_variables[2]
    concentration = stoichiometry_variables[3]
    delta_rho = stoichiometry_variables[4]
        
    print(fraction_d2o)
    print(partial_specific_volume)
    print(izero)
    print(concentration)
    print(delta_rho)
    print(str(delta_rho))

    return 


fraction_d2o = [u'0.99', u'0.12', u'0.41']  # 1 value for each contrast
partial_specific_volume = [u'0.745', u'0.903']  # 1 value for each component
izero = [u'11.8', u'0.6', u'0.17']  # 1 value for each contrast
concentration = [u'3.7', u'3.6', u'3.1']  # 1 value for each contrast
delta_rho = [[u'-3.2', u'-5.7'], [u'1.6', u'0.26'], [u'0.031', u'-1.74']]  # 2 values for each contrast since there are 2 components
#do we need a blank delta_rho for the case where it is read from file? This should be taken care of in the initialization of the array in read_contrast_output_files?  Will need to test.
stoichiometry_variables = [fraction_d2o,partial_specific_volume,izero,concentration,delta_rho]  

print('all variables: ', stoichiometry_variables)

test_izero = stoichiometry_variables[2]
print('izero variables: ', test_izero)
print('len(stoichiometry_variables): ', len(stoichiometry_variables))

process_input_variables(stoichiometry_variables) #add one by one
#print ?


