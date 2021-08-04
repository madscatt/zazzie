# -*- coding: utf-8 -*-

"""
Stoichiometry: Method to calculate the Mw of the components in a two-component system for which a contrast variation series of measurements has been performed. It should be used for systems that have multiple copies of the two components. While this module can be used for experiment planning, it is more likely that it would be used for data analysis. 

@author: susan.krueger
"""

#hard-wired for two components at this time; adding components would require changing the function definition on the fly, which could be done by going back to the basic equation, which is a "square of the sum"
#drho (contrast) is in units of 10**10 cm**-2
#conc is in units of 10^-3 g/cm**3 (mg/ml)
#izero is in cm^-1
#psv (partial specific volume) is in units of cm**3/g
#m1 and m2 are given in Da * 10^-6 and are converted to kDa

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import string
import locale
import numpy as np
from scipy.optimize import curve_fit


def get_mws(read_from_file,number_of_contrast_points,drho,psv,izero,conc,fname1,fname):

    def func(x, m1, m2):
        f = (x[0] * m1**2 + x[1] * m1*m2 + x[2] * m2**2 - x[3] * m1 - x[3] * m2)
        return f

    na = 6.023 #10**23
    
    v1 = locale.atof(psv[0])
    v2 = locale.atof(psv[1])
    print ('v1,v2: ', v1,v2)

    print ("number of contrast points: ", number_of_contrast_points)
#read the contrast values from file here.  They will be passed back as numbers. So, if we aren't reading the file, convert to numbers here and then change the definition of drho1 and drho2 below.

    outfile = io.open(fname,'w')
    if read_from_file == 1:
        outfile.write('input file: '+fname1+'\n')  #this is just to write the input filename into the output file so we know where the values came from

    x=np.zeros(shape = (number_of_contrast_points,4))

    for i in range(number_of_contrast_points):
        drho1 = locale.atof(drho[i][0])
        drho2 = locale.atof(drho[i][1])
        izero1 = locale.atof(izero[i])
        conc1 = locale.atof(conc[i])
        print('drho1, drho2, izero, conc: ', drho1, drho2, izero1, conc1)
        x0 = drho1**2*v1**2
        x1 = 2*drho1*drho2*v1*v2
        x2 = drho2**2*v2**2
        x3 = izero1*na/conc1
        x[i]=(x0,x1,x2,x3)
    
    print ('x: ', x)
    print ('length of x: ', len(x))

    y = np.zeros(len(x))
#    print ('y: ', y)

    x=x.T   #transpose x to have correct inputs for curve_fit

    popt, pcov = curve_fit(func, x, y)

    print ('calculated m1,m2: ', popt)
    outfile.write('calculated m1,m2: '+str(popt)+'\n')
    print ('calculated covariance: ', pcov)
    outfile.write('calculated covariance: '+str(pcov)+'\n')

    print ('-------------------------------')
    outfile.write('--------------------------------\n')
    print ('Final Results:')
    outfile.write('Final Results\n')

#convert molecular weights to kDa
    m1 = popt[0]*10**3
    m2 = popt[1]*10**3
    print ('m1, m2 (kDa): ', m1,m2)
    outfile.write('m1, m2 (kDa): '+str(m1)+'\t'+str(m2)+'\n')

    if np.all(pcov != np.inf):      #covariance matrix will be infinite if there are only two contrasts;  numpy.all tests whether all array elements along a given axis evaluate to True.

#rescale values in covariance matrix since molecular weights were converted to kDa
        rcov = pcov*10**6
        print ('covariance: ', rcov)
        outfile.write('covariance: ' + str(rcov)+'\n')
        sdm1 = np.sqrt(rcov[0][0])
        sdm2 = np.sqrt(rcov[1][1])
        print ('standard deviations, sdm1, sdm2: ', sdm1, sdm2)
        outfile.write('standard deviations, sdm1, sdm2: '+str(sdm1)+'\t'+str(sdm2)+'\n')
        corr = rcov[0][1]/(sdm1*sdm2)
        print ('correlation: ', corr)
        outfile.write('correlation: '+str(corr)+'\n')

    m = m1 + m2
    f1 = m1/m
    f2 = m2/m

#Do we want to output the covariance matrix?  How much meaning do the standard deviations have in this case?
#Do we want to calculate volume fractions here as well? Will they be needed?
#We will want to output all of the input values so we know what they were.

    print ('Mw (kDa), f1, f2: ', m, f1, f2)
    outfile.write('Mw (kDa), f1, f2: '+str(m)+'\t'+str(f1)+'\t'+str(f2))
    outfile.close()

if __name__ == "__main__":

#There are some inputs that will always be needed:
#   1) number of contrasts (integer >/= 2)
#   2) partial specific volumes (real, array)(We can make this part of the contrast calculator output; for now I am assuming the values will be input manually; we can provide default values for protein, RNA, DNA, but this will require similar input as for the contrast calculator where we ask whether the components are protein, DNA, RNA or other molecule)
#   3) I(0) values (real, array)
#   4) total concentration of complex (real array)
#Users will have a choice for other inputs.
#   1) Get contrasts (and perhaps partial specific volumes) from contrast calculator output file
#       a) %D2O (real, array) (for consistency with contrast calculator and SasCalc; converted to fraction D2O)
#       b) contrast calculator output filename  and read contrast values 
#   2) Input contrasts
#       a) contrast (real, array)
#       b) %D2O (real, array) (not really needed but maybe would help the user if they reattached to the job later)

#right now, we will just read the values from an input file that lists the fraction D2O and the contrasts in the same way as the contrast calculator output file OR input the everything by hand.

#SASSIE-web inputs 
    read_from_file = 0  #0 is no and 1 is yes to reading the contrasts from a contrast calculator output file
    number_of_contrast_points = 3   #must be >= number of components (2 in this case)
    number_of_components = 2
    fracd2o = [u'0.99', u'0.12', u'0.41'] #1 values for each contrast
    psv = [u'0.745', u'0.903']  #1 value for each component
    izero = [u'11.8', u'0.6', u'0.17']  #1 value for each contrast
    conc = [u'3.7',u'3.6', u'3.1']  #1 value for each contrast    
    path = ('./')
    fname = os.path.join(path,'99_12_41.out')    #output file (will always be specified)

    print (fname)

    
#Specify whether contrast values will be read from a file.
    if read_from_file == 1:  
        fname1 = os.path.join(path,'input_contrast.txt')   #input file (specify only if reading contrasts from file; default file can be '')
        drho = []
    else:
        fname1 = ''
        drho = [[u'-3.3',u'-5.7'],[u'1.6',u'0.2'],[u'0.033',u'-1.75']]  #2 values for each contrast point (since there are  2 components)

    print (fname1)

#        print ('drho: ', drho)
#        print (drho[0][0])
#        drho1 = locale.atof(drho[0][0])
#        print (drho1)

            
#An output filename won't need to be specified, as the output filename will be specified from the prefix given by the user (to be consistent with the contrast calculator).


    get_mws(read_from_file,number_of_contrast_points,drho,psv,izero,conc,fname1,fname)
