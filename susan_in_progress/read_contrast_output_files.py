# -*- coding: utf-8 -*-

"""
Read Contrast Output Files:  Reads the output files from the contrast calculator to obtain contrast, Mw, I(0) or other information for further calculation.  It can be used by both experiment planning and data analysis tools.

@author: susan.krueger
"""

#right now, I am only reading contrast values from the *_contrast.txt output files in order to obtain contrast values at a particular %D2O.  Values MUST BE READ the same order as the other variables are entered.  I am not sure how we are going to do this in SASSIE-web.  We need some way to know which partial specific volume, concentration, I(0), etc. goes with which contrast value.  I think we will need to parse the titles of the columns and ask the user which one is which.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import sys
import string
import locale
import numpy as np
from scipy import stats


def read_contrast_file(fname1,fracd2o,number_of_contrast_points,number_of_components):


    data_array = []
    if(number_of_contrast_points == len(fracd2o)):
        print ('number of contrasts: ',number_of_contrast_points)
    else:
        print ('number of contrasts is not equal to length of fracd2o')
        sys.exit()
        
    if(number_of_contrast_points < number_of_components):
        print ('number of contrast points must be >= number of components')
        sys.exit()

    with io.open(fname1, 'r') as infile:
        a = infile.readlines()
#    infile.close()


#Read the fracd2o values of interest and associate them with numpy array x
    x=np.zeros(shape = (number_of_contrast_points))    
    for i in range(number_of_contrast_points):
        x[i] = locale.atof(fracd2o[i])
#    print ('frac D2O: ', x)
    

    number_of_lines = len(a)
#    print ("number of lines ", number_of_lines)

#Find linear equation to calculate desired drho values. This equation can eventually be part of the contrast calculator output file so it can be read directly instead of having to re-calculate it.

#First, find out where the desired drho values are in the file.
    for i in range(number_of_lines):
#        print(i, a[i][0:4])
        if(a[i][0:4] == u'# fr'):
            starting_line = i
#            print('start here: ', i)
            
#    print(a[starting_line])
    l=string.split(a[starting_line],',')
#    print(l)
    
#Associate the contrast arrays with a title [protein, RNA, DNA, Molecule 1, etc.], as we may have to ask users which contrast goes first (to be in agreement with the way they are inputing I(0), concentration, etc.
    title = []
    title.append(u' frac D2O')
    for i in range(number_of_components+1):
        title.append(l[i+1])
#    print(title)
#    print(title[0])
#    print(title[1])
#    print(title[2])


#Read the lines and associate them with numpy array data [[fracd2o_from_file] [contrast1] [contrast2]]

    number_of_data_points = number_of_lines - starting_line - 1
    print('number of data points for linear fit: ', number_of_data_points)
    
    data = np.zeros(shape = (number_of_components + 1, number_of_data_points))
#    print('data: ', data)

    for i in range(starting_line+1,number_of_lines):
#        print('i: ', i)
        l=string.split(a[i],'\t')
#        print(l)
#        print(l[0])
#        print(l[1])
#        print(l[2])
        for j in range(number_of_components+1):
#            print('j,line(j): ',j, l[j])
            data[j][i-starting_line-1] = locale.atof(l[j])
#            print(data[j][i-starting_line-1])

#    print('data: ', data)
#    print('data[0]: ', data[0])
#    print('data[1]: ', data[1])
#    print('data[0][0]: ', data[0][0])


#Find the slope and intercept for each component and put them into arrays. (I'm doing it this way in case we end up just reading the values from an equation that is included in the contrast calculator output file.)
#Match points are printed as a check.
    slope_values = np.zeros(shape=number_of_components)
    intercept_values = np.zeros(shape=number_of_components)
#    print('slope_values: ', slope_values)
#    print('intercept_values: ', intercept_values)
    for i in range(number_of_components):
        slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
            data[0], data[i+1])
        x_intercept = -intercept / slope
#        print ('slope, intercept, x_intercept: ', slope, intercept, x_intercept)
        matchpoint = x_intercept * 100
        print('match point '+title[i+1]+': ', matchpoint)
        slope_values[i] = slope
        intercept_values[i] = intercept
#    print('slope_values: ', slope_values)
#    print('intercept_values: ', intercept_values)
    
#Find the drho values and put them into a drho array
#Converted the drho values into list of lists with unicode values to pass to other methods, i.e., drho = [[u'value', u'value'],[u'value',u'value'],[u'value',u'value]]
#The unicode conversion not be necessary for Python 3? Will it be unicode by default if just str is used?
    drho = [[] for n in range (number_of_contrast_points)]
#    print('drho: ', drho)
    for i in range(number_of_contrast_points):
        for j in range(number_of_components):
            test = '%.3f' % (slope_values[j]*x[i] + intercept_values[j])
            utest = unicode(test)
#            print(utest,type(utest))
            drho[i].append(utest)
#    print('drho: ', drho)

    return drho


if __name__ == "__main__":

    path = ('./')
    fname1 = os.path.join(path,'input_contrast.txt')
    print (fname1)
    number_of_contrast_points = 3
    number_of_components = 2 
    fracd2o = [u'0.99', u'0.12', u'0.41']

    drho = read_contrast_file(fname1,fracd2o,number_of_contrast_points,number_of_components)
    
    print('drho after return to main: ', drho)  
    
