# -*- coding: utf-8 -*-

#    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#       READ DATA FILE
#
#       3/2023   --  Python 3 coding (in Data Interpolation:  Joseph E. Curtis
#       5/2023   --  standalone helper program coding      :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Read Data File** contains the method to read SAS data files. It is called by the **Decomposition Analysis** method and can also be used by the **Data Interpolation** method.
    
    Calls **Read File**.
    
'''
import os
import io
import sys
import string
import locale
import numpy


def read_file(data_file):
    '''
    **Read File** is the method to read SAS data files consisting of q, I and I error.  
    
    Notes:
    
        If the I error column doesn't exist, the method currently returns an error defined as I * error_magnitude, where error_magnitude is currently set to 0.1.
    
        TODO:   Return the error message and error_magnitude back to the main program
                Determine how to handle error values of 0.00


    Parameters
    ----------
    
    data_file: (string)
        name of the data file
        
    Returns
    -------
    
    nval: (int)
        the number of data points in the data file
    x: (float)
        the q values read from the data file
    y: (float)
        the I values read from the data file
    z: (float)
        the I error values read from the data file

    '''
    #    print('in readfile')

    fake_error = False
    error_magnitude = 0.1

    x = []
    y = []
    z = []
    nval = 0

    with io.open(data_file, 'r') as infile:
        contents = infile.readlines()

#   following line replaces control M (\r) with \n

    contents = [x.replace("\r", "\n") for x in contents]
#    print(contents)

    for line in contents:
        #        print(line)
        this_line = line.split()
#        print(this_line)
        try:
            qval = locale.atof(this_line[0])
            ival = locale.atof(this_line[1])
#            ival=abs(locale.atof(this_line[1]))

            try:
                error_value = locale.atof(this_line[2])
            except:
                error_value = error_magnitude * ival
                fake_error = True
# TODO: figure out how to handle error = 0.  This will be important for model data where 0.00 is in the error column.  The above error handling will work if the model data contain only two columns.
#             if(error_value == 0.0):
#                error_value = error_value + 1E-5

            x.append(qval)
            y.append(ival)
            z.append(error_value)
            nval = nval + 1

        except:
            pass
# NOTE: we may want to return this flag or the error message?
    if fake_error:
        message = 'Error values in I(q) were set to be ' + str(
            100 * error_magnitude) + '% of I(0) values since appropriate values were not found'
        print(message)

#    print('nval: ', nval)
#    print('x: ', x)
#    print('y: ', y)
#    print('z: ', z)

    return nval, x, y, z


if __name__ == "__main__":

    path = ('./')
    number_of_contrast_points = 7
    data_file_name = [os.path.join((path), '0.dat'), os.path.join((path), '10.dat'), os.path.join((path), '20.dat'), os.path.join(
        (path), '40.dat'), os.path.join((path), '80.dat'), os.path.join((path), '90.dat'), os.path.join((path), '100.dat')]
    print('data_file_name: ', data_file_name)

    contrast_variation_data = []
    numdatpts = []

    for item in data_file_name:
        #        print('item: ', item)
        #        print('data file: ', item)
        q = []
        i = []
        ierr = []

        numpts, q, i, ierr = read_file(item)

        print('numpts: ', numpts)
#    print('q: ', q)
#    print('i: ', i)
#    print('ierr: ', ierr)

        data = numpy.zeros((numpts, 3))
        for j in range(numpts):
            data[j][0] = q[j]
            data[j][1] = i[j]
            data[j][2] = ierr[j]
#        print(data)
        numdatpts.append(numpts)
        contrast_variation_data.append(data)

#    print('contrast_variation_data: ', contrast_variation_data)
#    print('number of data points: ', number_of_data_points)

    contrast_variation_data = numpy.asarray(
        contrast_variation_data, dtype=float)

# NOTE: conversion to numpy array isn't necessary for the tests below.
#    print('contrast_variation_data after convert to nparray: ', contrast_variation_data)

    for i in range(1, number_of_contrast_points):
        if numdatpts[0] != numdatpts[i]:
            print("The number of data points, " + str(numdatpts[i]) + ", in " + data_file_name[i] + " differs from " + str(
                numdatpts[0]) + " in the reference data file, " + data_file_name[0] + ".\n")
            sys.exit()

# if number of data points are the same for all files, then set the scalar value
    number_of_data_points = numdatpts[0]

    for i in range(number_of_contrast_points):
        for j in range(number_of_data_points):
            if contrast_variation_data[0][j][0] != contrast_variation_data[i][j][0]:
                print("The binning for point " + str(j+1) + " in dataset " +
                      data_file_name[i] + " differs from the that in " + data_file_name[0] + ".\n")
                sys.exit()

    print('read_file_is done')
