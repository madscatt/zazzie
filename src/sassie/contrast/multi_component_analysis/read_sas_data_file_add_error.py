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
#       3/2023   --  Python 3 coding (in Data Interpolation) :  Joseph E. Curtis
#       5/2023   --  standalone helper program coding        :  Susan Krueger
#       1/2024   --  add random Gaussian noise to model data :  Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Read Data File** contains the method to read experimental or model SAS data files. It can be called by any method that reads SAS experimental or model data files.

    Calls **Read File** and **Get Random Gaussian Noise**

'''
import os
import io
import sys
import locale
import sassie.contrast.multi_component_analysis.add_noise_to_sas_data as add_noise
#import add_noise_to_sas_data as add_noise

# NOTE:  the passed variables can become mcavars since this version of read_file is only used by multi-component analysis; we don't want to have to provide the sn parameters when reading a data file in other modules, where experimental data is assumed.


def read_file(data_file, sn_amplitude, sn_mean, sn_stddev, sn_bgd):
    '''
    **Read File** is the method to read SAS data files consisting of q, I and I error.

    Notes:

        If the I error column doesn't exist or if all of the error values are 0.0, as is the case for model SAS data, random Gaussian noise is added at each data point and a descriptive message is returned.

        If one or more (but not all) of the error values are 0.0, the error for those values is set to error_offset, which is currently set to 1.OE-4, and a descriptive message is returned. This would most likely occur for experimental data containing an error value that, by chance, is equal to exactly 0.0.



    Parameters
    ----------

    data_file: string
        name of the data file
    sn_amplitude: float (default value can be changed as an advanced option)
        amplitude of the Gaussian function that describes the signal-to-noise as a function of q
    sn_mean: float (default value; currently can't be changed)
        mean of the Gaussian function that describes the signal-to-noise as a function of q
    sn_stddev: float (default value; currently can't be changed)
        standard deviation of the Gaussian function that describes the signal-to-noise as a function of q
    sn_bgd: float (default value; currently can't be changed)
        background term of the Gaussian function that describes the signal-to-noise as a function of q

    Returns
    -------

    nval: integer
        the number of data points in the data file
    x: float array (dimension = nval)
        the q values read from the data file
    y: float array (dimension = nval)
        the I values read from the data file
    z: float array (dimension = nval)
        the I error values, either read from the data file or generated as described above if no errors are present
    message: string
        a string containing information about the error values; if no changes were made to the error values, the message is blank

    '''
#    print('in read_file')

    no_error = False
    zero_error = False
    add_error_flag = False
    error_offset = 1.0E-4

    x = []
    y = []
    tempz = []
    nval = 0
    zero_error_number = 0

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
                if (error_value == 0.0):
                    # don't set a new error_value here since we want to test whether all error values are zero first
                    zero_error = True
                    zero_error_number = zero_error_number + 1
            except:
                error_value = 0.0  # need to set a value here so the "append" doesn't fail
                no_error = True

            x.append(qval)
            y.append(ival)
            tempz.append(error_value)
            nval = nval + 1

        except:
            #            print("Failed to read data from " + data_file + ". Does the file contain at least two columns of data?")
            #            print("Stopping here.")
            #            sys.exit()
            pass  # there is a check in the calling method to deal with the exception, i.e., if nval = 0, an error will be set in the multi-component analysis filter

#    print('no_error, zero_error, nval, zero_error_number: ', no_error, zero_error, nval, zero_error_number)

# We want to finalize error values here so that there are no issues after the call to read data file.
    z = []
    add_error_flag = False

    if no_error:
        message = "No error values were found in " + data_file + \
            ", so errors were generated by adding random Gaussian noise at each data point.\n"
        add_error_flag = True
    elif (zero_error and zero_error_number == nval):
        # if all error values are zero, treat the errors the same as if there was no error column in the data
        message = "All error values in " + data_file + \
            " were found to be 0.0, so errors were generated by adding random Gaussian noise at each data point.\n"
        add_error_flag = True
    elif zero_error:
        #  test to see which error values are zero and only replace those with a small offset
        for i in range(nval):
            if tempz[i] == 0.0:
                #                print('i, tempz[i]: ', i, tempz[i])
                z.append(tempz[i] + error_offset)
            else:
                z.append(tempz[i])
        message = "Error values of 0.0 in " + data_file + \
            " were set to " + str(error_offset) + ".\n"
    else:
        message = " "
        z = [1.0 * value for value in tempz]

# Add errors here if necessary

    if (add_error_flag == True):
        import numpy
# initialize the parameters of the gaussian that describes S/N vs q as numpy scalars (numpy.float64)
        amplitude = numpy.float64(sn_amplitude)
        mean = numpy.float64(sn_mean)
        stddev = numpy.float64(sn_stddev)
        bgd = numpy.float64(sn_bgd)

        noise = add_noise.get_random_gaussian_noise(
            nval, x, y, amplitude, mean, stddev, bgd)

#        print('after adding noise')
#        print('noise: ', noise)
#        print('type(noise): ', type(noise))

# convert numpy array to list to match x and y
        z = noise.tolist()

#    print('nval: ', nval)
#    print('x: ', x)
#    print('y: ', y)
#    print('z: ', z)
#    print('type x, y, z: ', type(x), type(y), type(z))

    return nval, message, x, y, z


if __name__ == "__main__":

    import matplotlib.pyplot as plt

#   input variables
    path = ('./')

#    data_file_name = [os.path.join((path), '0d_pai_vn.iq')]
#    data_file_name = [os.path.join((path), '0_single_zero_error.dat')]
#    data_file_name = [os.path.join((path), '0d_pai_vn.iq'), os.path.join((path), '100d_pai_vn.iq')]
#    data_file_name = [os.path.join((path), 'no_data.dat')]
#    data_file_name = [os.path.join((path), '0d_pai_vn.iq'), os.path.join((path), '0p.dat')]
    data_file_name = [os.path.join((path), '0d_pai_vn.iq'), os.path.join(
        (path), '20d_pai_vn.iq'), os.path.join((path), '85d_pai_vn.iq'), os.path.join((path), '100d_pai_vn.iq')]

#    output_file_name = [os.path.join((path), '0_single_zero_error.out')]
#    output_file_name = [os.path.join((path), 'no_data_test.out')]
#    output_file_name = [os.path.join((path), 'pai_vn_with_error_test.out')]
#    output_file_name = [os.path.join((path), '0d_pai_vn_with_error.iq'), os.path.join((path), '100d_pai_vn_with_error.iq')]
#    output_file_name = [os.path.join((path), '0d_pai_vn_with_error.iq'), os.path.join((path), '0d_pai_vn_expt_data.iq')]
#    output_file_name = [os.path.join((path), '0d_pai_vn_with_error.iq'), os.path.join((path), '20d_pai_vn_with_error.iq'), os.path.join((path), '85d_pai_vn_with_error.iq'), os.path.join((path), '100d_pai_vn_with_error.iq')]
    output_file_name = [os.path.join((path), '0d_pai_vn_with_error_test.iq'), os.path.join((path), '20d_pai_vn_with_error_test.iq'), os.path.join(
        (path), '85d_pai_vn_with_error_test.iq'), os.path.join((path), '100d_pai_vn_with_error_test.iq')]

#   advance option input variable (if errors need to be added)
#   scale factor to adjust S/N.  The higher the number, the less noisy the data. This would be useful for adding less noise to high contrast model data and more noise to low contrast model data.
#    sn_amplitude = [50.0]   #default value
#    sn_amplitude = [400.0]
#    sn_amplitude= [400.0, 100.0]
#    sn_amplitude= [400.0, 50.0]
    sn_amplitude = [50.0, 50.0, 50.0, 50.0]  # defaults CV series test
    sn_amplitude = [400.0, 300.0, 10.0, 100.0]  # pai-vn CV series

#   plot all data on the same plot for comparison
    plot_file_name = os.path.join((path), 'pai_vn_with_error_cv_series.png')
#    plot_file_name = os.path.join((path), '0d_pai_vn_compare.png')
#    plot_file_name = os.path.join((path), '0_single_error_test.png')
#    plot_file_name = os.path.join((path), 'test.png')


#   initialization
    number_of_contrast_points = len(data_file_name)
    if (len(output_file_name) != number_of_contrast_points):
        print("output_file_name must have a length of " +
              str(number_of_contrast_points) + ". Stopping here.")
        sys.exit()
    if (len(sn_amplitude) != number_of_contrast_points):
        print("sn_amplitude must have a length of " +
              str(number_of_contrast_points) + ". Stopping here.")
        sys.exit()
# TODO:  Do we want the user to be able to change mean, stddev and bgd of the gaussian that describes S/N as a function of q as an advanced option? The bgd value of 1 should be reasonable in most cases and the mean can be set somewhere between 0.03 and 0.04. #stddev is the width of the gaussian describing S/N; the smaller the stddev, the noisier the data at higher q values since S/N decreases faster as a function of q. A default value of 0.05 worked well for all tested data
    sn_mean = 0.035
    sn_bgd = 1.0
    sn_stddev = 0.05

    print('data_file_name: ', data_file_name)
    print('output_file_name: ', output_file_name)
    print('plot_file_name: ', plot_file_name)
    print('number of contrast points: ', number_of_contrast_points)


# set up the plot parameters
    plt.xscale("log", nonpositive='clip')
    plt.yscale("log", nonpositive='clip')
#    plt.title('model CV data with gaussian noise')
    plt.xlabel("q(1/A)")
    plt.ylabel("I(q)")


# loop to read multiple files for contrast variation data
    for i in range(number_of_contrast_points):
        #        print('i, data file, output file, plot file, sn_amplitude: ', i,
        #              data_file_name[i], output_file_name[i], plot_file_name, sn_amplitude[i])

        numpts, message, q, iq, ierr = read_file(
            data_file_name[i], sn_amplitude[i], sn_mean, sn_stddev, sn_bgd)

        if (numpts == 0):
            print("Failed to read data from " +
                  data_file_name[i] + ". Does the file contain at least two columns of data?")
            sys.exit()

#        print('numpts: ', numpts)
#        print('q: ', q)
#        print('I(q): ', iq)
#        print('I(q) error: ', ierr)
#        print('message: ', message)

        outfile = open(output_file_name[i], 'w')
        outfile.write("#Input file used : %s\n" % (data_file_name[i]))
        if (message != " "):
            outfile.write("#" + message)
            if (message[0:19] != "Error values of 0.0"):
                # maybe just write the amplitude to the output file since it is the only thing that can be changed?
                outfile.write("#Amplitude, mean, stddev of Gaussian and bgd values used to describe S/N vs q: %f\t%f\t%f\t%f\n" %
                              (sn_amplitude[i], sn_mean, sn_stddev, sn_bgd))
        outfile.write("#q, I, I error\n")
        for j in range(numpts):
            outfile.write('%f\t%f\t%f\n' % (q[j], iq[j], ierr[j]))
        outfile.close()

        plt.errorbar(q, iq, yerr=ierr, marker='o',
                     markersize=2, label=data_file_name[i])

    plt.legend(loc='best')
    plt.savefig(plot_file_name)
    plt.show()

    print('read_file_is done')
