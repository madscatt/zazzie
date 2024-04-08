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
#       ADD NOISE TO MODEL DATA
#
#       2/2024   --  initial coding                        :   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Get Gaussian Random Noise** contains the method to add Gaussian random noise to model SAS data. It is called by the **Read Data File** method.
    
    Calls **Gaussian Function**.
    
    Requires **numpy.random.default_rng().normal**.
    
'''
import numpy


def gaussian_function(x, amplitude, mean, stddev, bgd):
    r'''

    **Gaussian Function** is the function that describes the behavior of the signal-to-noise (S/N) of SAS data as a function of q, where S/N decreases (noise increases) as a function of q in a Gaussian-like manner. This function is used to determine the initial guess for the error on I(q) at each q value before random Gaussian noise is added.

    Called from **Get Random Gaussian Noise**. 

    Note
    ----

    Gaussian function is the equation:

    :math:`y = A e^{- \frac{(x - \mu )^{2}}{2 \sigma^{2}}} + B`, where

    :math:`A` is the amplitude

    :math:`\mu` is the mean

    :math:`\sigma` is the standard deviation and

    :math:`B` is the background 

    From S/N calculations using experimental SANS data, it was found that S/N can be approximated by a Gaussian function with amplitude depending on the contrast, concentration, etc. The amplitude of the Gaussian is ~ the maximum S/N value, which is not known *a priori* for model data. The background term is ~ 1.0, since S/N ~ 1.0 at high q. If the mean and standard deviation are fixed to values that approximate those found from the experimental SANS data, the amplitude can be changed to simulate more (lower amplitude) or less (higher amplitude) noise in the model data.  

    Parameters
    ----------

    x: float array (dimension = number_of_data_points)
        the array containing the q values in A\ :sup:`-1`\
    amplitude: numpy.float64
        the amplitude of the Gaussian function
    mean: numpy.float64
        the mean of the Gaussian function
    stddev: numpy.float64
        the standard deviation of the Gaussian function
    bgd: numpy.float64
        the background term in the Gaussian function

    Returns
    -------

    signal_to_noise: float array (dimension = number of data points)
        the Gaussian function describing the behavior of S/N vs q

    '''
# NOTE:  q is assumed to be in 1/Angstrom! This is important since the mean value of the Gaussian assumes units of 1/Angstrom.

    signal_to_noise = amplitude * \
        numpy.exp(-(x - mean)**2 / (2.0 * stddev**2)) + bgd
    return signal_to_noise


def get_random_gaussian_noise(number_of_data_points, x, y, amplitude, mean, stddev, bgd):
    r'''
    **Get Random Gaussian Noise** adds random Gaussian noise to model SAS data. 

    Calls **Gaussian Function**.

    Requires **numpy.random.default_rng().normal**.

    Note
    ----

    An initial guess for Ierror(q) at each q value is determined based on S/N vs q as defined in **Gaussian Function**.

    :math:`S/N(q) = \frac{I(q)}{Ierror(q)}`, thus

    :math:`Ierror(q)_{initial} = \frac{I(q)}{S/N(q)}`.

    Then, random Gaussian noise is generated from :math:`Ierror(q)_{initial}` at each q value using the numpy default normal random number generator. These randomized values define the error on each data point.

    **numpy.random.default_rng().normal** returns random samples from a Gaussian (normal) distribution. The probability density for the distribution is:

    :math:`p(x) = \frac{1}{\sqrt{2 \pi \sigma^{2}}} e^{- \frac{(x - \mu )^{2}}{2 \sigma^{2}}}`

    The peak of this function is at the mean, :math:`\mu`, and its spread increases with the standard deviation, :math:`\sigma`. The function reaches ~0.6 times its maximum at :math:`x + \sigma` and :math:`x - \sigma`. Thus, values lying close to the mean should be returned.

    Here, :math:`\mu = 0.0`, so the distribution has a mean value close to the initial guess and :math:`\sigma = Ierror(q)_{initial}` at each q value to generate values near that of the initial guess.

    Parameters
    ----------

    number_of_data_points: float
        the number of data points in the data file
    x: float array (dimension = number_of_data_points)
        the array containing the q values in A\ :sup:`-1`\
    y: float array (dimension = number_of_data_points)
        the array containing the I(q) values
    amplitude: numpy.float64
        the amplitude of the Gaussian function that describes the signal-to-noise as a function of q
    mean: numpy.float64
        the mean of the Gaussian function that describes the signal-to-noise as a function of q
    stddev: numpy.float64
        the standard deviation of the Gaussian function that describes the signal-to-noise as a function of q
    bgd: numpy.float64
        the background term in the Gaussian function that describes the signal-to-noise as a function of q

    Returns
    -------

    noise: float array (dimension = number_of_data_points)
        the random Gaussian error on each data point, given by the absolute value of the numpy normal random number generator result 

    '''

    print('in get_random_gaussian_noise')

    signal_to_noise = gaussian_function(x, amplitude, mean, stddev, bgd)
#    print('signal to noise: ', signal_to_noise)
#    print('length of signal to noise: ', len(signal_to_noise))

# find initial yerror based on the q-dependence of S/N
    initial_yerr = numpy.zeros(number_of_data_points)
    for i in range(number_of_data_points):
        initial_yerr[i] = y[i]/signal_to_noise[i]
#    print('initial yerr: ', initial_yerr)

# Add Gaussian noise to initial_yerr at each q value using numpy.random.default_rng().normal
    mu = 0.0
    sigma = initial_yerr
# NOTE: for debugging, a random seed can be included; it is set to 1260 in the example below
#    noise = abs(numpy.random.default_rng(1260).normal(
#        mu, sigma, number_of_data_points))
    noise = abs(numpy.random.default_rng().normal(
        mu, sigma, number_of_data_points))
#    print('noise: ', noise)

    return noise
