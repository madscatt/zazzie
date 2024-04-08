# -*- coding: utf-8 -*-
#
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
#       POLYNOMIAL FIT
#
#
#       2/2024       --      initial coding and testing:   Susan Krueger
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    **Polynomial Fit** contains the function that evaluates a polynomial of arbitrary order.

    INPUTS:
        x array, yarray, y error array, initial guess for polynomial coefficients
    
    OUTPUTS:
        optimized coefficients of the polynomial, covariance 2D array
        
    Requires **numpy.polynomial.polynomial.polyval** and **scipy.optimize.curve_fit**.
    
    Can be called by any method performing a fit to a polynomial.
    
'''

import numpy
import scipy.optimize
import chi_squared_correlation as chi_squared_correlation


# define the polynomial function
def polynomial_function(x, *params):
    '''
    **Polynomial Function** evaluates a polynomial of arbitrary order,

    :math:`y = c_{0} + c_{1}x + c_{2}x^2 ... + c_{n}x^n`, at points, :math:`x`,

    using **numpy.polynomial.polynomial.polyfit**. Note that the polynomial must have **order+1** coefficients. 

    The function is called using **scipy.optimize.curve_fit**, specifying sigma, the error in the y values, and with absolute_sigma = True, to perform a weighted fit: 

    optimized_coefficients, covariance = scipy.optimize.curve_fit(polynomial_function, x, y, sigma = y error, absolute_sigma = True, p0=[1.]*(order+1))

    An initial guess for the coefficients of the polynomial is **required**, as that is what defines the number of coefficients to be optimized.  **scipy.optimize.curve_fit** will fit data to a polynomial with order+1, if an initial guess, p0, of the same length is provided. Note that [1.0]*(order+1) defines an array of 1s of the correct length. A different value can be chosen as appropriate. 

    Parameters
    ----------
    x:  float array (dimension = number of contrast points)
        :math:`x` the x values at each contrast
    y:  float array (dimension = number of contrast points)
        :math:`y` the y values at each contrast
    y error: float array (dimension = number of contrast points)
        The errors on the :math:`y` values at each contrast
    p0:  float array (dimension = order + 1)
        Initial guess for the coefficients of the polynomial (default = 1.0 for all coefficients) 

    Returns
    -------
    optimized_coefficients: float array (dimension = order + 1)
        The optimized coefficients of the polynomial
    covariance:  2D float array (dimension = (order + 1) * (order + 1))
        The covariance matrix for the polynomial fit

    '''

    f = numpy.polynomial.polynomial.polyval(x, params)
    return f


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # test Stuhrmann equation
    # Values for Rg, Rg error and delta_rho are taken from PAI-VN data.

    output_file_name = 'test_polyfit.out'
    plot_file_name = 'test_polyfit.png'
    heatmap_file_name = 'test_heatmap_polyfit.png'

    # 3 contrast points
#    fraction_d2o = [0.0, 0.2, 0.85]
#    radius_of_gyration = [25.45, 24.95, 28.0]
#    radius_of_gyration_error = [0.07, 0.09, 3.0]
#    delta_rho = [[2.551, 5.104], [1.383, 3.928], [-2.415, 0.109]]

    # 4 contrast points
    fraction_d2o = [0.0, 0.2, 0.85, 1.0]
    radius_of_gyration = [25.45, 24.95, 28.0, 31.34]
    radius_of_gyration_error = [0.07, 0.09, 3.0, 0.4]
    delta_rho = [[2.551, 5.104], [1.383, 3.928],
                 [-2.415, 0.109], [-3.292, -0.773]]

    # variables for the two components
    number_of_components = 2  # hardwired at this time
    molecular_weight = [14.3, 44.2]  # kDa
    partial_specific_volume = [0.73, 0.73]

    # initial guess for the polynomial fit
    initial_guess = [1.5, 100.2, 10.0]
    initial_guess = [1., 1., 1.]  # default

    total_volume = 0.0
    volume_fraction = numpy.zeros(number_of_components, float)
    for i in range(number_of_components):
        total_volume = total_volume + \
            molecular_weight[i]*partial_specific_volume[i]

    for i in range(number_of_components):
        volume_fraction[i] = molecular_weight[i] * \
            partial_specific_volume[i]/total_volume

    number_of_data_points = len(fraction_d2o)
    print("number of data points: ", number_of_data_points)
    print('volume fraction: ', volume_fraction)
    print('fraction_d2o: ', fraction_d2o)
    print('radius_of_gyration: ', radius_of_gyration)
    print('radius_of_gyration_error: ', radius_of_gyration_error)
    print('delta_rho: ', delta_rho)

#    vector_shape = (number_of_data_points,1)
    delta_rho_inverse = numpy.zeros(number_of_data_points)
    rg_squared = numpy.zeros(number_of_data_points)
    rg_squared_error = numpy.zeros(number_of_data_points)
    rg_squared_calculated = numpy.zeros(number_of_data_points)
    diff = numpy.zeros(number_of_data_points)

    # STUHRMANN EQUATION:  Rg**2 = Rm**2 + alpha/delta_rho - beta/delta_rho**2
    # delta_rho_inverse are the x-coordinates (1/delta_rho) obtained from calculated delta_rho values
    # rg_squared are the y-coordinates (Rg**2) obtained from the experimental Rg values, and
    # rg_squared_error are the errors in the y-coordinates propagated from the experimental error in Rg: err(Rg**2) = 2*Rg*err(Rg)

    for i in range(number_of_data_points):
        delta_rho_inverse[i] = 1.0/(delta_rho[i][0] *
                                    volume_fraction[0] + delta_rho[i][1]*volume_fraction[1])
        rg_squared[i] = radius_of_gyration[i]*radius_of_gyration[i]
        rg_squared_error[i] = 2.0*radius_of_gyration[i] * \
            radius_of_gyration_error[i]
    print('delta_rho_inverse, rg_squared, rg_squared_error: ',
          delta_rho_inverse, rg_squared, rg_squared_error)


# fit a second-order polynomial to the data to get alpha, beta and Rm


# NOTE:  The errors on the parameters and chi-squared, as calculated here (as defined in Bevington), are the same as those calculated for a weighted polynomial fit using IGOR.

# The initial guess array is now an input that can be changed as an advanced option in case the fit fails with the default values of 1.  Since we know the degree of the polynomial we want to fit, we can allow a different number for each coefficient.
#    order = 2
#    initial_guess = [1.]*(order+1)

    print('initial_guess: ', initial_guess)
    fit, covariance = scipy.optimize.curve_fit(
        polynomial_function, delta_rho_inverse, rg_squared, sigma=rg_squared_error, absolute_sigma=True, p0=initial_guess)
    print('fit, covariance: ', fit, covariance)
    print('length fit: ', len(fit))

    standard_deviation, correlation = chi_squared_correlation.correlation_from_covariance(
        covariance)

    print('Rm**2, error: ', fit[0], standard_deviation[0])
    print('alpha, error: ', fit[1], standard_deviation[1])
    print('beta, error: ', - fit[2], standard_deviation[2])

    rg_squared_calculated = polynomial_function(delta_rho_inverse, *fit)
    print('rg_squared_calculated: ', rg_squared_calculated)

    reduced_chi_squared = chi_squared_correlation.get_reduced_chi_squared(
        rg_squared, rg_squared_error, rg_squared_calculated, number_of_data_points, 3)
    print('reduced chi-squared: ', reduced_chi_squared)

    outfile = open(output_file_name, 'w')
    outfile.write("#1/delta_rho, Rg**2, Rg**2 error, Rg**2 calc, diff\n")
    for i in range(number_of_data_points):
        outfile.write('%f\t%f\t%f\t%f\t%f\n' % (
            delta_rho_inverse[i], rg_squared[i], rg_squared_error[i], rg_squared_calculated[i], diff[i]))
    outfile.close()

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

#    ax0.plot(delta_rho_inverse, rg_squared,'bo', markersize = 3, label = "data")
    ax0.errorbar(delta_rho_inverse, rg_squared,
                 yerr=rg_squared_error, fmt='bo', markersize=3, label="data")
#    ax0.plot(delta_rho_inverse,rg_squared_calculated,'r--')
    ax0.plot(delta_rho_inverse, rg_squared_calculated,
             'ro', markersize=3, label="calc")
    ax0.set_title('Stuhrmann plot')
    ax0.set_ylabel('Rg**2')
    ax0.legend(loc="upper right")
#
    ax1.plot(delta_rho_inverse, diff, 'ko', markersize=3)
    ax1.set_title('residuals')
    ax1.set_xlabel('1/delta_rho')
    plt.savefig(plot_file_name)
    plt.show()


#   visualize correlation matrix as a heatmap
    plt.imshow(correlation, cmap='coolwarm', interpolation='none')
    plt.colorbar()
    plt.xticks(range(len(fit)), ['Rm**2', 'alpha', 'beta'])
    plt.yticks(range(len(fit)), ['Rm**2', 'alpha', 'beta'])
    plt.title('Correlation Heatmap')
    plt.savefig(heatmap_file_name)
    plt.show()

#
