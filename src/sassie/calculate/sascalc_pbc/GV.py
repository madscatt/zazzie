#
# Hailiang Zhang
# NIST & UTK
#

import sys
import numpy


class GV:
    '''
    Home of golden vectors
    '''

    def __init__(self, n_golden_vectors):
        '''
        Functionality: initialize the golden vectors
        Input: n_golden_vectors--number of golden vectors
        Return: None
        '''
        self.__num_golden_vectors = n_golden_vectors
        self.__golden_vectors = self.__setup_golden_vectors()

    def get_golden_vectors(self):
        '''
        Functionality: get the golden vectors
        Input: None
        Return: golden vectors
        '''
        return self.__golden_vectors

    def __setup_golden_vectors(self):
        '''
        Functionality: set up the golden vectors
        Input: None
        Return: golden vectors
        Note: this method is based on Max's igolden_vectors.py script
        '''
        Npoints = self.__num_golden_vectors
        golden_vectors = numpy.zeros(([Npoints, 3]))
        if(numpy.mod(Npoints, 2) == 1):
            pi = numpy.pi
            phi_inv = 2.0 / (1 + numpy.sqrt(5))  # 1/(Golden Ratio)
            N = (Npoints - 1) / 2
            i = numpy.arange(-N, N + 1)  # run from -N...0...N
            sin_theta = numpy.cos(numpy.arcsin(2.0 * i / (Npoints)))
            # sin_theta=numpy.sqrt(1-(2.0*i/(Npoints))**2) # mathematically
            # equivalent
            cos_theta = 2.0 * i / (Npoints)
            phi = 2 * pi * i * phi_inv
            golden_vectors[:, 0] = sin_theta * \
                numpy.cos(phi)  # x,y,z coordinates
            golden_vectors[:, 1] = sin_theta * numpy.sin(phi)
            golden_vectors[:, 2] = cos_theta
        else:
            print 'The number of points must be odd'
            sys.exit()
        return golden_vectors


"""
The following is for devoloper testing/debugging purpose
"""
if __name__ == '__main__':
    gv = GV(51)
    print "golden vectors:\n", gv.get_golden_vectors()
