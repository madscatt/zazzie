'''
Takes a pdb and dcd in the pdb_NAME and dcd_NAME variables.
parallalizes over frames via the process_frame function that acutally does the work.
parallel over the frames of a given dcd. Only takes one dcd/pdb input pair.
'''
from __future__ import division
from scipy import spatial
import numpy as np
import numexpr as ne
import sasmol.sasmol as sasmol


def pairwise_numpy(X):
    '''
    Can be done faster with scipy cdist.
    '''
    return np.sqrt(((X[:, None, :] - X) ** 2).sum(-1))


def pairwise_scipy(X):
    '''
    May be faster than numpy pairwise.
    '''

    return spatial.distance.cdist(X, X, 'euclidean')



def ne_sinc(x):
    '''
    sinc function that is faster than nummpy sinc.
    Not exactly equal but within np.isclose
    '''
    a = ne.evaluate("sin(x)/x")
    a[np.isnan(a)] = 1
    return a


def process_frame(frame):
    '''
    debye formula over Q range
    '''
    I = np.empty(len(q))

    # pw = pairwise_numpy(coor[frame])
    pw = pairwise_scipy(coor[frame])

    for i, Q in enumerate(q):
        I[i] = np.sum(ne_sinc(Q * pw))

    return I


def glatter_kratky(q, coor, frame):
    '''
    debye formula taken from Glatter & Kratky 1982
    '''
    n_q = len(q)

    pr = calc_pr(coor[frame], n_q)
    np.savetxt('pr_{}.dat'.format(frame), pr)

    I = np.empty(pr.shape)
    I[:, 0] = q

    debye_integral(pr, I)

    return I

def debye_integral(pr, I):
    '''
    Calculate the integral from the Debye formula using numpy
    matrix operations (~14x faster than explicit summation)

    >>> pr = np.empty((4, 2))
    >>> I = np.empty((5, 2))
    >>> pr[:, 0] = pr[:, 1] = np.arange(1, 5)
    >>> pr[:, 1] *=np.pi/4
    >>> I[:, 0] = np.arange(5)
    >>> debye_integral(pr, I)
    >>> print(I)
    [[  0.          98.69604401]
     [  1.          11.20284903]
     [  2.           4.25595936]
     [  3.          -0.8644126 ]
     [  4.          -1.46050525]]

    '''
    qr = np.outer(pr[:, 0], I[:, 0])
    sinc_qr = ne_sinc(qr)
    I[:, 1] = 4 * np.pi * pr[:, 1].dot(sinc_qr)


def debye_integral_as_sum(pr, I):
    '''
    Calculate the integral from the Debye formula
    as a summation (~14x slower than numpy method)

    >>> pr = np.empty((4, 2))
    >>> I = np.empty((5, 2))
    >>> pr[:, 0] = pr[:, 1] = np.arange(1, 5)
    >>> pr[:, 1] *=np.pi/4
    >>> I[:, 0] = np.arange(5)
    >>> debye_integral_as_sum(pr, I)
    >>> print(I)
    [[  0.          98.69604401]
     [  1.          11.20284903]
     [  2.           4.25595936]
     [  3.          -0.8644126 ]
     [  4.          -1.46050525]]

    '''
    I[:, 1] = 0
    for j in range(len(I)):
        for i in range(len(pr)):
            qr = pr[i, 0] * I[j, 0]
            sinc_qr = ne_sinc(qr)
            I[j, 1] += pr[i, 1] * sinc_qr

    I[:, 1] *= 4 * np.pi


def calc_pr(coor, n_q=100):

    pw = pairwise_scipy(coor)
    p, edges = np.histogram(pw, bins=n_q)
    r = (edges[1:] + edges[:-1])/2.0
    pr = np.concatenate((r.reshape(n_q,1), p.reshape(n_q, 1)), axis=1)

    return pr


# Ouput Results
def getOutputDir():
    '''
    Returns a path to a time stamped directory under the Outputs folder as a string
    '''
    import os
    from time import strftime
    directory = "Outputs/" + strftime("%Y-%m-%d_%H-%M")
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory + "/"




if __name__ == "__main__":
    import doctest
    doctest.testmod()
