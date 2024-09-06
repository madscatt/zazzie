from __future__ import division
from scipy import spatial
import numpy as np
import numexpr as ne


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


def simple_debye_as_sum(q, coor):
    '''
    debye formula over Q range
    '''
    I = np.empty((len(q), 2))

    pw = pairwise_scipy(coor)

    for i, Q in enumerate(q):
        I[i, :] = [np.sum(ne_sinc(Q * pw)), Q]

    n_atoms = len(coor)
    I[:, 1] /= n_atoms

    return I


def simple_debye(q, coor):
    '''
    debye formula over Q range using np matrix math
    '''
    I = np.empty(len(q))

    pw = pairwise_scipy(coor)
    qpw = 0
    sinc_qr = ne_sinc(qr)

    for i, Q in enumerate(q):
        I[i] = np.sum(ne_sinc(Q * pw))

    return I



def glatter_kratky(q, coor):
    '''
    debye formula taken from Glatter & Kratky 1982
    '''
    n_q = len(q)

    pr = calc_pr(coor, n_q)
    # np.savetxt('pr.dat', pr)

    I = np.empty(pr.shape)
    I[:, 0] = q

    debye_integral_as_sum_wo_r(pr, I)
    n_atoms = len(coor)
    I[:, 1] /= n_atoms

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


def debye_integral_as_sum_wo_r(pr, I):
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
    dr = pr[1:, 0] - pr[:-1, 0]
    assert np.allclose(dr/dr[0], np.ones(len(dr)))
    dr = dr[0]
    for j in range(len(I)):
        for i in range(len(pr)):
            qr = pr[i, 0] * I[j, 0]
            sinc_qr = ne_sinc(qr)
            I[j, 1] += pr[i, 1] * sinc_qr

    I[:, 1] *= 4 * np.pi * dr

def debye_integral_as_sum_w_r(pr, I):
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
    dr = pr[1:, 0] - pr[:-1, 0]
    assert np.allclose(dr/dr[0], np.ones(len(dr)))
    dr = dr[0]
    for j in range(len(I)):
        for i in range(len(pr)):
            q = I[j, 0]
            r = pr[i, 0]
            sinc_qr = ne_sinc(q*r)
            r2 = r ** 2
            I[j, 1] += pr[i, 1] * sinc_qr * r2

    I[:, 1] *= 4 * np.pi * dr

def calc_pr(coor, n_q=100):

    pw = pairwise_scipy(coor)
    p, edges = np.histogram(pw, bins=n_q)
    r = (edges[1:] + edges[:-1])/2.0
    pr = np.concatenate((r.reshape(n_q,1), p.reshape(n_q, 1)), axis=1)

    return pr


if __name__ == "__main__":
    import doctest
    doctest.testmod()
