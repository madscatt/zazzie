import string
import locale
from scipy.integrate import simps
import numpy
import matplotlib.pyplot as plt


def integrate_gr(r, gr, rho=.02, nQ=500, loQ=0, hiQ=15, cut=-1):
    from scipy.integrate import simps

    r = numpy.asanyarray(r)
    gr = numpy.asanyarray(gr)

    def integrand(q, r, gr):
        return (gr - 1.0) * r * numpy.sin(q * r)

    r = r[:cut]
    gr = gr[:cut]
    Q = numpy.linspace(loQ, hiQ, nQ)
    fact = 4 * numpy.pi * rho

    sq = numpy.zeros_like(Q)
    for i, q in enumerate(Q):
        sq[i] = 1 + fact / q * simps(integrand(q, r, gr), r)
    return (Q, sq)


def read_gr(filename):

    data = open(filename, 'r').readlines()
    nl = len(data)
    r = numpy.zeros(nl, numpy.float)
    gr = numpy.zeros(nl, numpy.float)

    for i in xrange(nl):
        lin = string.split(data[i])
        tr = locale.atof(lin[0])
        tgr = locale.atof(lin[1])
        r[i] = tr
        gr[i] = tgr

    return r, gr

if __name__ == '__main__':

    filename = 'gr_argon.dat'

    r, gr = numpy.loadtxt(filename, unpack=True, usecols=[0, 2])
    q, sq = integrate_gr(r, gr)

    # Save output
    numpy.savez('output', q=q, sq=sq)

    # plot output
    plt.style.use('ggplot')
    plt.xlabel('Q')
    plt.ylabel('s(Q)')
    plt.plot(q, sq)
    plt.show()
