from scipy.integrate import simps
import numpy

def integrate_gr(r, gr, rho=.0213, nQ=1000, loQ=0, hiQ=15, cut=-1):
    '''
    default rho taken from Yarnell et al., DOI: 10.1103/PhysRevA.7.2130
    '''
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


if __name__ == '__main__':
    '''
    reproduces integration from g(r) to S(q) published by Yarnell et al.,
    DOI: 10.1103/PhysRevA.7.2130
    '''
    do_plot = True

    filename = 'argon_85K_gr.dat'

    r, gr = numpy.genfromtxt(filename, unpack=True)
    q, sq = integrate_gr(r, gr)

    # Save output
    n_q = len(q)
    s_of_q = numpy.concatenate((q.reshape(n_q,1), sq.reshape(n_q, 1)), axis=1)
    numpy.savetxt('ian_sq_from_gr.dat', s_of_q)

    if do_plot:
        import matplotlib.pyplot as plt
        show = False

        plt.style.use('ggplot')
        plt.xlabel('Q')
        plt.ylabel('s(Q)')
        plt.plot(q, sq,label='Integrated',lw=3)
        q_file, sq_file = numpy.genfromtxt('argon_85K_sq.dat',unpack=True)
        plt.plot(q_file,sq_file,label='argon_85K_sq.dat')
        plt.legend()
        save_name = 'ian_sq_from_gr'
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        if show:
            plt.show()

    print('\m/ >.< \m/')