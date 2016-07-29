import string,locale
from scipy.integrate import simps
from scipy.integrate import trapz
import numpy


def integrand(i,q,r,gr):

	value = (gr[i] - 1.0) * r[i] * numpy.sin(q * r[i])
	
	return value

def read_gr(filename):

	data=open(filename,'r').readlines()
	nl=len(data)
	r=numpy.zeros(nl,numpy.float)
	gr=numpy.zeros(nl,numpy.float)

	for i in xrange(nl):
		lin=string.split(data[i])
		tr=locale.atof(lin[0])
		tgr=locale.atof(lin[1])
		r[i]=tr ; gr[i]=tgr

	return r,gr

if __name__=='__main__': 

	pi = numpy.pi

	filename = 'argon_gr.txt'

	r,gr = read_gr(filename)

	highq = 12.0
	nsteps = len(r)
	stepq = (highq)/nsteps

	stepq = 0.0294
	nsteps = 399

	lowq = stepq
	
	outfile=open('gr.int','w')

#	no = 1.0/864.0
#	fact = 4 * pi * no

	correction = 2.7013/2.986804
	print 'correction = ',correction

	no = 0.02125
	fact = no * correction

	for s in xrange(nsteps+1):
		this_q = lowq + stepq*s
		value = numpy.zeros(len(r),numpy.float)
		for i in xrange(len(r)):
			value[i] = integrand(i,this_q,r,gr)

		sq = (fact*simps(value)/this_q) + 1
		sq2 = (fact*trapz(value)/this_q) + 1

		outfile.write('%f\t%f\t%f\n' % (this_q,sq,sq2))

	outfile.close()

