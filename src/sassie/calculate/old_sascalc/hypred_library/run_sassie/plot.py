import locale
import sys
from matplotlib import pylab as plt
from scipy.interpolate import interp1d
import numpy

def parse_pRDF(f):
	pRDFs={}
	count=0
	for line in open(f).readlines():
		words=line.split()
		if words[0]=="dstep":
			dstep=locale.atof(words[1])
			continue
		elif words[0]=="nstep":
			nstep=locale.atof(words[1])
			continue
		else:
			atom_type = words[0]
			pRDF=[]
			for word in words[1:]:
				pRDF.append(locale.atof(word))
			pRDFs[atom_type]=pRDF
	return (pRDFs,dstep)


f=sys.argv[1]
(pRDFs, dstep)=parse_pRDF(f)
colors={"C":"k", "HC":"r", "N":"b", "HN":"g", "O":"m", "HO":"y", "S":"c"}
X = numpy.linspace(0, 10+dstep, int(10/dstep)+1)
scale=3
Xf = numpy.linspace(0, 10+dstep, (int(10/dstep)+1)*scale)
for atom_type, pRDF in pRDFs.iteritems():
	Y=numpy.array(pRDF)
	#f = interp1d(X,Y,kind='cubic')
	plt.plot(X[:-1],Y[:-1],colors[atom_type], label=atom_type)
	#plt.plot(Xf[:-scale],f(Xf)[:-scale],colors[atom_type], label=atom_type)
	plt.hold(True)
plt.legend()
plt.xlabel("R (A)")
plt.ylabel("protons/A3")
plt.show()
