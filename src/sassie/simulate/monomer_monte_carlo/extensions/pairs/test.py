import sys,numpy
import sasmol.sasmol as s
sys.path.append('./')
import pairs


m = s.SasMol(0)
m.read_pdb('min3.pdb')
m.set_average_vdw()

na = m.natoms()
npairs = na*(na-1)/2.0

cutoffarray = numpy.zeros(npairs,numpy.float)

#cutoffarray = pairs.pairs(m.atom_vdw(),npairs)
pairs.pairs(m.atom_vdw(),cutoffarray)

print 'ca[0]  = ',cutoffarray[0]
print 'ca[-1]  = ',cutoffarray[-1]

#subroutine pairs(avdw,natoms,npairs,cutoffarray)
