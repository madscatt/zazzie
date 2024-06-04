import sys,numpy
import sasmol.system as system
sys.path.append('./')
sys.path.append('./build/lib.linux-x86_64-3.9/')

import pairs


m = system.Molecule(0)
m.read_pdb('min3.pdb')
m.set_average_vdw()

na = m.natoms()
npairs = int(na*(na-1)/2.0)

cutoffarray = numpy.zeros(npairs,numpy.float32)

#cutoffarray = pairs.pairs(m.atom_vdw(),npairs)
pairs.pairs(m.atom_vdw(),cutoffarray)

print('ca[0]  = ',cutoffarray[0])
print('ca[-1]  = ',cutoffarray[-1])

#subroutine pairs(avdw,natoms,npairs,cutoffarray)
