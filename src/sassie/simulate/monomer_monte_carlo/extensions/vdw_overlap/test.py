import sys,numpy
import sasmol.system as system
sys.path.append('../pairs/build/lib.linux-x86_64-3.9/')
sys.path.append('./build/lib.linux-x86_64-3.9/')
import pairs
import vdw_overlap

m = system.Molecule(0)
m.read_pdb('min3.pdb')
m.set_average_vdw()

na = m.natoms()
npairs = int(na*(na-1)/2.0)

cutoffarray = numpy.zeros(npairs,float)

#cutoffarray = pairs.pairs(m.atom_vdw(),npairs)
pairs.pairs(m.atom_vdw(),cutoffarray)

print('ca[0]  = ',cutoffarray[0])
print('ca[-1]  = ',cutoffarray[-1])

coor = m.coor()[0]

print('calling overlap')

fact = 0.31

check = vdw_overlap.overlap(coor,cutoffarray,fact)

print('check = ',check)

#subroutine overlap(coor1,natoms1,cutoffarray,npairs,check)
#subroutine pairs(avdw,natoms,npairs,cutoffarray)