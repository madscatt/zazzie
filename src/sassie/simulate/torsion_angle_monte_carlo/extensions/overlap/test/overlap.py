import sys,os,random
import sasmol.sasmol as sasmol
#import sassie.sasmol.sasmol as sasmol
#import sassie.simulate.monte_carlo.overlap as overlap
import locale
import numpy
#sys.path.append('./')
#import overlap

def getBondList(file_psf):
    count = 0
    qstart = False
    nbonds = 99999999
    for line in open(file_psf).readlines():
        if 'bonds' in line:
            nbonds = locale.atoi(line.split()[0])
            bond_list = numpy.empty(shape=[nbonds,2],dtype=numpy.int32)
            qstart = True
            continue
        if qstart:
            words = line.split()
            for i in xrange(len(words)/2):
                a = locale.atoi(words[2*i])-1
                b = locale.atoi(words[2*i+1])-1
                bond_list[count][0] = a
                bond_list[count][1] = b
                count += 1
                if count >= nbonds: return numpy.array(bond_list)


R_VDW={'H':1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'S':1.8, 'P':1.8}

m = sasmol.SasMol(0)
m.read_pdb("100_mer_min0.pdb")
coor = m.coor()[0]

radius = numpy.zeros(m.natoms())
for i,atom in enumerate(m.name()):
    radius[i] = R_VDW[atom[0]]

bond_list = getBondList("100_mer.psf")

#print coor
#print radius
#print bond_list

scale = 0.5
overlap.overlap(coor,radius,m.natoms(),scale,bond_list,len(bond_list))
#overlap.overlap(coor,radius,m.natoms(),scale)
