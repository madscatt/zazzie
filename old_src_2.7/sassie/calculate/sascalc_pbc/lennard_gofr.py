from __future__ import division
import numpy as np
import sys
import os
import sasmol.sasmol as sasmol
import compiledUtils.dna_overlap as dna_overlap


def deduceBoxSize(xCoor, yCoor, zCoor):
    # first check if there is CRYST1 thing
    xLen = xCoor.max() - xCoor.min()
    yLen = yCoor.max() - yCoor.min()
    zLen = zCoor.max() - zCoor.min()
#    print(xLen)
#    print(yLen)
#    print(zLen)
#    print('lens ^^^^^^^^^^^^^^^^^^^')
    return max(map(round, [xLen, yLen, zLen]))


'''
for this_segname in system.segnames():
    basis = 'segname[i] == "' + this_segname + '"'
    error, this_mask = system.get_subset_mask(basis)
    if len(error) > 0:
        print("ERROR: "+str(error))
        sys.exit()
    print(this_mask)
    print(type(this_mask))
    x = system.coor()[0][:,0]
    y = system.coor()[0][:,1]
    z = system.coor()[0][:,2]
    xCoor[int(this_segname)-1]=np.ma.masked_array(x,this_mask).compressed()
    yCoor[int(this_segname)-1]=np.ma.masked_array(y,this_mask).compressed()
    zCoor[int(this_segname)-1]=np.ma.masked_array(z,this_mask).compressed()
'''


class gofr_calc:
    ngr = 0
    g = 0
    length = 0
    delg = 0
    coors = 0
    npart = 0

    def __init__(self, mol, box_length, nhis=600):
        self.ngr = 0
        self.g = np.zeros(nhis)
        self.nhis = nhis
        self.coors = mol.coor()
        xCoor = self.coors[0][:, 0]
        yCoor = self.coors[0][:, 1]
        zCoor = self.coors[0][:, 2]
        self.length = deduceBoxSize(xCoor, yCoor, zCoor)
        self.delg = self.length / (2 * self.nhis)
        self.npart = len(xCoor)

    def g_hist(self, frame=0):
        '''
        takes in a sasmol object+ frame number
        outputs the g histogram for that frame
        '''
        xCoor = self.coors[frame][:, 0]
        yCoor = self.coors[frame][:, 1]
        zCoor = self.coors[frame][:, 2]
        npart = len(xCoor)
        self.ngr += 1
        self.g = dna_overlap.gr(
            1, self.coors[frame], self.g, self.length, self.delg, 0)
        '''
        dist = np.zeros((self.npart,self.npart))
        for part1 in xrange(self.npart-1):
            for part2 in xrange(part1+1,self.npart):
               # ngr += 1
                #dx = xCoor[part1] - xCoor[part2]
                #dy = yCoor[part1] - yCoor[part2]
                #dz = zCoor[part1] - zCoor[part2]

                #dx = dx - self.length*int(dx/self.length)
                #dy = dy - self.length*int(dy/self.length)
                #dz = dz - self.length*int(dz/self.length)

                #dr = np.sqrt(dx**2+dy**2+dz**2)
                dr = dist[part1,part2]
                if(dr<self.length/2): # can extend this to use the corners
                    ig =int(dr/self.delg)#int(dr/delg)
                    self.g[ig] += 2
        '''

    def g_of_r(self, sigma=3.405):
        for i in range(self.nhis):
            r = self.delg * (i + .5)
            vb = ((i + 1)**3 - i**3) * self.delg**3
            rho = self.npart / self.length**3
            nid = (4 / 3) * np.pi * vb * rho
            self.g[i] = self.g[i] / (self.npart * nid * self.ngr)
        x = np.linspace(0, self.length / 2, len(self.g)) * sigma
        return (x, self.g)


if __name__ == '__main__':

    path = '/home/schowell/ellipsoids_simulation/simulations/LJ_sphere_monomer'
    pdb_file = os.path.join(path, 'run_0.pdb')
    dcd_file = os.path.join(path, 'run_1.dcd')
    box_length_file = os.path.join(path, 'box_length.txt')

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file)
    mol.read_dcd(dcd_file)
    box_length = np.loadtxt(box_length_file)[:, 1]

    gc = gofr_calc(mol, box_length, nhis=200)
    gc.g_hist(200)

    r, g_of_r = gc.g_of_r()

    if True:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(r, g_of_r)
        plt.savefig('g_of_r.png', dpi=400)

    print('\m/ >.< \m/')
l