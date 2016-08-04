"""
A collection of useful function for extracting info from a pdb file
"""
from __future__ import division


def writePDB(xCoor, yCoor, zCoor, outFile='converted.pdb'):
    fo = open(outFile, 'w')
    for i in range(len(xCoor)):
        x = xCoor[i]
        y = yCoor[i]
        z = zCoor[i]

        name = 'O'

        fo.write('{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:<1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6}{:>6}          {:>2}{:>2}\n'.format(
            'ATOM', i, name, ' ', 'xxx', 'A', '1', '', x, y, z, 1, 0, name, ''))


def writeXYZ(x, y, z, names=None, outFile='converted.xyz'):
    fo = open(outFile, 'w')
    fo.write(str(len(x)) + '\n')
    fo.write('convertedFile\n')
    for i in range(len(x)):
        if names is None:
            name = 'C'
        else:
            name = names[i]
        fo.write('{:3>}\t{:>5}\t{:>5}\t{:>5}\n'.format(name, x[i], y[i], z[i]))
    fo.close()


def deduceBoxLen_mol(mol, frame=0, verbose=False):
    coors = mol.coor()[frame]
    x = coors[:, 0]
    y = coors[:, 1]
    z = coors[:, 2]
    return deduceBoxLen(x, y, z, verbose)


def deduceBoxLen(xCoor, yCoor, zCoor, verbose=False):
    from numpy import ceil
    # first check if there is CRYST1 thing
    xLen = xCoor.max() - xCoor.min()
    yLen = yCoor.max() - yCoor.min()
    zLen = zCoor.max() - zCoor.min()
    if(verbose):
        print(xLen)
        print(yLen)
        print(zLen)
        print('lens ^^^^^^^^^^^^^^^^^^^')
        print('minimum volume: ' + str(xLen * yLen * zLen))
    return max(ceil([xLen, yLen, zLen]))


def deduceBoxVol(xCoor, yCoor, zCoor):
    '''sloppy repeptition of code iwth above method'''
    xLen = xCoor.max() - xCoor.min()
    yLen = yCoor.max() - yCoor.min()
    zLen = zCoor.max() - zCoor.min()
    print(xLen)
    print(yLen)
    print(zLen)
    print('lens ^^^^^^^^^^^^^^^^^^^')
    print('minimum volume: ' + str(xLen * yLen * zLen))
    return xLen * yLen * zLen


def density(mol):
    '''mol is a sasmol object'''
    mol.mass()
    mass = mol.calcmass()
    print('mass: ' + str(mass))
    coors = mol.coor()[0]
    x = coors[:, 0]
    y = coors[:, 1]
    z = coors[:, 2]
    vol = deduceBoxVol(x, y, z)
    return mass / vol


def maxDist(coor):
    from scipy.spatial.distance import cdist
    return cdist(coor, coor)


def maxDist_mol(mol, frame=0):
    coors = mol.coor()[frame]
    x = coors[:, 0]
    y = coors[:, 1]
    z = coors[:, 2]
    return maxDist(x, y, z)
