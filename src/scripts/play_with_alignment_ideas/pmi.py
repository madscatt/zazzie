# -*- coding: utf-8 -*-
import sasmol.sasmol as sasmol
import numpy, math
import sys 

def make_and_write_pdb(coor_array, pdbfilename):

    m1 = sasmol.SasMol(0)
    natoms = len(coor_array)
    m1.setNatoms = natoms

    atom = []
    index = []
    name = []
    loc = []
    resname = []
    chain = []
    resid = []
    rescode = []
    x = []
    y = []
    z = []
    occupancy = []
    beta = []
    segname = []
    element = []
    charge = []
    moltype = []

    coor = numpy.zeros((1, natoms, 3), numpy.float32)

    resid_count = -999
    resid_count = 1
    for i in xrange(natoms):
        atom.append('ATOM')
        if(i < 99998):
            index.append(i + 1)
        else:
            index.append(99999)

        name.append('CA')
        loc.append(' ')
        resname.append('GLY')
        chain.append(' ')
        if(resid_count < 9999):
            resid.append(resid_count)
            resid_count += 1
        else:
            resid.append(9999)
        rescode.append(' ')
        x.append(coor_array[i][0])
        y.append(coor_array[i][1])
        z.append(coor_array[i][2])
        occupancy.append("  0.00")
        beta.append("  0.00")
        element.append('C')
        segname.append('DUM')
        charge.append(' ')
        moltype.append('protein')

    m1.setAtom(atom)
    m1.setIndex(numpy.array(index, numpy.int))
    m1.setName(name)
    m1.setLoc(loc)
    m1.setResname(resname)
    m1.setChain(chain)
    m1.setResid(numpy.array(resid, numpy.int))
    m1.setRescode(rescode)
    m1.setOccupancy(occupancy)
    m1.setBeta(beta)
    m1.setElement(element)
    m1.setSegname(segname)
    m1.setCharge(charge)
    m1.setMoltype(moltype)

    x = numpy.array(x, numpy.float32)
    y = numpy.array(y, numpy.float32)
    z = numpy.array(z, numpy.float32)

    coor[0, :, 0] = x
    coor[0, :, 1] = y
    coor[0, :, 2] = z

    m1.setCoor(coor)

    m1.write_pdb(pdbfilename, 0, 'w')

    return m1


def align_on_pmi(input_filename):

    frame = 0

    m1 = sasmol.SasMol(0)
    m1.read_pdb(input_filename)
    m1.center(frame)

    uk, ak, I = m1.calcpmi(frame)

    ak[0] = -1 * ak[0]
    print 
    print 'initial eigenvectors:'
    print 'uk = ', uk
    print 'ak = ', ak
    print
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]
    print
    print 'I = ', I
    print 'det ak = ', numpy.linalg.det([ak])

    print
    print "first rotation"
    # first rotation
    print 'vector 1 = ', ak[2]
    print 'vector 2 = ', [0,0,1.0]

    rotvec = numpy.cross(ak[2],numpy.array([0,0,1.0]))
    sine = numpy.linalg.norm(rotvec)
    rotvec = rotvec/numpy.linalg.norm(rotvec)
    cosine = numpy.dot(ak[2], numpy.array([0,0,1.0]))
    try:
        angle = 1.0*math.atan(sine/cosine)
    except:
        print 'cosine = 0\nstopping here\n\n\n'
        sys.exit()
    print 'angle = ',angle
    print 'rotvec= ',rotvec
    
    #print 'angle should be = ', 1.570457424779631
    #print 'rotvec should be = ', '-0.9999999425 2.113354e-7 0.0'
    
    
    #m1.general_axis_rotate(frame,theta,ux,uy,uz)
    m1.center(frame)
    m1.general_axis_rotate(frame,angle,rotvec[0],rotvec[1],rotvec[2])

    uk, ak, I = m1.calcpmi(frame)
    print 'and the new eigenvectors are:'
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]

    ak[0] = -1 * ak[0]
    #ak = ak/numpy.linalg.norm(ak)
    print
    print 'det ak = ', numpy.linalg.det([ak])

    print "second rotation"
    print 'vector 1 = ', ak[1]
    print 'vector 2 = ', [0,1.0,0]
    # second rotation
    com = m1.calccom(frame)
    rotvec = 1.0*numpy.cross(ak[1],numpy.array([0,1.0,0]))
    sine = numpy.linalg.norm(rotvec)
    rotvec = rotvec/numpy.linalg.norm(rotvec)
    cosine = numpy.dot(ak[1], numpy.array([0,1.0,0]))
    try:
        angle = 1.0*math.atan(sine/cosine)
    except:
        print 'cosine = 0\nstopping here\n\n\n'
        sys.exit()
    print 'angle = ',angle
    print 'rotvec= ',rotvec
    #print 'angle should be = ', 1.5707962131085322
    #print 'rotvec should be = ', '5.959831601542385e-8 0.0 -0.9999999999999918'
    
    m1.center(frame)
    m1.general_axis_rotate(frame,angle,rotvec[0],rotvec[1],rotvec[2])

    uk, ak, I = m1.calcpmi(frame)
    ak[0] = -1 * ak[0]
    
    print 'and the new eigenvectors are:'
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]

    m1.write_pdb('pmi_' + input_filename, frame, 'w')

if __name__ == '__main__':

    input_filename = 'test.pdb'
    input_filename = 'tor3.pdb'
    input_filename = 'fc_rough.pdb'
    input_filename = 'lysozyme.pdb'
    input_filename = 'min3.pdb'

    align_on_pmi(input_filename)
