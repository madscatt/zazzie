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


def calc_pmi(input_filename):

    frame = 0

    m1 = sasmol.SasMol(0)
    m1.read_pdb(input_filename)
    m1.center(frame)

    uk, ak, I = m1.calcpmi(frame)

    idx = uk.argsort()[::-1]   
    uk = uk[idx]
    ak = -1.0*ak[:,idx]

    nak1=numpy.array([ak[0][0],ak[1][0],ak[2][0]])
    nak2=numpy.array([ak[0][1],ak[1][1],ak[2][1]])
    nak3=numpy.array([ak[0][2],ak[1][2],ak[2][2]])

    ak=numpy.array([nak1,nak2,nak3])
    

    print 'uk = ', uk
    print 'ak = ', ak
    print
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]
    print
    print 'I = ', I

    print
    print "first rotation"
    # first rotation
    print 'vector 1 = ', ak[2]
    print 'vector 2 = ', [0,0,1.0]

    rotvec = numpy.cross(ak[2],numpy.array([0,0,1.0]))
    sine = numpy.linalg.norm(rotvec)
    cosine = numpy.dot(ak[2], numpy.array([0,0,1.0]))
    try:
        angle = 1.0*math.atan(sine/cosine)
    except:
        print 'cosine = 0\nstopping here\n\n\n'
        sys.exit()
    print 'angle = ',angle
    print 'rotvec= ',rotvec
    print 'angle should be = ', 1.570457424779631
    print 'rotvec should be = ', '-0.9999999425 2.113354e-7 0.0'
    #m1.general_axis_rotate(frame,theta,ux,uy,uz)
    m1.center(frame)
    m1.general_axis_rotate(frame,angle,rotvec[0],rotvec[1],rotvec[2])

    uk, ak, I = m1.calcpmi(frame)
    idx = uk.argsort()[::-1]   
    uk = uk[idx]
    ak = -1.0*ak[:,idx]

    nak1=numpy.array([ak[0][0],ak[1][0],ak[2][0]])
    nak2=numpy.array([ak[0][1],ak[1][1],ak[2][1]])
    nak3=numpy.array([ak[0][2],ak[1][2],ak[2][2]])

    ak=numpy.array([nak1,nak2,nak3])
    
    #ak=-1.0*ak
    print 'and the new eigenvectors are:'
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]
    

    
    print
    print "second rotation"
    print 'vector 1 = ', ak[1]
    print 'vector 2 = ', [0,1.0,0]
    # second rotation
    com = m1.calccom(frame)
    rotvec = 1.0*numpy.cross(ak[1],numpy.array([0,1.0,0]))
    sine = numpy.linalg.norm(rotvec)
    cosine = numpy.dot(ak[1], numpy.array([0,1.0,0]))
    try:
        angle = 1.0*math.atan(sine/cosine)
    except:
        print 'cosine = 0\nstopping here\n\n\n'
        sys.exit()
    print 'angle = ',angle
    print 'rotvec= ',rotvec
    print 'angle should be = ', 1.5707962131085322
    print 'rotvec should be = ', '5.959831601542385e-8 0.0 -0.9999999999999918'
    
    #m1.general_axis_rotate(frame,theta,ux,uy,uz)
    m1.center(frame)
    m1.general_axis_rotate(frame,angle,rotvec[0],rotvec[1],rotvec[2])

    uk, ak, I = m1.calcpmi(frame)
    #ak=-1.0*ak
    print 'and the new eigenvectors are:'
    print 'ak[0] = ', ak[0]
    print 'ak[1] = ', ak[1]
    print 'ak[2] = ', ak[2]

    '''
    coor_array = ak
    sub_m1 = make_and_write_pdb(coor_array, 'pmi.pdb')
    coor_sub_m1 = sub_m1.coor()[0]

    unit_array = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    sub_m2 = make_and_write_pdb(unit_array, 'ua.pdb')
    coor_sub_m2 = sub_m2.coor()[0]

    com_sub_m1 = numpy.zeros(3, numpy.float32)
    com_sub_m2 = numpy.zeros(3, numpy.float32)
    print 'before calling align'
    print m1.coor()[0][0]

    # m1.align(frame,coor_sub_m1,com_sub_m1,coor_sub_m2,com_sub_m2)
    # sub_m1.align(frame,coor_sub_m1,com_sub_m1,coor_sub_m2,com_sub_m2)
    sub_m1.align(frame, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)

    coor_sub_m1 = sub_m1.coor()[0]
    m1.align(frame, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)
    sub_m1.align(frame, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)
    print 'after calling align'
    print m1.coor()[0][0]

    #uk,ak,I = sub_m1.calcpmi(frame)
    # print 'ak[0] = ',ak[0]
    # print 'ak[1] = ',ak[1]
    # print 'ak[2] = ',ak[2]
    sub_m1.write_pdb('sub_m1_pmi_' + input_filename, frame, 'w')
    '''

    m1.write_pdb('pmi_' + input_filename, frame, 'w')


if __name__ == '__main__':

    input_filename = 'test.pdb'
    input_filename = 'tor3.pdb'

    calc_pmi(input_filename)
