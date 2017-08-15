# -*- coding: utf-8 -*-
import sasmol.sasmol as sasmol
import numpy, math
import sys 

def main():
    m1=sasmol.SasMol(0)

    starting_vector = [0.0, 0.0, 1.0]
    starting_vector = [-0.08711655, -0.512547, 0.85422847]
    m1.setCoor(numpy.array([starting_vector]))

    print m1.coor()[:]
  
    goal_vector = [0.0, 0.0, 1.0]

    gvector = numpy.array(goal_vector)
    rotvec = numpy.cross(m1.coor()[0],gvector)
    sine = numpy.linalg.norm(rotvec)
    rotvec = rotvec/numpy.linalg.norm(rotvec)
    cosine = numpy.dot(m1.coor()[0], gvector)
    try:
        angle = 1.0*math.atan(sine/cosine)
    except:
        print 'cosine = 0\nstopping here\n\n\n'
        sys.exit()
    print 'angle = ',angle
    print 'rotvec= ',rotvec
    ux=rotvec[0] ; uy=rotvec[1] ; uz=rotvec[2]
    #m1.general_axis_rotate(frame,theta,ux,uy,uz)
    m1.general_axis_rotate(0,angle,ux,uy,uz)
    #m1.rotate(0,'z',angle)
    print m1.coor()[:]

if __name__ == '__main__':
    main()
