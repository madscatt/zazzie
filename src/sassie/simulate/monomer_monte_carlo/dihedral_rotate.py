'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import sys
import numpy
import math
import sasmol.linear_algebra as linear_algebra
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.monomer_monte_carlo.vdw_overlap as vdw_overlap
import sassie.simulate.monomer_monte_carlo.overlap as overlap

#	rotate dihedral (protein)
#
#	12/17/04	--	adapted from trehalose project for gag modeling : jc
#	12/20/04	--	align N-terminal of CA in 1L6N and 1E6J		: jc
#	01/04/05	--	align 1L6N to complete bottom piece (1E6J+..) 	: jc
#	03/01/05	--	calculate dihedrals for 1L6N linker		: jc
#	09/25/05	--	rotate a dihedral angle				: jc
#	01/12/11	--	added sasmol support				: jc
#	12/01/11	--	added nucleic acid move set			: jc
#


def measure(coor, indices, an, this_mask, q0, first_last_resid, molecule_type):

    ind = numpy.nonzero(this_mask * numpy.arange(1, len(this_mask) + 1))[0]
    this_frame_coor = coor[0, :, :]

    lcoor = numpy.take(this_frame_coor[:, :], ind, 0)

    error = []

    if(molecule_type == 'protein'):

        if(an == 'phi'):
            angle = linear_algebra.dihedral_angle(
                lcoor[0, :], lcoor[1, :], lcoor[2, :], lcoor[3, :])
        elif(an == 'psi'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[0, :], lcoor[1, :], lcoor[2, :], lcoor[3, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[1, :], lcoor[2, :], lcoor[3, :], lcoor[4, :])
        else:
            angle = 0.0
            message = 'error (in rotate.measure): no phi/psi angle specified'
            error.append(message)
            sys.exit()
#
# OPEN	Stub for error handling
#
    elif(molecule_type == 'rna'):

        if(an == 'alpha'):
            angle = linear_algebra.dihedral_angle(
                lcoor[0, :], lcoor[1, :], lcoor[2, :], lcoor[3, :])
        elif(an == 'beta'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[0, :], lcoor[1, :], lcoor[2, :], lcoor[3, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[1, :], lcoor[2, :], lcoor[3, :], lcoor[4, :])
        elif(an == 'gamma'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[1, :], lcoor[2, :], lcoor[3, :], lcoor[4, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[2, :], lcoor[3, :], lcoor[4, :], lcoor[5, :])
        elif(an == 'delta'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[2, :], lcoor[3, :], lcoor[4, :], lcoor[5, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[3, :], lcoor[4, :], lcoor[5, :], lcoor[6, :])
        elif(an == 'epsilon'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[3, :], lcoor[4, :], lcoor[5, :], lcoor[6, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[4, :], lcoor[5, :], lcoor[6, :], lcoor[7, :])
        elif(an == 'eta'):
            if(q0 == first_last_resid[0]):
                angle = linear_algebra.dihedral_angle(
                    lcoor[4, :], lcoor[5, :], lcoor[6, :], lcoor[7, :])
            else:
                angle = linear_algebra.dihedral_angle(
                    lcoor[5, :], lcoor[6, :], lcoor[7, :], lcoor[8, :])
        else:
            angle = 0.0
            message = 'error (in rotate.measure): no alpha/beta/delta/epsilon/eta angle specified'
            error.append(message)
            sys.exit()

    return angle


def rotate_dihedral(coor, m1, frame, q0, itheta, an, indices, this_mask, first_last_resid, molecule_type):

    result = 1

    theta = itheta * (math.pi / 180.0)

    this_frame_coor = coor[0, :, :]
    lcoor = numpy.take(this_frame_coor[:, :], indices, 0)

    v = numpy.zeros(3, numpy.float)
    tee = numpy.identity(4, numpy.float)

    if(molecule_type == 'protein'):

        if(q0 == first_last_resid[0]):
            n1 = lcoor[0, :]
            ca1 = lcoor[1, :]
            c1 = lcoor[2, :]
            n2 = lcoor[3, :]
        elif(q0 == first_last_resid[1]):
            c0 = lcoor[0, :]
            n1 = lcoor[1, :]
            ca1 = lcoor[2, :]
            c1 = lcoor[3, :]
        else:
            c0 = lcoor[0, :]
            n1 = lcoor[1, :]
            ca1 = lcoor[2, :]
            c1 = lcoor[3, :]
            n2 = lcoor[4, :]

        if(an == 'phi'):
            v[0] = ca1[0] - n1[0]
            v[1] = ca1[1] - n1[1]
            v[2] = ca1[2] - n1[2]
            tee[0][3] = n1[0]
            tee[1][3] = n1[1]
            tee[2][3] = n1[2]
            watch_value = indices[2] - 1

        elif(an == 'psi'):
            v[0] = c1[0] - ca1[0]
            v[1] = c1[1] - ca1[1]
            v[2] = c1[2] - ca1[2]
            tee[0][3] = ca1[0]
            tee[1][3] = ca1[1]
            tee[2][3] = ca1[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[2] - 1
            else:
                watch_value = indices[3] - 1

        try:
            lv = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            if(lv > 0.0):
                nv = v / lv
            else:
                result = 0
                print('v = ', v)
                print('lv = ', lv)
                return result
        except:
            result = 0
            print('v = ', v)
            print('lv = ', lv)
            return result

    elif(molecule_type == 'rna'):

        if(q0 == first_last_resid[0]):
            p1 = lcoor[0, :]
            o5p = lcoor[1, :]
            c5p = lcoor[2, :]
            c4p = lcoor[3, :]
            c3p = lcoor[4, :]
            o3p = lcoor[5, :]
            p2 = lcoor[6, :]
            o5p2 = lcoor[7, :]
        elif(q0 == first_last_resid[1]):
            o3p0 = lcoor[0, :]
            p1 = lcoor[1, :]
            o5p = lcoor[2, :]
            c5p = lcoor[3, :]
            c4p = lcoor[4, :]
            c3p = lcoor[5, :]
            o3p = lcoor[6, :]
        else:

            o3p0 = lcoor[0, :]
            p1 = lcoor[1, :]
            o5p = lcoor[2, :]
            c5p = lcoor[3, :]
            c4p = lcoor[4, :]
            c3p = lcoor[5, :]
            o3p = lcoor[6, :]
            p2 = lcoor[7, :]
            o5p2 = lcoor[8, :]

        if(an == 'alpha'):
            v[0] = o5p[0] - p1[0]
            v[1] = o5p[1] - p1[1]
            v[2] = o5p[2] - p1[2]
            tee[0][3] = p1[0]
            tee[1][3] = p1[1]
            tee[2][3] = p1[2]
            watch_value = indices[2] - 1
        elif(an == 'beta'):
            v[0] = c5p[0] - o5p[0]
            v[1] = c5p[1] - o5p[1]
            v[2] = c5p[2] - o5p[2]
            tee[0][3] = o5p[0]
            tee[1][3] = o5p[1]
            tee[2][3] = o5p[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[2] - 1
            else:
                watch_value = indices[3] - 1
        elif(an == 'gamma'):
            v[0] = c4p[0] - c5p[0]
            v[1] = c4p[1] - c5p[1]
            v[2] = c4p[2] - c5p[2]
            tee[0][3] = c5p[0]
            tee[1][3] = c5p[1]
            tee[2][3] = c5p[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[3] - 1
            else:
                watch_value = indices[4] - 1
        elif(an == 'delta'):
            v[0] = c3p[0] - c4p[0]
            v[1] = c3p[1] - c4p[1]
            v[2] = c3p[2] - c4p[2]
            tee[0][3] = c4p[0]
            tee[1][3] = c4p[1]
            tee[2][3] = c4p[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[4] - 1
            else:
                watch_value = indices[5] - 1
        elif(an == 'epsilon'):
            v[0] = o3p[0] - c3p[0]
            v[1] = o3p[1] - c3p[1]
            v[2] = o3p[2] - c3p[2]
            tee[0][3] = c3p[0]
            tee[1][3] = c3p[1]
            tee[2][3] = c3p[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[5] - 1
            else:
                watch_value = indices[6] - 1
        elif(an == 'eta'):
            v[0] = p2[0] - o3p[0]
            v[1] = p2[1] - o3p[1]
            v[2] = p2[2] - o3p[2]
            tee[0][3] = o3p[0]
            tee[1][3] = o3p[1]
            tee[2][3] = o3p[2]
            if(q0 == first_last_resid[0]):
                watch_value = indices[6] - 1
            else:
                watch_value = indices[7] - 1

        try:
            lv = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            if(lv > 0.0):
                nv = v / lv
            else:
                result = 0
                print('>>> ')
                print('v = ', v)
                print('lv = ', lv)
                print('>>> ')
                print('>>> ')
                print('an = ', an)
                print('q0 = ', q0)
                print('>>> ')
                print('p1 = ', p1)
                print('o5p = ', o5p)
                print('c5p = ', c5p)
                print('c4p = ', c4p)
                print('c3p = ', c3p)
                print('o3p = ', o3p)
                print('p2 = ', p2)
                print('o5p2 = ', o5p2)
                print('>>> ')
                print('>>> ')
                return result
        except:
            result = 0
            print('v = ', v)
            print('lv = ', lv)

            print('>>> ')
            print('an = ', an)
            print('q0 = ', q0)
            print('>>> ')
            return result

    vx = nv[0]
    vy = nv[1]
    vz = nv[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    cot = math.cos(theta)
    sit = math.sin(theta)
    r = numpy.zeros((4, 4), numpy.float)

    r[0][0] = cot + (1.0 - cot) * vx2
    r[0][1] = (1.0 - cot) * vx * vy - sit * vz
    r[0][2] = (1.0 - cot) * vx * vz + sit * vy
    r[0][3] = 0.0
    r[1][0] = (1.0 - cot) * vx * vy + sit * vz
    r[1][1] = cot + (1.0 - cot) * vy2
    r[1][2] = (1.0 - cot) * vy * vz - sit * vx
    r[1][3] = 0.0
    r[2][0] = (1.0 - cot) * vx * vz - sit * vy
    r[2][1] = (1.0 - cot) * vy * vz + sit * vx
    r[2][2] = cot + (1.0 - cot) * vz2
    r[2][3] = 0.0
    r[3][0] = 0.0
    r[3][1] = 0.0
    r[3][2] = 0.0
    r[3][3] = 1.0

    itee = numpy.matrix(tee).I
    ir = numpy.matrix(r) * itee
    br = tee * ir

    natoms = m1.natoms()

    qsel = numpy.arange(watch_value + 1, natoms)

    newa = []
    for i in range(natoms):
        if(i > watch_value):
            ta = coor[frame][i]
            far = ta.tolist()
            far.append(1.0)
            nta = numpy.array(far)
            nta.shape = (-1, 4)
            p = br * numpy.matrix(numpy.transpose(nta))
            p = numpy.array(p)
            newa.append(numpy.transpose(p))
    ncoords = numpy.array(newa)
    ncoords.shape = (-1, 4)
    temp = numpy.take(ncoords, (0, 1, 2), -1)

    totsnz1 = numpy.zeros(len(temp) * 3, numpy.int32)
    lts = 0
    for ts in range(len(temp)):
        totsnz1[lts] = (qsel[ts] * 3)
        totsnz1[lts + 1] = (qsel[ts] * 3) + 1
        totsnz1[lts + 2] = (qsel[ts] * 3) + 2
        lts = lts + 3

    numpy.put(coor[frame], totsnz1, temp)

    return result


def rotate(coor, m1, q0, th, an, cut, lowrg, highrg, re, taccepted, zflag, zval, cflag, dcdoutfile, indices, this_mask, basis_mask, sub_m2, align_mask, coor_sub_m1, com_sub_m1, mask_a_array, mask_b_array, distance_array, type_array, first_last_resid, molecule_type, cutoff_array, vdw_factor):

    over = 0
    badrg = 0
    accepted = 0
    arg = 0.0
    lowestrg = re[5]
    hrg = re[6]
    frame = 0
    badz = 0
    badc = 0

    check = 0
    minmax = []

    result = rotate_dihedral(coor, m1, frame, q0, th, an,
                             indices, this_mask, first_last_resid, molecule_type)

    if result == 1:
        error, basis_coor = m1.get_coor_using_mask(frame, basis_mask)
        print('type(basis_coor) = ', type(basis_coor))
        print('type(basis_coor[0]) = ', type(basis_coor[0]))
        print('type(cut) = ', type(cut))
        print('basis_coor[0] = ', basis_coor[0])
        print('basis_coor[0][0] = ', basis_coor[0][0])
        print('basis_coor[0][1] = ', basis_coor[0][1])
        print('cut = ', cut)
        
        if(len(cutoff_array) > 0):
            check = vdw_overlap.overlap(
                basis_coor[0], cutoff_array, vdw_factor)
        else:
            check = overlap.overlap(basis_coor[0], float(cut))
        filename = ''
    else:
        filename = ''
        check = 1

    #thisrg = m1.calcrg(frame)
    thisrg = m1.calculate_radius_of_gyration(frame)
    if(thisrg > hrg):
        hrg = thisrg
    if(thisrg < lowestrg):
        lowestrg = thisrg

    if(check == 0):
        if(thisrg > lowrg and thisrg < highrg):
            filename = 'winner'

            m1.center(frame)

            error, sub_m2.coor = m1.get_coor_using_mask(frame, align_mask)

            sub_m2.setCoor(sub_m2.coor)

            com_sub_m2 = sub_m2.calccom(0)
            sub_m2.center(0)
            coor_sub_m2 = sub_m2.coor[0]

            m1.align(frame, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)

            if(zflag == 1):
                error, sub_m2.coor = m1.get_coor_using_mask(frame, basis_mask)
                sub_m2.setCoor(sub_m2.coor)
                zee = sub_m2.coor[0, :, 2]
                zcheck = numpy.alltrue(numpy.greater_equal(zee, zval))
                if(zcheck == 0):
                    check = 1
                    badz = 1

            if(check == 0 and cflag == 1):
                check = constraints.check_constraints(
                    m1, mask_a_array, mask_b_array, distance_array, type_array)
                if(check == 1):
                    badc = 1

            if(check == 0):
                m1.write_dcd_step(dcdoutfile, 0, taccepted + 1)
                minmax = m1.calcminmax_frame(0)
                accepted = 1
                arg = thisrg
        else:
            badrg = 1

    else:
        over = 1

    re[0] = accepted
    re[1] = over
    re[2] = badrg
    re[3] = thisrg
    re[4] = arg
    re[5] = lowestrg
    re[6] = hrg
    re[7] = badz
    re[8] = badc
    re[9] = minmax

    return filename
