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
import sasmol.sasmol as sasmol
import sasmol.sasmath as sasmath
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.monomer_monte_carlo.dihedral_rotate as nrotate
from sassie.simulate.monomer_monte_carlo.dihedral_rotate import overlap
#from sassie.simulate.complex_monte_carlo.nmer_overlap_check import moloverlap

#	NMER-ROTATE
#
#	12/17/04	--	adapted from trehalose project for gag modeling : jc
#	12/20/04	--	align N-terminal of CA in 1L6N and 1E6J		: jc
#	01/04/05	--	align 1L6N to complete bottom piece (1E6J+..) 	: jc
#	03/01/05	--	calculate dihedrals for 1L6N linker		: jc
#	09/25/05	--	rotate a dihedral angle				: jc
#	07/01/09	--	added code for complex dihedral			: jc
#	12/05/11	--	added sasmol support				: jc
#


def overlap_check(m1, frame, cut, asegs, abasis, all_segment_basis_full_mask, interatom, interres):

    # print 'ZHL overlap ',interatom
    # Check for each individual segments
    for i in range(len(asegs)):
        mask_seg = all_segment_basis_full_mask[i]
        #import pprint; pprint.pprint(mask_seg.tolist()); exit(0)
        error, coor_tmp = m1.get_coor_using_mask(frame, mask_seg)
        check = overlap.overlap(coor_tmp[0], float(cut))
        if check:
            return check
        # print 'ZHL check ', check
    # Check between segments
    for i in range(len(asegs) - 1):
        mask_seg = all_segment_basis_full_mask[i]
        error, coor_1 = m1.get_coor_using_mask(frame, mask_seg)
        for j in range(i + 1, len(asegs)):
            mask_seg = all_segment_basis_full_mask[j]
            error, coor_2 = m1.get_coor_using_mask(frame, mask_seg)
            # atomlist=moloverlap(coor_1,coor_2,float(cut),interres)
            # if (len(atomlist)):
            #	check = 1
            check = overlap.moloverlap(
                coor_1[0], coor_2[0], float(cut), numpy.array(interres))
            if check:
                return check
    '''
	for i in range(len(asegs)):
		tmask_seg = 'segname[i] == "'+asegs[i]+'" and name[i] == "'+abasis[i].strip()+'"'
		error,mask_seg = m1.get_subset_mask(tmask_seg)
		error,coor_tmp = m1.get_coor_using_mask(frame,mask_seg)
		check = overlap.overlap(coor_tmp[0],float(cut))
		#print 'ZHL check ', check
	# Check between segments
	for i in range(len(asegs)-1):
		tmask_seg = 'segname[i] == "'+asegs[i]+'" and name[i] == "'+abasis[i].strip()+'"'
		error,mask_seg = m1.get_subset_mask(tmask_seg)
		error,coor_1 = m1.get_coor_using_mask(frame,mask_seg)
		for j in range(i+1,len(asegs)):
			tmask_seg = 'segname[i] == "'+asegs[j]+'" and name[i] == "'+abasis[j].strip()+'"'
			error,mask_seg = m1.get_subset_mask(tmask_seg)
			error,coor_2 = m1.get_coor_using_mask(frame,mask_seg)
			#atomlist=moloverlap(coor_1,coor_2,float(cut),interres)
			#if (len(atomlist)):
			#	check = 1
			check=overlap.moloverlap(coor_1[0],coor_2[0],float(cut),numpy.array(interres))
	'''
    return check


def rotate(coor, m1, q0, th, an, cut, lowrg, highrg, re, taccepted, zflag, zval, cflag, dcdoutfile, indices, this_mask, basis_mask, sub_m2, align_mask, coor_sub_m1, com_sub_m1, mask_a_array, mask_b_array, distance_array, type_array, first_last_resid, molecule_type, segment_mask, segment_full_mask, all_segment_basis_full_mask, basis_full_mask, segment_mol, asegs, abasis, interatom, interres):

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

    error, new_coor = segment_mol.get_coor_using_mask(frame, segment_mask)

    segment_mol.setCoor(new_coor)

    result = nrotate.rotate_dihedral(
        new_coor, segment_mol, frame, q0, th, an, indices, this_mask, first_last_resid, molecule_type)
    # nrotate.rotate_dihedral(coor,segment_mol,frame,q0,th,an,indices,this_mask,first_last_resid,molecule_type)
    # nrotate.rotate_dihedral(coor,m1,frame,q0,th,an,indices,this_mask,first_last_resid,molecule_type)

    #error = segment_mol.set_coor_using_mask(segment_mol,frame,segment_basis_mask)

    if(error != []):
        print 'error = ', error
        sys.exit()

    #error,coor = m1.get_coor_using_mask(frame,segment_basis_mask)

#	error,basis_coor = m1.get_coor_using_mask(frame,basis_mask)
#	check=overlap.overlap(basis_coor[0],float(cut))

    thisrg = m1.calcrg(frame)

    if(thisrg > hrg):
        hrg = thisrg
    if(thisrg < lowestrg):
        lowestrg = thisrg
    filename = ''

    '''
	tmask = ''
	for i in range(len(asegs)):
		tmask+='segname[i] == "'+asegs[i]+'" and name[i] =="'+abasis[i].strip()+'"'
		if i!=len(asegs)-1:
			tmask+=' or '
	error,basis_full_mask= m1.get_subset_mask(tmask)
	'''
    if(result == 0):
        check = 1

    if(check == 0):
        if(thisrg > lowrg and thisrg < highrg):
            filename = 'winner'

            # m1.center(frame)

            # print 'sum(align_mask)= ',numpy.sum(align_mask)

            '''
			error,sub_m2.coor = m1.get_coor_using_mask(frame,align_mask)
			sub_m2.setCoor(sub_m2.coor)
			com_sub_m2 = sub_m2.calccom(0)
			sub_m2.center(0)
			coor_sub_m2 = sub_m2.coor[0]
			'''
            error, coor_tmp = segment_mol.get_coor_using_mask(
                frame, align_mask)
            sub_m2.setCoor(coor_tmp)
            com_sub_m2 = sub_m2.calccom(0)
            sub_m2.center(0)
            coor_sub_m2 = coor_tmp[0]

            '''
			print '\n'
			print 'com_sub_m1 = ',com_sub_m1
			print 'com_sub_m2 = ',com_sub_m2,'\n'

			print 'first_coor: m1.coor()[0][0] = ',m1.coor()[0][0]
			print 'last_coor: m1.coor()[0][-1] = ',m1.coor()[0][-1]
			print '>>>>> CALLING ALIGN <<<<<<'
			'''

            segment_mol.align(frame, coor_sub_m2, com_sub_m2,
                              coor_sub_m1, com_sub_m1)
            error = m1.set_coor_using_mask(
                segment_mol, frame, segment_full_mask)

            '''
			print 'first_coor: m1.coor()[0][0] = ',m1.coor()[0][0]
			print 'last_coor: m1.coor()[0][-1] = ',m1.coor()[0][-1]
			'''
            check = overlap_check(
                m1, frame, cut, asegs, abasis, all_segment_basis_full_mask, interatom, interres)

            if(check == 0 and zflag == 1):
                zee = m1.coor()[0, :, 2]
                zcheck = numpy.alltrue(numpy.greater_equal(zee, zval))
                if( not zcheck ):
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
                over = 1

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

    # align_segment(a,md1,st,tseg,abasis[tsegn],seglow[tsegn],seghigh[tsegn],zrefa1[tsegn],rcm1[tsegn],segatm_first[tsegn],segatm_last[tsegn])

    # check=nmer_filter.noverlap(a,md1,st,cut,npairs,interpairs,tseg,tsegn,asegs,abasis)

#
