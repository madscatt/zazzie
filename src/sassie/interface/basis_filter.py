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
import os
import sys
import string
import locale
import numpy
import sasmol.sasmol as sasmol


def check_basis(pdbfile, basis_string, residue, num_basis):
    error = []
    mask_array = []

    supported_basis = {}
    supported_basis['calpha'] = ('CA',)
    supported_basis['backbone'] = ('N', 'CA', 'C')
    supported_basis['all'] = ('ALLATOMS',)

# parse the supplied basis which turns it into an array

    basis = string.split(basis_string, ',')

# check and see if the number of values found matches the number of basis
# definitions that are needed

    if(len(basis) != num_basis):
        error.append(
            'incorrect number of basis values entered : ' + basis_string)
        return error, mask_array

    try:
        m1 = sasmol.SasMol(0)
        m1.read_pdb(pdbfile)

    except:
        error.append('unable to read file : ' +
                     pdbfile + ' while checking basis')
        return error, mask_array

    segmentlist = []
    segments = m1.segname()
    for seg in segments:
        if seg not in segmentlist:
            segmentlist.append(seg)

    print('segmentlist= ', segmentlist)

# check and see if a keyword was supplied, if this passes then you get the
# atom array for that segment

    segment_basis = []

    for key in supported_basis:
        print('key = ', key)
        for i in range(len(segmentlist)):
            if key == basis[i]:
                segment_basis.append(supported_basis[key])

    print('segment_basis = ', segment_basis)

    print('len(segment_basis) = ', len(segment_basis))

    if(len(segment_basis) > 0 and len(segment_basis) != num_basis):
        error.append(
            'can not use a mixture of keyword basis definitons and atom names')
        return error, mask_array

    elif(len(segment_basis) == 0):
        for i in range(len(segmentlist)):
            segment_basis.append((basis[i],))

        print('atom based basis = ', segment_basis)
        print('keyword based basis = ', segment_basis)


# check and see if the atoms requested for each basis exist over the
# ranges requested

    if(len(residue) != num_basis):
        error.append(
            'number of alignment ranges does not match the number of basis ranges')
        return error, mask_array

    resid = m1.resid()
    name = m1.name()
    natoms = m1.natoms()
    print('natoms = ', natoms)
    goal_count = 0
    if(basis[0] == 'all'):
        goal_count = 'all'
    elif(basis[0] == 'backbone'):
        goal_count = 3
    elif(basis[0] == 'calpha'):
        goal_count = 1

    total_atoms_in_basis = 0

    '''
#
# the following section goes through each basis residue range set and searches for the atoms
# that should be included in the basis (first check and see if the rediue is in the residue
# list then check if that atom is in the list.  If so, count the atom and then set a counter
# flag to 1 in the "preliminary_mask_array" for that index (i.e. if it set to 1 then keep that
# atom later on in the module.  Each preliminary_mask_array is then written into the overall
# numpy mask_array to hold the mask information for each segment.  
#
	'''

    mask_array = numpy.zeros((num_basis, natoms), numpy.int32)

    segname = m1.segname()

    unique_segment_list = []
    natoms_per_segment = []
    natoms_last_segment = 0
    for i in range(natoms):
        this_segment = segname[i]
        if this_segment not in unique_segment_list:
            unique_segment_list.append(this_segment)
            if(i != 0):
                natoms_per_segment.append((i + 1) - natoms_last_segment)
                natoms_last_segment = i + 1
    natoms_per_segment.append((i + 1) - natoms_last_segment)

    print('unique_segment_list = ', unique_segment_list)
    print('natoms_per_segment = ', natoms_per_segment)

    if(num_basis != len(unique_segment_list)):
        error.append(
            'number of basis ranges does not equal the number of unique segment names in pdb file')
        return error, mask_array

    for i in range(num_basis):
        thislow = residue[i][0]
        thishigh = residue[i][1]
        print('this low = ', thislow)
        print('this high = ', thishigh)
        newresid = 1
        preliminary_mask_array = []
        for j in range(natoms_per_segment[i]):
            thisres = resid[j]
            thisatom = name[j]
            if (j == 0):
                lastres = thisres - 1
                if thisres != 1:
                    error.append(
                        'the first residue needs to be 1 for each segment in PDB file')
                    return error, mask_array

            if (thisres >= thislow and thisres <= thishigh):
                if(thisres != lastres):
                    local_count = 0
                    lastres = thisres
                if(goal_count == 'all'):
                    local_count += 1
                    total_atoms_in_basis += 1
                    preliminary_mask_array.append(1)
                elif(thisatom in segment_basis[i]):
                    local_count += 1
                    total_atoms_in_basis += 1
                    preliminary_mask_array.append(1)

                else:
                    preliminary_mask_array.append(0)
            else:
                preliminary_mask_array.append(0)

        if(numpy.sum(preliminary_mask_array) < 3):
            error.append(
                'you need to specify at least 3 atoms to form a basis for segment = ' + unique_segment_list[i])
            return array, mask_array

        print('sum pa = ', numpy.sum(preliminary_mask_array))

        mask_array[i] = preliminary_mask_array

    print('mask_array = ', mask_array)
    print('len(mask_array[0]) = ', len(mask_array[0]))

    print('total_atoms_in_basis = ', total_atoms_in_basis)

    coor = m1.coor()

    print('coor[0,0,0] = ', coor[0, 0, 0])
    print('coor[0,0] = ', coor[0, 0])
    print('coor[0] = ', coor[0])

    print('len(residue) = ', len(residue))

    for i in range(num_basis):
        pass

    # clearly this is not done

    #	check the atoms in each range, create a numpy array for each
    #	return an error if zero atoms are found per segment (< 3 atoms??)
    #	the numpy array needs to have the right dimensions to use as take/put
    #	for use later.  also make sure that this array gets packaged and
    #	sent back to the module
    #
    #	test this with protein, dimer, mcm, integrase, complex, complex+dna, zn
    #

    return error, mask_array

if __name__ == '__main__':

    pdbfile = '/home/curtisj/svn/smol/coord_files/pdbfiles/min3.pdb'
    pdbfile = '/home/curtisj/svn/sassie_version_1_svn/test_sasmol/protein0.pdb'
    basis_string = 'calpha'
    basis_string = 'backbone'
    basis_string = 'all,all,all,all,all,all'
    residue = [[1, 48]]
    residue = [[1, 48], [1, 48], [1, 48], [1, 48], [1, 48], [1, 48]]
    num_basis = 6
    error, mask_array = check_basis(pdbfile, basis_string, residue, num_basis)
    print(error, mask_array)
