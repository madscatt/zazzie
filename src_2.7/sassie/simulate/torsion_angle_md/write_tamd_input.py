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
import string
import sys

import numpy


#       WRITE_TAMD_INPUT

#

#       10/01/2007       --      initial coding                   :       jc

#       12/09/2013       --      adapted for tamd            :       jc

#

# LC      1         2         3         4         5         6         7

# LC4567890123456789012345678901234567890123456789012345678901234567890123456789

#                                                                      *      **

'''

        WRITE_TAMD_INPUT is the function to write a generic input

    file for the program TAMD run through CHARMM



    This module is called by TAMD.PY





    REFERENCE:

    

    #### NOTE that this will not work for PROLINE in setting

    ####    cluster definitions ... also need to check other

    ####    amino acid and patch definitions

    ####

'''


def get_rigid_regions(nregions, reslow, numcont, residues):

    # print '!!',nregions,",",reslow,",",numcont,"\n"

    rigid_regions = []

    flexible_regions = []

    count = 0

    for i in xrange(nregions):

        this_region = [i for i in range(reslow[i], reslow[i] + numcont[i] + 1)]

        flexible_regions.append(this_region)

    # there has to be at least ONE flexible region

    for i in xrange(nregions):

        print '!!!', i, nregions

        this_rigid_region = []

        # check and see if the last flexible region goes to the end of the
        # protein

        if((i == nregions - 1) and (flexible_regions[-1][-1] == residues[-1]) and (nregions > 1)):

            print 'last region is flexible 1'

            for residue in xrange(flexible_regions[i - 1][-1] + 1, flexible_regions[i][0]):

                print 'residue=', residue

                this_rigid_region.append(residue)

            rigid_regions.append(this_rigid_region)

            # for residue in
            # xrange(flexible_regions[i][-1]+1,flexible_regions[i+1][0]):

            #    this_rigid_region.append(residue)

            # rigid_regions.append(this_rigid_region)

        # if the first flexible region does not include the first residue

        # therefore it is a rigid region

        elif((i == 0) and (flexible_regions[0][0] != residues[0])):

            print 'first residue is rigid: will accumulate first rigid region from first residue 2'

            for residue in xrange(residues[0], flexible_regions[0][0]):

                this_rigid_region.append(residue)

            rigid_regions.append(this_rigid_region)

            if (nregions == 1):

                this_rigid_region = []

                for residue in xrange(flexible_regions[-1][-1] + 1, residues[-1] + 1):

                    this_rigid_region.append(residue)

                rigid_regions.append(this_rigid_region)

        elif(i == 0):

            print 'first residue is flexible: first rigid region starts past the first flexible region 3'

            pass

        elif((i == nregions - 1)):

            print 'last residue is rigid 4'

            for residue in xrange(flexible_regions[i - 1][-1] + 1, flexible_regions[i][0]):

                this_rigid_region.append(residue)

            rigid_regions.append(this_rigid_region)

            this_rigid_region = []

            for residue in xrange(flexible_regions[i][-1] + 1, residues[-1] + 1):

                this_rigid_region.append(residue)

            rigid_regions.append(this_rigid_region)

        else:

            print '5'

            for residue in xrange(flexible_regions[i - 1][-1] + 1, flexible_regions[i][0]):

                this_rigid_region.append(residue)

            rigid_regions.append(this_rigid_region)

        print '>>>flexible regions:', flexible_regions

        print '>>>rigid regions:', rigid_regions

    return flexible_regions, rigid_regions


def get_protein_cluster_strings(all_strings, segname, numranges, reslow, numcont, residue):
    '''

    This method takes the input regions of flexiblity and constructs the

    tamd cluster defintions that are defined by the non-flexible regions.

    '''

    low_residue = residue[0]
    high_residue = residue[-1]

    residues = [i for i in range(low_residue, high_residue + 1)]


# assumes residue numbering increases

    flexible_regions, rigid_regions = get_rigid_regions(
        numranges, reslow, numcont, residues)

    if (rigid_regions[0][0] == residues[0]):

        print '>>> first amino acid is rigid'

    else:

        print '>>> first amino acid is flexible'

    if (rigid_regions[-1][-1] == residues[-1]):

        print '>>> last amino acid is rigid'

    else:

        print '>>> last amino acid is flexible'

    # for i in xrange(len(rigid_regions)):

    #    print 'rigid regions : ',rigid_regions[i][0],':',rigid_regions[i][-1]

    # for i in xrange(numranges):

    # print 'flexible regions :
    # ',flexible_regions[i][0],':',flexible_regions[i][-1]

    # find last region type

    if (rigid_regions[-1][-1] < flexible_regions[-1][-1]):

        last = "flexible"

    else:

        last = "rigid"

    cluster_string = []

    for i in xrange(len(rigid_regions)):

        print 'len(rigid_regions)', len(rigid_regions), '\n'

        if((i == 0) and (rigid_regions[0][0] == residues[0])):

            # X ===== rigid ===== Y ----- flexible -----

            # X is the first residue

            # cluster select resid X:Y .or. -

            #    ( resid X+1 .and. ( type N .or. type H ) ) end

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            st = '\tcluster select segid ' + segname + \
                ' .and. resid ' + str(X) + ':' + str(Y) + ' .or. -\n'

            st = st + '\t\t(segid ' + segname + ' .and. resid ' + \
                str(Y + 1) + ' .and. ( type N .or. type HN ) ) end\n'

        #    print st

            cluster_string.append(st)

            all_strings.append(st)

        elif((last == "rigid") and (i == len(rigid_regions) - 1)):

            # last == rigid

            # flexible -------- X ===== rigid ===== Y

            #  cluster select resid X:Y .or. -

            #      ( resid X-1 .and. ( type C .or. type O ) ) end

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            st = '\tcluster select segid ' + segname + \
                ' .and. resid ' + str(X) + ':' + str(Y) + ' .or. -\n'

            st = st + '\t\t(segid ' + segname + ' .and. resid ' + \
                str(X - 1) + ' .and. ( type C .or. type O ) ) end\n'

        #    print st

            cluster_string.append(st)

            all_strings.append(st)

        else:

            # ----- flexible ----- X ===== rigid ===== Y ----- flexible -----

            # cluster select resid X:Y .or. -

            #    ( resid X-1 .and. ( type C .or. type O ) ) .or. -

            #    ( resid Y+1 .and. (type N .or. type H ) ) end

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            st = '\tcluster select segid ' + segname + \
                ' .and. resid ' + str(X) + ':' + str(Y) + ' .or. -\n'

            st = st + '\t\t(segid ' + segname + ' .and. resid ' + \
                str(X - 1) + ' .and. ( type C .or. type O ) ) .or. -\n'

            st = st + '\t\t(segid ' + segname + ' .and. resid ' + \
                str(Y + 1) + ' .and. ( type N .or. type HN ) ) end\n'

        #    print st

            cluster_string.append(st)

            all_strings.append(st)

    '''

#### if first residue is flexible then the first cluster begins one residue after

####    the first flexible region



    # case 1: the first flexible region does not go to the end of the protein



        # ----- flexible ----- X ===== rigid ===== Y ----- flexible ----- Z ===== rigid =====

        #

        # reslow + numcont takes you to X 



        # cluster select resid X:Y .or. -

        # ( resid X-1 .and. ( type C .or. type O ) ) .or. -

        # ( resid Y+1 .and. (type N .or. type H ) ) end



    # case 2: the first flexible region does go to the end of the protein

    

        # there are no clusters!



#### else if first residue is rigid then the first cluster begins at the first residue

####    to the first flexible region



    # case 1: the first flexible region does not go to the end of the protein

        # NOTE: if the rigid region went to the end there would not be any

        #     flexible regions.

        #    

        # X ===== rigid ===== Y ----- flexible -----

        # X is the first residue

        # cluster select resid X:Y .or. -

            #    ( resid Y+1 .and. ( type N .or. type H ) ) end



    # case 2: subsequent (middle regions and possible end) region(s) are defined as



        # ----- flexible ----- X ===== rigid ===== Y ----- flexible -----

        # cluster select resid X:Y .or. -

            #    ( resid X-1 .and. ( type C .or. type O ) ) .or. -

        #    ( resid Y+1 .and. (type N .or. type H ) ) end





    # case 3: if final region is not flexible 



        # flexible -------- X ===== rigid ===== Y

        # reslow + numcont takes you to X 

        #  cluster select resid X:Y .or. -

        #      ( resid X-1 .and. ( type C .or. type O ) ) end



#cluster select resid 1:121 .or. -

#             ( resid 122 .and. ( type N .or. type H ) ) end

#

#  cluster select resid 145:276 .or. -

#               ( resid 144 .and. ( type C .or. type O ) ) .or. -

#               ( resid 277 .and. ( type N .or. type H ) ) end

#

#  cluster select resid 284:353 .or. -

#               ( resid 283 .and. ( type C .or. type O ) ) .or. -

#               ( resid 354 .and. ( type N .or. type H ) ) end

#

#  cluster select resid 390:405 .or. -

#               ( resid 389 .and. ( type C .or. type O ) ) .or. -

#               ( resid 406 .and. ( type N .or. type H ) ) end

#

#  cluster select resid 412:431 .or. -

#               ( resid 411 .and. ( type C .or. type O ) ) end

    '''

    return


def get_cluster_strings_dna(cluster_file_name, flexible_regions, rigid_regions, resnum, resname, segres):

    # output dna clustering file

    # assuming dna pdb file is properly  renamed to match the residue name as
    # in top/par files

    clusterfile = open(cluster_file_name, 'w')

    cluster_strings = []

    cluster_strings.append("*string file clustering PDB file")

    cluster_strings.append("*\n")

    cluster_strings.append("TAMD\n")

    cluster_strings.append("! rigid regions\n")

    for i in xrange(len(rigid_regions)):

        if rigid_regions[i]:

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            # if the rigid region is across different segment, cluster as two
            # rigid bodies.

            for this_segres in segres:

                if ((X < this_segres) and (Y >= this_segres)):

                    st = '\tcluster select resid ' + \
                        str(X) + ':' + str(this_segres - 1) + ' end\n'

                    st = st + '\tcluster select resid ' + \
                        str(this_segres) + ':' + str(Y) + ' end\n'

                else:

                    st = '\tcluster select resid ' + \
                        str(X) + ':' + str(Y) + ' end'

            cluster_strings.append(st)

    cluster_strings.append("\n! flexible regions\n")

    for i in xrange(len(flexible_regions)):

        if flexible_regions[i]:

            for flexible_resnum in flexible_regions[i]:

                # find residue name

                flexible_resname = resname[resnum.index(flexible_resnum)]

                st = "\t!" + str(flexible_resname) + "," + \
                    str(flexible_resnum) + ":\n"

                st = st + "\t!sugar\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C1' .or. type H1' ) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C2' .or. type H2'' .or. type H2') end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C3' .or. type H3' ) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O3' .or. type H3T) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C4' .or. type H4' ) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C5' .or. type H5' .or. type H5'' ) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O5' .or. type H5T) end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + " .and.  type O4' end\n"

                st = st + "\t!phosphate group\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + " .and. type P end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + " .and. type O1P end\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + " .and. type O2P end\n"

                st = st + "\t!base\n"

                st = st + "\t cluster select resid " + \
                    str(flexible_resnum) + " .and. -\n"

                print "!!!!!test:", flexible_resname, "\n"

                if (flexible_resname == "ADE"):

                    st = st + \
                        "\t\t(type N1  .or. type C2 .or. type H2  .or. -\n"

                    st = st + "\t\t type N3  .or. type C4 .or. type C5  .or. -\n"

                    st = st + "\t\t type C6  .or. type N6 .or. type H61 .or. -\n"

                    st = st + "\t\t type H62 .or. type N7 .or. type C8  .or. -\n"

                    st = st + "\t\t type H8  .or. type N9 ) end\n"

                elif (flexible_resname == "GUA"):

                    st = st + \
                        "\t\t(type H1 .or. type N1  .or. type C2  .or. -\n"

                    st = st + "\t\t type N2 .or. type H21 .or. type H22 .or. -\n"

                    st = st + "\t\t type N3 .or. type C4  .or. type C5  .or. -\n"

                    st = st + "\t\t type C6 .or. type O6  .or. type N7  .or. -\n"

                    st = st + "\t\t type C8 .or. type H8  .or. type N9) end\n"

                elif (flexible_resname == "THY"):

                    st = st + "\t\t(type N1 .or. type C2 .or. type O2 .or. -\n"

                    st = st + "\t\t type N3 .or. type H3 .or. type C4 .or. -\n"

                    st = st + "\t\t type O4 .or. type C5 .or. type C5M .or. -\n"

                    st = st + "\t\t type H51 .or. type H52 .or. type H53 .or. -\n"

                    st = st + "\t\t type C6 .or. type H6 ) end\n"

                elif (flexible_resname == "CYT"):

                    st = st + \
                        "\t\t(type N1  .or. type C2  .or. type O2 .or. -\n"

                    st = st + "\t\t type N3  .or. type C4  .or. type N4 .or. -\n"

                    st = st + "\t\t type H41 .or. type H42 .or. type C5 .or. -\n"

                    st = st + "\t\t type H5  .or. type C6  .or. type H6 ) end\n"

                cluster_strings.append(st)

    cluster_strings.append("END\n")

    for i in xrange(len(cluster_strings)):

        clusterfile.write("%s\n" % cluster_strings[i])

    clusterfile.close()

    return


def get_rna_cluster_strings(all_strings, segname, numranges, reslow, numcont, residue, resnames):

    cluster_file_name = 'cluster_' + segname + '.str'

    clusterfile = open(cluster_file_name, 'w')

    low_residue = residue[0]
    high_residue = residue[-1]

    residues = [i for i in range(low_residue, high_residue + 1)]

    flexible_regions, rigid_regions = get_rigid_regions(
        numranges, reslow, numcont, residues)

    if (rigid_regions[0][0] == residues[0]):

        print '>>> first residue is rigid'

    else:

        print '>>> first residue is flexible'

    if (rigid_regions[-1][-1] == residues[-1]):

        print '>>> last residue is rigid'

    else:

        print '>>> last residue is flexible'

    if (rigid_regions[-1][-1] < flexible_regions[-1][-1]):

        last = "flexible"

    else:

        last = "rigid"

    cluster_string = []

    cluster_string.append("*string file for clustering")

    cluster_string.append("*\n")

    cluster_string.append("! rigid regions\n")

    for i in xrange(len(rigid_regions)):

        if rigid_regions[i]:

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            st = '\tcluster select segid ' + segname + \
                ' .and. resid ' + str(X) + ':' + str(Y) + ' end'

            cluster_string.append(st)

    cluster_string.append("\n! flexible regions\n")

    for i in xrange(len(flexible_regions)):

        if flexible_regions[i]:

            for flexible_resnum in flexible_regions[i]:

                # find residue name

                flexible_resname = resnames[flexible_resnum - 1]

                st = "\t!" + str(flexible_resname) + "," + \
                    str(flexible_resnum) + ":\n"

                st = st + "\t!sugar\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C1' .or. type H1' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C2' .or. type H2'') end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O2' .or. type H2') end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C3' .or. type H3' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O3' .or. type H3T) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C4' .or. type H4' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C5' .or. type H5' .or. type H5'' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O5' .or. type H5T) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and.  type O4' end\n"

                st = st + "\t!phosphate group\n"

                st = st + "\t cluster select segid " + segname + \
                    " .and. resid " + str(flexible_resnum) + \
                    " .and. type P end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and. type O1P end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and. type O2P end\n"

                st = st + "\t!base\n"

                st = st + "\t cluster select segid " + segname + \
                    " .and. resid " + str(flexible_resnum) + " .and. -\n"

                if (flexible_resname == "ADE"):

                    st = st + \
                        "\t\t(type N1  .or. type C2 .or. type H2  .or. -\n"

                    st = st + "\t\t type N3  .or. type C4 .or. type C5  .or. -\n"

                    st = st + "\t\t type C6  .or. type N6 .or. type H61 .or. -\n"

                    st = st + "\t\t type H62 .or. type N7 .or. type C8  .or. -\n"

                    st = st + "\t\t type H8  .or. type N9 ) end\n"

                elif (flexible_resname == "GUA"):

                    st = st + \
                        "\t\t(type H1 .or. type N1  .or. type C2  .or. -\n"

                    st = st + "\t\t type N2 .or. type H21 .or. type H22 .or. -\n"

                    st = st + "\t\t type N3 .or. type C4  .or. type C5  .or. -\n"

                    st = st + "\t\t type C6 .or. type O6  .or. type N7  .or. -\n"

                    st = st + "\t\t type C8 .or. type H8  .or. type N9) end\n"

                elif (flexible_resname == "URA"):

                    st = st + "\t\t(type N1 .or. type C2 .or. type O2 .or. -\n"

                    st = st + "\t\t type N3 .or. type H3 .or. type C4 .or. -\n"

                    st = st + "\t\t type O4 .or. type C5 .or. type H5 .or. -\n"

                    st = st + "\t\t type C6 .or. type H6 ) end\n"

                elif (flexible_resname == "CYT"):

                    st = st + \
                        "\t\t(type N1  .or. type C2  .or. type O2 .or. -\n"

                    st = st + "\t\t type N3  .or. type C4  .or. type N4 .or. -\n"

                    st = st + "\t\t type H41 .or. type H42 .or. type C5 .or. -\n"

                    st = st + "\t\t type H5  .or. type C6  .or. type H6 ) end\n"

                cluster_string.append(st)

    cluster_string.append("return\n")

    for i in xrange(len(cluster_string)):

        clusterfile.write("%s\n" % cluster_string[i])

    clusterfile.close()

    all_strings.append('\tstream "' + cluster_file_name + '"\n')

    return


def get_dna_cluster_strings(all_strings, segname, numranges, reslow, numcont, residue, resnames):

    cluster_file_name = 'cluster_' + segname + '.str'

    clusterfile = open(cluster_file_name, 'w')

    low_residue = residue[0]
    high_residue = residue[-1]

    residues = [i for i in range(low_residue, high_residue + 1)]

    flexible_regions, rigid_regions = get_rigid_regions(
        numranges, reslow, numcont, residues)

    if (rigid_regions[0][0] == residues[0]):

        print '>>> first residue is rigid'

    else:

        print '>>> first residue is flexible'

    if (rigid_regions[-1][-1] == residues[-1]):

        print '>>> last residue is rigid'

    else:

        print '>>> last residue is flexible'

    if (rigid_regions[-1][-1] < flexible_regions[-1][-1]):

        last = "flexible"

    else:

        last = "rigid"

    cluster_string = []

    cluster_string.append("*string file for clustering")

    cluster_string.append("*\n")

    cluster_string.append("! rigid regions\n")

    for i in xrange(len(rigid_regions)):

        if rigid_regions[i]:

            X = rigid_regions[i][0]
            Y = rigid_regions[i][-1]

            st = '\tcluster select segid ' + segname + \
                ' .and. resid ' + str(X) + ':' + str(Y) + ' end'

            cluster_string.append(st)

    cluster_string.append("\n! flexible regions\n")

    for i in xrange(len(flexible_regions)):

        if flexible_regions[i]:

            for flexible_resnum in flexible_regions[i]:

                # find residue name

                flexible_resname = resnames[flexible_resnum - 1]

                st = "\t!" + str(flexible_resname) + "," + \
                    str(flexible_resnum) + ":\n"

                st = st + "\t!sugar\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C1' .or. type H1' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C2' .or. type H2'' .or. type H2') end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C3' .or. type H3' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O3' .or. type H3T) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C4' .or. type H4' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type C5' .or. type H5' .or. type H5'' ) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + \
                    " .and. (type O5' .or. type H5T) end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and.  type O4' end\n"

                st = st + "\t!phosphate group\n"

                st = st + "\t cluster select segid " + segname + \
                    " .and. resid " + str(flexible_resnum) + \
                    " .and. type P end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and. type O1P end\n"

                st = st + "\t cluster select segid " + segname + " .and. resid " + \
                    str(flexible_resnum) + " .and. type O2P end\n"

                st = st + "\t!base\n"

                st = st + "\t cluster select segid " + segname + \
                    " .and. resid " + str(flexible_resnum) + " .and. -\n"

                if (flexible_resname == "ADE"):

                    st = st + \
                        "\t\t(type N1  .or. type C2 .or. type H2  .or. -\n"

                    st = st + "\t\t type N3  .or. type C4 .or. type C5  .or. -\n"

                    st = st + "\t\t type C6  .or. type N6 .or. type H61 .or. -\n"

                    st = st + "\t\t type H62 .or. type N7 .or. type C8  .or. -\n"

                    st = st + "\t\t type H8  .or. type N9 ) end\n"

                elif (flexible_resname == "GUA"):

                    st = st + \
                        "\t\t(type H1 .or. type N1  .or. type C2  .or. -\n"

                    st = st + "\t\t type N2 .or. type H21 .or. type H22 .or. -\n"

                    st = st + "\t\t type N3 .or. type C4  .or. type C5  .or. -\n"

                    st = st + "\t\t type C6 .or. type O6  .or. type N7  .or. -\n"

                    st = st + "\t\t type C8 .or. type H8  .or. type N9) end\n"

                elif (flexible_resname == "THY"):

                    st = st + "\t\t(type N1 .or. type C2 .or. type O2 .or. -\n"

                    st = st + "\t\t type N3 .or. type H3 .or. type C4 .or. -\n"

                    st = st + "\t\t type O4 .or. type C5 .or. type C5M .or. -\n"

                    st = st + "\t\t type H51 .or. type H52 .or. type H53 .or. -\n"

                    st = st + "\t\t type C6 .or. type H6 ) end\n"

                elif (flexible_resname == "CYT"):

                    st = st + \
                        "\t\t(type N1  .or. type C2  .or. type O2 .or. -\n"

                    st = st + "\t\t type N3  .or. type C4  .or. type N4 .or. -\n"

                    st = st + "\t\t type H41 .or. type H42 .or. type C5 .or. -\n"

                    st = st + "\t\t type H5  .or. type C6  .or. type H6 ) end\n"

                cluster_string.append(st)

    cluster_string.append("return\n")

    for i in xrange(len(cluster_string)):

        clusterfile.write("%s\n" % cluster_string[i])

    clusterfile.close()

    all_strings.append('\tstream "' + cluster_file_name + '"\n')

    return


def get_residue_sequence(mol):
    '''

    NOTE: numpy SORTS the array when unique is called : so out of order

    sequences are going to be wrong

    '''

    resid = mol.resid()

    resname = mol.resname()

    resids, unique_indices = numpy.unique(resid, return_index=True)

    resnames = numpy.take(resname, unique_indices)

    return resids, resnames


def add_dna_patches(segname_molecule, input_strings):

    this_segname = segname_molecule.segname()[0]

    resname = segname_molecule.resname()

    resids, resnames = get_residue_sequence(segname_molecule)

    for i in xrange(len(resnames)):

        this_resid = resids[i]

        if(resnames[i] == "CYT" or resnames[i] == "THY"):

            input_strings.append(
                "patch deo1 " + this_segname + " " + str(this_resid))

        elif(resnames[i] == "ADE" or resnames[i] == "GUA"):

            input_strings.append(
                "patch deo2 " + this_segname + " " + str(this_resid))

    return


def get_dna_rna_setup_strings(all_strings, segname, flex_variables, residue):

    numranges = flex_variables[0]

    reslow = flex_variables[1]

    numcont = flex_variables[2]

    flexible_regions, rigid_regions = get_rigid_regions(
        numranges, reslow, numcont, residue)

    noest = []

    st = []

    noest = 'set dist = 1.52\n'

    noest = noest + 'NOE\n'

    for i in xrange(len(flexible_regions)):

        if flexible_regions[i]:

            for flexible_resnum in flexible_regions[i]:

                st = 'delete bond select segid ' + segname + ' -\n'

                st = st + '\t.and. resid ' + \
                    str(flexible_resnum) + ' .and. type C3\' end -\n'

                st = st + '\tselect segid ' + segname + ' -\n'

                st = st + '\t.and. resid ' + \
                    str(flexible_resnum) + ' .and. type C2\' end\n\n'

                all_strings.append(st)

                noest = noest + '\tassign select segid ' + segname + ' -\n'

                noest = noest + '\t.and. resid ' + \
                    str(flexible_resnum) + ' .and. type C3\' end -\n'

                noest = noest + '\tselect segid ' + segname + ' -\n'

                noest = noest + '\t.and. resid ' + \
                    str(flexible_resnum) + ' .and. type C2\' end -\n'

                noest = noest + '\tkmax 100 fmax 150 rmax @dist rmin @dist\n\n'

    noest = noest + 'END\n'

    all_strings.append(noest)

    return


def get_all_cluster_strings(all_strings, segname_molecules, flexible_segment_variables):

    print '\n\n'

    print '>>>>>>>>>>>>>>>>>>>> cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    print '>>>>>>>>>>>>>>>>>>>> cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    print '>>>>>>>>>>>>>>>>>>>> cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    print '\n\n'

    for i in xrange(len(segname_molecules)):

        this_moltype = segname_molecules[i].moltype()[0]

        this_segname = segname_molecules[i].segname()[0]

        resid, resnames = get_residue_sequence(segname_molecules[i])

        if(this_moltype == 'protein' and this_segname in flexible_segment_variables):

            print '\n>> assigning protein cluster strings\n'

            print '>>>> this segname = ', this_segname, ' <<<<\n'

            flex_variables = flexible_segment_variables[this_segname]

            residue = segname_molecules[i].resid()

            numranges = flex_variables[0]

            reslow = flex_variables[1]

            numcont = flex_variables[2]

            print 'reslow=', reslow, ' numcont=', numcont, ' numranges=', numranges, '\n'

            get_protein_cluster_strings(
                all_strings, this_segname, numranges, reslow, numcont, residue)

        elif(this_moltype == 'rna' and this_segname in flexible_segment_variables):

            print '\n>> assigning rna cluster strings\n'

            print '>>>> this segname = ', this_segname, ' <<<<\n'

            flex_variables = flexible_segment_variables[this_segname]

            residue = segname_molecules[i].resid()

            numranges = flex_variables[0]

            reslow = flex_variables[1]

            numcont = flex_variables[2]

            # default RNA clustering for now --> HACK for debugging

            #all_strings.append('\tcluster select segid'+this_segname+' end\n')

            get_rna_cluster_strings(
                all_strings, this_segname, numranges, reslow, numcont, residue, resnames)

        elif(this_moltype == 'dna' and this_segname in flexible_segment_variables):

            print '\n>> assigning dna cluster strings\n'

            print '>>>> this segname = ', this_segname, ' <<<<\n'

            flex_variables = flexible_segment_variables[this_segname]

            residue = segname_molecules[i].resid()

            numranges = flex_variables[0]

            reslow = flex_variables[1]

            numcont = flex_variables[2]

            #all_strings.append('\tcluster select segid '+this_segname+' end\n')

            get_dna_cluster_strings(
                all_strings, this_segname, numranges, reslow, numcont, residue, resnames)

        else:

            all_strings.append('\tcluster select segid ' +
                               this_segname + ' end\n')

    print '\n\n'

    print '>>>>>>>>>>>>>>>>>>>> end cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<'

    print '>>>>>>>>>>>>>>>>>>>> end cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<'

    print '>>>>>>>>>>>>>>>>>>>> end cluster string debugging <<<<<<<<<<<<<<<<<<<<<<<<'

    print '\n\n'

    return


def write_tamd_input(m1, segname_molecules, input_file_name, number_of_steps, dcdfreq, temp_pdb_files, topology_file_name, parameter_file_name, output_dcd_file_name, temperature, rgforce, rgvalue, residue, flexible_segment_variables, pretamd_min_steps):
    '''

    WRITE_TAMD_INPUT is the function to write a generic input

    file for the program TAMD run through CHARMM





    INPUT:  variable descriptions:





            number_of_steps:    number of dynamics steps

            dcdfreq:        frequency to write dcd file to disk



            temp_pdb_files:        name of pdb files for base structure to perform dynamics_name

            topology_file_name:    path and name of topology file

            parameter_file_name:    path and name of parameter file

            output_dcd_file_name:    name of final structure dcd file



            pretamd_min_steps:    number of energy minization steps to carry out before starting tamd    



            temperature:        temperature





    OUTPUT:

            input_file_name:          name of input file that function will write





    NOTE:  this has been tested for single segment protein molecules only



    '''

    outfile = open(input_file_name, 'w')

    auto_directed_rg = False

    directed_rg = False

    if rgforce > 0.0:

        if rgvalue > 0.0:

            directed_rg = True

        else:

            auto_directed_rg = True

    all_strings = []

    all_strings.append("open read card unit 1 name " +
                       topology_file_name + "\n")

    all_strings.append("read rtf card unit 1\nclose unit 1\n")

    all_strings.append("open read card unit 2 name " +
                       parameter_file_name + "\n")

    all_strings.append("read para card unit 2\nclose unit 2")

    # generate each segment separately

    base_unit = 10

    dnarnaflag = 0

    for i in xrange(len(temp_pdb_files)):

        pdb_file_name = temp_pdb_files[i]

        this_segname = segname_molecules[i].segname()[0]

        this_moltype = segname_molecules[i].moltype()[0]

        residue = segname_molecules[i].resid()

        all_strings.append("read sequence pdb name " +
                           pdb_file_name + " unit " + str(base_unit) + "\n")

        if(this_moltype == 'protein'):

            all_strings.append("generate " + this_segname + " setup warn\n")

        elif(this_moltype == 'rna'):

            if (this_segname in flexible_segment_variables):

                dnarnaflag = 1

            all_strings.append("generate " + this_segname +
                               " setup warn first 5TER last 3TER\n")

        elif(this_moltype == 'dna'):

            if (this_segname in flexible_segment_variables):

                dnarnaflag = 1

            all_strings.append("generate " + this_segname +
                               " setup warn first 5TER last 3TER\n")

            add_dna_patches(segname_molecules[i], all_strings)

        all_strings.append("rewind unit " + str(base_unit) + "\n")

        if(i == 0):

            all_strings.append("read coor pdb name " + pdb_file_name +
                               " sele segid " + this_segname + " END\n")

        else:

            all_strings.append("read coor pdb name " + pdb_file_name +
                               " sele segid " + this_segname + " END APPEND\n")

        all_strings.append("close unit " + str(base_unit) + "\n")

        if(this_moltype == 'dna' or this_moltype == 'rna' and this_segname in flexible_segment_variables):

            get_dna_rna_setup_strings(
                all_strings, this_segname, flexible_segment_variables[this_segname], residue)

        base_unit += 1

        #all_strings.append("auto angle dihe")

    if (dnarnaflag == 1):

        all_strings.append("block")

        all_strings.append("\tclear")

        all_strings.append("end")

        all_strings.append("block 2")

        all_strings.append("\tcall 2 select type C3' .or. type C2' end")

        all_strings.append("\tcoef 1 1 1.0")

        all_strings.append("\tcoef 2 2 0.0")

        all_strings.append("\tcoef 1 2 1.0")

        all_strings.append("end\n")

    all_strings.append("ic param\n")

    all_strings.append("set tmpNIC ?NIC\n")

    all_strings.append(
        "coor copy comp\nic build comp\ncoor copy select .not hydrogen end")

    all_strings.append(
        "hbuild atom cdie eps 80 cutnb 10.0 ctofnb 7.5 ctonnb 6.5 shift vshift bygr\n")

    all_strings.append("faster on\n")

    all_strings.append(
        "update atom rdie eps 4 cutnb 20 ctofnb 18 ctonnb 16 shift vshift bygr\n")

    all_strings.append("nbond elec switch rdie\n")

    all_strings.append('\nopen unit 1 write form name tamd_output.psf\n')

    all_strings.append('write psf card unit 1\n* test\n*\nclose unit 1\n')

    all_strings.append('\nopen unit 2 write form name tamd_output.pdb\n')

    all_strings.append('write coor pdb unit 2\n* test\n*\nclose unit 2\n')

    all_strings.append("energy\n")

    all_strings.append("\ncoor stat mass sele all end\n")

    all_strings.append("set xcm = ?XAVE\n")

    all_strings.append("set ycm = ?YAVE\n")

    all_strings.append("set zcm = ?ZAVE\n")

    all_strings.append("MMFP\n")

    all_strings.append(
        "\tGEO sphere rcm force 5.0 droff 0.0 xref @xcm yref @ycm zref @zcm -\nsele all end\n")

    all_strings.append("END\n")

    if auto_directed_rg:
        all_strings.append('set rgforce = ' + str(rgforce) + '\n')
        all_strings.append('coor rgyr select all end\n')
        all_strings.append('set rgyr = ?rgyr\n')
        all_strings.append('rgyr force @rgforce refe @rgyr select all end\n')
    elif directed_rg:
        all_strings.append('set rgforce = ' + str(rgforce) + '\n')
        all_strings.append('set rgyr = ' + str(rgvalue) + '\n')
        all_strings.append('rgyr force @rgforce refe @rgyr select all end\n')

    all_strings.append("\ntamd\n")

    all_strings.append("\treset\n")

    all_strings.append("end\n")

    all_strings.append("mini sd   nstep " +
                       pretamd_min_steps + " nprint 10 step 0.01\n")

    all_strings.append("mini abnr nstep " +
                       pretamd_min_steps + " nprint 10 step 0.01\n")

    all_strings.append("tamd\n")

    all_strings.append("\tbomlev -1")

    get_all_cluster_strings(
        all_strings, segname_molecules, flexible_segment_variables)

    #all_strings.append("\ttree setup topv 19\n")

    all_strings.append("\ttree setup")

    all_strings.append("\tbomlev 0\n")

    all_strings.append("\ttree print\n")

    all_strings.append("\ttree check\n")

    all_strings.append("\topen write unit 101 card name tamd.tree")

    all_strings.append("\ttree write unit 101\n\tclose unit 101\n")

    all_strings.append("\topen write unit 131 file name " +
                       output_dcd_file_name + "\n")

    all_strings.append("\topen write unit 130 card name tamd_loops.rst\n")


#    all_strings.append("\tmini nstep "+pretamd_min_steps+"\n")

#    all_strings.append("\tmini nstep 2000\n")

    #all_strings.append("\tmini abnr nstep 500 nprint 50 tolg 0.01\n")

    all_strings.append("\tcoor rms\n")


#    open unit 2 write form name min_c.pdb

#    write coor pdb unit 2

#    * test

#    *

#    close unit 2

    all_strings.append("\tdyna start echeck 200 -")

    all_strings.append("\t\tnstep " + str(number_of_steps) + " timestep 0.002 qref 20 tref " +
                       str(temperature) + " first " + str(temperature) + " -")

    #all_strings.append("\t\tnsavc 500 nprint 10000 iprfrq 100000 nsavv 0 isvfrq 2000 -")

    all_strings.append("\t\tnsavc " + str(dcdfreq) +
                       " nprint 10000 iprfrq 100000 nsavv 0 isvfrq 2000 -")

    all_strings.append("\t\tiunrea -29 iunwri 130 iuncrd 131 iunvel -1 -")

    all_strings.append("\t\tntrfrq 5000 iasors 1")

    all_strings.append("end\nstop")

    for i in xrange(len(all_strings)):

        outfile.write("%s\n" % all_strings[i])

    outfile.close()

    return


if __name__ == '__main__':

    print 'running as a main process'

    write_tamd_input()

else:

    print 'running as a spawned process'

