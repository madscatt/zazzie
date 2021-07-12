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
import bisect
import random
import time
import platform
import numpy
import random
import Gnuplot
import Gnuplot.PlotItems
import Gnuplot.funcutils
import sassie.sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.monte_carlo.monomer.dihedral_monte_carlo as dihedral_monte_carlo
import sassie.simulate.monte_carlo.monomer.dihedral_rotate as dihedral_rotate
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monte_carlo.monomer.pairs as pairs
import sassie.simulate.monte_carlo.monomer.step as step

#import sassie.simulate.monte_carlo.complex.nmer_overlap_check as nmer_overlap_check
#import sassie.simulate.monte_carlo.complex.nmer_nrotate as nmer_nrotate
import nmer_overlap_check
import nmer_nrotate

#       NMER_DIHEDRAL
#
#	09/26/05	--    gag-dihedral search		:     jc
#	11/19/05	--	gag-dimer dihedral search	:	jc
#	06/29/09	--	generalized to nmer		:	jc/sr
#	11/17/11	--	added sasmol support		:	jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	NMR_DIHEDRAL is the module that contains the functions
	that are used to generate ensembles of structures by varying
	protein dihedral angles.  This particular version allows multiple
	flexible proteins in the presence of non-flexible proteins and
	nucleic acids.

	This module calls to C / Python extension modules to speed up
	calculations.
'''


def unpack_variables(variables):

    runname = variables['runname'][0]
    dcdfile = variables['dcdfile'][0]
    path = variables['path'][0]
    pdbfile = variables['pdbfile'][0]
    trials = variables['trials'][0]
    goback = variables['goback'][0]
    temp = variables['temp'][0]
    nsegments = variables['nsegments'][0]
    segbasis = variables['segbasis'][0]
    npsegments = variables['npsegments'][0]
    flpsegname = variables['flpsegname'][0]
    sseglow = variables['seglow'][0]
    sseghigh = variables['seghigh'][0]
    #cutoff		= variables['cutoff'][0]
    lowrg = variables['lowrg'][0]
    highrg = variables['highrg'][0]
    zflag = variables['zflag'][0]
    zcutoff = variables['zcutoff'][0]
    cflag = variables['cflag'][0]
    confile = variables['confile'][0]
    plotflag = variables['plotflag'][0]
    directedmc = variables['directedmc'][0]
    seed = variables['seed'][0]

    return runname, dcdfile, path, pdbfile, trials, goback, temp, nsegments, segbasis, npsegments, flpsegname, sseglow, sseghigh, lowrg, highrg, zflag, zcutoff, cflag, confile, plotflag, directedmc, seed


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def wait(sti=None, prompt='Plot will clear in 2 seconds ...\n'):
    '''
    WAIT is the function to prompt the user to clear a plot on a screen
    '''
    if sti is not None:
        print sti
    try:
        if(platform.system() == "Linux"):
            import curses
            stdscr = curses.initscr()
            stdscr.addstr('press a key to continue')
            c = stdscr.getch()
            curses.endwin()
    except:
        time.sleep(1)


def alignment_initialization(all_segment_mol, asegs, abasis, flexible_segments, seglow, seghigh):

    all_flexible_align_mask = []

    all_flexible_coor_sub_m1 = []

    all_flexible_com_sub_m1 = []

    all_flexible_sub_m2 = []

    for i in xrange(len(flexible_segments)):
        this_segment = flexible_segments[i]
        idx = asegs.index(this_segment)

        m1 = all_segment_mol[idx]

        if(m1.moltype()[0] == 'protein'):
            this_basis = 'CA'
        elif(m1.moltype()[0] == 'rna' or m1.moltype()[0] == 'dna'):
            this_basis = 'P'
        else:
            print 'NO ALIGNMENT BASIS ATOM DEFINED FOR SEGNAME'
            sys.exit()

        # TODO need to handle the exception in complex_filter.py
        # ONLY protein and RNA need this alignment

    # get alignment sub molecule

        align_filter = 'name[i] == "' + this_basis + '" and (segname[i] == "' + this_segment + '") and (resid[i] >= ' + str(
            seglow[i]) + ' and resid[i] <= ' + str(seghigh[i]) + ')'
        error, align_mask = m1.get_subset_mask(align_filter)

        all_flexible_align_mask.append(align_mask)

        sub_m1 = sasmol.SasMol(2)
        error = m1.copy_molecule_using_mask(sub_m1, align_mask, 0)
        com_sub_m1 = sub_m1.calccom(0)
        sub_m1.center(0)
        coor_sub_m1 = sub_m1.coor()[0]

        all_flexible_coor_sub_m1.append(coor_sub_m1)
        all_flexible_com_sub_m1.append(com_sub_m1)

        sub_m2 = sasmol.SasMol(4)
        error = m1.copy_molecule_using_mask(sub_m2, align_mask, 0)

        all_flexible_sub_m2.append(sub_m2)

    return all_flexible_align_mask, all_flexible_coor_sub_m1, all_flexible_com_sub_m1, all_flexible_sub_m2


def run_file_utilities(runname, pdbpath, pdbfile, dcdfile):

    direxist = os.path.exists(runname)
    if(direxist == 0):
        os.system('mkdir -p ' + runname + '/')
#
#       global run administration
#
    genpath = runname + '/complex_monte_carlo'
    genpaths = genpath + '/'
    direxist = os.path.exists(genpath)
    if(direxist == 0):
        os.system('mkdir -p ' + genpath)

    cpst = 'cp ' + pdbpath + '/' + pdbfile + ' ' + genpaths
    os.system(cpst)
#
#       write global run name, pdb, and dcd filenames to .last_sas
#
    fileexist = os.path.exists('.last_sas')
    if(fileexist == 1):
        os.system('mv -f .last_sas .last_sas_bu')
    lastsasfile = open('./.last_sas', 'w')
    lastsasfile.write('run_name\t' + runname + '\n')
    lastsasfile.write('pdb_name\t' + pdbfile + '\n')
    lastsasfile.write('dcd_name\t' + dcdfile + '\n')

    return lastsasfile, genpaths


def process_input_variables(psegvariables, segbasis, sseglow, sseghigh, flpsegname):

    allsith = []
    allsnumranges = []
    allsrlow = []
    allsrnum = []
    allmoltype = []

    for i in range(len(psegvariables)):
        allsnumranges.append(psegvariables[i][0])
        allsith.append(psegvariables[i][1])
        allsrlow.append(psegvariables[i][2])
        allsrnum.append(psegvariables[i][3])
        allmoltype.append(psegvariables[i][4])

    # abasis=string.split(segbasis,',')
    abasis = [item.strip() for item in string.split(segbasis, ',')]

    # seglow=string.split(sseglow,',')
    # seghigh=string.split(sseghigh,',')

    aith = []
    anumranges = []
    arlow = []
    arnum = []
    amoltype = []
    for i in range(len(allsith)):
        linith = string.split(allsith[i], ',')
        locith = []
        for i in range(len(linith)):
            tith = linith[i]
            fith = locale.atof(tith)
            if(fith > 180.0):
                fith = 180.0
            elif(fith < 0.0):
                fith = 0.0
            locith.append(fith)
        aith.append(locith)

    for i in range(len(allsnumranges)):
        nr = locale.atoi(allsnumranges[i])
        anumranges.append(nr)

    for i in range(len(allsrlow)):
        linrlow = string.split(allsrlow[i], ',')
        linrnum = string.split(allsrnum[i], ',')
        rlow = []
        rnum = []
        for k in range(len(linrlow)):
            trlow = locale.atoi(linrlow[k])
            trnum = locale.atoi(linrnum[k])
            rlow.append(trlow)
            rnum.append(trnum)
        # print 'rlow = ',rlow
        # print 'rnum = ',rnum
        arlow.append(rlow)
        arnum.append(rnum)

    for i in range(len(psegvariables)):
        moltype = allmoltype[i].strip()
        amoltype.append(moltype)

    '''
	print 'anumranges = ',anumranges
	print 'aith = ',aith
	print 'arlow = ',arlow
	print 'arnum = ',arnum
	'''

    raw_flexible_segments = string.split(flpsegname, ",")

    flexible_segments = []
    for fp in raw_flexible_segments:
        flexible_segments.append(fp.strip())

    # print 'flexible_segments = ',flexible_segments

    return amoltype, allsith, allsnumranges, allsrlow, allsrnum, abasis, sseglow, sseghigh, anumranges, aith, arlow, arnum, flexible_segments


def initialize_segments(m1, flexible_segments, nsegments, abasis):

    segname = m1.segname()

    asegs = []
    for tseg in segname:
        if(tseg not in asegs):
            asegs.append(tseg)
    numsegs = len(asegs)
    print 'found ', numsegs, ' segment names'

    first_last_resid = []

    all_segment_mask = []
    all_segment_full_mask = []
    all_segment_basis_full_mask = []

    all_segment_mol = []

    tmask = ''

    keyword_basis = False
    if(len(abasis) == 1):
        basis = abasis[0].strip()
        if(basis.lower() == 'all' or basis.lower() == 'heavy' or basis.lower() == 'backbone'):
            keyword_basis = True

    for i in xrange(numsegs):
        segmol = sasmol.SasMol(0)
        error, segment_full_mask = m1.get_subset_mask(
            'segname[i] == "' + asegs[i] + '"')
        m1.copy_molecule_using_mask(segmol, segment_full_mask, 0)
        this_resid = segmol.resid()

        first_last_resid.append([this_resid[0], this_resid[-1]])
        all_segment_full_mask.append(segment_full_mask)
        all_segment_mol.append(segmol)

    # this is where abasis is used  --> and this is where it matters!
        if keyword_basis:
            if(basis.lower() == 'all'):
                # print 'setting up all atom overlap arrays'
                segmol.set_average_vdw()
                npairs = segmol.natoms() * (segmol.natoms() - 1) / 2
                cutoff_array = numpy.zeros(npairs, numpy.float)
                pairs.pairs(segmol.atom_vdw(), cutoff_array)
                keyword_basis_filter = 'segname[i] == "' + \
                    asegs[i] + '" and (not name[i] == "None") '
            elif(basis.lower() == 'backbone'):
                this_moltype = segmol.moltype()[0]
                # print 'this_moltype = ',this_moltype  ### check this

                if(segmol.moltype()[0] == 'protein'):
                    keyword_basis_filter = 'segname[i] == "' + asegs[
                        i] + '" and (name[i] == "N" or name[i] == "CA" or name[i] == "C") '
                elif(segmol.moltype()[0] == 'rna' or segmol.moltype()[0] == 'dna'):
                    keyword_basis_filter = 'segname[i] == "' + asegs[
                        i] + '" and (name[i] == "P" or name[i] == "O5\'" or name[i] == "C5\'" or name[i] == "C4\'" or name[i] == "C3\'" or name[i] == "O3\'") '
                else:
                    keyword_basis_filter = 'segname[i] == "' + \
                        asegs[i] + '" and (not name[i][0] == "H") '
                # TODO --> add to complex_filter so the following hack is not
                # needed
            elif(basis.lower() == 'heavy'):
                keyword_basis_filter = 'segname[i] == "' + \
                    asegs[i] + '" and (not name[i][0] == "H") '

            error, segment_basis_mask = m1.get_subset_mask(
                keyword_basis_filter)
        else:
            error, segment_basis_mask = m1.get_subset_mask(
                'segname[i] == "' + asegs[i] + '" and name[i] =="' + abasis[i].strip() + '"')
        all_segment_basis_full_mask.append(segment_basis_mask)

        error, segment_mask = all_segment_mol[i].get_subset_mask(
            'segname[i] == "' + asegs[i] + '"')
        all_segment_mask.append(segment_mask)

    # TODO ... this is probably why flexible segments need to be first!!
    # should just take the NAMES of the flexible segnames to make this
    ###
    # this is also where abasis is used  --> but basis_full_mask is ONLY used for zcut
    # checking: abasis itself is passed to check_overlap in nmer_nrotate
    ###
    # OPTIONS: use moltype()[0] for each asegs[i] to set the basis (CA--> protein, P --> RNA)
    # or better yet, use not hydrogen instead ... as this is ONLY used for z-cut check
    ###
        tmask += '(segname[i] == "' + asegs[i] + \
            '" and (not name[i][0] == "H")) '
        #tmask+='segname[i] == "'+asegs[i]+'" and name[i] =="'+abasis[i].strip()+'"'
        if i != len(flexible_segments) - 1:
            tmask += ' or '

    error, basis_full_mask = m1.get_subset_mask(tmask)

    # print 'first_last_resid = ',first_last_resid

    return asegs, first_last_resid, all_segment_mask, all_segment_full_mask, all_segment_basis_full_mask, basis_full_mask, all_segment_mol, keyword_basis


def initialize_interaction_regions(m1, interpairs, npairs, cutoff, sseglow, asegs, abasis):

    if(len(interpairs) > 0):
        print 'pair distances < cut == ', cutoff, ' angstroms between segments have been found'
        print 'these distances will be ignorned in overlap check'
        print 'interpairs = ', interpairs
    else:
        print 'all distances between segments are greater than cut == ', cutoff
        print 'normal overlap checking will be used'
    print 'npairs = ', npairs

    ### initialize interaction regions in each segment ###

    interres = []
    interatom = []
    for i in range(len(interpairs)):
        segnum_1 = interpairs[i][0][0]
        segnum_2 = interpairs[i][0][1]
        for j in range(len(interpairs[i][1])):
            resnum_1 = interpairs[i][1][j][0]
            resnum_2 = interpairs[i][1][j][1]

        # TODO --> need to match basis here as well
        # TODO --> need to match basis here as well
        # TODO --> need to match basis here as well

            basis_segment_1 = '(segname[i] == "' + asegs[segnum_1] + \
                '" and name[i] =="' + abasis[segnum_1].strip() + '")'
            error, basis_mask_segment_1 = m1.get_subset_mask(basis_segment_1)
            # idx_1 = numpy.where(basis_mask_segment_1==1.0)[0][resnum_1] # a
            # ugly numpy function
            idx_1 = filter(lambda x: basis_mask_segment_1[x] == 1.0, range(
                len(basis_mask_segment_1)))[resnum_1]
            basis_segment_2 = '(segname[i] == "' + asegs[segnum_2] + \
                '" and name[i] =="' + abasis[segnum_2].strip() + '")'
            error, basis_mask_segment_2 = m1.get_subset_mask(basis_segment_2)
            # idx_2 = numpy.where(basis_mask_segment_2==1.0)[0][resnum_2] # a
            # ugly numpy function
            idx_2 = filter(lambda x: basis_mask_segment_2[x] == 1.0, range(
                len(basis_mask_segment_2)))[resnum_2]
            interres.append([resnum_1, resnum_2])
            interatom.append([idx_1, idx_2])

    print 'interres  = ', interres
    print 'interatom = ', interatom

    return interatom, interres


def set_up_dihedral_arrays(all_segment_mol, asegs, abasis, amoltype, first_last_resid, flexible_segments, anumranges, arlow, arnum, keyword_basis, txtOutput):

    flexible_dihedral_parameters = []

    all_flexible_basis_mask = []

    for i in xrange(len(flexible_segments)):

        this_segname = flexible_segments[i]
        idx = asegs.index(this_segname)
        m1 = all_segment_mol[idx]

        # TODO --> need to deal with specific basis here
        # TODO --> need to deal with specific basis here
        # TODO --> need to deal with specific basis here

        if(keyword_basis):

            if amoltype[i] == 'protein':
                basis_atom = "CA"

            elif amoltype[i] == 'rna':
                #basis_atom = "P"
                basis_atom = "O5\'"

            basis_filter = 'name[i] == "' + basis_atom + \
                '" and segname[i] == "' + this_segname + '"'

        else:
            basis_filter = 'name[i] == "' + abasis[idx] + \
                '" and segname[i] == "' + this_segname + '"'

        error, basis_mask = m1.get_subset_mask(basis_filter)
        all_flexible_basis_mask.append(basis_mask)

        basis_m1 = sasmol.SasMol(1)
        error = m1.copy_molecule_using_mask(basis_m1, basis_mask, 0)
        basis_resname = basis_m1.resname()
        basis_resid = basis_m1.resid()

        arespsi = []
        aresphi = []
        numranges = anumranges[i]
        reslow = arlow[i]
        numcont = arnum[i]

        if amoltype[i] == 'protein':
            respsi = []
            resphi = []
            dihedral_energy.protein_initialization(
                respsi, resphi, basis_resid, basis_resname, numranges, reslow, numcont, first_last_resid[idx], txtOutput)
            flexible_dihedral_parameters.append([respsi, resphi])

        elif amoltype[i] == 'rna':
            resalpha = []
            resbeta = []
            resgamma = []
            resdelta = []
            resepsilon = []
            reseta = []
            dihedral_energy.rna_initialization(resalpha, resbeta, resgamma, resdelta, resepsilon, reseta,
                                               basis_resid, basis_resname, numranges, reslow, numcont, first_last_resid[idx], txtOutput)
            flexible_dihedral_parameters.append(
                [resalpha, resbeta, resgamma, resdelta, resepsilon, reseta])

    return flexible_dihedral_parameters, all_flexible_basis_mask


def set_up_constraints(m1, cflag, confile):

    if(cflag == 1):

        filter_flag = 0
        error, constraint_basis1_array, constraint_basis2_array, distance_array, type_array = constraints.read_constraints(
            m1, confile, filter_flag)

        mask_a_array = []
        mask_b_array = []

        for i in xrange(len(distance_array)):
            print constraint_basis1_array[i]
            print constraint_basis2_array[i]
            print distance_array[i]
            print type_array[i]

            error, local_mask_a_array = m1.get_subset_mask(
                constraint_basis1_array[i])
            error, local_mask_b_array = m1.get_subset_mask(
                constraint_basis2_array[i])

            mask_a_array.append(local_mask_a_array)
            mask_b_array.append(local_mask_b_array)

    else:
        mask_a_array = []
        mask_b_array = []
        distance_array = []
        type_array = []

    return mask_a_array, mask_b_array, distance_array, type_array


def setup_flexible_residue_mask_arrays(m1, flexible_segments, anumranges, arlow, arnum, amoltype, txtOutput):

    all_flexible_residues = []
    all_flexible_residue_rotation_indices = []
    all_flexible_residue_rotation_mask = []

    for i in xrange(len(flexible_segments)):

        numranges = anumranges[i]
        reslow = arlow[i]
        numcont = arnum[i]
        flexible_residues = dihedral_monte_carlo.get_flexible_residues(
            numranges, reslow, numcont)
        all_flexible_residues.append(flexible_residues)

        segment_filter = 'segname[i] == "' + flexible_segments[i] + '"'
        error, segment_mask = m1.get_subset_mask(segment_filter)
        # print 'segment_filter = ',segment_filter
        # print 'error = ',error
        segment_m1 = sasmol.SasMol(98)
        error = m1.copy_molecule_using_mask(segment_m1, segment_mask, 0)

        molecule_type = amoltype[i]

        residue_rotation_indices, residue_rotation_mask = dihedral_monte_carlo.get_rotation_indices(
            segment_m1, molecule_type, flexible_residues, txtOutput)
        all_flexible_residue_rotation_indices.append(residue_rotation_indices)
        all_flexible_residue_rotation_mask.append(residue_rotation_mask)

    return all_flexible_residues, all_flexible_residue_rotation_indices, all_flexible_residue_rotation_mask


def evaluate_rg(rg_difference_list, directed_rg_list, accepted_rg_list, this_rg_difference, this_rg, accepted):

    maximum_value = max(rg_difference_list)

    if(maximum_value > this_rg_difference):
        index = rg_difference_list.index(maximum_value)
        rg_difference_list[index] = this_rg_difference
        directed_rg_list[index] = this_rg
        accepted_rg_list[index] = accepted

    return

### main method ###


def dihedralgenerate(variables, psegvariables, txtOutput):

        # amoltype=['protein','protein']
        # amoltype=['rna','protein']
        # amoltype=['protein']

        # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))

    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    # unpack variables
    runname, dcdfile, path, pdbfile, trials, goback, temp, nsegments, segbasis, npsegments, flpsegname, sseglow, sseghigh, lowrg, highrg, zflag, zcutoff, cflag, confile, plotflag, directedmc, seed = unpack_variables(
        variables)

    segbasis.strip()

    # process variables
    amoltype, allsith, allsnumranges, allsrlow, allsrnum, abasis, seglow, seghigh, anumranges, aith, arlow, arnum, flexible_segments = process_input_variables(
        psegvariables, segbasis, sseglow, sseghigh, flpsegname)
    import pprint
    fout = open('a.txt', 'w')
    pprint.pprint(variables, fout)
    pprint.pprint(psegvariables, fout)
    pprint.pprint(segbasis, fout)
    pprint.pprint(seglow, fout)
    pprint.pprint(seghigh, fout)
    pprint.pprint(flpsegname, fout)
    fout.close()

    # set up run file I/O
    lastsasfile, genpaths = run_file_utilities(runname, path, pdbfile, dcdfile)

    kb = 1.380658E-23  # J/K
    beta = 1.0 / (temp * kb)

    m1 = sasmol.SasMol(0)
    m1.read_pdb(path + pdbfile)

    nf1 = m1.number_of_frames()
    # print 'nf1 = %d\n' % nf1

    dcdoutfile = m1.open_dcd_write(genpaths + dcdfile)

    # set up segment arrays
    asegs, first_last_resid, all_segment_mask, all_segment_full_mask, all_segment_basis_full_mask, basis_full_mask, all_segment_mol, keyword_basis = initialize_segments(
        m1, flexible_segments, nsegments, abasis)

    # set up constraints variables
    mask_a_array, mask_b_array, distance_array, type_array = set_up_constraints(
        m1, cflag, confile)

    # set up segment alignment coordinates and com arrays
    all_flexible_align_mask, all_flexible_coor_sub_m1, all_flexible_com_sub_m1, all_flexible_sub_m2 = alignment_initialization(
        all_segment_mol, asegs, abasis, flexible_segments, seglow, seghigh)

    if(keyword_basis):
        if(segbasis.lower() == 'all'):
            cutoff = 0.8
        elif(segbasis.lower() == 'heavy' or segbasis.lower() == 'backbone'):
            cutoff = 0.8
    else:
        cutoff = 2.0

    print 'cutoff = ', cutoff

    check_initial_interactions = False
    if(check_initial_interactions):

        # survey interaction between segments
        interpairs, npairs = nmer_overlap_check.nmer_overlap_check(
            m1, path, pdbfile, cutoff, abasis, keyword_basis)

        interatom, interres = initialize_interaction_regions(
            m1, interpairs, npairs, cutoff, sseglow, asegs, abasis)

    else:

        interpairs = []
        npairs = 0
        interatom = []
        interres = []

    # set up dihedral parameters for each flexible segment
    flexible_dihedral_parameters, all_flexible_basis_mask = set_up_dihedral_arrays(
        all_segment_mol, asegs, abasis, amoltype, first_last_resid, flexible_segments, anumranges, arlow, arnum, keyword_basis, txtOutput)

#	if(segbasis.lower() == 'all' or segbasis.lower() == 'heavy' or segbasis.lower() == 'backbone'):
#		print 'segbasis = ',segbasis,' so I should stop for now\n'
#		#sys.exit()
#	else:
#		print 'segbasis = ',segbasis,' so I should continue\n'

    # set up flexible residue rotation mask arrays
    all_flexible_residues, all_flexible_residue_rotation_indices, all_flexible_residue_rotation_mask = setup_flexible_residue_mask_arrays(
        m1, flexible_segments, anumranges, arlow, arnum, amoltype, txtOutput)

    step_parameters = step.Setup()

    hrg = 0.0
    lowestrg = 1000.0
    an = 'psi'
    accepted = 0
    over = 0
    badrg = 0
    badz = 0
    badc = 0
    nsteps = 0
    arg = 0.0
    trg = 0.0

    coor = m1.coor()
    frame = 0

    # MAIN LOOP

    q0 = 1
    th = 1.0
    seg = asegs[0]
    pairdat = [an, q0, th, seg]

    all_rg_tally = []
    accepted_rg_tally = []
    phi_tally = []
    aphi_tally = []
    psi_tally = []
    apsi_tally = []
    atpsi_tally = []
    atphi_tally = []
    atphipsi_tally = []

    if(plotflag == 1):
        graph = Gnuplot.Gnuplot(debug=1)
        graph.clear()
        graph('set title "Rg Results"')
        graph.xlabel('Structure Number')
        graph.ylabel('Rg (Angstrom^2)')

    nonbondflag = 0

    if(seed[0] == 1):
        from numpy.random import RandomState
        seed_object = RandomState(seed[1])
    else:
        seed_object = -1

    failtally = 0
    acc = 0
    afile = ''
    accfile = []

    minx = []
    miny = []
    minz = []
    maxx = []
    maxy = []
    maxz = []

    if(directedmc > 0):
        rg_difference_list = []
        directed_rg_list = []
        accepted_rg_list = []
        rg_list_length = 10  # hardwired

    for i in range(trials):

        if(seed[0] == 1):
            ran_num = seed_object.rand()
            tflexsegn = int(len(flexible_segments) * ran_num)
            tsegn = asegs.index(flexible_segments[tflexsegn])
        else:
            tflexsegn = int(len(flexible_segments) * random.random())
            tsegn = asegs.index(flexible_segments[tflexsegn])

        tseg = asegs[tsegn]

        molecule_type = amoltype[tflexsegn]

        dtheta = aith[tflexsegn]
        numranges = anumranges[tflexsegn]
        reslow = arlow[tflexsegn]
        numcont = arnum[tflexsegn]

        segment_full_mask = all_segment_full_mask[tsegn]

        error, new_coor = m1.get_coor_using_mask(frame, segment_full_mask)

        segment_mol = all_segment_mol[tsegn]

        segment_mol.setCoor(new_coor)

        '''
		if(i<10):
			print 'segment_mol.coor()[0,0,0] = ',segment_mol.coor()[0,0,0]

		else:
			sys.exit()
		'''

        vdi, vdf, indices, this_mask = step_parameters.chooser(new_coor, segment_mol, pairdat, dtheta, numranges, reslow, numcont, flexible_dihedral_parameters[
                                                               tflexsegn], beta, all_flexible_residue_rotation_indices[tflexsegn], all_flexible_residue_rotation_mask[tflexsegn], nonbondflag, first_last_resid[tsegn], molecule_type, seed_object)

        '''
		print 'len(indices) = ',len(indices)
		print 'indices[0] = ',indices[0]
		print 'indices[-1] = ',indices[-1]
		print 'tsegn = ',tsegn
		'''

        pairdat[3] = tseg
        an = pairdat[0]
        q0 = pairdat[1]
        th = pairdat[2]
        seg = pairdat[3]
        nsteps += 1
        re = [0, 0, 0, 0.0, 0.0, lowestrg, hrg, 0, 0, []]

        newafile = nmer_nrotate.rotate(coor, m1, q0, th, an, cutoff, lowrg, highrg, re, accepted, zflag, zcutoff, cflag, dcdoutfile, indices, this_mask, all_flexible_basis_mask[tflexsegn], all_flexible_sub_m2[tflexsegn], all_flexible_align_mask[tflexsegn], all_flexible_coor_sub_m1[
                                       tflexsegn], all_flexible_com_sub_m1[tflexsegn], mask_a_array, mask_b_array, distance_array, type_array, first_last_resid[tsegn], molecule_type, all_segment_mask[tsegn], segment_full_mask, all_segment_basis_full_mask, basis_full_mask, all_segment_mol[tsegn], asegs, abasis, interatom, interres)

        print '.',
        sys.stdout.flush()

        accepted = accepted + re[0]
        over = over + re[1]
        badrg = badrg + re[2]
        rg_value = re[3]
        trg = trg + re[3]
        arg = arg + re[4]
        lowestrg = re[5]
        hrg = re[6]
        badz = badz + re[7]
        badc = badc + re[8]

        if(len(re[9]) > 0):
            minmax = re[9]

            minx.append(minmax[0][0])
            miny.append(minmax[0][1])
            minz.append(minmax[0][2])
            maxx.append(minmax[1][0])
            maxy.append(minmax[1][1])
            maxz.append(minmax[1][2])

        all_rg_tally.append([i, rg_value])

        if(re[0] == 1):
            accepted_rg_tally.append([i, accepted, rg_value])
            if(directedmc > 0):
                if(len(rg_difference_list) <= rg_list_length):
                    this_rg_difference = abs(rg_value - directedmc)
                    rg_difference_list.append(this_rg_difference)
                    directed_rg_list.append(rg_value)
                    accepted_rg_list.append(accepted)
                else:
                    this_rg_difference = abs(rg_value - directedmc)
                    evaluate_rg(rg_difference_list, directed_rg_list,
                                accepted_rg_list, this_rg_difference, rg_value, accepted)

        if(re[0] == 0):
            if(failtally == goback):
                failtally = 0
                if(accepted > 0):
                    if(seed[0] == 1):
                        ran_num = seed_object.rand()
                        dum = int(accepted * ran_num) - 1

                    elif(directedmc > 0):
                        local_rg_list_length = len(directed_rg_list)
                        ran_num = random.randrange(0, local_rg_list_length)
                        dum = accepted_rg_list[ran_num]

                    else:
                        dum = int(accepted * random.random()) - 1
                    if(dum == -1):
                        print '\nreloading coordinates from original starting structure'
                        m1.read_pdb(path + pdbfile, fastread=True,
                                    saspdbrx_topology=True)
                        coor = m1.coor()
                    else:
                        print '\nreloading coordinates from a previously accepted structure'

                        m1.read_single_dcd_step(genpaths + dcdfile, dum + 1)
                        # m1.read_single_dcd_step(genpaths+dcdfile,dum)
                        coor = m1.coor()
                else:
                    print '\n>>>>>reloading coordinates from original starting structure'
                    m1.read_pdb(path + pdbfile, fastread=True,
                                saspdbrx_topology=True)
                    coor = m1.coor()
            else:
                failtally = failtally + 1

        if(((i + 1) % (float(trials) / 100.0) == 0 or (trials < 10))):
            fraction_done = (float(i + 1) / float(trials))
            progress_string = '\nCOMPLETED ' + \
                str(i + 1) + ' of ' + str(trials) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            print accepted, ' configurations accepted out of ', nsteps, (float(accepted) / nsteps) * 100.0, ' %\n\n'
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

        if(i > 9):
            if((i + 1) % (trials / 10) == 0 and accepted > 0 and i + 1 > 10):
                if(plotflag == 1):
                    graph.plot(Gnuplot.Data(all_rg_tally, using='1:2 w p ps 4', title='all Rg'), Gnuplot.Data(
                        accepted_rg_tally, using='1:3 w lp pt 5 ps 2', title='accepted'))
                fraction_done = (float(i + 1) / float(trials))
                report_string = 'STATUS\t' + str(fraction_done)
                txtOutput.put(report_string)

        elif(accepted > 0):
            if(plotflag == 1):
                graph.plot(Gnuplot.Data(all_rg_tally, using='1:2 w p ps 4', title='all Rg'), Gnuplot.Data(
                    accepted_rg_tally, using='1:3 w lp pt 5 ps 2', title='accepted'))

            fraction_done = (float(i + 1) / float(trials))
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

    m1.close_dcd_write(dcdoutfile)

    rgplot = open('./' + runname + '/complex_monte_carlo/' +
                  dcdfile + '.all_rg_results_data.txt', 'w')
    rgplot.write('# structure number (structure 1 = 1; not 0), Rg (all)\n')
    for ii in range(len(all_rg_tally)):
        rgplot.write('%i\t%f\n' %
                     (all_rg_tally[ii][0] + 1, all_rg_tally[ii][1]))
    rgplot.close()
    rgplot = open('./' + runname + '/complex_monte_carlo/' +
                  dcdfile + '.accepted_rg_results_data.txt', 'w')
    rgplot.write(
        '# structure number (structure 1 = 1; not 0), Rg (accepted)\n')
    for ii in range(len(accepted_rg_tally)):
        rgplot.write('%i\t%f\t%i\n' % (accepted_rg_tally[ii][
                     1] - 1, accepted_rg_tally[ii][2], accepted_rg_tally[ii][0] + 1))
    rgplot.close()

    '''
	outfile2=open(genpaths+dcdfile+'.phi','w')
	outfile3=open(genpaths+dcdfile+'.psi','w')
	outfile5=open(genpaths+dcdfile+'.aphi','w')
	outfile6=open(genpaths+dcdfile+'.apsi','w')
	outfile7=open(genpaths+dcdfile+'.aphivsapsi','w')

	outfile7.write('# ACCEPTED STRUCTURES\n')
	outfile7.write('# AA phi psi\n')

	for i in range(len(phi_tally)):
		outfile2.write('%i\t%f\n' % (phi_tally[i][0],phi_tally[i][1]))
	for i in range(len(psi_tally)):
		outfile3.write('%i\t%f\n' % (psi_tally[i][0],psi_tally[i][1]))
	for i in range(len(aphi_tally)):
		outfile5.write('%i\t%f\n' % (aphi_tally[i][0],aphi_tally[i][1]))
	for i in range(len(apsi_tally)):
		outfile6.write('%i\t%f\n' % (apsi_tally[i][0],apsi_tally[i][1]))
	for i in range(len(atphipsi_tally)):
		outfile7.write('%i\t%f\t%f\n' % (atphipsi_tally[i][0],atphipsi_tally[i][1],atphipsi_tally[i][2]))

	outfile2.close()
	outfile3.close()
	outfile5.close()
	outfile6.close()
	outfile7.close()
	'''

    ttxt = time.ctime()

    st = ''.join(['=' for x in xrange(60)])

    if(accepted > 0):
        txtOutput.put("Average accepted rg2 = %lf\n" % (arg / (accepted)))
        txtOutput.put(
            "Configurations and statistics saved in %s directory\n" % ('./' + genpaths))
    else:
        txtOutput.put("Average accepted rg2 = %lf\n" % (0.0))
        txtOutput.put(
            "\n NO ACCEPTED MOVES\n\n Statistics saved in %s directory\n" % (genpaths))

    outfile7 = open(genpaths + dcdfile + '.stats', 'w')
    outfile7.write('%s\t%f\t%s\t%f\n' %
                   ('lowest Rg = ', lowestrg, 'highest Rg = ', hrg))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('accepted ', accepted,
                                                   ' out of ', nsteps, ' moves : ', (accepted / float(nsteps)) * 100.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('overlapped ', over,
                                                   ' out of ', nsteps, ' moves : ', (over / float(nsteps)) * 100.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad rg2 ', badrg,
                                                   ' out of ', nsteps, ' moves : ', (badrg / float(nsteps)) * 100.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad z-filter ', badz,
                                                   ' out of ', nsteps, ' moves : ', (badz / float(nsteps)) * 100.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad constaints ', badc,
                                                   ' out of ', nsteps, ' moves : ', (badc / float(nsteps)) * 100.0, ' %'))

    if(accepted > 0):
        outfile7.write('%s\t%f\n' %
                       ('average accepted rg2 = ', arg / (accepted)))
    else:
        outfile7.write('%s\t%f\n' % ('average accepted rg2 = ', 0.0))
    outfile7.write('%s\t%f\n' %
                   ('average total rg2 of ensemble = ', trg / (nsteps)))

    print '\n\nDCD data were written to %s\n' % ('./' + genpaths + dcdfile)
    txtOutput.put("\nDCD data were written to %s\n\n" %
                  ('./' + genpaths + dcdfile))
    txtOutput.put("lowest Rg = %lf\t highest Rg = %lf\n" % (lowestrg, hrg))
    txtOutput.put("accepted %d out of %d : %lf percent\n" %
                  (accepted, nsteps, (accepted / float(nsteps)) * 100.0))
    txtOutput.put("overlapped %d out of %d moves : %lf percent\n" %
                  (over, nsteps, (float(over) / float(nsteps)) * 100.0))
    txtOutput.put("bad rg2 %d out of %d moves : %lf percent\n" %
                  (badrg, nsteps, (float(badrg) / float(nsteps)) * 100.0))
    if(zflag == 1):
        txtOutput.put("bad zcut %d out of %d moves : %lf percent\n\n\n" % (
            badz, nsteps, (float(badz) / float(nsteps)) * 100.0))
    if(cflag == 1):
        txtOutput.put("constraint filter rejected %d out of %d moves : %lf percent\n\n\n" % (
            badc, nsteps, (float(badc) / float(nsteps)) * 100.0))

    if(len(minx) > 0 and len(miny) > 0 and len(minz) > 0 and len(maxx) > 0 and len(maxy) > 0 and len(maxz) > 0):
        min_x = numpy.min(minx)
        min_y = numpy.min(miny)
        min_z = numpy.min(minz)
        max_x = numpy.max(maxx)
        max_y = numpy.max(maxy)
        max_z = numpy.max(maxz)

        txtOutput.put("\nminimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                      (min_x, max_x, (max_x - min_x)))

        txtOutput.put("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                      (min_y, max_y, (max_y - min_y)))
        txtOutput.put("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" %
                      (min_z, max_z, (max_z - min_z)))

        outfile7.write("\nminimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                       (min_x, max_x, (max_x - min_x)))

        outfile7.write("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                       (min_y, max_y, (max_y - min_y)))

        outfile7.write("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" %
                       (min_z, max_z, (max_z - min_z)))

        outfile7.close()
    else:
        outfile7.close()

    txtOutput.put("\n%s \n" % (st))

    lastsasfile.close()
    print 'COMPLEX DIHEDRAL IS DONE'
    time.sleep(1.5)

    if(plotflag == 1):
        wait('\n')

    return()


if __name__ == '__main__':

    runname = 'run_0'
    dcdfile = 'run_0.dcd'
    path = './'
    pdbfile = 'fram601.pdb'
    trials = '50'
    goback = '50'
    nsegments = '2'
    npsegments = '2'
    flpsegname = 'ENDA,ENDB'
    segbasis = 'CA, CA'
    #segbasis = 'all'
    #segbasis = 'heavy'
    #segbasis = 'backbone'
    seglow = '95, 95'
    seghigh = '110, 110'
    temp = '300.0'
    lowrg = '20.0'
    highrg = '185.0'
    zflag = '0'
    zcutoff = '0.0'
    cflag = '0'
    confile = 'constraints.txt'
    directedmc = '0'
    psffilepath = './'
    psffilename = 'refgag.psf'
    import sassie.sasconfig as sasconfig
    parmfilepath = sasconfig._bin_path + 'toppar'
    parmfilename = 'par_all27_prot_na.inp'
    plotflag = '1'
    seed = '0, 123'

    svariables = {}

    svariables['cflag'] = (cflag, 'int')
    svariables['confile'] = (confile, 'string')
    svariables['dcdfile'] = (dcdfile, 'string')
    svariables['directedmc'] = (directedmc, 'float')
    svariables['flpsegname'] = (flpsegname, 'string')
    svariables['goback'] = (goback, 'int')
    svariables['highrg'] = (highrg, 'float')
    svariables['lowrg'] = (lowrg, 'float')
    svariables['npsegments'] = (npsegments, 'int')
    svariables['nsegments'] = (nsegments, 'int')
    svariables['parmfilename'] = (parmfilename, 'string')
    svariables['path'] = (path, 'string')
    svariables['pdbfile'] = (pdbfile, 'string')
    svariables['plotflag'] = (plotflag, 'int')
    svariables['psffilename'] = (psffilename, 'string')
    svariables['runname'] = (runname, 'string')
    svariables['seed'] = (seed, 'int_array')
    svariables['segbasis'] = (segbasis, 'string')
    svariables['seghigh'] = (seghigh, 'int_array')
    svariables['seglow'] = (seglow, 'int_array')
    svariables['temp'] = (temp, 'float')
    svariables['trials'] = (trials, 'int')
    svariables['zcutoff'] = (zcutoff, 'float')
    svariables['zflag'] = (zflag, 'int')

    psegvariables = [['1', '30', '2', '30', 'protein'],
                     ['1', '30', '2', '30', 'protein']]

    import sassie.interface.input_filter as input_filter
    error, variables = input_filter.type_check_and_convert(svariables)

    # error=generate_filter.check_protein(variables,eflag,monflag)

    if(len(error) > 0):
        print 'error = ', error
        sys.exit()

    runname = variables['runname'][0]

    import multiprocessing

    txtQueue = multiprocessing.JoinableQueue()
    dihedralgenerate(variables, psegvariables, txtQueue)
