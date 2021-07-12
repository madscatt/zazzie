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
import time
import platform
import numpy
import random
import Gnuplot
import Gnuplot.PlotItems
import Gnuplot.funcutils
import sassie.sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monte_carlo.monomer.pairs as pairs
import sassie.simulate.monte_carlo.monomer.dihedral_rotate as dihedral_rotate
import sassie.simulate.monte_carlo.monomer.step as step

#       DIHEDRAL
#
#	09/26/2005 	--	gag-dihedral search		:	jc
#	07/14/2008	--	file management changes		:	jc
#	01/12/2011	--	added sasmol support		:	jc
#	01/15/2011	--	generalized rotation basis	:	jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        DIHEDRAL is the module that contains the functions
        that are used to generate ensembles of structures by varying
	protein dihedral angles.

        This module is called from Protein Monomer Dihedral Generation from the main
        GUI through the graphical_generate.py script.

        This module calls to C / Python extension modules to speed up
        calculations (see: overlap.c).

	REFERENCE:

    	J. A. D. MacKerell et al.
    	Journal of Physical Chemistry B,  102  3586-3616  (1998)

    	B. R. Brooks et al.
    	Journal of Computational Chemistry  4  187--217  (1983)

'''


def unpack_variables(variables):

    runname = variables['runname'][0]
    dcdfile = variables['dcdfile'][0]
    path = variables['path'][0]
    pdbfile = variables['pdbfile'][0]
    trials = variables['trials'][0]
    goback = variables['goback'][0]
    temp = variables['temp'][0]
    molecule_type = variables['moltype'][0]
    numranges = variables['numranges'][0]
    dtheta = variables['dtheta'][0]
    reslow = variables['reslow'][0]
    numcont = variables['numcont'][0]
    lowres1 = variables['lowres1'][0]
    highres1 = variables['highres1'][0]
    basis = variables['basis'][0]
    cutoff = variables['cutoff'][0]
    lowrg = variables['lowrg'][0]
    highrg = variables['highrg'][0]
    zflag = variables['zflag'][0]
    zcutoff = variables['zcutoff'][0]
    cflag = variables['cflag'][0]
    confile = variables['confile'][0]
    nonbondflag = variables['nonbondflag'][0]
    nonbondscale = variables['nonbondscale'][0]
    psffilepath = variables['psffilepath'][0]
    psffilename = variables['psffilename'][0]
    parmfilepath = variables['parmfilepath'][0]
    parmfilename = variables['parmfilename'][0]
    plotflag = variables['plotflag'][0]
    directedmc = variables['directedmc'][0]
    seed = variables['seed'][0]

    return runname, dcdfile, path, pdbfile, trials, goback, temp, molecule_type, numranges, dtheta, reslow, numcont, lowres1, highres1, basis, cutoff, lowrg, highrg, zflag, zcutoff, cflag, confile, nonbondflag, nonbondscale, psffilepath, psffilename, parmfilepath, parmfilename, plotflag, directedmc, seed


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


def get_flexible_residues(numranges, reslow, numcont):
    '''
    Method to determine residues that are going to be rotated.
    '''
    flexible_residues = []
    templist = []
    for i in xrange(numranges):
        thisreslow = reslow[i]
        thisnumcont = numcont[i]
        thisrange = numpy.arange(thisnumcont) + thisreslow
        templist.append(thisrange.tolist())

    flexible_residues = [item for sublist in templist for item in sublist]

    print 'flexible_resiudes = ', flexible_residues

    return flexible_residues


def get_rotation_indices(m1, molecule_type, flexible_residues, txtOutput):

    print 'getting rotation indices for molecule'

    residue_rotation_indices = {}
    residue_rotation_mask = {}

    if(molecule_type == 'protein'):

        mtype = 0
        mask = m1.get_dihedral_subset_mask(flexible_residues, mtype)

        for i in xrange(len(mask)):
            this_mask = mask[i][:]
            q0 = flexible_residues[i]
            residue_rotation_mask[q0] = this_mask.tolist()
            indices = m1.get_indices_from_mask(this_mask)
            residue_rotation_indices[q0] = indices.tolist()

        # print 'residue_rotation_indices = ',residue_rotation_indices
        # print 'residue_rotation_mask = ',residue_rotation_mask

    elif(molecule_type == 'rna'):

        mtype = 1
        mask = m1.get_dihedral_subset_mask(flexible_residues, mtype)

        for i in xrange(len(mask)):
            this_mask = mask[i][:]
            q0 = flexible_residues[i]
            residue_rotation_mask[q0] = this_mask.tolist()
            indices = m1.get_indices_from_mask(this_mask)
            residue_rotation_indices[q0] = indices.tolist()

#		print 'residue_rotation_indices = ',residue_rotation_indices
#		print 'residue_rotation_mask = ',residue_rotation_mask

    elif(molecule_type == 'dna'):
        message = 'rotation basis set not defined for molecule type = ' + molecule_type
        pass

    else:
        message = 'rotation basis set not defined for molecule type = ' + molecule_type
        print_failure(message, txtOutput)

    print 'done getting rotation indices'
    sys.stdout.flush()

    return residue_rotation_indices, residue_rotation_mask


def evaluate_rg(rg_difference_list, directed_rg_list, accepted_rg_list, this_rg_difference, this_rg, accepted):

    maximum_value = max(rg_difference_list)

    if(maximum_value > this_rg_difference):
        index = rg_difference_list.index(maximum_value)
        rg_difference_list[index] = this_rg_difference
        directed_rg_list[index] = this_rg
        accepted_rg_list[index] = accepted

    return


def dihedralgenerate(variables, txtOutput):

    runname, dcdfile, path, pdbfile, trials, goback, temp, molecule_type, numranges, dtheta, reslow, numcont, lowres1, highres1, basis, cutoff, lowrg, highrg, zflag, zcutoff, cflag, confile, nonbondflag, nonbondscale, psffilepath, psffilename, parmfilepath, parmfilename, plotflag, directedmc, seed = unpack_variables(
        variables)

    alignflag = '1'
    saveflag = '0'

    pdbfilename = pdbfile

    if(runname[-1] == '/'):
        lin = len(runname)
        runname = runname[:lin - 1]

    direxist = os.path.exists(runname)
    if(direxist == 0):
        os.system('mkdir -p ' + runname)

    genpath = runname + '/monomer_monte_carlo'
    genpaths = genpath + '/'
    direxist = os.path.exists(genpath)
    if(direxist == 0):
        os.system('mkdir -p ' + genpath)

    if(path != ""):
        if path[-1] != "/":
            path += '/'

    cpfile = os.path.join(path, pdbfilename)
    cpst = 'cp ' + cpfile + ' ' + genpaths
    os.system(cpst)

    parmfilepath += '/'

    fileexist = os.path.exists('.last_sas')
    if(fileexist == 1):
        os.system('mv -f .last_sas .last_sas_bu')
    lastsasfile = open('./.last_sas', 'w')
    lastsasfile.write('run_name\t' + runname + '\n')
    lastsasfile.write('pdb_name\t' + pdbfile + '\n')
    lastsasfile.write('dcd_name\t' + dcdfile + '\n')

    kb = 1.380658E-23  # J/K
    beta = 1.0 / (temp * kb)

    m1 = sasmol.SasMol(0)
    m2 = sasmol.SasMol(1)

    m1.read_pdb(path + pdbfile, saspdbrx_topology=True)
    nf1 = m1.number_of_frames()
    m1.calccom(0)

    dcdoutfile = m1.open_dcd_write(genpaths + dcdfile)

    resid = m1.resid()
    first_last_resid = [resid[0], resid[-1]]

    print 'first_last_resid = ', first_last_resid

    cutoff_array = []
    vdw_factor = 0.3

    if(molecule_type == 'protein'):

        if(basis.lower() == 'all'):
            print 'setting up all atom overlap arrays'
            m1.set_average_vdw()
            npairs = m1.natoms() * (m1.natoms() - 1) / 2
            cutoff_array = numpy.zeros(npairs, numpy.float)
            print 'cutoff_array: ', cutoff_array
            pairs.pairs(m1.atom_vdw(), cutoff_array)
            basis_filter = 'not name[i] == None'

        elif(basis.lower() == 'backbone'):
            basis_filter = 'name[i] == "N" or name[i] == "CA" or name[i] == "C" '
        elif(basis.lower() == 'heavy'):
            basis_filter = 'not name[i][0] == "H" '
        else:
            string_basis = string.split(basis, ",")
            number_basis_types = len(string_basis)
            for i in xrange(number_basis_types):
                if(i == 0):
                    basis_filter = 'name[i] == "' + string_basis[i] + '" '
                else:
                    basis_filter = basis_filter + \
                        ' or name[i] == "' + string_basis[i] + '" '

        print 'basis_filter = ', basis_filter

        error, basis_mask = m1.get_subset_mask(basis_filter)

        basis_atom = "CA"
        setup_basis_filter = 'name[i] == "' + basis_atom + '"'
        error, setup_basis_mask = m1.get_subset_mask(setup_basis_filter)

        basis_m1 = sasmol.SasMol(1)
        error = m1.copy_molecule_using_mask(basis_m1, setup_basis_mask, 0)
        basis_resname = basis_m1.resname()
        basis_resid = basis_m1.resid()

        respsi = []
        resphi = []
        dihedral_energy.protein_initialization(
            respsi, resphi, basis_resid, basis_resname, numranges, reslow, numcont, first_last_resid, txtOutput)
        dihedral_parameters = [respsi, resphi]

    elif(molecule_type == 'rna'):

        #	ATOM N1   NN2    -0.13  !                            \
        #	ATOM C6   CN3     0.05  !        O1P    H5' H4'  O4'  \
        #	ATOM H6   HN3     0.17  !         |      |    \ /   \  \
        #	ATOM C5   CN3    -0.13  !        -P-O5'-C5'---C4'    C1'
        #	ATOM H5   HN3     0.07  !         |      |     \     / \
        #	ATOM C2   CN1     0.52  !        O2P    H5''   C3'--C2' H1'
        #	ATOM O2   ON1C   -0.49  !                     / \   / \
        #	ATOM N3   NN3    -0.66  !                  O3' H3' O2' H2''
        #	ATOM C4   CN2     0.65  !                   |       |
        #	ATOM N4   NN1    -0.75  !                          H2'

        if(basis.lower() == 'all'):
            print 'setting up all atom overlap arrays'
            m1.set_average_vdw()
            npairs = m1.natoms() * (m1.natoms() - 1) / 2
            cutoff_array = numpy.zeros(npairs, numpy.float)
            pairs.pairs(m1.atom_vdw(), cutoff_array)
            basis_filter = 'not name[i] == None'

        elif(basis.lower() == 'backbone'):
            basis_filter = 'name[i] == "P" or name[i] == "O5\'" or name[i] == "C5\'" or name[i] == "C4\'" or name[i] == "C3\'" or name[i] == "O3\'" '
        elif(basis.lower() == 'heavy'):
            basis_filter = 'not name[i][0] == "H" '
        else:
            string_basis = string.split(basis, ",")
            number_basis_types = len(string_basis)
            for i in xrange(number_basis_types):
                if(i == 0):
                    basis_filter = 'name[i] == "' + string_basis[i] + '" '
                else:
                    basis_filter = basis_filter + \
                        ' or name[i] == "' + string_basis[i] + '" '

        error, basis_mask = m1.get_subset_mask(basis_filter)

        basis_atom = "P"
        setup_basis_filter = 'name[i] == "' + basis_atom + '"'

        error, setup_basis_mask = m1.get_subset_mask(setup_basis_filter)
        basis_m1 = sasmol.SasMol(1)
        error = m1.copy_molecule_using_mask(basis_m1, setup_basis_mask, 0)
        basis_resname = basis_m1.resname()
        basis_resid = basis_m1.resid()

        resalpha = []
        resbeta = []
        resgamma = []
        resdelta = []
        resepsilon = []
        reseta = []
        dihedral_energy.rna_initialization(resalpha, resbeta, resgamma, resdelta, resepsilon, reseta,
                                           basis_resid, basis_resname, numranges, reslow, numcont, first_last_resid, txtOutput)
        dihedral_parameters = [resalpha, resbeta,
                               resgamma, resdelta, resepsilon, reseta]
    else:
        print 'molecule_type = ', molecule_type
        print 'm1.moltype() = ', m1.moltype()
        message = 'rotation basis is not setup for molecule type = ' + \
            m1.moltype()[0]
        message += ' : stopping here'
        print_failure(message, txtOutput)
        return

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

    if(nonbondflag == 1):
        message = '> energy implementation is not implemented in this version '
        message += ' : stopping here'
        print_failure(message, txtOutput)
        return

    else:
        charge = []
        exclusionlist = []
        onefourlist = []
        vdwparam = []
        angleparam = []
        dihedralparam = []

        # angleoutfile=open(genpaths+runname+'.angle_energies.txt','a')

    hrg = 0.0
    lowestrg = 1000.0
    accepted = 0
    over = 0
    badrg = 0
    nsteps = 0
    arg = 0.0
    trg = 0.0
    badz = 0
    badc = 0

    pairdat = ['dum', 1, 1.0]
    all_rg_tally = []
    accepted_rg_tally = []

    if(plotflag == 1):
        graph = Gnuplot.Gnuplot(debug=1)
        graph.clear()
        graph('set title "Rg Results"')
        graph.xlabel('Structure Number')
        graph.ylabel('Rg (Angstrom^2)')

    failtally = 0
    acc = 0
    oofile = 'objects.txt'

    lineintxtOutput = ''.join(['=' for x in xrange(60)])
    # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))
    txtOutput.put("\n%s \n" % (lineintxtOutput))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    coor = m1.coor()

    frame = 0

    flexible_residues = get_flexible_residues(numranges, reslow, numcont)

    print 'molecule_type = ', molecule_type

    residue_rotation_indices, residue_rotation_mask = get_rotation_indices(
        m1, molecule_type, flexible_residues, txtOutput)
    name = m1.name()

    step_parameters = step.Setup()

    # get alignment sub molecule

    align_filter = 'name[i] == "' + basis_atom + '" and (resid[i] >= ' + str(
        lowres1) + ' and resid[i] <= ' + str(highres1) + ')'
    error, align_mask = m1.get_subset_mask(align_filter)

    sub_m1 = sasmol.SasMol(2)
    error = m1.copy_molecule_using_mask(sub_m1, align_mask, 0)

    if len(error) > 0:
        message = 'unable to create subset molecule for align'
        print_failure(message, txtOutput)

    com_sub_m1 = sub_m1.calccom(0)
    sub_m1.center(0)
    coor_sub_m1 = sub_m1.coor()[0]

    sub_m2 = sasmol.SasMol(4)
    error = m1.copy_molecule_using_mask(sub_m2, align_mask, 0)

    # MAIN LOOP

    if(seed[0] == 1):
        from numpy.random import RandomState
        seed_object = RandomState(seed[1])
    else:
        seed_object = -1

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
        print '.',
        sys.stdout.flush()
        # print i ; sys.stdout.flush()

        vdi, vdf, indices, this_mask = step_parameters.chooser(coor, m1, pairdat, dtheta, numranges, reslow, numcont, dihedral_parameters,
                                                               beta, residue_rotation_indices, residue_rotation_mask, nonbondflag, first_last_resid, molecule_type, seed_object)

        an = pairdat[0]
        q0 = pairdat[1]
        th = pairdat[2]
        nsteps += 1
        re = [0, 0, 0, 0.0, 0.0, lowestrg, hrg, 0, 0, []]

        this_first_coor_x = m1.coor()[0, 0, 0]
        this_first_coor_y = m1.coor()[0, 0, 1]
        this_first_coor_z = m1.coor()[0, 0, 2]
        this_final_coor_x = m1.coor()[0, -1, 0]
        this_final_coor_y = m1.coor()[0, -1, 1]
        this_final_coor_z = m1.coor()[0, -1, 2]

        if(this_first_coor_x == this_final_coor_x and this_first_coor_y == this_final_coor_y and this_first_coor_z == this_final_coor_z):
            print 'this_first_coor_x = ', this_first_coor_x
            print 'this_first_coor_y = ', this_first_coor_y
            print 'this_first_coor_z = ', this_first_coor_z
            print 'this_final_coor_x = ', this_final_coor_x
            print 'this_final_coor_y = ', this_final_coor_y
            print 'this_final_coor_z = ', this_final_coor_z

            if(i > 0):
                print 'previous_first_coor = ', previous_first_coor
                print 'previous_final_coor = ', previous_final_coor

            print 'failtally = ', failtally
            print 'goback = ', goback

            print '>>> STOPPING NOW'

            sys.stdout.flush()

            sys.exit()

        else:
            previous_first_coor = m1.coor()[0, 0, :]
            previous_final_coor = m1.coor()[0, -1, :]

        newafile = dihedral_rotate.rotate(coor, m1, q0, th, an, cutoff, lowrg, highrg, re, accepted, zflag, zcutoff, cflag, dcdoutfile, indices, this_mask, basis_mask,
                                          sub_m2, align_mask, coor_sub_m1, com_sub_m1, mask_a_array, mask_b_array, distance_array, type_array, first_last_resid, molecule_type, cutoff_array, vdw_factor)

#		print 'back from dihedral_rotate'; sys.stdout.flush()

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
                        print 'accepted = ', accepted
                        print 'dum = ', dum
                    else:
                        print '\nreloading coordinates from a previously accepted structure'
                        # print 'accepted = ',accepted
                        #dum = int(accepted*1)-1
                        # print 'dum = ',dum
                        # print 'dum+1 = ',dum+1

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

        # print 'finished loop for i = ',i,'\n' ; sys.stdout.flush()

    m1.close_dcd_write(dcdoutfile)

    rgplot = open('./' + runname + '/monomer_monte_carlo/' +
                  dcdfile + '.all_rg_results_data.txt', 'a')
    rgplot.write('# structure number (structure 1 = 1; not 0), Rg (all)\n')
    for ii in range(len(all_rg_tally)):
        rgplot.write('%i\t%f\n' %
                     (all_rg_tally[ii][0] + 1, all_rg_tally[ii][1]))
    rgplot.close()

    rgplot = open('./' + runname + '/monomer_monte_carlo/' +
                  dcdfile + '.accepted_rg_results_data.txt', 'a')
    rgplot.write(
        '# structure number (structure 1 = 1; not 0), Rg (accepted), trial number\n')
    for ii in range(len(accepted_rg_tally)):
        rgplot.write('%i\t%f\t%i\n' % (accepted_rg_tally[ii][
                     1] + 0, accepted_rg_tally[ii][2], accepted_rg_tally[ii][0] + 1))
    rgplot.close()

    outfile7 = open(genpaths + dcdfile + '.stats', 'a')
    outfile7.write('%s\t%f\t%s\t%f\n' %
                   ('lowest Rg = ', lowestrg, 'highest Rg = ', hrg))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('accepted ', accepted,
                                                   ' out of ', nsteps, ' moves : ', (accepted / float(nsteps)) * 99.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('overlapped ', over,
                                                   ' out of ', nsteps, ' moves : ', (over / float(nsteps)) * 100.0, ' %'))
    outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad rg2 ', badrg,
                                                   ' out of ', nsteps, ' moves : ', (badrg / float(nsteps)) * 100.0, ' %'))
    if(zflag == 1):
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % (
            'bad zcut ', badz, ' out of ', nsteps, ' moves : ', (badz / float(nsteps)) * 100.0, ' %'))
    if(cflag == 1):
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('constraint filter rejected ',
                                                       badc, ' out of ', nsteps, ' moves : ', (badz / float(nsteps)) * 100.0, ' %'))
    if(accepted > 0):
        outfile7.write('%s\t%f\n' %
                       ('average accepted rg2 = ', arg / (accepted)))
    else:
        outfile7.write('%s\t%f\n' % ('average accepted rg2 = ', 0.0))
    outfile7.write('%s\t%f\n' %
                   ('average total rg2 of ensemble = ', trg / (nsteps)))

    lastsasfile.close()
    # angleoutfile.close()

    if(accepted > 0):
        txtOutput.put("Average accepted rg2 = %lf\n" % (arg / (accepted)))
        txtOutput.put(
            "\nConfigurations and statistics saved in %s directory\n\n" % ('./' + genpaths))
    else:
        txtOutput.put(
            "\n NO ACCEPTED MOVES\n\n Statistics saved in %s directory\n\n" % (genpaths))

    txtOutput.put("lowest Rg = %lf\t highest Rg = %lf\n" % (lowestrg, hrg))
    txtOutput.put("accepted %d out of %d : %lf percent\n" %
                  (accepted, nsteps, (accepted / float(nsteps)) * 100.0))
    txtOutput.put("overlap check discarded %d out of %d moves : %lf percent\n" % (
        over, nsteps, (float(over) / float(nsteps)) * 100.0))
    txtOutput.put("Rg cutoffs discarded %d out of %d moves : %lf percent\n\n" % (
        badrg, nsteps, (float(badrg) / float(nsteps)) * 100.0))
    if(zflag == 1):
        txtOutput.put("Z coordinate filter discarded %d out of %d moves : %lf percent\n\n" % (
            badz, nsteps, (float(badz) / float(nsteps)) * 100.0))
    if(cflag == 1):
        txtOutput.put("constraint filter(s) discarded %d out of %d moves : %lf percent\n\n" % (
            badc, nsteps, (float(badc) / float(nsteps)) * 100.0))

    if(len(minx) > 0 and len(miny) > 0 and len(minz) > 0 and len(maxx) > 0 and len(maxy) > 0 and len(maxz) > 0):
        min_x = numpy.min(minx)
        min_y = numpy.min(miny)
        min_z = numpy.min(minz)
        max_x = numpy.max(maxx)
        max_y = numpy.max(maxy)
        max_z = numpy.max(maxz)

        txtOutput.put("minimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                      (min_x, max_x, (max_x - min_x)))

        txtOutput.put("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                      (min_y, max_y, (max_y - min_y)))
        txtOutput.put("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n" %
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

    txtOutput.put("\n%s \n" % (lineintxtOutput))
    print 'DIHEDRAL IS DONE'
    time.sleep(1.5)

    if(plotflag == 1):
        wait('\n')

    return()


if __name__ == '__main__':

    runname = 'run_0'
    dcdfile = 'run_0.dcd'
    moltype = 'protein'
    path = ''
    pdbfile = 'min3.pdb'
    trials = '50'
    goback = '10'
    temp = '300.0'
    numranges = '5'
    dtheta = '30.0,30.0,30.0,30.0,30.0'
    reslow = '123,278,354,378,408'
    numcont = '21,5,24,11,4'
    lowres1 = '284'
    highres1 = '350'
    basis = 'CA'
    cutoff = '3.0'
    lowrg = '0.0'
    highrg = '400.0'
    zflag = '0'
    zcutoff = '0.0'
    nbflag = '0'
    nbscale = '1.0'
    cflag = '1'
    confile = 'constraints.txt'
    directedmc = '0'
    psffilepath = './'
    psffilename = 'refgag.psf'
    import sassie.sasconfig as sasconfig
    parmfilepath = sasconfig._bin_path + 'toppar'
    parmfilename = 'par_all27_prot_na.inp'
    plotflag = '0'
    seed = '0,123'

    svariables = {}

    svariables['runname'] = (runname, 'string')
    svariables['dcdfile'] = (dcdfile, 'string')
    svariables['moltype'] = (moltype, 'string')
    svariables['path'] = (path, 'string')
    svariables['pdbfile'] = (pdbfile, 'string')
    svariables['trials'] = (trials, 'int')
    svariables['goback'] = (goback, 'int')
    svariables['temp'] = (temp, 'float')
    svariables['numranges'] = (numranges, 'int')
    svariables['dtheta'] = (dtheta, 'float_array')
    svariables['reslow'] = (reslow, 'int_array')
    svariables['numcont'] = (numcont, 'int_array')
    svariables['lowres1'] = (lowres1, 'int')
    svariables['highres1'] = (highres1, 'int')
    svariables['basis'] = (basis, 'string')
    svariables['cutoff'] = (cutoff, 'float')
    svariables['lowrg'] = (lowrg, 'float')
    svariables['highrg'] = (highrg, 'float')
    svariables['zflag'] = (zflag, 'int')
    svariables['zcutoff'] = (zcutoff, 'float')
    svariables['cflag'] = (cflag, 'int')
    svariables['confile'] = (confile, 'string')
    svariables['directedmc'] = (directedmc, 'float')
    svariables['nonbondflag'] = (nbflag, 'int')
    svariables['nonbondscale'] = (nbscale, 'float')
    svariables['psffilepath'] = (psffilepath, 'string')
    svariables['psffilename'] = (psffilename, 'string')
    svariables['parmfilepath'] = (parmfilepath, 'string')
    svariables['parmfilename'] = (parmfilename, 'string')
    svariables['plotflag'] = (plotflag, 'int')
    svariables['seed'] = (seed, 'int_array')

    import sassie.simulate.monte_carlo.monomer.dihedral_monte_carlo as dihedral
    import sassie.interface.input_filter as input_filter
    #import sassie.interface.generate_filter as generate_filt
    error, variables = input_filter.type_check_and_convert(svariables)

    #eflag=0 ; monflag=1

    # error=generate_filter.check_protein(variables,eflag,monflag)

    if(len(error) > 0):
        print 'error = ', error
        sys.exit()

    runname = variables['runname'][0]

    import multiprocessing

    txtQueue = multiprocessing.JoinableQueue()
    dihedralgenerate(variables, txtQueue)
    # process=multiprocessing.Process(target=dihedral.dihedralgenerate,args=(variables,txtQueue))
    # process.start()
