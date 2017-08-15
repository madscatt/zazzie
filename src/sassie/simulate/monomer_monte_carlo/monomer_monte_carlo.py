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
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monomer_monte_carlo.pairs as pairs
import sassie.simulate.monomer_monte_carlo.dihedral_rotate as dihedral_rotate
import sassie.simulate.monomer_monte_carlo.step as step


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
        MONOMER MONTE CARLO is the module that contains the functions
        that are used to generate ensembles of structures by varying
	    protein dihedral angles.

        This module calls to C / Python extension modules to speed up
        calculations (see: overlap.c).

	    REFERENCE:

    	    J. A. D. MacKerell et al.
    	    Journal of Physical Chemistry B,  102  3586-3616  (1998)

        	B. R. Brooks et al.
    	J   ournal of Computational Chemistry  4  187--217  (1983)

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'monomer_monte_carlo'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class monomer_monte_carlo_input_variables():

    def __init__(self, parent=None):
        pass


class monomer_monte_carlo():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = monomer_monte_carlo_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization(txtOutput)      #txtOutput is needed for call to dihedral_energy

        self.dihedralgenerate()

        self.epilogue()

        return


    def unpack_variables(self,variables):

        '''
        method to extract variables into system wise class instance
        '''
        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.dcdfile = variables['dcdfile'][0]
        mvars.path = variables['path'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.trials = variables['trials'][0]
        mvars.goback = variables['goback'][0]
        mvars.temp = variables['temp'][0]
        mvars.molecule_type = variables['moltype'][0]
        mvars.numranges = variables['numranges'][0]
        mvars.dtheta = variables['dtheta'][0]
        mvars.reslow = variables['reslow'][0]
        mvars.numcont = variables['numcont'][0]
        mvars.lowres1 = variables['lowres1'][0]
        mvars.highres1 = variables['highres1'][0]
        mvars.basis = variables['basis'][0]
        mvars.cutoff = variables['cutoff'][0]
        mvars.lowrg = variables['lowrg'][0]
        mvars.highrg = variables['highrg'][0]
        mvars.zflag = variables['zflag'][0]
        mvars.zcutoff = variables['zcutoff'][0]
        mvars.cflag = variables['cflag'][0]
        mvars.confile = variables['confile'][0]
        mvars.nonbondflag = variables['nonbondflag'][0]
        mvars.nonbondscale = variables['nonbondscale'][0]
        mvars.psffilepath = variables['psffilepath'][0]
        mvars.psffilename = variables['psffilename'][0]
        mvars.parmfilepath = variables['parmfilepath'][0]
        mvars.parmfilename = variables['parmfilename'][0]
        mvars.plotflag = variables['plotflag'][0]
        mvars.directedmc = variables['directedmc'][0]
        mvars.seed = variables['seed'][0]

        log.debug(vars(mvars))

        return 

#mvars: runname, dcdfile, path, pdbfile, trials, goback, temp, molecule_type, numranges, dtheta, reslow, numcont, lowres1, highres1, basis, cutoff, lowrg, highrg, zflag, zcutoff, cflag, confile, nonbondflag, nonbondscale, psffilepath, psffilename, parmfilepath, parmfilename, plotflag, directedmc, seed
#avars: m1, flexible_residues, lastsasfile, dcdoutfile, first_last_resid, cutoff_array, vdw_factor, basis_atom, pairdat, all_rg_tally, accepted_rg_tally, graph, hrg, lowestrg, accepted, over, badrg, nsteps, arg, trg, badz, badc, mask_a_array, mask_b_array, beta, dihedral_parameters, basis_mask, distance_array, type_array, genpaths, minx, miny, minz, maxx, maxy, maxz



#    pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return


    def wait(self,sti=None, prompt='Plot will clear in 2 seconds ...\n'):
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


    def get_flexible_residues(self):
        '''
        Method to determine residues that are going to be rotated.
        '''
        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        flexible_residues = []
        templist = []
        for i in xrange(mvars.numranges):
            thisreslow = mvars.reslow[i]
            thisnumcont = mvars.numcont[i]
            thisrange = numpy.arange(thisnumcont) + thisreslow
            templist.append(thisrange.tolist())

        flexible_residues = [item for sublist in templist for item in sublist]

        message = 'flexible_residues = '+str(flexible_residues)
        log.debug(message)

        return flexible_residues


    def get_rotation_indices(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui

        log.debug('getting rotation indices for molecule')

        residue_rotation_indices = {}
        residue_rotation_mask = {}

        if(mvars.molecule_type == 'protein'):

            mtype = 0
            mask = avars.m1.get_dihedral_subset_mask(avars.flexible_residues, mtype)

            for i in xrange(len(mask)):
                this_mask = mask[i][:]
                q0 = avars.flexible_residues[i]
                residue_rotation_mask[q0] = this_mask.tolist()
                indices = avars.m1.get_indices_from_mask(this_mask)
                residue_rotation_indices[q0] = indices.tolist()

#            print 'residue_rotation_indices = ',residue_rotation_indices
#            print 'residue_rotation_mask = ',residue_rotation_mask

        elif(mvars.molecule_type == 'rna'):

            mtype = 1
            mask = avars.m1.get_dihedral_subset_mask(avars.flexible_residues, mtype)

            for i in xrange(len(mask)):
                this_mask = mask[i][:]
                q0 = avars.flexible_residues[i]
                residue_rotation_mask[q0] = this_mask.tolist()
                indices = avars.m1.get_indices_from_mask(this_mask)
                residue_rotation_indices[q0] = indices.tolist()

#            print 'residue_rotation_indices = ',residue_rotation_indices
#            print 'residue_rotation_mask = ',residue_rotation_mask

        elif(mvars.molecule_type == 'dna'):
            message = 'rotation basis set not defined for molecule type = ' + mvars.molecule_type
            pgui(message)
            pass

        else:
            message = 'rotation basis set not defined for molecule type = ' + mvars.molecule_type
            pgui(message)

        log.debug('done getting rotation indices')

        return residue_rotation_indices, residue_rotation_mask


    def evaluate_rg(self,rg_difference_list, directed_rg_list, accepted_rg_list, this_rg_difference, this_rg, accepted):

        maximum_value = max(rg_difference_list)

        if(maximum_value > this_rg_difference):
            index = rg_difference_list.index(maximum_value)
            rg_difference_list[index] = this_rg_difference
            directed_rg_list[index] = this_rg
            accepted_rg_list[index] = accepted

        return


    def initialization(self,txtOutput):
        '''
        method to prepare for dihedral generate
        '''
#   txtOutput is needed for call to dihedral_energy

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        pdbfilename = mvars.pdbfile

        if(mvars.runname[-1] == '/'):
            lin = len(mvars.runname)
            mvars.runname = mvars.runname[:lin - 1]


        direxist = os.path.exists(mvars.runname)
        if(direxist == 0):
            os.system('mkdir -p ' + mvars.runname)

        genpath = mvars.runname + '/monomer_monte_carlo'
        avars.genpaths = genpath + '/'
        direxist = os.path.exists(genpath)
        if(direxist == 0):
            os.system('mkdir -p ' + genpath)

        if(mvars.path != ""):
            if mvars.path[-1] != "/":
                mvars.path += '/'

        cpfile = os.path.join(mvars.path, pdbfilename)
        cpst = 'cp ' + cpfile + ' ' + avars.genpaths
        os.system(cpst)

        mvars.parmfilepath += '/'

        fileexist = os.path.exists('.last_sas')
        if(fileexist == 1):
            os.system('mv -f .last_sas .last_sas_bu')
        avars.lastsasfile = open('./.last_sas', 'w')
        avars.lastsasfile.write('run_name\t' + mvars.runname + '\n')
        avars.lastsasfile.write('pdb_name\t' + mvars.pdbfile + '\n')
        avars.lastsasfile.write('dcd_name\t' + mvars.dcdfile + '\n')

        kb = 1.380658E-23  # J/K
        avars.beta = 1.0 / (mvars.temp * kb)

        avars.m1 = sasmol.SasMol(0)
        m2 = sasmol.SasMol(1)       #not used?

#       The 2.0 version of read_pdb with the saspdbrx_topology=True option gives different result than 1.0 version with same option
#       The last few atoms of the last residue aren't included in the coordinates if this option is True
#       But they are read in the 1.0 version of read_pdb.  It is necessary to disable this option to pass the 1.0 tests.
#       THIS ISSUE NEEDS TO BE ADDRESSED BEFORE FURTHER 2.0 DEVELOPMENT! 
#        avars.m1.read_pdb(mvars.path + mvars.pdbfile, saspdbrx_topology=True)  
        avars.m1.read_pdb(mvars.path + mvars.pdbfile)
        nf1 = avars.m1.number_of_frames()       #not used?
        avars.m1.calccom(0)                     #not used?
#        print avars.m1.calccom(0)
#        print avars.m1.coor()

        avars.dcdoutfile = avars.m1.open_dcd_write(avars.genpaths + mvars.dcdfile)

        resid = avars.m1.resid()
        avars.first_last_resid = [resid[0], resid[-1]]

        message = 'first_last_resid = '+str(avars.first_last_resid)
        log.debug(message)

        avars.cutoff_array = []
        avars.vdw_factor = 0.3

        if(mvars.molecule_type == 'protein'):

            if(mvars.basis.lower() == 'all'):
                log.debug('setting up all atom overlap arrays')
                avars.m1.set_average_vdw()
                npairs = avars.m1.natoms() * (avars.m1.natoms() - 1) / 2
                avars.cutoff_array = numpy.zeros(npairs, numpy.float)
                message = 'cutoff_array: '+str(avars.cutoff_array)
                log.debug(message) 
                pairs.pairs(avars.m1.atom_vdw(), avars.cutoff_array)
                basis_filter = 'not name[i] == None'

            elif(mvars.basis.lower() == 'backbone'):
                basis_filter = 'name[i] == "N" or name[i] == "CA" or name[i] == "C" '
            elif(mvars.basis.lower() == 'heavy'):
                basis_filter = 'not name[i][0] == "H" '
            else:
                string_basis = string.split(mvars.basis, ",")
                number_basis_types = len(string_basis)
                for i in xrange(number_basis_types):
                    if(i == 0):
                        basis_filter = 'name[i] == "' + string_basis[i] + '" '
                    else:
                        basis_filter = basis_filter + \
                            ' or name[i] == "' + string_basis[i] + '" '

            log.debug('basis_filter = %s\n' %(basis_filter))

            error, avars.basis_mask = avars.m1.get_subset_mask(basis_filter)

            avars.basis_atom = "CA"
            setup_basis_filter = 'name[i] == "' + avars.basis_atom + '"'
            error, setup_basis_mask = avars.m1.get_subset_mask(setup_basis_filter)

            basis_m1 = sasmol.SasMol(1)
            error = avars.m1.copy_molecule_using_mask(basis_m1, setup_basis_mask, 0)
            basis_resname = basis_m1.resname()
            basis_resid = basis_m1.resid()

            respsi = []
            resphi = []
            dihedral_energy.protein_initialization(
                respsi, resphi, basis_resid, basis_resname, mvars.numranges, mvars.reslow, mvars.numcont, avars.first_last_resid, txtOutput)
            avars.dihedral_parameters = [respsi, resphi]

        elif(mvars.molecule_type == 'rna'):

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

            if(mvars.basis.lower() == 'all'):
                log.debug('setting up all atom overlap arrays')
                avars.m1.set_average_vdw()
                npairs = avars.m1.natoms() * (avars.m1.natoms() - 1) / 2
                avars.cutoff_array = numpy.zeros(npairs, numpy.float)
                pairs.pairs(avars.m1.atom_vdw(), avars.cutoff_array)
                basis_filter = 'not name[i] == None'

            elif(mvars.basis.lower() == 'backbone'):
                basis_filter = 'name[i] == "P" or name[i] == "O5\'" or name[i] == "C5\'" or name[i] == "C4\'" or name[i] == "C3\'" or name[i] == "O3\'" '
            elif(mvars.basis.lower() == 'heavy'):
                basis_filter = 'not name[i][0] == "H" '
            else:
                string_basis = string.split(mvars.basis, ",")
                number_basis_types = len(string_basis)
                for i in xrange(number_basis_types):
                    if(i == 0):
                        basis_filter = 'name[i] == "' + string_basis[i] + '" '
                    else:
                        basis_filter = basis_filter + \
                            ' or name[i] == "' + string_basis[i] + '" '

            error, avars.basis_mask = avars.m1.get_subset_mask(basis_filter)

            avars.basis_atom = "P"
            setup_basis_filter = 'name[i] == "' + avars.basis_atom + '"'

            error, setup_basis_mask = avars.m1.get_subset_mask(setup_basis_filter)
            basis_m1 = sasmol.SasMol(1)
            error = avars.m1.copy_molecule_using_mask(basis_m1, setup_basis_mask, 0)
            basis_resname = basis_m1.resname()
            basis_resid = basis_m1.resid()

            resalpha = []
            resbeta = []
            resgamma = []
            resdelta = []
            resepsilon = []
            reseta = []
            dihedral_energy.rna_initialization(resalpha, resbeta, resgamma, resdelta, resepsilon, reseta,
                                           basis_resid, basis_resname, mvars.numranges, mvars.reslow, mvars.numcont, avars.first_last_resid, txtOutput)
            avars.dihedral_parameters = [resalpha, resbeta,
                               resgamma, resdelta, resepsilon, reseta]
        else:
            pgui('molecule_type = %s\n' %(mvars.molecule_type))
            pgui('m1.moltype() = %s\n' %(avars.m1.moltype()))
            message = 'rotation basis is not setup for molecule type = ' + \
                avars.m1.moltype()[0]
            message += ' : stopping here'
            pgui(message)
            return

        if(mvars.cflag == 1):

            filter_flag = 0
            error, constraint_basis1_array, constraint_basis2_array, avars.distance_array, avars.type_array = constraints.read_constraints(
                avars.m1, mvars.confile, filter_flag)

            avars.mask_a_array = []
            avars.mask_b_array = []

            for i in xrange(len(avars.distance_array)):
#                print constraint_basis1_array[i]
#                print constraint_basis2_array[i]
#                print avars.distance_array[i]
#                print avars.type_array[i]

                error, local_mask_a_array = avars.m1.get_subset_mask(
                    constraint_basis1_array[i])
                error, local_mask_b_array = avars.m1.get_subset_mask(
                    constraint_basis2_array[i])

                avars.mask_a_array.append(local_mask_a_array)
                avars.mask_b_array.append(local_mask_b_array)

        else:
            avars.mask_a_array = []
            avars.mask_b_array = []
            avars.distance_array = []
            avars.type_array = []

        if(mvars.nonbondflag == 1):
            message = '> energy implementation is not implemented in this version '
            message += ' : stopping here'
            pgui(message)
            return

        else:
            charge = []
            exclusionlist = []
            onefourlist = []
            vdwparam = []
            angleparam = []
            dihedralparam = []

        avars.hrg = 0.0
        avars.lowestrg = 1000.0
        avars.accepted = 0
        avars.over = 0
        avars.badrg = 0
        avars.nsteps = 0
        avars.arg = 0.0
        avars.trg = 0.0
        avars.badz = 0
        avars.badc = 0

        avars.pairdat = ['dum', 1, 1.0]
        avars.all_rg_tally = []
        avars.accepted_rg_tally = []

        if(mvars.plotflag == 1):
            avars.graph = Gnuplot.Gnuplot(debug=1)
            avars.graph.clear()
            avars.graph('set title "Rg Results"')
            avars.graph.xlabel('Structure Number')
            avars.graph.ylabel('Rg (Angstrom^2)')


        log.debug(vars(mvars))
        log.debug(vars(avars))

        return

    def dihedralgenerate(self):

        '''
        DIHEDRAL generates ensembles of structures by varying
        protein dihedral angles.

        INPUT:	variable descriptions

        runname:        string      project name                          
        path:           string      input file path                 
        dcdfile:        string      name of output dcd file containing accepted structures       
        moltype:        string      molecule type ('protein' or 'rna')
        pdbfile:        string      name of input pdb file containing intial structure
        trials:         integer     number of Monte Carlo move attempts
        goback:         integer     number of failed Monte Carlo attempts before returning to previously accepted structure
        temp            float       run temperature (K)
        numranges       integer     number of flexible regions
        dtheta          float_array maximum angle that torsion can sample (in each flexible region)
        reslow          int_array   low residue number for each flexible region
        numcont         int_array   number of contiguous residues per flexible region (not enetered directly; parsed from entered residue range in GenApp)
        lowres1         integer     low residue for (non-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
        highres1        integer     high residue for (no-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
        basis           string      type of basis for overlap check ("all", "heavy", "backbone" or specific atom name, i.e., "CA")
        cutoff          float       overlap cutoff distance (all=0.8, heavy=0.8, backbone=1.0, specific atom=user's choice)
        lowrg           float       low Rg cutoff value if Advanced Input is chosen
        highrg          float       high Rg cutoff value if Advanced Input is chosen
        zflag           integer     enable zcutoff flag (0=no, 1=yes)
        zcutoff         float       zcutoff value (discard structures with any z-axis coordinates less than this value)
        cflag           integer     enable atomic constraint flag (0=no, 1=yes)
        confile         string      name of file describing additional constraints to check before accepting a structure
        directedmc      float       non-zero Rg value to guide Monte Carlo run; 0=no directed Monte Carlo (used if Advanced Input is chosen)
        
        Input options not implemented:

        nonbondflag     integer     flag for nonbonded option
        nonbondedscale  float       nonbondedscale value
        psffilepath     string      path to psf file
        psffilename     string      psf file name
        parmfilepath    string      path to CHARMM parameter file
        parmfilename    string      name of CHARMM parameter file
        plotflag        integer     option to plot structure number vs Rg

        OUTPUT: files stored in "runname"/monomer_monte_carlo directory:

        original PDB file with molecular structure data
        DCD file containing accepted structures
        file containing Rg values for all trial structures
        file containing Rg value for accepted structures
        file containing run statistics
        '''

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in dihedralgenerate')

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))


        alignflag = '1'     #not used?
        saveflag = '0'      #not used?

        failtally = 0
        acc = 0
        oofile = 'objects.txt'

        lineintxtOutput = ''.join(['=' for x in xrange(60)])
#        ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        pgui("\n%s \n" % (lineintxtOutput))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        coor = avars.m1.coor()

#        print 'coor: ', avars.m1.coor()

        frame = 0

        avars.flexible_residues = self.get_flexible_residues()


        residue_rotation_indices, residue_rotation_mask = self.get_rotation_indices()
#        print 'residue_rotation_indices after call = ',residue_rotation_indices
#        print 'residue_rotation_mask after call = ',residue_rotation_mask        
        name = avars.m1.name()

        step_parameters = step.Setup()

#       get alignment sub molecule

        align_filter = 'name[i] == "' + avars.basis_atom + '" and (resid[i] >= ' + str(
            mvars.lowres1) + ' and resid[i] <= ' + str(mvars.highres1) + ')'
        error, align_mask = avars.m1.get_subset_mask(align_filter)

        sub_m1 = sasmol.SasMol(2)
        error = avars.m1.copy_molecule_using_mask(sub_m1, align_mask, 0)

        if len(error) > 0:
            message = 'unable to create subset molecule for align'
            pgui(message)

        com_sub_m1 = sub_m1.calccom(0)
        sub_m1.center(0)
        coor_sub_m1 = sub_m1.coor()[0]

        sub_m2 = sasmol.SasMol(4)
        error = avars.m1.copy_molecule_using_mask(sub_m2, align_mask, 0)

#       MAIN LOOP

        if(mvars.seed[0] == 1):
            from numpy.random import RandomState
            seed_object = RandomState(mvars.seed[1])
        else:
            seed_object = -1

        avars.minx = []
        avars.miny = []
        avars.minz = []
        avars.maxx = []
        avars.maxy = []
        avars.maxz = []

        if(mvars.directedmc > 0):
            rg_difference_list = []
            directed_rg_list = []
            accepted_rg_list = []
            rg_list_length = 10  # hardwired

        for i in range(mvars.trials):
            print '.',
            sys.stdout.flush()
#           print i ; sys.stdout.flush()

            vdi, vdf, indices, this_mask = step_parameters.chooser(coor, avars.m1, avars.pairdat, mvars.dtheta, mvars.numranges, mvars.reslow, mvars.numcont, avars.dihedral_parameters, avars.beta, residue_rotation_indices, residue_rotation_mask, mvars.nonbondflag, avars.first_last_resid, mvars.molecule_type, seed_object)

            an = avars.pairdat[0]
            q0 = avars.pairdat[1]
            th = avars.pairdat[2]
            avars.nsteps += 1
            re = [0, 0, 0, 0.0, 0.0, avars.lowestrg, avars.hrg, 0, 0, []]

            this_first_coor_x = avars.m1.coor()[0, 0, 0]
            this_first_coor_y = avars.m1.coor()[0, 0, 1]
            this_first_coor_z = avars.m1.coor()[0, 0, 2]
            this_final_coor_x = avars.m1.coor()[0, -1, 0]
            this_final_coor_y = avars.m1.coor()[0, -1, 1]
            this_final_coor_z = avars.m1.coor()[0, -1, 2]

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
                print 'goback = ', mvars.goback

                print '>>> STOPPING NOW'

                sys.stdout.flush()

                sys.exit()

            else:
                previous_first_coor = avars.m1.coor()[0, 0, :]
                previous_final_coor = avars.m1.coor()[0, -1, :]

#                print 'previous_first_coor = ', previous_first_coor
#                print 'previous_final_coor = ', previous_final_coor                

            newafile = dihedral_rotate.rotate(coor, avars.m1, q0, th, an, mvars.cutoff, mvars.lowrg, mvars.highrg, re, avars.accepted, mvars.zflag, mvars.zcutoff, mvars.cflag, avars.dcdoutfile, indices, this_mask, avars.basis_mask, sub_m2, align_mask, coor_sub_m1, com_sub_m1, avars.mask_a_array, avars.mask_b_array, avars.distance_array, avars.type_array, avars.first_last_resid, mvars.molecule_type, avars.cutoff_array, avars.vdw_factor)

#            print 'back from dihedral_rotate'; sys.stdout.flush()

            avars.accepted = avars.accepted + re[0]
            avars.over = avars.over + re[1]
            avars.badrg = avars.badrg + re[2]
            rg_value = re[3]
            avars.trg = avars.trg + re[3]
            avars.arg = avars.arg + re[4]
            avars.lowestrg = re[5]
            avars.hrg = re[6]
            avars.badz = avars.badz + re[7]
            avars.badc = avars.badc + re[8]
            if(len(re[9]) > 0):
                minmax = re[9]

                avars.minx.append(minmax[0][0])
                avars.miny.append(minmax[0][1])
                avars.minz.append(minmax[0][2])
                avars.maxx.append(minmax[1][0])
                avars.maxy.append(minmax[1][1])
                avars.maxz.append(minmax[1][2])

            avars.all_rg_tally.append([i, rg_value])

            if(re[0] == 1):
                avars.accepted_rg_tally.append([i, avars.accepted, rg_value])
                if(mvars.directedmc > 0):
                    if(len(rg_difference_list) <= rg_list_length):
                        this_rg_difference = abs(rg_value - mvars.directedmc)
                        rg_difference_list.append(this_rg_difference)
                        directed_rg_list.append(rg_value)
                        accepted_rg_list.append(avars.accepted)
                    else:
                        this_rg_difference = abs(rg_value - mvars.directedmc)
                        self.evaluate_rg(rg_difference_list, directed_rg_list,
                                accepted_rg_list, this_rg_difference, rg_value, avars.accepted)


            if(re[0] == 0):
                if(failtally == mvars.goback):
                    failtally = 0
                    if(avars.accepted > 0):
                        if(mvars.seed[0] == 1):
                            ran_num = seed_object.rand()
                            log.debug('ran_num: %f\n' % (ran_num))
                            dum = int(avars.accepted * ran_num) - 1

                        elif(mvars.directedmc > 0):
                            local_rg_list_length = len(directed_rg_list)
                            ran_num = random.randrange(0, local_rg_list_length)
                            dum = accepted_rg_list[ran_num]

                        else:
                            dum = int(avars.accepted * random.random()) - 1
                        if(dum == -1):
                            pgui('\nreloading coordinates from original starting structure')
#                            avars.m1.read_pdb(mvars.path + mvars.pdbfile, fastread=True,saspdbrx_topology=True)  #to pass 1.0 tests
                            avars.m1.read_pdb(mvars.path + mvars.pdbfile, fastread=True)
                            coor = avars.m1.coor()
                            pgui('accepted = %i\n' % (avars.accepted))
                            pgui('dum = %i\n' % (dum))
                        else:
                            pgui('\nreloading coordinates from a previously accepted structure')
                            # print 'accepted = ',avars.accepted
                            #dum = int(avars.accepted*1)-1
                            # print 'dum = ',dum
                            # print 'dum+1 = ',dum+1

                            avars.m1.read_single_dcd_step(avars.genpaths + mvars.dcdfile, dum + 1)
#                           avars.m1.read_single_dcd_step(avars.genpaths+mvars.dcdfile,dum)
                            coor = avars.m1.coor()
                    else:
                        pgui('\n>>>>>reloading coordinates from original starting structure')
#                        avars.m1.read_pdb(mvars.path + mvars.pdbfile, fastread=True,saspdbrx_topology=True)    #to pass 1.0 tests
                        avars.m1.read_pdb(mvars.path + mvars.pdbfile, fastread=True)
                        coor = avars.m1.coor()
                else:
                    failtally = failtally + 1

            if(((i + 1) % (float(mvars.trials) / 100.0) == 0 or (mvars.trials < 10))):
                fraction_done = (float(i + 1) / float(mvars.trials))
                progress_string = '\nCOMPLETED ' + \
                    str(i + 1) + ' of ' + str(mvars.trials) + ' : ' + \
                    str(fraction_done * 100.0) + ' % done'
                pgui('%s\n' % progress_string)
                accepted_percent = (float(avars.accepted) / avars.nsteps) * 100.0 
                accepted_string = str(avars.accepted) + ' configurations accepted out of ' + str(avars.nsteps) +': ' +str(accepted_percent) + ' %\n\n'
                pgui(accepted_string)
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

            if(i > 9):
                if((i + 1) % (mvars.trials / 10) == 0 and avars.accepted > 0 and i + 1 > 10):
                    if(mvars.plotflag == 1):
                        avars.graph.plot(Gnuplot.Data(avars.all_rg_tally, using='1:2 w p ps 4', title='all Rg'), Gnuplot.Data(
                            avars.accepted_rg_tally, using='1:3 w lp pt 5 ps 2', title='accepted'))
                fraction_done = (float(i + 1) / float(mvars.trials))
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

            elif(avars.accepted > 0):
                if(mvars.plotflag == 1):
                    avars.graph.plot(Gnuplot.Data(avars.all_rg_tally, using='1:2 w p ps 4', title='all Rg'), Gnuplot.Data(
                        avars.accepted_rg_tally, using='1:3 w lp pt 5 ps 2', title='accepted'))

                fraction_done = (float(i + 1) / float(mvars.trials))
                report_string = 'STATUS\t' + str(fraction_done)
                pgui(report_string)

#        print 'finished loop for i = ',i,'\n' ; sys.stdout.flush()

        avars.m1.close_dcd_write(avars.dcdoutfile)

        avars.lastsasfile.close()


        return


    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug('in epilogue')

        log.debug(vars(mvars))
        log.debug(vars(avars))

        rgplot = open('./' + mvars.runname + '/monomer_monte_carlo/' +
                  mvars.dcdfile + '.all_rg_results_data.txt', 'a')
        rgplot.write('# structure number (structure 1 = 1; not 0), Rg (all)\n')
        for ii in range(len(avars.all_rg_tally)):
            rgplot.write('%i\t%f\n' %
                     (avars.all_rg_tally[ii][0] + 1, avars.all_rg_tally[ii][1]))
        rgplot.close()

        rgplot = open('./' + mvars.runname + '/monomer_monte_carlo/' +
                  mvars.dcdfile + '.accepted_rg_results_data.txt', 'a')
        rgplot.write(
            '# structure number (structure 1 = 1; not 0), Rg (accepted), trial number\n')
        for ii in range(len(avars.accepted_rg_tally)):
            rgplot.write('%i\t%f\t%i\n' % (avars.accepted_rg_tally[ii][
                     1] + 0, avars.accepted_rg_tally[ii][2], avars.accepted_rg_tally[ii][0] + 1))
        rgplot.close()

        outfile7 = open(avars.genpaths + mvars.dcdfile + '.stats', 'a')
        outfile7.write('%s\t%f\t%s\t%f\n' %
                   ('lowest Rg = ', avars.lowestrg, 'highest Rg = ', avars.hrg))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('accepted ', avars.accepted,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.accepted / float(avars.nsteps)) * 99.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('overlapped ', avars.over,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.over / float(avars.nsteps)) * 100.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad rg2 ', avars.badrg,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.badrg / float(avars.nsteps)) * 100.0, ' %'))
        if(mvars.zflag == 1):
            outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % (
                'bad zcut ', avars.badz, ' out of ', avars.nsteps, ' moves : ', (avars.badz / float(avars.nsteps)) * 100.0, ' %'))
        if(mvars.cflag == 1):
            outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('constraint filter rejected ',
                                                       avars.badc, ' out of ', avars.nsteps, ' moves : ', (avars.badz / float(avars.nsteps)) * 100.0, ' %'))
        if(avars.accepted > 0):
            outfile7.write('%s\t%f\n' %
                       ('average accepted rg2 = ', avars.arg / (avars.accepted)))
        else:
            outfile7.write('%s\t%f\n' % ('average accepted rg2 = ', 0.0))
        outfile7.write('%s\t%f\n' %
                   ('average total rg2 of ensemble = ', avars.trg / (avars.nsteps)))

        if(avars.accepted > 0):
            pgui("Average accepted rg2 = %lf\n" % (avars.arg / (avars.accepted)))
            pgui(
                "\nConfigurations and statistics saved in %s directory\n\n" % ('./' + avars.genpaths))
        else:
            pgui(
                "\n NO ACCEPTED MOVES\n\n Statistics saved in %s directory\n\n" % (avars.genpaths))

        pgui("lowest Rg = %lf\t highest Rg = %lf\n" % (avars.lowestrg, avars.hrg))
        pgui("accepted %d out of %d : %lf percent\n" %
                (avars.accepted, avars.nsteps, (avars.accepted / float(avars.nsteps)) * 100.0))
        pgui("overlap check discarded %d out of %d moves : %lf percent\n" % (
                avars.over, avars.nsteps, (float(avars.over) / float(avars.nsteps)) * 100.0))
        pgui("Rg cutoffs discarded %d out of %d moves : %lf percent\n\n" % (
                avars.badrg, avars.nsteps, (float(avars.badrg) / float(avars.nsteps)) * 100.0))
        if(mvars.zflag == 1):
            pgui("Z coordinate filter discarded %d out of %d moves : %lf percent\n\n" % (
                avars.badz, avars.nsteps, (float(avars.badz) / float(avars.nsteps)) * 100.0))
        if(mvars.cflag == 1):
            pgui("constraint filter(s) discarded %d out of %d moves : %lf percent\n\n" % (
                avars.badc, avars.nsteps, (float(avars.badc) / float(avars.nsteps)) * 100.0))

        if(len(avars.minx) > 0 and len(avars.miny) > 0 and len(avars.minz) > 0 and len(avars.maxx) > 0 and len(avars.maxy) > 0 and len(avars.maxz) > 0):
            min_x = numpy.min(avars.minx)
            min_y = numpy.min(avars.miny)
            min_z = numpy.min(avars.minz)
            max_x = numpy.max(avars.maxx)
            max_y = numpy.max(avars.maxy)
            max_z = numpy.max(avars.maxz)

            pgui("minimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                      (min_x, max_x, (max_x - min_x)))

            pgui("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                      (min_y, max_y, (max_y - min_y)))
            pgui("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n" %
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

        self.run_utils.clean_up(log)

        lineintxtOutput = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" % (lineintxtOutput))
        pgui('DIHEDRAL IS DONE')
        time.sleep(1.5)

        if(mvars.plotflag == 1):
            self.wait('\n')

        return




