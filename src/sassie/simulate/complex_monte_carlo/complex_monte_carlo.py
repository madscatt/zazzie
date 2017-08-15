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
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monomer_monte_carlo.pairs as pairs
import sassie.simulate.monomer_monte_carlo.step as step
import sassie.simulate.complex_monte_carlo.nmer_overlap_check as nmer_overlap_check
import sassie.simulate.complex_monte_carlo.nmer_nrotate as nmer_nrotate


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
	COMPLEX MONTE CARLO is the module that contains the functions
	that are used to generate ensembles of structures by varying
	protein dihedral angles.  This particular version allows multiple
	flexible proteins in the presence of non-flexible proteins and
	nucleic acids.

	This module calls to C / Python extension modules to speed up
	calculations.
'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'complex_monte_carlo'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class complex_monte_carlo_input_variables():

    def __init__(self, parent=None):
        pass


class complex_monte_carlo():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, psegvariables, txtOutput):

        self.mvars = module_variables()

        self.avars = complex_monte_carlo_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.process_input_variables(psegvariables)

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
        mvars.nsegments = variables['nsegments'][0]
        mvars.segbasis = variables['segbasis'][0]
        mvars.npsegments = variables['npsegments'][0]
        mvars.flpsegname = variables['flpsegname'][0]
        mvars.sseglow = variables['seglow'][0]
        mvars.sseghigh = variables['seghigh'][0]
        mvars.lowrg = variables['lowrg'][0]
        mvars.highrg = variables['highrg'][0]
        mvars.zflag = variables['zflag'][0]
        mvars.zcutoff = variables['zcutoff'][0]
        mvars.cflag = variables['cflag'][0]
        mvars.confile = variables['confile'][0]
        mvars.plotflag = variables['plotflag'][0]
        mvars.directedmc = variables['directedmc'][0]
        mvars.seed = variables['seed'][0]

        log.debug(vars(mvars))

        return 

#mvars: runname, dcdfile, path, pdbfile, trials, goback, temp, nsegments, segbasis, npsegments, flpsegname, sseglow, sseghigh, lowrg, highrg, zflag, zcutoff, cflag, confile, plotflag, directedmc, seed, allsith, allsnumranges, allsrlow, allsrnum

#avars: amoltype, abasis, anumranges, aith, arlow, arnum, flexible_segments, lastsasfile, genpaths, beta, m1, dcdoutfile, asegs, first_last_resid, all_segment_mask, all_segment_full_mask, all_segment_basis_full_mask, basis_full_mask, all_segment_mol, keyword_basis, mask_a_array, mask_b_array, distance_array, type_array, all_flexible_align_mask, all_flexible_coor_sub_m1, all_flexible_com_sub_m1, all_flexible_sub_m2, cutoff, interpairs, npairs, interatom, interres, flexible_dihedral_parameters, all_flexible_basis_mask, all_flexible_residues, all_flexible_residue_rotation_indices, all_flexible_residue_rotation_mask, graph, accepted, nsteps, all_rg_tally, accepted_rg_tally, rg_difference_list, directed_rg_list, accepted_rg_list, this_rg_difference, rg_value, arg, minx, miny, minz, maxx, maxy, maxz, trg, lowestrg, hrg, over, badrg, badz, badc

    def process_input_variables(self,psegvariables):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in process_input_variables')

        log.debug('psegvariables: %s' %(psegvariables))    
    
        mvars.allsith = []
        mvars.allsnumranges = []
        mvars.allsrlow = []
        mvars.allsrnum = []
        mvars.allmoltype = []

        log.debug('len(psegvariables): %i' % len(psegvariables))
        log.debug('npsegments: %i' % (mvars.npsegments))

        for i in range(len(psegvariables)):
            mvars.allsnumranges.append(psegvariables[i][0])
            mvars.allsith.append(psegvariables[i][1])
            mvars.allsrlow.append(psegvariables[i][2])
            mvars.allsrnum.append(psegvariables[i][3])
            mvars.allmoltype.append(psegvariables[i][4])

        mvars.segbasis.strip()
        avars.abasis = [item.strip() for item in string.split(mvars.segbasis, ',')]

        avars.aith = []
        avars.anumranges = []
        avars.arlow = []
        avars.arnum = []
        avars.amoltype = []

        for i in range(len(mvars.allsith)):
            linith = string.split(mvars.allsith[i], ',')
            locith = []
            for i in range(len(linith)):
                tith = linith[i]
                fith = locale.atof(tith)
                if(fith > 180.0):
                    fith = 180.0
                elif(fith < 0.0):
                    fith = 0.0
                locith.append(fith)
            avars.aith.append(locith)

        for i in range(len(mvars.allsnumranges)):
            nr = locale.atoi(mvars.allsnumranges[i])
            avars.anumranges.append(nr)

        for i in range(len(mvars.allsrlow)):
            linrlow = string.split(mvars.allsrlow[i], ',')
            linrnum = string.split(mvars.allsrnum[i], ',')
            rlow = []
            rnum = []
            for k in range(len(linrlow)):
                trlow = locale.atoi(linrlow[k])
                trnum = locale.atoi(linrnum[k])
                rlow.append(trlow)
                rnum.append(trnum)
#        print 'rlow = ',rlow
#        print 'rnum = ',rnum
            avars.arlow.append(rlow)
            avars.arnum.append(rnum)

        for i in range(len(psegvariables)):
            moltype = mvars.allmoltype[i].strip()
            avars.amoltype.append(moltype)

        '''
        	print 'anumranges = ',avars.anumranges
	    print 'aith = ',avars.aith
	    print 'arlow = ',avars.arlow
	    print 'arnum = ',avars.arnum
        	'''

        raw_flexible_segments = string.split(mvars.flpsegname, ",")

        avars.flexible_segments = []
        for fp in raw_flexible_segments:
            avars.flexible_segments.append(fp.strip())

#        print 'flexible_segments = ',avars.flexible_segments

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return

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


    def alignment_initialization(self):

        log = self.log
        log.debug('in alignment_initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars

        avars.all_flexible_align_mask = []

        avars.all_flexible_coor_sub_m1 = []

        avars.all_flexible_com_sub_m1 = []

        avars.all_flexible_sub_m2 = []

        for i in xrange(len(avars.flexible_segments)):
            this_segment = avars.flexible_segments[i]
            idx = avars.asegs.index(this_segment)

            local_m1 = avars.all_segment_mol[idx]

            if(local_m1.moltype()[0] == 'protein'):
                this_basis = 'CA'
            elif(local_m1.moltype()[0] == 'rna' or local_m1.moltype()[0] == 'dna'):
                this_basis = 'P'
            else:
                pgui('NO ALIGNMENT BASIS ATOM DEFINED FOR SEGNAME')
                sys.exit()

            # TODO need to handle the exception in complex_filter.py
            # ONLY protein and RNA need this alignment

        # get alignment sub molecule

            align_filter = 'name[i] == "' + this_basis + '" and (segname[i] == "' + this_segment + '") and (resid[i] >= ' + str(
                avars.seglow[i]) + ' and resid[i] <= ' + str(avars.seghigh[i]) + ')'
            log.debug('align filter = %s\n' % (align_filter))
            error, align_mask = local_m1.get_subset_mask(align_filter)

            avars.all_flexible_align_mask.append(align_mask)

            sub_m1 = sasmol.SasMol(2)
            error = local_m1.copy_molecule_using_mask(sub_m1, align_mask, 0)
            com_sub_m1 = sub_m1.calccom(0)
            sub_m1.center(0)
            coor_sub_m1 = sub_m1.coor()[0]

            avars.all_flexible_coor_sub_m1.append(coor_sub_m1)
            avars.all_flexible_com_sub_m1.append(com_sub_m1)

            sub_m2 = sasmol.SasMol(4)
            error = local_m1.copy_molecule_using_mask(sub_m2, align_mask, 0)

            avars.all_flexible_sub_m2.append(sub_m2)

        return


    def run_file_utilities(self):

        log = self.log
        log.debug('in run_file_utilities')
        mvars = self.mvars
        avars = self.avars

        pdbfilename = mvars.pdbfile

        if(mvars.runname[-1] == '/'):
            lin = len(mvars.runname)
            mvars.runname = mvars.runname[:lin - 1]

        direxist = os.path.exists(mvars.runname)
        if(direxist == 0):
            os.system('mkdir -p ' + mvars.runname + '/')
#
#           global run administration
#
        genpath = mvars.runname + '/complex_monte_carlo'
        avars.genpaths = genpath + '/'
        direxist = os.path.exists(genpath)
        if(direxist == 0):
            os.system('mkdir -p ' + genpath)

        cpfile = os.path.join(mvars.path, pdbfilename)
        cpst = 'cp ' + cpfile + ' ' + avars.genpaths
        os.system(cpst)
#
#           write global run name, pdb, and dcd filenames to .last_sas
#
        fileexist = os.path.exists('.last_sas')
        if(fileexist == 1):
            os.system('mv -f .last_sas .last_sas_bu')
        avars.lastsasfile = open('./.last_sas', 'w')
        avars.lastsasfile.write('run_name\t' + mvars.runname + '\n')
        avars.lastsasfile.write('pdb_name\t' + mvars.pdbfile + '\n')
        avars.lastsasfile.write('dcd_name\t' + mvars.dcdfile + '\n')

        return


    def initialize_segments(self):

        log = self.log
        log.debug('in initialize_segments')
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars

        segname = avars.m1.segname()

        avars.asegs = []
        for tseg in segname:
            if(tseg not in avars.asegs):
                avars.asegs.append(tseg)
        numsegs = len(avars.asegs)
        message = 'found '+str(numsegs)+' segment names'
        log.debug(message)
        pgui(message)

        avars.first_last_resid = []

        avars.all_segment_mask = []
        avars.all_segment_full_mask = []
        avars.all_segment_basis_full_mask = []

        avars.all_segment_mol = []

        tmask = ''

        avars.keyword_basis = False
        if(len(avars.abasis) == 1):
            basis = avars.abasis[0].strip()
            log.debug('basis = %s\n' % (basis))
            if(basis.lower() == 'all' or basis.lower() == 'heavy' or basis.lower() == 'backbone'):
                avars.keyword_basis = True

        for i in xrange(numsegs):
            segmol = sasmol.SasMol(0)
            error, segment_full_mask = avars.m1.get_subset_mask(
                'segname[i] == "' + avars.asegs[i] + '"')
            avars.m1.copy_molecule_using_mask(segmol, segment_full_mask, 0)
            this_resid = segmol.resid()

            avars.first_last_resid.append([this_resid[0], this_resid[-1]])
            avars.all_segment_full_mask.append(segment_full_mask)
            avars.all_segment_mol.append(segmol)

        # this is where abasis is used  --> and this is where it matters!
            if avars.keyword_basis:
                if(basis.lower() == 'all'):
                    log.debug('setting up all atom overlap arrays')
                    segmol.set_average_vdw()
                    npairs = segmol.natoms() * (segmol.natoms() - 1) / 2
                    cutoff_array = numpy.zeros(npairs, numpy.float)
                    pairs.pairs(segmol.atom_vdw(), cutoff_array)
                    keyword_basis_filter = 'segname[i] == "' + \
                        avars.asegs[i] + '" and (not name[i] == "None") '
                elif(basis.lower() == 'backbone'):
                    this_moltype = segmol.moltype()[0]
                    log.debug('this_moltype: %s\n' % (this_moltype))  ### check this

                    if(segmol.moltype()[0] == 'protein'):
                        keyword_basis_filter = 'segname[i] == "' + avars.asegs[
                            i] + '" and (name[i] == "N" or name[i] == "CA" or name[i] == "C") '
                    elif(segmol.moltype()[0] == 'rna' or segmol.moltype()[0] == 'dna'):
                        keyword_basis_filter = 'segname[i] == "' + avars.asegs[
                            i] + '" and (name[i] == "P" or name[i] == "O5\'" or name[i] == "C5\'" or name[i] == "C4\'" or name[i] == "C3\'" or name[i] == "O3\'") '
                    else:
                        keyword_basis_filter = 'segname[i] == "' + \
                            avars.asegs[i] + '" and (not name[i][0] == "H") '
                    # TODO --> add to complex_filter so the following hack is not
                    # needed
                elif(basis.lower() == 'heavy'):
                    keyword_basis_filter = 'segname[i] == "' + \
                        avars.asegs[i] + '" and (not name[i][0] == "H") '

                error, segment_basis_mask = avars.m1.get_subset_mask(
                    keyword_basis_filter)
            else:
                error, segment_basis_mask = avars.m1.get_subset_mask(
                    'segname[i] == "' + avars.asegs[i] + '" and name[i] =="' + avars.abasis[i].strip() + '"')
            avars.all_segment_basis_full_mask.append(segment_basis_mask)

            error, segment_mask = avars.all_segment_mol[i].get_subset_mask(
                'segname[i] == "' + avars.asegs[i] + '"')
            avars.all_segment_mask.append(segment_mask)

        # TODO ... this is probably why flexible segments need to be first!!
        # should just take the NAMES of the flexible segnames to make this
        ###
        # this is also where abasis is used  --> but basis_full_mask is ONLY used for zcut
        # checking: abasis itself is passed to check_overlap in nmer_nrotate
        ###
        # OPTIONS: use moltype()[0] for each asegs[i] to set the basis (CA--> protein, P --> RNA)
        # or better yet, use not hydrogen instead ... as this is ONLY used for z-cut check
        ###
            tmask += '(segname[i] == "' + avars.asegs[i] + \
                '" and (not name[i][0] == "H")) '
            if i != len(avars.flexible_segments) - 1:
                tmask += ' or '

        error, avars.basis_full_mask = avars.m1.get_subset_mask(tmask)

        message = 'first_last_resid = '+str(avars.first_last_resid)
        log.debug(message)

        return

    def initialize_interaction_regions(self):

        log = self.log
        log.debug('in initialize_interactions_regions')
        mvars = self.mvars
        avars = self.avars

        if(len(avars.interpairs) > 0):
            message = 'pair distances < cut == ', avars.cutoff, ' angstroms between segments have been found'
            pgui(message)
            pgui('these distances will be ignorned in overlap check')
            pgui('interpairs = %f\n' % (avars.interpairs))
        else:
            pgui('all distances between segments are greater than cut == %f\n' % (avars.cutoff))
            pgui('normal overlap checking will be used')

        log.debug('npairs: %i\n' % (avars.npairs))
    
        ### initialize interaction regions in each segment ###

        avars.interres = []
        avars.interatom = []
        for i in range(len(avars.interpairs)):
            segnum_1 = avars.interpairs[i][0][0]
            segnum_2 = avars.interpairs[i][0][1]
            for j in range(len(avars.interpairs[i][1])):
                resnum_1 = avars.interpairs[i][1][j][0]
                resnum_2 = avars.interpairs[i][1][j][1]

            # TODO --> need to match basis here as well
            # TODO --> need to match basis here as well
            # TODO --> need to match basis here as well

                basis_segment_1 = '(segname[i] == "' + avars.asegs[segnum_1] + \
                    '" and name[i] =="' + avars.abasis[segnum_1].strip() + '")'
                error, basis_mask_segment_1 = avars.m1.get_subset_mask(basis_segment_1)
                # idx_1 = numpy.where(basis_mask_segment_1==1.0)[0][resnum_1] # a
                # ugly numpy function
                idx_1 = filter(lambda x: basis_mask_segment_1[x] == 1.0, range(
                    len(basis_mask_segment_1)))[resnum_1]
                basis_segment_2 = '(segname[i] == "' + avars.asegs[segnum_2] + \
                    '" and name[i] =="' + avars.abasis[segnum_2].strip() + '")'
                error, basis_mask_segment_2 = avars.m1.get_subset_mask(basis_segment_2)
                # idx_2 = numpy.where(basis_mask_segment_2==1.0)[0][resnum_2] # a
                # ugly numpy function
                idx_2 = filter(lambda x: basis_mask_segment_2[x] == 1.0, range(
                    len(basis_mask_segment_2)))[resnum_2]
                avars.interres.append([resnum_1, resnum_2])
                avars.interatom.append([idx_1, idx_2])

        log.debug('interres: %i\n' %(avars.interres))
        log.debug('interatom: %i\n' %(avars.interatom))

        return


    def set_up_dihedral_arrays(self, txtOutput):

        log = self.log
        log.debug('in set_up_dihedral_arrays')
        mvars = self.mvars
        avars = self.avars

        avars.flexible_dihedral_parameters = []

        avars.all_flexible_basis_mask = []

        for i in xrange(len(avars.flexible_segments)):

            this_segname = avars.flexible_segments[i]
            idx = avars.asegs.index(this_segname)
            local_m1 = avars.all_segment_mol[idx]

            # TODO --> need to deal with specific basis here
            # TODO --> need to deal with specific basis here
            # TODO --> need to deal with specific basis here

            if(avars.keyword_basis):

                if avars.amoltype[i] == 'protein':
                    basis_atom = "CA"

                elif avars.amoltype[i] == 'rna':
                    #basis_atom = "P"
                    basis_atom = "O5\'"

                basis_filter = 'name[i] == "' + basis_atom + \
                    '" and segname[i] == "' + this_segname + '"'

            else:
                basis_filter = 'name[i] == "' + avars.abasis[idx] + \
                    '" and segname[i] == "' + this_segname + '"'

            log.debug('basis filter = %s\n' % (basis_filter))
            error, basis_mask = local_m1.get_subset_mask(basis_filter)
            avars.all_flexible_basis_mask.append(basis_mask)

            basis_m1 = sasmol.SasMol(1)
            error = local_m1.copy_molecule_using_mask(basis_m1, basis_mask, 0)
            basis_resname = basis_m1.resname()
            basis_resid = basis_m1.resid()

            arespsi = []
            aresphi = []
            numranges = avars.anumranges[i]
            reslow = avars.arlow[i]
            numcont = avars.arnum[i]

            if avars.amoltype[i] == 'protein':
                respsi = []
                resphi = []
                dihedral_energy.protein_initialization(
                    respsi, resphi, basis_resid, basis_resname, numranges, reslow, numcont, avars.first_last_resid[idx], txtOutput)
                avars.flexible_dihedral_parameters.append([respsi, resphi])

            elif avars.amoltype[i] == 'rna':
                resalpha = []
                resbeta = []
                resgamma = []
                resdelta = []
                resepsilon = []
                reseta = []
                dihedral_energy.rna_initialization(resalpha, resbeta, resgamma, resdelta, resepsilon, reseta,
                                                basis_resid, basis_resname, numranges, reslow, numcont, avars.first_last_resid[idx], txtOutput)
                avars.flexible_dihedral_parameters.append(
                    [resalpha, resbeta, resgamma, resdelta, resepsilon, reseta])

        return


    def set_up_constraints(self):

        log = self.log
        log.debug('in set_up_constraints')
        mvars = self.mvars
        avars = self.avars

        if(mvars.cflag == 1):

            filter_flag = 0
            error, constraint_basis1_array, constraint_basis2_array, avars.distance_array, avars.type_array = constraints.read_constraints(
                avars.m1, mvars.confile, filter_flag)

            avars.mask_a_array = []
            avars.mask_b_array = []

            for i in xrange(len(avars.distance_array)):
                log.debug(constraint_basis1_array[i])
                log.debug(constraint_basis2_array[i])
                log.debug(avars.distance_array[i])
                log.debug(avars.type_array[i])

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

        return 


    def setup_flexible_residue_mask_arrays(self):

        log = self.log
        log.debug('in setup_flexible_residue_mask_arrays')
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars

        avars.all_flexible_residues = []
        avars.all_flexible_residue_rotation_indices = []
        avars.all_flexible_residue_rotation_mask = []

        for i in xrange(len(avars.flexible_segments)):

            avars.bnumranges = avars.anumranges[i]
            avars.breslow = avars.arlow[i]
            avars.bnumcont = avars.arnum[i]
            log.debug('avars.bnumranges: %i\n' % (avars.bnumranges))
            log.debug('avars.breslow: %s\n' % (avars.breslow))
            log.debug('avars.bnumcont: %s\n' % (avars.bnumcont))
            avars.flexible_residues = self.get_flexible_residues()
            avars.all_flexible_residues.append(avars.flexible_residues)

            segment_filter = 'segname[i] == "' + avars.flexible_segments[i] + '"'
            error, segment_mask = avars.m1.get_subset_mask(segment_filter)
            log.debug('segment_filter = %s\n' % (segment_filter))
            log.debug('error = %s\n' % (error))
            avars.segment_m1 = sasmol.SasMol(98)
            error = avars.m1.copy_molecule_using_mask(avars.segment_m1, segment_mask, 0)

            avars.bmolecule_type = avars.amoltype[i]

            residue_rotation_indices, residue_rotation_mask = self.get_rotation_indices()
            avars.all_flexible_residue_rotation_indices.append(residue_rotation_indices)
            avars.all_flexible_residue_rotation_mask.append(residue_rotation_mask)

        return 

    def get_flexible_residues(self):
        '''
        Method to determine residues that are going to be rotated.
        '''
        log = self.log
        log.debug('in get_flexible_residues')
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui

        flexible_residues = []
        templist = []
        for i in xrange(avars.bnumranges):
            thisreslow = avars.breslow[i]
            thisnumcont = avars.bnumcont[i]
            thisrange = numpy.arange(thisnumcont) + thisreslow
            templist.append(thisrange.tolist())

        flexible_residues = [item for sublist in templist for item in sublist]

        message = 'flexible_residues = '+str(flexible_residues)
        log.debug(message)
        pgui(message)

        return flexible_residues


    def get_rotation_indices(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui

        log.debug('getting rotation indices for molecule')

        residue_rotation_indices = {}
        residue_rotation_mask = {}

        log.debug('avars.bmolecule_type = %s\n' %(avars.bmolecule_type))
        if(avars.bmolecule_type == 'protein'):

            mtype = 0
            mask = avars.segment_m1.get_dihedral_subset_mask(avars.flexible_residues, mtype)

            for i in xrange(len(mask)):
                this_mask = mask[i][:]
                q0 = avars.flexible_residues[i]
                residue_rotation_mask[q0] = this_mask.tolist()
                indices = avars.segment_m1.get_indices_from_mask(this_mask)
                residue_rotation_indices[q0] = indices.tolist()

#            print 'residue_rotation_indices = ',residue_rotation_indices
#            print 'residue_rotation_mask = ',residue_rotation_mask

        elif(avars.bmolecule_type == 'rna'):

            mtype = 1
            mask = avars.segment_m1.get_dihedral_subset_mask(avars.flexible_residues, mtype)

            for i in xrange(len(mask)):
                this_mask = mask[i][:]
                q0 = avars.flexible_residues[i]
                residue_rotation_mask[q0] = this_mask.tolist()
                indices = avars.segment_m1.get_indices_from_mask(this_mask)
                residue_rotation_indices[q0] = indices.tolist()

#            print 'residue_rotation_indices = ',residue_rotation_indices
#            print 'residue_rotation_mask = ',residue_rotation_mask

        elif(avars.bmolecule_type == 'dna'):
            message = 'rotation basis set not defined for molecule type = ' + avars.bmolecule_type
            pgui(message)
            pass

        else:
            message = 'rotation basis set not defined for molecule type = ' + avars.bmolecule_type
            pgui(message)

        log.debug('done getting rotation indices')

        return residue_rotation_indices, residue_rotation_mask


    def evaluate_rg(self):

        avars = self.avars

        maximum_value = max(avars.rg_difference_list)

        if(maximum_value > avars.this_rg_difference):
            index = avars.rg_difference_list.index(maximum_value)
            avars.rg_difference_list[index] = avars.this_rg_difference
            avars.directed_rg_list[index] = avars.rg_value
            avars.accepted_rg_list[index] = avars.accepted

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

        # set up run file I/O
        self.run_file_utilities()

        kb = 1.380658E-23  # J/K
        avars.beta = 1.0 / (mvars.temp * kb)

        avars.m1 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.path + mvars.pdbfile)

        nf1 = avars.m1.number_of_frames()       #not used?       

        avars.dcdoutfile = avars.m1.open_dcd_write(avars.genpaths + mvars.dcdfile)       

        # set up segment arrays
        self.initialize_segments()

        # set up constraints variables
        self.set_up_constraints()

        avars.seglow = mvars.sseglow
        avars.seghigh = mvars.sseghigh

        # set up segment alignment coordinates and com arrays
        self.alignment_initialization()
        
        if(avars.keyword_basis):
            if(mvars.segbasis.lower() == 'all'):
                avars.cutoff = 0.8
            elif(mvars.segbasis.lower() == 'heavy' or mvars.segbasis.lower() == 'backbone'):
                avars.cutoff = 0.8
        else:
            avars.cutoff = 2.0

        message = 'cutoff = '+str(avars.cutoff)
        pgui(message)
    
        check_initial_interactions = False
        if(check_initial_interactions):

            # survey interaction between segments
            avars.interpairs, avars.npairs = nmer_overlap_check.nmer_overlap_check(
                avars.m1, mvars.path, mvars.pdbfile, avars.cutoff, avars.abasis, avars.keyword_basis)

            self.initialize_interaction_regions()

        else:

            avars.interpairs = []
            avars.npairs = 0
            avars.interatom = []
            avars.interres = []

        # set up dihedral parameters for each flexible segment
        self.set_up_dihedral_arrays(txtOutput)

        # set up flexible residue rotation mask arrays
        self.setup_flexible_residue_mask_arrays()

        if(mvars.plotflag == 1):
            avars.graph = Gnuplot.Gnuplot(debug=1)
            avars.graph.clear()
            avars.graph('set title "Rg Results"')
            avars.graph.xlabel('Structure Number')
            avars.graph.ylabel('Rg (Angstrom^2)')

        avars.hrg = 0.0
        avars.lowestrg = 1000.0
        avars.accepted = 0
        avars.over = 0
        avars.badrg = 0
        avars.badz = 0
        avars.badc = 0
        avars.nsteps = 0
        avars.arg = 0.0
        avars.trg = 0.0

#        log.debug(vars(mvars))
#        log.debug(vars(avars))
        
        return


    def dihedralgenerate(self):

        '''
        DIHEDRAL ensembles of structures by varying
        protein dihedral angles. 

        INPUT:  variable descriptions

        runname:        string      project name                          
        path:           string      input file path                 
        dcdfile:        string      name of output dcd file containing accepted structures       
        pdbfile:        string      name of input pdb file containing intial structure
        trials:         integer     number of Monte Carlo move attempts
        goback:         integer     number of failed Monte Carlo attempts before returning to previously accepted structure
        temp:           float       run temperature (K)
        nsegments:      integer     total number of segments   
        npsegments:     integer     number of segments containing flexible regions
        flpsegname:     string      names of segments with flexible regions (separated by commas if more than one)
        segbasis:       string      type of basis for overlap check ("all", "heavy", "backbone" or specific atom name, i.e., "CA")     
        seglow:         integer     low residue for (non-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
        seghigh:        integer     high residue for (no-flexible) structure alignment region (not entered directly; parsed from entered alignment range in GenApp)
        lowrg:          float       low Rg cutoff value if Advanced Input is chosen
        highrg:         float       high Rg cutoff value if Advanced Input is chosen
        zflag:          integer     enable zcutoff flag (0=no, 1=yes)
        zcutoff:        float       zcutoff value (discard structures with any z-axis coordinates less than this value)
        cflag:          integer     enable atomic constraint flag (0=no, 1=yes)
        confile:        string      name of file describing additional constraints to check before accepting a structure
        directedmc:     float       non-zero Rg value to guide Monte Carlo run; 0=no directed Monte Carlo (used if Advanced Input is chosen)

        psegvariables:              flexible segment variables
                        integer     number of flexible regions
                        float_array maximum angle that torsion can sample (in each flexible region)
                        int_array   low residue number for each flexible region
                        int_array   number of contiguous residues per flexible region (not enetered directly; parsed from entered residue range in GenApp)
                        string      molecule type ('protein' or 'rna')                                       


        Input options not implemented:

        psffilepath     string      path to psf file
        psffilename     string      psf file name
        parmfilepath    string      path to CHARMM parameter file
        parmfilename    string      name of CHARMM parameter file
        plotflag        integer     option to plot structure number vs Rg

        OUTPUT:  files stored in "runname"/complex_monte_carlo directory:

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

#        log.debug(vars(mvars))
#        log.debug(vars(avars))


        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        step_parameters = step.Setup()
        coor = avars.m1.coor()
        frame = 0

        # MAIN LOOP

        an = 'psi'
        q0 = 1
        th = 1.0
        seg = avars.asegs[0]
        pairdat = [an, q0, th, seg]

        avars.all_rg_tally = []
        avars.accepted_rg_tally = []
        phi_tally = []
        aphi_tally = []
        psi_tally = []
        apsi_tally = []
        atpsi_tally = []
        atphi_tally = []
        atphipsi_tally = []

        nonbondflag = 0

        if(mvars.seed[0] == 1):
            from numpy.random import RandomState
            seed_object = RandomState(mvars.seed[1])
        else:
            seed_object = -1

        failtally = 0
        acc = 0
        afile = ''
        accfile = []

        avars.minx = []
        avars.miny = []
        avars.minz = []
        avars.maxx = []
        avars.maxy = []
        avars.maxz = []

        if(mvars.directedmc > 0):
            avars.rg_difference_list = []
            avars.directed_rg_list = []
            avars.accepted_rg_list = []
            rg_list_length = 10  # hardwired

        for i in range(mvars.trials):

            if(mvars.seed[0] == 1):
                ran_num = seed_object.rand()
                tflexsegn = int(len(avars.flexible_segments) * ran_num)
                tsegn = avars.asegs.index(avars.flexible_segments[tflexsegn])
            else:
                tflexsegn = int(len(avars.flexible_segments) * random.random())
                tsegn = avars.asegs.index(avars.flexible_segments[tflexsegn])

            tseg = avars.asegs[tsegn]

            molecule_type = avars.amoltype[tflexsegn]

            dtheta = avars.aith[tflexsegn]
            numranges = avars.anumranges[tflexsegn]
            reslow = avars.arlow[tflexsegn]
            numcont = avars.arnum[tflexsegn]

            segment_full_mask = avars.all_segment_full_mask[tsegn]

            error, new_coor = avars.m1.get_coor_using_mask(frame, segment_full_mask)

            segment_mol = avars.all_segment_mol[tsegn]

            segment_mol.setCoor(new_coor)

            '''
		    if(i<10):
			    print 'segment_mol.coor()[0,0,0] = ',segment_mol.coor()[0,0,0]

		    else:
			    sys.exit()
		    '''

            vdi, vdf, indices, this_mask = step_parameters.chooser(new_coor, segment_mol, pairdat, dtheta, numranges, reslow, numcont, avars.flexible_dihedral_parameters[
                                                               tflexsegn], avars.beta, avars.all_flexible_residue_rotation_indices[tflexsegn], avars.all_flexible_residue_rotation_mask[tflexsegn], nonbondflag, avars.first_last_resid[tsegn], molecule_type, seed_object)

            
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
            avars.nsteps += 1
            re = [0, 0, 0, 0.0, 0.0, avars.lowestrg, avars.hrg, 0, 0, []]

            newafile = nmer_nrotate.rotate(coor, avars.m1, q0, th, an, avars.cutoff, mvars.lowrg, mvars.highrg, re, avars.accepted, mvars.zflag, mvars.zcutoff, mvars.cflag,  avars.dcdoutfile, indices, this_mask, avars.all_flexible_basis_mask[tflexsegn], avars.all_flexible_sub_m2[tflexsegn], avars.all_flexible_align_mask[tflexsegn],     avars.all_flexible_coor_sub_m1[tflexsegn], avars.all_flexible_com_sub_m1[tflexsegn], avars.mask_a_array, avars.mask_b_array, avars.distance_array, avars.type_array, avars.first_last_resid[tsegn], molecule_type, avars.all_segment_mask[tsegn], segment_full_mask, avars.all_segment_basis_full_mask, avars.basis_full_mask, avars.all_segment_mol[tsegn], avars.asegs, avars.abasis, avars.interatom, avars.interres)

        
            print '.',
            sys.stdout.flush()

            avars.accepted = avars.accepted + re[0]
            avars.over = avars.over + re[1]
            avars.badrg = avars.badrg + re[2]
            avars.rg_value = re[3]
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

            avars.all_rg_tally.append([i, avars.rg_value])

            if(re[0] == 1):
                avars.accepted_rg_tally.append([i, avars.accepted, avars.rg_value])
                if(mvars.directedmc > 0):
                    if(len(avars.rg_difference_list) <= rg_list_length):
                        avars.this_rg_difference = abs(avars.rg_value - mvars.directedmc)
                        avars.rg_difference_list.append(avars.this_rg_difference)
                        avars.directed_rg_list.append(avars.rg_value)
                        avars.accepted_rg_list.append(avars.accepted)
                    else:
                        avars.this_rg_difference = abs(avars.rg_value - mvars.directedmc)
                        self.evaluate_rg()

            if(re[0] == 0):
                if(failtally == mvars.goback):
                    failtally = 0
                    if(avars.accepted > 0):
                        if(mvars.seed[0] == 1):
                            ran_num = seed_object.rand()
                            dum = int(avars.accepted * ran_num) - 1

                        elif(mvars.directedmc > 0):
                            local_rg_list_length = len(avars.directed_rg_list)
                            ran_num = random.randrange(0, local_rg_list_length)
                            dum = avars.accepted_rg_list[ran_num]

                        else:
                            dum = int(avars.accepted * random.random()) - 1
                        if(dum == -1):
                            pgui('\nreloading coordinates from original starting structure')
#           The 2.0 version of read_pdb with the saspdbrx_topology=True option gives different result than 1.0 version with same option
#           The last few atoms of the last residue aren't included in the coordinates if this option is True
#           But they are read in the 1.0 version of read_pdb.  It is necessary to disable this option to pass the 1.0 tests.
#           THIS ISSUE NEEDS TO BE ADDRESSED BEFORE FURTHER 2.0 DEVELOPMENT! 
#                           avars.m1.read_pdb(mvars.path + mvars.pdbfile, saspdbrx_topology=True)  
                            avars.m1.read_pdb(mvars.path + mvars.pdbfile)
                            coor = avars.m1.coor()
                        else:
                            pgui('\nreloading coordinates from a previously accepted structure')

                            avars.m1.read_single_dcd_step(avars.genpaths + mvars.dcdfile, dum + 1)
                            coor = avars.m1.coor()
                    else:
                        pgui('\n>>>>>reloading coordinates from original starting structure')
#                        avars.m1.read_pdb(mvars.path + mvars.pdbfile, fastread=True,saspdbrx_topology=True)  #to pass 1.0 tests
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

#        log.debug(vars(mvars))
#        log.debug(vars(avars))

        rgplot = open('./' + mvars.runname + '/complex_monte_carlo/' +
                  mvars.dcdfile + '.all_rg_results_data.txt', 'w')
        rgplot.write('# structure number (structure 1 = 1; not 0), Rg (all)\n')
        for ii in range(len(avars.all_rg_tally)):
            rgplot.write('%i\t%f\n' %
                     (avars.all_rg_tally[ii][0] + 1, avars.all_rg_tally[ii][1]))
        rgplot.close()
        rgplot = open('./' + mvars.runname + '/complex_monte_carlo/' +
                  mvars.dcdfile + '.accepted_rg_results_data.txt', 'w')
        rgplot.write(
            '# structure number (structure 1 = 1; not 0), Rg (accepted)\n')
        for ii in range(len(avars.accepted_rg_tally)):
            rgplot.write('%i\t%f\t%i\n' % (avars.accepted_rg_tally[ii][
                     1] - 1, avars.accepted_rg_tally[ii][2], avars.accepted_rg_tally[ii][0] + 1))
        rgplot.close()


        if(avars.accepted > 0):
            pgui("Average accepted rg2 = %lf\n" % (avars.arg / (avars.accepted)))
            pgui(
                "Configurations and statistics saved in %s directory\n" % ('./' + avars.genpaths))
        else:
            pgui("Average accepted rg2 = %lf\n" % (0.0))
            pgui(
                "\n NO ACCEPTED MOVES\n\n Statistics saved in %s directory\n" % (avars.genpaths))

        outfile7 = open(avars.genpaths + mvars.dcdfile + '.stats', 'w')
        outfile7.write('%s\t%f\t%s\t%f\n' %
                   ('lowest Rg = ', avars.lowestrg, 'highest Rg = ', avars.hrg))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('accepted ', avars.accepted,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.accepted / float(avars.nsteps)) * 100.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('overlapped ', avars.over,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.over / float(avars.nsteps)) * 100.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad rg2 ', avars.badrg,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.badrg / float(avars.nsteps)) * 100.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad z-filter ', avars.badz,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.badz / float(avars.nsteps)) * 100.0, ' %'))
        outfile7.write('%s\t%i\t%s\t%i\t%s\t%f%s\n' % ('bad constaints ', avars.badc,
                                                   ' out of ', avars.nsteps, ' moves : ', (avars.badc / float(avars.nsteps)) * 100.0, ' %'))

        if(avars.accepted > 0):
            outfile7.write('%s\t%f\n' %
                       ('average accepted rg2 = ', avars.arg / (avars.accepted)))
        else:
            outfile7.write('%s\t%f\n' % ('average accepted rg2 = ', 0.0))
        outfile7.write('%s\t%f\n' %
                   ('average total rg2 of ensemble = ', avars.trg / (avars.nsteps)))

        pgui("\nDCD data were written to %s\n\n" %
                  ('./' + avars.genpaths + mvars.dcdfile))
        pgui("lowest Rg = %lf\t highest Rg = %lf\n" % (avars.lowestrg, avars.hrg))
        pgui("accepted %d out of %d : %lf percent\n" %
                  (avars.accepted, avars.nsteps, (avars.accepted / float(avars.nsteps)) * 100.0))
        pgui("overlapped %d out of %d moves : %lf percent\n" %
                  (avars.over, avars.nsteps, (float(avars.over) / float(avars.nsteps)) * 100.0))
        pgui("bad rg2 %d out of %d moves : %lf percent\n" %
                  (avars.badrg, avars.nsteps, (float(avars.badrg) / float(avars.nsteps)) * 100.0))
        if(mvars.zflag == 1):
            pgui("bad zcut %d out of %d moves : %lf percent\n\n\n" % (
                avars.badz, avars.nsteps, (float(avars.badz) / float(avars.nsteps)) * 100.0))
        if(mvars.cflag == 1):
            pgui("constraint filter rejected %d out of %d moves : %lf percent\n\n\n" % (
            avars.badc, avars.nsteps, (float(avars.badc) / float(avars.nsteps)) * 100.0))

        if(len(avars.minx) > 0 and len(avars.miny) > 0 and len(avars.minz) > 0 and len(avars.maxx) > 0 and len(avars.maxy) > 0 and len(avars.maxz) > 0):
            min_x = numpy.min(avars.minx)
            min_y = numpy.min(avars.miny)
            min_z = numpy.min(avars.minz)
            max_x = numpy.max(avars.maxx)
            max_y = numpy.max(avars.maxy)
            max_z = numpy.max(avars.maxz)

            pgui("\nminimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                      (min_x, max_x, (max_x - min_x)))

            pgui("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                      (min_y, max_y, (max_y - min_y)))
            pgui("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" %
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
        pgui('COMPLEX DIHEDRAL IS DONE')
        time.sleep(1.5)

        if(mvars.plotflag == 1):
            self.wait('\n')

