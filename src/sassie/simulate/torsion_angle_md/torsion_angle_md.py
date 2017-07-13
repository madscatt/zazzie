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
import locale
import time
import numpy
import glob
import subprocess
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
from sassie.simulate.torsion_angle_md.write_tamd_input import *
#from write_tamd_input import *
import sasmol.sasmol as sasmol

#       TAMD
#
#       12/05/2004       --      initial coding            :       jc
#       01/02/2011       --      added sasmol support :       jc
#       08/26/2011       --      adapted for mdx :       jc
#       06/16/2012       --      adapted for namd v. 2.9 	:       jc
#       12/08/2013       --      adapted for tamd :       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	TAMD is the module that contains the functions
    	that are used to run a series of tamd dynamics calculations
	on a set of structures in a supplied pdb/dcd file.

        REFERENCES:

        J. Chen et al.
        Journal of Computational Chemistry  26  1565-1578  (2005)

        W. Zhang et al.
        Journal of Molecular Graphics and Modeling  73  179-190  (2017)        
'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'torsion_angle_md'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class torsion_angle_md_input_variables():

    def __init__(self, parent=None):
        pass


class torsion_angle_md():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, psegvariables, txtOutput):

        self.mvars = module_variables()

        self.avars = torsion_angle_md_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization(psegvariables)      #psegvariables are processed here with a call to process_input_variables

        self.tamd()

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
        mvars.infile = variables['infile'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.outfile = variables['outfile'][0]
        mvars.nsteps = variables['nsteps'][0]

        mvars.topfile = variables['topfile'][0]
        mvars.parmfile = variables['parmfile'][0]
        mvars.keepout = variables['keepout'][0]
        mvars.dcdfreq = variables['dcdfreq'][0]
        mvars.charmmexe = variables['charmmexe'][0]
        mvars.temperature = variables['temperature'][0]
        mvars.rgforce = variables['rgforce'][0]
        mvars.rgvalue = variables['rgvalue'][0]

        mvars.dna_segnames = variables['dna_segnames'][0]

        mvars.number_flexible_segments = variables['number_flexible_segments'][0]

        mvars.pretamd_min_steps = variables['pretamd_min_steps'][0]
        mvars.poll_frequency = variables['poll_frequency'][0]

        log.debug(vars(mvars))

        return

#mvars:  runname, infile, pdbfile, outfile, nsteps, topfile, parmfile, keepout, dcdfreq, charmmexe, temperature, rgforce, rgvalue, dna_segnames, number_flexible_segments, pretamd_min_steps, poll_frequency
#avars:  segname_molecules, segname_masks, flexible_segment_variables, temp_pdb_files, path, m1, nf, residue, infiletype, dcdfile


    def process_input_variables(self,psegvariables):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in process_input_variables')

        log.debug('psegvariables: %s' %(psegvariables))               

        all_flexible_segnames = []
        all_snumranges = []
        all_srlow = []
        all_srnum = []
        all_moltype = []
        all_segnames = []
        all_resid = []


        for mol in avars.segname_molecules:
            all_segnames.append(mol.segname()[0])
            all_resid.append(mol.resid()[0])

        offset = 0

        for i in range(len(psegvariables)):
            this_flexible_segname = psegvariables[i][0]
            for j in xrange(len(all_segnames)):
                if this_flexible_segname == all_segnames[j]:
                    offset = - all_resid[j] + 1

            all_flexible_segnames.append(psegvariables[i][0])
            all_snumranges.append(psegvariables[i][1])

            this_srlow = psegvariables[i][2]
            tmp_list = [int(val.strip()) +
                        offset for val in this_srlow.split(',')]
            new_list = ', '.join(str(x) for x in tmp_list)
            all_srlow.append(new_list)
            all_srnum.append(psegvariables[i][3])
            all_moltype.append(psegvariables[i][4])

        all_numranges = []

        for i in range(len(all_snumranges)):
            nr = locale.atoi(all_snumranges[i])
            all_numranges.append(nr)

        all_rlow = []
        all_rnum = []

        for i in range(len(all_srlow)):
            linrlow = string.split(all_srlow[i], ',')
            linrnum = string.split(all_srnum[i], ',')
            rlow = []
            rnum = []
            for k in range(len(linrlow)):
                trlow = locale.atoi(linrlow[k])
                trnum = locale.atoi(linrnum[k])
                rlow.append(trlow)
                rnum.append(trnum)
            all_rlow.append(rlow)
            all_rnum.append(rnum)

        log.debug('all_flexible_segnames = %s' % (all_flexible_segnames))
        log.debug('all_numranges = %s' % str(all_numranges))
        log.debug('all_rlow = %s' % str(all_rlow))
        log.debug('all_rnum = %s' % str(all_rnum))
        log.debug('all_moltype = %s' % (all_moltype))

        avars.flexible_segment_variables = {}

        for i in xrange(len(all_flexible_segnames)):

            avars.flexible_segment_variables[all_flexible_segnames[i]] = [
                all_numranges[i], all_rlow[i], all_rnum[i], all_moltype[i]]

        pgui('flexible_segment_variables: %s' % (avars.flexible_segment_variables))

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return 

#pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return

    def get_segname_molecules(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in get_segname_molecules')

        frame = 0

        avars.m1.initialize_children()

        segnames = avars.m1.segnames()
        avars.segname_masks = avars.m1.segnames_mask()

        number_of_segnames = len(segnames)

        avars.segname_molecules = []

        for i in xrange(number_of_segnames):
            this_mol = sasmol.SasMol(0)
            error = avars.m1.copy_molecule_using_mask(this_mol, avars.segname_masks[i], frame)
            avars.segname_molecules.append(this_mol)

        pgui('segnames = %s' % (segnames))

        return


    def renumber_residues(self, mol):

        '''
        method to renumber resid to start at 1
        '''

        log = self.log
        log.debug('in renumber_residues')

        resid = mol.resid()

        number = []
        resid_array = []
        count = 1
        for i in xrange(len(resid)):
            this_resid = resid[i]
            if(i == 0):
                last_resid = this_resid
            else:
                if(this_resid != last_resid):
                    count += 1
            resid_array.append(count)
            last_resid = this_resid
            number.append(i + 1)

        mol.setResid(resid_array)
        mol.setIndex(number)

        return


    def update_segname_molecules(self, frame):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        log.debug('in update_segname_molecules')

        avars.temp_pdb_files = []

        for i in xrange(len(avars.segname_molecules)):
            error, coor = avars.m1.get_coor_using_mask(frame, avars.segname_masks[i])
            avars.segname_molecules[i].setCoor(coor)

            self.renumber_residues(avars.segname_molecules[i])

            temp_pdb_name = avars.path + 'temp_' + str(i) + '.pdb'
            avars.segname_molecules[i].write_pdb(temp_pdb_name, 0, 'w')
            avars.temp_pdb_files.append(temp_pdb_name)

        return 


    def fix_moltypes(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in fix_moltypes')        

        mvars.dna_segnames = mvars.dna_segnames.split(',')
        pgui('dna_segnames = %s' % (mvars.dna_segnames))

        if len(mvars.dna_segnames) > 0:
            for i in xrange(len(avars.segname_molecules)):
                pgui('segname: %s' % (avars.segname_molecules[i].segname()[0]))
                pgui('natoms: %i' % (avars.segname_molecules[i].natoms()))
                if(avars.segname_molecules[i].segname()[0] in mvars.dna_segnames):
                    moltype = ['dna'] * avars.segname_molecules[i].natoms()
                    avars.segname_molecules[i].setMoltype(moltype)
        else:
            pgui('>> no DNA containing segments in system')

        return


    def initialization(self, psegvariables):
        '''
        method to prepare for tamd
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        avars.path = mvars.runname + '/torsion_angle_md/'
        log.debug('path = %s' % (avars.path))
        mvars.outfile = avars.path + mvars.outfile

        direxist = os.path.exists(avars.path)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + avars.path)
                os.system('cp ' + mvars.pdbfile + ' ' + avars.path)
            except:
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pgui(message)
            if(result != 0):
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pgui(message)

        avars.m1 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.pdbfile)

        avars.residue = avars.m1.resid()
        number_of_residues = avars.residue[-1] - avars.residue[0] + 1
        log.debug('number_of_residues = %i' % (number_of_residues))

        log.debug('infile = %s' % (mvars.infile))
        avars.infiletype = os.path.splitext(mvars.infile)[1][1:]
        log.debug('infiletype = %s' % (avars.infiletype))

        if(avars.infiletype == 'dcd'):
            avars.dcdfile = avars.m1.open_dcd_read(mvars.infile)
            avars.nf = avars.dcdfile[2]
        elif(avars.infiletype == 'pdb'):
            avars.m1.read_pdb(mvars.infile)
            avars.nf = avars.m1.coor()[:, 0, 0].shape[0]
        else:
            st = 'only .dcd or .pdb suffix allowed for infile'
            pgui(st)
            sys.exit()

        log.debug('number of frames = %i' % (avars.nf))

        self.get_segname_molecules()

        self.process_input_variables(psegvariables)

        self.fix_moltypes()

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return


    def tamd(self):
        '''
	    TAMD is the module that contains the functions
    	    that are used to run a series of tamd dynamics calculations
	    on a set of structures in a supplied pdb/dcd file.
      

        INPUT:  variable descriptions:

 	    runname:                        string      run name 
        infile:                         string      input pdb or dcd file name
        pdbfile:                        string      input (reference) pdb file name
        outfile:                        string      output dcd file name
        nsteps:                         integer     number of TAMD steps
        topfile:                        string      path and name of topology file
        parmfile:                       string      path and name of parameter file
        keepout:                        integer     keep output files (0==no, 1==yes)
        dcdfreq:                        integer     save individual dcd frequency
        charmexe                        string      path and name of charmm executable file
        temperature                     float       temperature
        rgforce                         float       rg force
        rgvalue                         float       rg value (Angstroms: ignored if rg force = 0.0)
        dna_segnames                    string      names of dsDNA segments
        number_of_flexible_segments     integer     number of flexible segments
        pretamd_min_steps               string      number of pre-TAMD minimization steps
        poll_frequency                  float       time used in time.sleep command (not input by user)    

        psegvariables:                              flexible segment variables
                                        string      flexible segment name for each flexible segment
                                        integer     number of flexible regions for each flexible segment
                                        int_array   low residue number for each flexible region
                                        int_array   number of contiguous residues per flexible regions
                                        string      molecule type ('protein', 'rna' or 'dna') for each flexible segment    


        OUTPUT:
                
        files stored in ~/run_name/torsion_angle_md directory:
                
        original PDB file with molecular structure data
        tamd_output.pdb (PDB file created by TAMD)
        tamd_output.psf (PSF file created by TAMD)
        TAMD tree file
        TAMD restart file
        DCD file containing TAMD trajectories (size depends on DCD write frequency)
        tamd_dyn_00001.dcd file (trajectory files -- if keep output files is chosen)
        min_00001.out (output files) 
        input file for TAMD run (temp.inp)
        temp_0.pdb (temporary PDB files)
        

        '''

        log = self.log
        log.debug('in tamd')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))

        vers = 'version 0.1 : 12/08/13 : jc'
        #ttxt = time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        dcdlist = []
        coorlist = []
        for i in range(avars.nf):
            message ='\nTAMD dynamics on frame '+ str(i + 1) + ' of ' + str(avars.nf)
            pgui(message)
            pgui('writing TAMD PDB file')

            if(avars.infiletype == 'dcd'):
                avars.m1.read_dcd_step(avars.dcdfile, i)
                self.update_segname_molecules(0)
            else:
                self.update_segname_molecules(i)

            pgui('writing TAMD input file')
            istr = str(i + 1).zfill(5)  # 99999 maximum number of frames

            thisdcd = avars.path + 'tamd_dyn_' + istr + '.dcd'
            dcdlist.append(thisdcd)
            if (mvars.rgvalue == "-1"):                     #not used; negative Rg is not allowed in filter
                my_rg = str(avars.m1.calcrg(0))
            else:
                my_rg = mvars.rgvalue

            write_tamd_input(avars.m1, avars.segname_molecules, 'temp.inp', mvars.nsteps, mvars.dcdfreq, avars.temp_pdb_files, mvars.topfile, mvars.parmfile,
                            thisdcd, mvars.temperature, mvars.rgforce, my_rg, avars.residue, avars.flexible_segment_variables, mvars.pretamd_min_steps)

            message = 'starting TAMD dynamics on frame '+ str(i + 1) + ' of ' + str(avars.nf)
            pgui(message)
            ttime = time.ctime()
            runstring = vers + ' : ' + mvars.outfile + ' run stated at : ' + ttime
            pgui(runstring)

            nst = mvars.charmmexe + ' < temp.inp >& junk.out &'
            log.debug('> nst = %s' % (nst))

            p = subprocess.Popen(nst, shell=True, executable='/bin/bash')
            sts = os.waitpid(p.pid, 0)[1]
            pgui('p.pid = %s' % (str(p.pid)))
            thisjob = str(int(p.pid) + 1)

            run = 1
            esteps = 0
            run_fail = False
            fail_check = None

            while(run == 1):
                time.sleep(mvars.poll_frequency)
                lsst = 'ls junk.out | grep -c "junk.out" '
                lsfile = os.popen(lsst, 'r').readlines()
                stls = string.split(lsfile[0])
                nstls = locale.atoi(stls[0])
                if(nstls > 0):
                    tout2 = os.popen(
                        'tail -15 junk.out | grep "NORMAL TERMINATION"', 'r').readlines()

                    fail_check = os.popen(
                        'tail -15 junk.out | grep "ABNORMAL TERMINATION"', 'r').readlines()

                if(len(tout2) > 0):
                    pgui('finished tamd run')
                    pgui(tout2)
                    run = 0

                if(len(fail_check) > 0):
                    pgui("%s \n" % ('*' * 60))
                    pgui('CHARMM execution has failed: check output file junk.out for details')
                    pgui("%s \n" % ('*' * 60))
                    run_fail = True
                    run = 0

            if not run_fail:
                fraction_done = (float(i + 1) / float(avars.nf))
                progress_string = 'COMPLETED ' + \
                    str(i + 1) + ' of ' + str(avars.nf) + ' : ' + \
                    str(fraction_done * 100.0) + ' % done'
                pgui('%s\n' % (progress_string))
                report_string = 'STATUS\t' + str(fraction_done)
                sys.stdout.flush()
                pgui(report_string)

            else:
                os.system('mv junk.out ' + avars.path + 'min_' + istr + '.out')
                report_string = 'STATUS\t' + str(100.0)
                pgui(report_string)
                pgui("%s \n" % ('*' * 60))
                pgui("RUN FAILED : read contents of min_%s.out file in : %s\n" % (istr, avars.path))
                pgui("Files from the partial run saved to : %s\n" % ('./' + avars.path))
                pgui("%s \n" % ('*' * 60))
                time.sleep(1)
                break

            if(mvars.keepout == 1):
                os.system('mv junk.out ' + avars.path + 'min_' + istr + '.out')
            elif(i == 0):
                os.system('mv junk.out ' + avars.path + 'min_' + istr + '.out')
            else:
                os.system('rm -f junk.out')

            pgui('output trajectory file: %s' % (thisdcd))
            temp_mol = sasmol.SasMol(0)
            temp_mol.read_pdb('tamd_output.pdb', fastread=True)

            header = temp_mol.open_dcd_read(thisdcd)
            temp_mol.close_dcd_read(header[0])

            ndcdfiles = header[2]

            # get the last configuration from the dcd file
            if ndcdfiles > 1:
                if(mvars.keepout != 1):
                    pgui('keeping only the final frame from this trajectory')
                    temp_mol.read_dcd(thisdcd)
                    nframes = temp_mol.number_of_frames()
                    natoms = temp_mol.natoms()
                    coor = numpy.zeros((1, natoms, 3), numpy.float32)
                    coor[0, :, :] = temp_mol.coor()[nframes - 1]
                    temp_mol.setCoor(coor)
                    os.system('rm -f ' + thisdcd)
                    temp_mol.write_dcd(thisdcd)
            if ndcdfiles < 1:
                pgui('number of frames = %i' % (ndcdfiles))
                message = 'Did not save any trajectory files.  Decrease dcd write frequency?'
                message + ' :  stopping here'
                pgui(message)
                sys.exit()

        if not run_fail:
            pgui('\n> finished dynamics all frames\n')
            pgui('writing final data to tamd directory')

            if(avars.infiletype == 'dcd'):
                avars.m1.close_dcd_read(avars.dcdfile[0])

            final_mol = sasmol.SasMol(0)
            final_molw = sasmol.SasMol(1)
            final_mol.read_pdb("tamd_output.pdb", fastread=True)
            final_molw.read_pdb("tamd_output.pdb", fastread=True)

            finaldcdfile = final_molw.open_dcd_write(mvars.outfile)

            for i in range(len(dcdlist)):
                log.debug('i = %i' % (i))
                sys.stdout.flush()
                final_mol.read_dcd(dcdlist[i])

                nframes = final_mol.number_of_frames()
                natoms = final_mol.natoms()
                coor = numpy.zeros((1, natoms, 3), numpy.float32)
                coor[0, :, :] = final_mol.coor()[nframes - 1]
                final_molw.setCoor(coor)
                final_molw.write_dcd_step(finaldcdfile, 0, i + 1)

            final_molw.close_dcd_write(finaldcdfile)

            if(mvars.keepout != 1):
                rmcmd = 'rm -f '
                for i in range(len(dcdlist)):
                    rmcmd = rmcmd + dcdlist[i] + ' '
                log.debug('rmcmd = %s' % (rmcmd))
                os.system(rmcmd)

            pgui("Total number of frames = %d\n\n" % (avars.nf))
            pgui("Structures saved to : %s\n" % ('./' + avars.path))

        return


    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in epilogue')

        for file in glob.glob("cluster_*.str"):
            mvst = 'mv '+ file +' ' + avars.path 
            log.debug(mvst)      
            os.system(mvst)               

        os.system('mv temp.inp ' + avars.path)
        os.system('mv tamd.tree ' + avars.path)
        os.system('mv tamd_loops.rst ' + avars.path)
        os.system('mv tamd_output.psf ' + avars.path)
        os.system('mv tamd_output.pdb ' + avars.path)

        self.run_utils.clean_up(log)

        pgui("%s \n" % ('=' * 60))
        pgui('TAMD IS DONE')
        time.sleep(1.5)

        return        

