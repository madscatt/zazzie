import os
import sys
import locale
import string
import time
import subprocess
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
import prody as prody
import numpy

'''
	PRODY is the module that contains functions to setup and run normal mode
	analysis, using ProDy, based on the structure in a supplied
	pdb file.
'''
if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'prody'

prody_exe = os.path.join(os.path.dirname(sys.executable), 'prody')

class module_variables():

    def __init__(self, parent=None):
        self.app = app


class prody_input_variables():

    def __init__(self, parent=None):
        pass


class prody_anm():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = prody_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.prody_driver()

        self.epilogue()

        return

    def unpack_variables(self, variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.number_modes = variables['number_modes'][0]
        mvars.number_conformations_samp = variables[
            'number_conformations_samp'][0]
        mvars.number_steps_traverse = variables['number_steps_traverse'][0]
        mvars.rmsd_conformations_samp = variables['rmsd_conformations_samp'][0]
        mvars.rmsd_traverse = variables['rmsd_traverse'][0]
        mvars.advanced_usage = variables['advanced_usage'][0]
        mvars.advanced_usage_cmd = variables['advanced_usage_cmd'][0]

        log.debug(vars(mvars))

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

    def initialization(self):

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in initialization')

        avars.path = os.path.join(mvars.runname, 'prody')
        direxist = os.path.exists(avars.path)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + avars.path)
            except:
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pqui(message)
            if(result != 0):
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pgui(message)

        return


    def prody_driver(self):
        '''
        PRODY DRIVER is the function to read in variables from GUI input and
        used to run a prody normal mode calculation using the anisotropic network model
        (ANM) on a structure provide in a pdb file.

        INPUT:  variable descriptions:

        pdbfile:          input pdb file (reference)

        OUTPUT:
                        model_anm_extended_bb.nmd
                        model_traverse.dcd
                        model_samp.dcd
                        model_samp.pdb
                        model_anm_sqflucts.txt
                        model_anm_kirchhoff.txt
                        model_anm_hessian.txt
                        model_anm_cross-correlations.hm
                        model_anm_cross-correlations.txt
                        model_anm_covariance.txt
                        model_anm_beta.txt
                        model_anm_evalues.txt
                        model_anm_evectors.txt
                        model_anm_extended_all.nmd
                        model_anm.nmd

        txtOutput:        TK handler for output to GUI textbox

        files stored in ~/runname/prody directory:
        outfile:          output filename

        '''
        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in prody_anm')

        # start gui output
        pgui("\n%s \n" % ('=' * 60))
        pgui("DATA FROM RUN: %s \n\n" % time.asctime( time.gmtime( time.time() ) ))

        #prody or evol commands work from the command line (with prody installed, of course), but not from within this program.
        #if we do implement this in the future, we should bypass filter checks except for checking the command so that the filter isn't checking for the
        #existence of mvars.pdbfile, etc.
        if mvars.advanced_usage == 1:
            run_cmd = prody_exe + mvars.advanced_usage_cmd
            os.system(run_cmd)
            run_cmd = 'mv *.nmd *.txt *.hm prody'
            os.system(run_cmd)
            exit()

        # display progress
        fraction_done = (0 + 1) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        prody.confProDy(verbosity='none')  #option to set silent verbosity
        model = mvars.pdbfile[0:len(mvars.pdbfile) - 4]
        run_cmd = prody_exe + ' anm ' + \
            mvars.pdbfile + ' -t all -n ' + str(mvars.number_modes) + ' -a'
#        print 'run_cmd: ', run_cmd
        log.info('starting prody_exe %s' % run_cmd)
        prody_run = subprocess.Popen(run_cmd,shell=True,executable='/bin/bash')
        prody_run.wait() 
        #prody.confProDy(verbosity='none')  #option to set silent verbosity
        file_anm = model + '_anm_extended_all.nmd'

        # display progress
        fraction_done = (1 + 1) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        # parse nmd file with resuts extended to all atoms
        log.info('starting prody.parseNMD %s' % file_anm)
        mod, ag = prody.parseNMD(file_anm, type=None)
        allatoms = ag.copy()
        # set up to randomly sample number_conformations_samp modes
        log.info('starting prody.sampleModes')
        ensemble = prody.sampleModes(mod[:mvars.number_modes],
                                     ag,
                                     n_confs=mvars.number_conformations_samp,
                                     rmsd=mvars.rmsd_conformations_samp)
        ensemble
        log.info('starting prody ensemble and writing pdb/dcd files')
        allatoms.addCoordset(ensemble)
        prody.writePDB('model_samp.pdb', allatoms)
        prody.writeDCD('model_samp.dcd', allatoms)
        trajectory_names = []
        
        # display progress
        fraction_done = (1 + 2) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        log.info('starting prody traverse')
        for i in xrange(0, mvars.number_modes):
            #print i
            # setup to tranverse slowest mode
            traverse = prody.traverseMode(
                mod[i],
                allatoms,
                n_steps=mvars.number_steps_traverse,
                rmsd=mvars.rmsd_traverse)
            traverse
            prody.writeDCD('traverse.dcd', traverse)
            this_dcd = str(os.path.join(avars.path, 'traverse_' + str(i) + '.dcd'))
            cmd = 'mv traverse.dcd ' + this_dcd
            os.system(cmd)
            trajectory_names.append(this_dcd)

        # display progress
        fraction_done = (1 + 7) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        m1 = sasmol.SasMol(0)
        m2 = sasmol.SasMol(0)
        m1.read_pdb(mvars.pdbfile)
        m2.read_pdb(mvars.pdbfile,fastread=True)

        mvars.dcdfile = mvars.runname + '.dcd'
        log.info('opening new dcd file to store trajectory: %s' %
                 os.path.join(self.runpath, mvars.dcdfile))

        outfile_name = str(os.path.join(avars.path, mvars.dcdfile))
        dcdoutfile = m2.open_dcd_write(outfile_name)
        count = 0
        coor = numpy.zeros((1,m2.natoms(),3),numpy.float32)
        for this_trajectory_name in trajectory_names:

            dcdfile = m1.open_dcd_read(this_trajectory_name)
            number_of_frames = dcdfile[2]

            for j in xrange(number_of_frames):
                m1.read_dcd_step(dcdfile,j)
                coor[0,:,:] = m1.coor()[0]
                m2.setCoor(coor)
                m2.write_dcd_step(dcdoutfile,0, count + 1)
                count += 1

        m2.close_dcd_write(dcdoutfile)

        log.info('moving files to runname / prody')

        file_anm = model + '_anm.nmd'
        mod, ag = prody.parseNMD(file_anm, type=None)
        mod1 = prody.parsePDB(mvars.pdbfile)
        calphas = mod1.select('calpha')
        bb_anm, bb_atoms = prody.extendModel(mod, calphas, mod1.select(
            'backbone'))  # extend model to backbone atoms
        prody.writeNMD('model_anm_extended_bb.nmd', bb_anm, bb_atoms)

        cmd = 'mv model_samp.pdb ' + avars.path + os.sep + os.path.basename(model) + '_samp.pdb'
        os.system(cmd)

        cmd = 'mv model_samp.dcd ' + avars.path + os.sep + os.path.basename(model) + '_samp.dcd'
        os.system(cmd)

        cmd = 'mv model_anm_extended_bb.nmd ' + \
            model + '_anm_extended_bb.nmd'
        os.system(cmd)

        cmd = 'mv *.hm *.nmd *.txt ' + avars.path + os.sep
        os.system(cmd)
        
        # display progress
        fraction_done = (1 + 9) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)

        return

    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        pgui('Configurations and statistics saved in %s directory\n\n' % (self.runpath+os.path.sep))

        pgui('DCD data were written to %s\n\n' % os.path.join(self.runpath, self.mvars.dcdfile))

        pgui("\n" + "=" * 60 + " \n")
        time.sleep(0.1)

        return
