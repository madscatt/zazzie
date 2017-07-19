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
import subprocess
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
#import sasconfig as sasconfig
import sassie.util.folder_management as folder_management
from sassie.simulate.energy_minimization.write_namd_input import *
from sassie.simulate.energy_minimization.prepend_namd_input import *
#from write_namd_input import *
#from prepend_namd_input import *

#       NAMD_MINIMIZE
#
#       12/05/2004       --      initial coding                	:       jc
#       01/02/2011       --      added sasmol support 		:       jc
#       08/26/2011       --      adapted for mdx 		:       jc
#       06/16/2012       --      adapted for namd v. 2.9 	:       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        ENERGY MINIMIZATION is the module that contains the functions
        that are used to run a series of energy minimization calculations
        on a set of structures in a supplied pdb/dcd file.

        This module is called from Structure Minization in the main 
        GUI through the graphical_minimize.py script.

        REFERENCE:

        J. C. Phillips et al.
        Journal of Computational Chemistry  26  1781-1802  (2005)

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'energy_minimization'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class energy_minimization_input_variables():

    def __init__(self, parent=None):
        pass


class energy_minimization():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = energy_minimization_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()      

        self.minimize()

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
        mvars.resparmfile = variables['resparmfile'][0]
        mvars.psffile = variables['psffile'][0]
        mvars.ncpu = variables['ncpu'][0]
        mvars.keepout = variables['keepout'][0]
        mvars.dcdfreq = variables['dcdfreq'][0]
        mvars.infiletype = variables['infiletype'][0]
        mvars.md = variables['md'][0]
        mvars.mdsteps = variables['mdsteps'][0]
        mvars.dielect = variables['dielect'][0]
        mvars.temperature = variables['temperature'][0]
        mvars.use_external_input_file = variables['use_external_input_file'][0]
        mvars.external_input_file = variables['external_input_file'][0]
        mvars.velocity_restart_file = variables['velocity_restart_file'][0]
        mvars.extended_system_restart_file = variables['extended_system_restart_file'][0]

        log.debug(vars(mvars))

        return 


#    pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(message)
#        st = ''.join(['=' for x in xrange(60)])
#        txtOutput.put("\n%s \n" % (st))
#        time.sleep(1.5)
#
#        return


    def initialization(self):
        '''
        method to prepare for minimize
        '''
        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        avars.path = mvars.runname + '/energy_minimization/'
        pgui('\npath = '+avars.path)
        pgui('infile = '+mvars.infile)

        avars.vers = 'version 0.7 : 06/16/12 : jc'
        direxist = os.path.exists(avars.path)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + avars.path)
            except:
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pgui(message)
                sys.exit()
            if(result != 0):
                message = 'can not create project directory: ' + avars.path
                message += '\nstopping here\n'
                pgui(message)
                sys.exit()

        pgui('psffile = '+ mvars.psffile)
        pgui('pdbfile = '+ mvars.pdbfile)
        os.system('cp ' + mvars.psffile + ' ' + avars.path)
        os.system('cp ' + mvars.pdbfile + ' ' + avars.path)

        avars.m1 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.pdbfile)

        try:
            if(mvars.infile[-3:] == 'dcd'):
                mvars.infiletype = 'dcd'
            elif(mvars.infile[-3:] == 'pdb'):
                mvars.infiletype = 'pdb'
        except:
            message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message += ' :  stopping here'
            pgui(message)
            sys.exit()

        pgui('\ninfiletype = '+mvars.infiletype)
        
        if(mvars.infiletype == 'dcd'):
            avars.dcdfile = avars.m1.open_dcd_read(mvars.infile)
            avars.nf = avars.dcdfile[2]
        else:
            avars.m1.read_pdb(mvars.infile)
            avars.nf = avars.m1.coor()[:, 0, 0].shape[0]

        log.debug('\nnumber of frames = '+str(avars.nf))

        return


    def minimize(self):
        '''
        MINIMIZE is the module that contains the functions
        that are used to run a series of energy minimization calculations
	    on a set of structures in a supplied pdb/dcd file.

        This module is called from Structure Minization in the main 
        GUI through the graphical_minimize.py script.

        REFERENCE:

        J. C. Phillips et al.
        Journal of Computational Chemistry  26  1781-1802  (2005)

        INPUT:  variable descriptions:
              
        reference PDB file
        input PDB or DCD file
        number of minimization steps
        path and name of topology file
        input PSF file
        output (DCD) file
        number of CPUs
        keep output files (0==no, 1==yes)
        MD flag (0=min, 1=min+md, 2=min+md+min)
        number of MD steps (if md = 1 or 2)
        solvent dielectric constant
        temperature (K)
        frequency to save individual DCD files
        flag to use external input file (True or False)
        external input file name
        velocity restart file name
        extended system restart file name


        Advanced input:

        name of user-supplied CHARMM parameter file

        
        OUTPUT:
                
        files stored in ~/run_name/energy_minimization directory:
                
        original PDB file with molecular structure data
        original PSF file
        DCD file containing energy minimized (and md, if applicable) trajectories (size depends on DCD write frequency)
        PDB file containing final energy minimized (and md, if applicable) structure
        output file (if keep output files is chosen)
        input file for energy minimization 
        '''


        log = self.log
        log.debug('in minimize')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))



#        ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        avars.dcdlist = []
        coorlist = []
        for i in range(avars.nf):
            pgui('\nminimizing frame '+str(i + 1)+ ' of '+str(avars.nf))
            log.debug('writing temporary PDB file')

            if(mvars.infiletype == 'dcd'):
                avars.m1.read_dcd_step(avars.dcdfile, i)
                avars.m1.write_pdb(avars.path + 'junk.pdb', 0, 'w')
            else:
                avars.m1.write_pdb(avars.path + 'junk.pdb', i, 'w')

            log.debug('writing temporary NAMD input file')
            if(i < 9):
                istr = '0000' + str(i + 1)
            elif(i < 99):
                istr = '000' + str(i + 1)
            elif(i < 999):
                istr = '00' + str(i + 1)
            elif(i < 9999):
                istr = '0' + str(i + 1)
            elif(i < 99999):
                istr = str(i + 1)
            else:
                print 'wow, man!'
                istr = str(i + 1)

            thisdcd = avars.path + 'min_' + istr + '.dcd'
            avars.dcdlist.append(thisdcd)

            if not mvars.use_external_input_file:
                write_namd_input('temp.inp', str(mvars.nsteps), str(
                    mvars.dcdfreq), avars.path + 'junk.pdb', mvars.psffile, thisdcd, mvars.resparmfile, mvars.md, mvars.mdsteps,
                    mvars.dielect,mvars.temperature)
            else:
                prepend_namd_input('temp.inp', avars.path + 'junk.pdb', mvars.psffile, thisdcd, mvars.resparmfile,
                   mvars.external_input_file, mvars.velocity_restart_file, mvars.extended_system_restart_file)

            print_string = 'starting minimization (nfiles = '+str(avars.nf)+')'
            pgui(print_string)
            ttime = time.ctime()
            runstring = avars.vers + ' : ' + mvars.outfile + ' run started at : ' + ttime
            pgui(runstring)
            bin_path = sasconfig.__bin_path__
            arch = sasconfig.__arch__
            if(mvars.ncpu == 1):
                #nst='/usr/local/bin/namd/namd2 temp.inp >& junk.out &'
                nst = bin_path + '/namd/namd2 temp.inp >& junk.out &'
                log.debug('\n> nst = ' + nst+'\n')
            else:
                pgui('\n> using ' + str(mvars.ncpu) + ' cpus\n')
                #nst='/usr/local/bin/namd/charmrun +p'+str(mvars.ncpu)+' /usr/local/bin/namd/namd2 temp.inp >& junk.out &'
                if arch == "cluster":
                    nst = bin_path + '/namd/charmrun +p' + \
                        str(mvars.ncpu) + ' ++remote-shell ssh ' + bin_path + \
                        '/namd/namd2 temp.inp >& junk.out &'
                    log.debug('\n> nst = ' + nst+'\n')
                else:
                    nst = bin_path + '/namd/charmrun +p' + \
                        str(mvars.ncpu) + ' ' + bin_path + \
                        '/namd/namd2 temp.inp >& junk.out &'
                    log.debug('\n> nst = ' + nst+'\n')

            p = subprocess.Popen(nst, shell=True, executable='/bin/bash')
            sts = os.waitpid(p.pid, 0)[1]
            print_string = 'process id = '+str(p.pid)
            pgui('%s\n' %print_string) 
            thisjob = str(int(p.pid) + 1)

            run = 1
            esteps = 0
            while(run == 1):
                time.sleep(10)
                lsst = 'ls junk.out | grep -c "junk.out" '
                lsfile = os.popen(lsst, 'r').readlines()
                stls = string.split(lsfile[0])
                nstls = locale.atoi(stls[0])
                if(nstls > 0):
                    if(mvars.ncpu == 1 or arch != "cluster"):
                        tout2 = os.popen(
                            'tail -15 junk.out | grep "Program finished"', 'r').readlines()
                    else:
                        tout2 = os.popen(
                            'tail -15 junk.out | grep "WallClock:"', 'r').readlines()

                    fail_check = os.popen(
                        'tail -200 junk.out | grep "Abort"', 'r').readlines()
                    if len(fail_check) == 0:
                        fail_check = os.popen(
                            'tail -200 junk.out | grep "FATAL"', 'r').readlines()
                    if len(fail_check) == 0:
                        fail_check = os.popen(
                            'tail -200 junk.out | grep "ERROR: Exiting prematurely"', 'r').readlines()
                    if(mvars.ncpu > 1 and arch == "cluster"):
                        tout2 = os.popen(
                            'tail -15 junk.out | grep "WallClock:"', 'r').readlines()

                if(len(tout2) > 0):
                    pgui('finished minimization')
                    run = 0
                if(len(fail_check) > 0):
                    message = '\n>>>> investigate error in junk.out file <<<<\n\n'
                    message += "".join(fail_check) + '\n'
                    pgui(message)
                    avars.dcdlist=[]    #to prevent error in epilogue
                    avars.nf = 0
                    return

            time.sleep(10)
            fraction_done = (float(i + 1) / float(avars.nf))
            progress_string = 'COMPLETED ' + \
                str(i + 1) + ' of ' + str(avars.nf) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            pgui('%s\n' % progress_string)
            pgui('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            pgui(report_string)

            if(mvars.keepout == 1):
                os.system('mv junk.out ' + avars.path + 'min_' + istr + '.out')
            else:
                os.system('rm -f junk.out')
            os.system('rm -f junk.coor junk.xs* junk.vel ')
            try:
                os.system('rm -f ' + avars.path + 'junk.pdb ')
            except:
                pgui('\n> could not move junk.pdb')

            os.system('rm -f ' + avars.path + '*.BAK')

            print_string = 'thisdcd = '+thisdcd
            log.debug('%s\n' % print_string)
            temp_mol = sasmol.SasMol(0)
            temp_mol.read_pdb(mvars.pdbfile, fastread=True)

            header = temp_mol.open_dcd_read(thisdcd)
            temp_mol.close_dcd_read(header[0])

            ndcdfiles = header[2]

            # get the last configuration from the dcd file
            if ndcdfiles > 1:
                temp_mol.read_dcd(thisdcd)
                nframes = temp_mol.number_of_frames()
                natoms = temp_mol.natoms()
                coor = numpy.zeros((1, natoms, 3), numpy.float32)
                coor[0, :, :] = temp_mol.coor()[nframes - 1]
                os.system('rm -f ' + thisdcd)
                temp_mol.setCoor(coor)
                temp_mol.write_dcd(thisdcd)
                # temp_mol.write_pdb(thisdcd+'.pdb',0,'w')
            if ndcdfiles < 1:
                print_string = 'ndcdfiles = '+str(ndcdfiles)
                pgui('%s\n' % print_string)
                message = 'Did not save any dcd files.  Decrease dcd write frequency?'
                message + ' :  stopping here'
                pgui(message)
                sys.exit()

        pgui('\n> finished minimizing all frames\n')

        if(mvars.infiletype == 'dcd'):
            avars.m1.close_dcd_read(avars.dcdfile[0])

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

        final_mol = sasmol.SasMol(0)
        final_molw = sasmol.SasMol(1)
        final_mol.read_pdb(mvars.pdbfile, fastread=True)
        final_molw.read_pdb(mvars.pdbfile, fastread=True)

        finaldcdfile = final_molw.open_dcd_write(avars.path + mvars.outfile)

        for i in range(len(avars.dcdlist)):
#            print 'i = ', i
            sys.stdout.flush()
            final_mol.read_dcd(avars.dcdlist[i])
            final_molw.setCoor(final_mol.coor())
            final_molw.write_dcd_step(finaldcdfile, 0, i + 1)

        final_molw.write_pdb(avars.path + mvars.outfile + '.pdb', 0, 'w')

        final_molw.close_dcd_write(finaldcdfile)

        rmcmd = 'rm -f '
        for i in range(len(avars.dcdlist)):
            rmcmd = rmcmd + avars.dcdlist[i] + ' '
#        print 'rmcmd = ', rmcmd
        os.system(rmcmd)
        os.system('mv temp.inp ' + avars.path)

        self.run_utils.clean_up(log)

        sys.stdout.flush()
        st = ''.join(['=' for x in xrange(60)])
        pgui("Total number of frames = %d\n\n" % (avars.nf))
        pgui("Minimized structures saved to : %s\n" % ('./' + avars.path))
        pgui("\n%s \n" % (st))
        time.sleep(0.5)
        pgui('NAMD MINIMIZATION IS DONE')

        return


