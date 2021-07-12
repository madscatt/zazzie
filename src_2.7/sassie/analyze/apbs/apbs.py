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
import string
import time
import subprocess
from sassie.analyze.apbs.write_apbs_input import *
#from write_apbs_input import *
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig

#       APBS
#
#       12/05/2004       --      initial coding                    :    jc
#       01/02/2011       --      added sasmol support         :    jc
#       08/26/2011       --      adapted for mdx             :    jc
#       06/16/2012       --      adapted for namd v. 2.9         :    jc
#       09/10/2012       --      adapted for apbs            :    jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
    APBS is the module that contains the functions
    that are used to run a series of electrostatic calculations
    on a set of structures in a supplied pdb/dcd file.

    REFERENCES:

    Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of
    nanosystems: application to microtubules and the ribosome.
    Proc. Natl. Acad. Sci. USA 98, 10037-10041 2001.

    M. Holst and F. Saied, Multigrid solution of the Poisson-Boltzmann equation.
    J. Comput. Chem. 14, 105-113, 1993.

    M. Holst and F. Saied, Numerical solution of the nonlinear Poisson-Boltzmann
    equation: Developing more robust and efficient methods.
    J. Comput. Chem. 16, 337-364, 1995.

    M. Holst, Adaptive numerical treatment of elliptic systems on manifolds.
    Advances in Computational Mathematics 15, 139-191, 2001.

    R. Bank and M. Holst, A New Paradigm for Parallel Adaptive Meshing Algorithms.
    SIAM Review 45, 291-323, 2003.

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'apbs'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class apbs_input_variables():

    def __init__(self, parent=None):
        pass


class apbs():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = apbs_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.apbs()

        self.epilogue()

        return

    def unpack_variables(self, variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.infile = variables['infile'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.temperature = variables['temperature'][0]
        mvars.ph = variables['ph'][0]
        mvars.ion_charge = variables['ion_charge'][0]
        mvars.ion_conc = variables['ion_conc'][0]
        mvars.ion_radius = variables['ion_radius'][0]
        mvars.manual_flag = variables['manual_flag'][0]  # not yet implemented
        mvars.manual_file = variables['manual_file'][0]  # not yet implemented

        log.debug(vars(mvars))

        return

# mvars: runname, infile, pdbfile, temperature, ph, ion_charge, ion_radius, manual_flag, manual_file
# avars: path, nf, m1, m2, maximum_dimensions, infiletype, dcdfile


#   pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return

    def rename_his(self, m):

        natoms = m.natoms()
        resname = m.resname()
        new_resname = []
        for i in xrange(natoms):
            this_resname = resname[i]
            if(this_resname == 'HSE' or this_resname == 'HSD' or this_resname == 'HSP'):
                new_resname.append('HIS')
            else:
                new_resname.append(this_resname)

        m.setResname(new_resname)

        return

    def initialization(self):
        '''
        method to prepare for apbs
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        avars.path = mvars.runname + '/apbs/'
        print 'path = ', avars.path
        print 'infile = ', mvars.infile

        vers = 'version 0.1 : 09/10/12 : jc'
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

        avars.m1 = sasmol.SasMol(0)
        avars.m2 = sasmol.SasMol(0)
        avars.m1.read_pdb(mvars.pdbfile)
        avars.m2.read_pdb(mvars.pdbfile, fastread=True)

        self.rename_his(avars.m1)
        self.rename_his(avars.m2)

        try:
            if(mvars.infile[-3:] == 'dcd'):
                avars.infiletype = 'dcd'
            elif(mvars.infile[-3:] == 'pdb'):
                avars.infiletype = 'pdb'

        except:
            message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message += ' :  stopping here'
            pgui(message)
            sys.exit()

        pgui('infiletype = %s' % (avars.infiletype))

        if(avars.infiletype == 'dcd'):
            min_max = avars.m2.calc_minmax_all_steps(mvars.infile)
            avars.dcdfile = avars.m1.open_dcd_read(mvars.infile)
            avars.nf = avars.dcdfile[2]
        else:
            avars.m1.read_pdb(mvars.infile)
            avars.nf = avars.m1.coor()[:, 0, 0].shape[0]
            min_max = avars.m2.calc_minmax_all_steps(mvars.infile, pdb='pdb')

        pgui('number of frames = %i' % (avars.nf))

        pgui('min_max = %s' % (str(min_max)))

        avars.maximum_dimensions=[min_max[1][0] - min_max[0][0],
                          min_max[1][1] - min_max[0][1], min_max[1][2] - min_max[0][2]]
        pgui('maximum_dimensions = %s' % (str(avars.maximum_dimensions)))

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return

    def apbs(self):
        '''
        APBS is the function to run a series of apbs calculations
        on a set of structures in a supplied pdb/dcd file.

        INPUT:  variable descriptions:

            runname:        project name
            pdbfile:        reference pdb name
            infile:         input trajectory filename (pdb or dcd)
            ph:             pH value
            temperature:    temperature value (K)
            ion_conc:       concentration of solute ion (M)
            ion_radius:     radius of solute ion (angstroms)


        OUTPUT:

            files stored in ~/runname/apbs directory:

            apbs_00001_io.mc            MC-shell I/O capture file
            apbs_00001_pdb2pqr.dat      Output captured from PDB2PQR program
            apbs_00001_pot.dx.mc        File containing electrostatic potential information
            apbs_00001.in               Inputs for APBS calculation
            apbs_00001.out              Output from APBS calculation
            apbs_00001.pdb              PDB file with APBS outputs
            apbs_00001.pqr              PDB file generated by PDB2PQR program
            apbs_00001.propka           Output from protein PKA predictor program
            .
            .
            .
            (depending on number of frames)

        '''

        log=self.log
        pgui=self.run_utils.print_gui
        mvars=self.mvars
        avars=self.avars
        log.debug('in apbs')


        ttxt=time.asctime(time.gmtime(time.time()))
        st=''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        final_energy=[]
        coorlist=[]
        for i in range(avars.nf):
            pgui('apbs calculation for frame %s of %i' % (str(i+1),avars.nf))
            pgui('writing temporary PDB file')

            if(avars.infiletype == 'dcd'):
                avars.m1.read_dcd_step(avars.dcdfile, i)
                avars.m1.write_pdb(avars.path + 'junk.pdb', 0, 'w')
            else:
                avars.m1.write_pdb(avars.path + 'junk.pdb', i, 'w')

            pgui('writing temporary APBS input file')
            if(i < 9):
                istr='0000' + str(i + 1)
            elif(i < 99):
                istr='000' + str(i + 1)
            elif(i < 999):
                istr='00' + str(i + 1)
            elif(i < 9999):
                istr='0' + str(i + 1)
            elif(i < 99999):
                istr=str(i + 1)
            else:
                print 'wow, man!'
                istr=str(i + 1)

            thisdcd=avars.path + 'min_' + istr + '.dcd'

            if(mvars.manual_flag == 0):
                inputfilename='junk.in'
                write_apbs_input(avars.maximum_dimensions, mvars.temperature,
                             inputfilename, mvars.ion_charge, mvars.ion_conc, mvars.ion_radius)
            else:
                inputfilename=mvars.manual_file

            pgui('starting apbs calculation (nfiles = %i)' % (avars.nf))
            ttime=time.ctime()
            ncpu=1
            bin_path=sasconfig.__bin_path__
            if(ncpu == 1):
                pgui('starting pdb2pqr calculation number: %s' % (istr))
#                runstring = 'python ' + bin_path + '/pdb2pqr.py --ff=charmm --with-ph=' + \
#                    str(mvars.ph) + ' -v ' + avars.path + \
#                        'junk.pdb junk.pqr >& pdb2pqr.out'
#                print 'runstring = ', runstring
#                run_pdb2pqr = 'python ' + bin_path + '/pdb2pqr.py --ff=charmm --with-ph=' + \
#                    str(mvars.ph) + ' -v ' + avars.path + \
#                        'junk.pdb junk.pqr >& pdb2pqr.out'
                run_pdb2pqr = bin_path + '/pdb2pqr.py --ff=charmm --with-ph=' + \
                    str(mvars.ph) + ' -v ' + avars.path + \
                        'junk.pdb junk.pqr >& pdb2pqr.out'
                        
                os.system(run_pdb2pqr)

                pgui('starting apbs calculation number: %s' % (istr))
                nst=bin_path + '/apbs junk.in >& junk.out &'

                p=subprocess.Popen(nst, shell=True, executable='/bin/bash')
                sts=os.waitpid(p.pid, 0)[1]
                pgui('p.pid = %s' % (str(p.pid)))
                thisjob=str(int(p.pid) + 1)

            run=1
            esteps=0
            while(run == 1):
                # time.sleep(5)
                lsst='ls junk.out | grep -c "junk.out" '
                lsfile=os.popen(lsst, 'r').readlines()
                stls=string.split(lsfile[0])
                nstls=locale.atoi(stls[0])
                if(nstls > 0):
                    tout2=os.popen(
                        'tail -15 junk.out | grep "Thanks for using"', 'r').readlines()

                if(len(tout2) > 0):
                    pgui('finished apbs calculation')
                    run=0

            fraction_done=(float(i + 1) / float(avars.nf))
            progress_string='COMPLETED ' + \
                str(i + 1) + ' of ' + str(avars.nf) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            pgui('%s\n' % progress_string)
            pgui('%s\n' % progress_string)
            report_string='STATUS\t' + str(fraction_done)
            pgui(report_string)

            pgui('finished run')

            mvst = 'mv io.mc ' + avars.path + 'apbs_' + istr + '_io.mc'
            print mvst
            os.system(mvst)
            mvst = 'mv pot.dx ' + avars.path + 'apbs_' + istr + '_pot.dx.mc'
            print mvst
            os.system(mvst)
            mvst = 'mv pdb2pqr.out ' + avars.path + 'apbs_' + istr + '_pdb2pqr.dat'
            print mvst
            os.system(mvst)
            mvst = 'mv ' + avars.path + 'junk.pdb ' + avars.path + 'apbs_' + istr + '.pdb'
            print mvst
            os.system(mvst)
            mvst = 'mv junk.out ' + avars.path + 'apbs_' + istr + '.out'
            print mvst
            os.system(mvst)
            mvst = 'mv junk.pqr ' + avars.path + 'apbs_' + istr + '.pqr'
            print mvst
            os.system(mvst)
            mvst = 'mv junk.propka ' + avars.path + 'apbs_' + istr + '.propka'
            print mvst
            os.system(mvst)
            mvst = 'mv junk.in ' + avars.path + 'apbs_' + istr + '.in'
            print mvst
            os.system(mvst)


        if(avars.infiletype == 'dcd'):
            avars.m1.close_dcd_read(avars.dcdfile[0])

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log=self.log
        pgui=self.run_utils.print_gui
        mvars=self.mvars
        avars=self.avars

        log.debug('in epilogue')

        pgui("Total number of frames = %d\n\n" % (avars.nf))
        pgui("output energies saved to : %s\n" % ('./' + avars.path))

        self.run_utils.clean_up(log)

        st=''.join(['=' for x in xrange(60)])
        pgui("\n%s \n\n" % (st))
        pgui('APBS IS DONE')

        time.sleep(0.5)

        return
