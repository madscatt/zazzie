# -*- coding: utf-8 -*-
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

import logging

import sassie.util.sasconfig as sasconfig
import sassie.util.communication as communication
import sassie.util.module_utilities as module_utilities

import sassie.build.pdbscan.pdbscan as pdbscan
import sassie.build.pdbscan.pdbscan.report as report
import sassie.build.pdbrx.pdbrx as pdbrx

import os
import sys
import time

if sasconfig.__level__ == "DEBUG": DEBUG = True
	
app = 'pdbrx'

class module_variables():
    def __init__(self,parent = None):
        self.app = app

class PDBRx():

    def __init__(self, parent = None):
        pass
    
    def main(self, input_variables, txtOutput, **kwargs):
      
        self.json_variables = False 
        
        try: 
            if kwargs:
                self.json_variables = {}
                self.json_variables['_uuid'] = kwargs['_uuid']
                self.json_variables['_tcphost'] = kwargs['_tcphost']
                self.json_variables['_tcpport'] = kwargs['_tcpport']
        except:
            pass
             
        self.mvars = module_variables()
        
        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)
        
        self.log.debug('in main')
        
        self.unpack_variables(input_variables)
        
        self.run_utils.general_setup(self)
        
        self.run_scan()
       
        self.epilogue()

        return
        
    def unpack_variables(self,variables):
        
        mvars = self.mvars
        
        self.log.debug('in unpack_variables')
        
        mvars.runname = variables['runname'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.topfile = variables['topfile'][0]
        mvars.use_defaults = variables['defaults'][0]
        mvars.user_interface = variables['user_interface'][0]

        # TODO: Think about topology file
        
        return

    def run_scan(self):

        mvars = self.mvars
        pgui = self.run_utils.print_gui
        log = self.log

        log.debug('in run_scan')

        pgui("\n"+"="*60+" \n")
        pgui("DATA FROM RUN: %s \n\n" %(time.asctime( time.gmtime( time.time() ) ) ))
 
        mol = pdbscan.SasMolScan()
        mol.read_pdb(mvars.pdbfile)
        ###TODO: evaluate whether this should be assigned or not to enable an error
        ###      condition in pdbscan
        #mol.setPdbname(pdbfile)

        pgui('Initiating PDB scan\n')

        fraction_done = 0.01
        pgui('STATUS\t'+str(fraction_done))

        try:
            mol.run_scan()
        except:
            pgui('FATAL ERROR running header scan\n')
            pgui('HALTING RUN\n')
            self.epilogue()
            return

        # if mvars.use_defaults and not mol.any_charmm_ready_segments():
        #
        #     pgui('-' * 50)
        #     pgui('Run terminated: No processable segments found.')
        #     pgui('-' * 50)
        #     sys.exit(0)

        mol.copy_biomt_segments()

        fraction_done = 0.25
        pgui('STATUS\t'+str(fraction_done))

        if not mvars.use_defaults:

            pdbscan_report = report.generate_simulation_prep_report(mol)

            pgui('processing user input for segment(s), sequence(s) and biomt record(s)\n')

            try:
                preprocessor = pdbrx.preprocessor.user_input(self, mol, pdbscan_report)
            except:
                pgui('FATAL ERROR running preprocessor\n')
                pgui('HALTING RUN\n')
                self.epilogue()
                return

        fraction_done = 0.5
        pgui('STATUS\t'+str(fraction_done))
        
        pgui('building scaffold structure\n')

        try:

            scaffold_builder = pdbrx.scaffold_builder.ScaffoldBuilder(mol=mol,
                                                                  default_subs=True,
                                                                  ui=mvars.user_interface)
        except:
            pgui('FATAL ERROR running scaffold builder\n')
            pgui('HALTING RUN\n')
            self.epilogue()
            return

        if mvars.use_defaults:
            
            try:
                scaffold_builder.create_default_scaffold()

            except:
                pgui('FATAL ERROR running scaffold builder: create_default_scaffold\n')
                pgui('HALTING RUN\n')
                self.epilogue()
                return
        
        else:

            try:

                scaffold_builder.user_system_selection(self)

            except:
                pgui('FATAL ERROR running scaffold builder user selection\n')
                pgui('HALTING RUN\n')
                self.epilogue()
                return
 
        tmp_struct_path = self.runpath + os.sep + 'tmp_struct'
        if not os.path.isdir(tmp_struct_path):
            os.mkdir(tmp_struct_path)

        pgui('building structure\n')

        try:
            structure_builder = pdbrx.structure_builder.StructureBuilder(scaffold_builder.scaffold_model,
                                                 tmp_struct_path)

        except:
            pgui('FATAL ERROR running structure builder\n')
            pgui('HALTING RUN\n')
            self.epilogue()
            return

        fraction_done = 0.75
        pgui('STATUS\t'+str(fraction_done))
        
        #try:
        completed_mol = structure_builder.complete_structure()

        #except:
        #    pgui('FATAL ERROR running structure builder: complete_structure\n')
        #    pgui('HALTING RUN\n')
        #    self.epilogue()
        #    return

        out_prefix = os.path.splitext(os.path.basename(mvars.pdbfile))[0] + '_charmm'

        segname_info = scaffold_builder.scaffold_model.segname_info

        psfgen = pdbrx.apply_psfgen.PsfgenDriver(completed_mol, segname_info, mvars.topfile, self.runpath, out_prefix)

        try:

            psfgen.run_psfgen()
       
        except:
            pgui('FATAL ERROR running psfgen\n')
            pgui('HALTING RUN\n')
            self.epilogue()
            return

        pgui('\nfinal structure saved as:: ' + out_prefix + '.pdb\n')
        pgui('\npsf saved as:' + out_prefix + '.psf\n')
        pgui('\nxplor formatted psf saved as: ' + out_prefix + '_xplor.psf\n\n')
 
        fraction_done = 0.9
        pgui('STATUS\t'+str(fraction_done))

    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        fraction_done = 1.0
        pgui('STATUS\t'+str(fraction_done))
        pgui("\n"+"="*60+" \n")

        time.sleep(0.1)
