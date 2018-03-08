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
        
        self.run_utils = module_utilities.run_utils(app,txtOutput)

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
        mvars.gui = variables['gui'][0]

        # TODO: Think about topology file
        
        return

    def run_scan(self):

        mvars = self.mvars
        pgui = self.run_utils.print_gui
        log = self.log

        log.debug('in run_scan')
    
        pgui('-'*50)
        pgui('PDB Rx')
        pgui('-'*50)

        mol = pdbscan.SasMolScan()
        mol.read_pdb(mvars.pdbfile)

        pgui('Initiating scan')
        mol.run_scan()

        # if mvars.use_defaults and not mol.any_charmm_ready_segments():
        #
        #     pgui('-' * 50)
        #     pgui('Run terminated: No processable segments found.')
        #     pgui('-' * 50)
        #     sys.exit(0)

        mol.copy_biomt_segments()

        if not mvars.use_defaults:

            pgui('Preprocessing')
            #preprocessor = pdbrx.preprocessor.PreProcessor(mol=mol,default_subs=True,ui=mvars.gui)
            #preprocessor = pdbrx.preprocessor.PreProcessor(mol=mol,default_subs=True,ui=mvars.gui,logger=log, json=self.json_variables)
            preprocessor = pdbrx.preprocessor.user_input(mvars, log, pgui, mol, self.json_variables)

            pdbscan_report = report.generate_simulation_prep_report(mol)

            pgui('Printing pdbscan_report')
            for line in pdbscan_report:

                pgui(line)

            preprocessor.user_edit_options(pdbscan_report)

        pgui('Build scaffold structure')

        scaffold_builder = pdbrx.scaffold_builder.ScaffoldBuilder(mol=mol,
                                                                  default_subs=True)

        if mvars.use_defaults:

            scaffold_builder.create_default_scaffold()

        else:

            scaffold_builder.user_system_selection()

        
        tmp_struct_path = self.runpath + os.sep + 'tmp_struct'
        if not os.path.isdir(tmp_struct_path):
            os.mkdir(tmp_struct_path)

        pgui('Start structure completion')
        structure_builder = pdbrx.structure_builder.StructureBuilder(scaffold_builder.scaffold_model,
                                                 tmp_struct_path)

        completed_mol = structure_builder.complete_structure()

        out_prefix = os.path.splitext(os.path.basename(mvars.pdbfile))[0] + '_charmm'

        segname_info = scaffold_builder.scaffold_model.segname_info

        psfgen = pdbrx.apply_psfgen.PsfgenDriver(completed_mol, segname_info, mvars.topfile, self.runpath, out_prefix)

        psfgen.run_psfgen()


    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        #pgui('\n%s IS DONE\n' % app)
        #pgui("\n"+"="*60+" \n")
        time.sleep(0.1)
