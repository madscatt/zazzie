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
import sassie.util.module_utilities as module_utilities

import sassie.build.pdbscan.pdbscan as pdbscan
import sassie.build.pdbscan.pdbscan.report as report
import sassie.build.pdbrx.pdbrx as pdbrx

import os
import sys

if sasconfig.__level__ == "DEBUG": DEBUG = True
	
app = 'pdbrx'

class module_variables():
    def __init__(self,parent = None):
        self.app = app

class PDBRx():

    def __init__(self, parent = None):
        pass
    
    def main(self, input_variables, txtOutput):
        
        self.mvars = module_variables()
        
        self.run_utils = module_utilities.run_utils(app,txtOutput)

        self.run_utils.setup_logging(self)
        
        self.log.debug('in main')
        
        self.unpack_variables(input_variables)
        
        self.run_utils.general_setup(self)
        
        self.run_scan()
        
        return
        
    def unpack_variables(self,variables):
        
        mvars = self.mvars
        
        self.log.debug('in unpack_variables')
        
        mvars.runname = variables['runname'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.topfile = variables['topfile'][0]
        mvars.use_defaults = variables['defaults'][0]

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

        if mvars.use_defaults and not mol.any_charmm_ready_segments():

            pgui('-' * 50)
            pgui('Run terminated: No processable segments found.')
            pgui('-' * 50)
            sys.exit(0)

        mol.copy_biomt_segments()

        if not mvars.use_defaults:

            pgui('Preprocessing starts here')
            preprocessor = pdbrx.preprocessor.PreProcessor(mol=mol)

            for line in report.generate_simulation_prep_report(mol):

                pgui(line)

            preprocessor.user_edit_options()

        pgui('Build scaffold structure')
        scaffold_builder = pdbrx.scaffold_builder.ScaffoldBuilder(mol=mol)

        if mvars.use_defaults:

            scaffold_builder.create_default_scaffold()

        else:

            scaffold_builder.user_system_selection()

        if not os.path.isdir('./tmp_struct'):
            os.mkdir('./tmp_struct')

        pgui('Start structure completion')
        structure_builder = pdbrx.structure_builder.StructureBuilder(scaffold_builder.scaffold_model, './tmp_struct')

        completed_mol = structure_builder.complete_structure()

        # completed_mol.write_pdb('wonder.pdb',0,'w')

        out_prefix = os.path.splitext(os.path.basename(mvars.pdbfile))[0] + '_charmm'

        # top_file_path = os.path.join(sasconfig.__bin_path__,'toppar','top_all27_prot_na.inp')

        segname_info = scaffold_builder.scaffold_model.segname_info

        # psfgen = pdbrx.apply_psfgen.PsfgenDriver(completed_mol, segname_info, top_file_path, self.runpath, out_prefix)
        psfgen = pdbrx.apply_psfgen.PsfgenDriver(completed_mol, segname_info, mvars.topfile, self.runpath, out_prefix)

        psfgen.run_psfgen()

