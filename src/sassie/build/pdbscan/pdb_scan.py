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

from __future__ import print_function

import logging
import os
import subprocess

import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities

import sassie.build.pdbscan.pdbscan as pdbscan
import sassie.build.pdbscan.pdbscan.report as report
import sassie.build.pdbscan.pdbscan.pdbscan_utils as pdbscan_utils

if sasconfig.__level__ == "DEBUG": DEBUG = True
	
app = 'pdbscan'

class module_variables():
    def __init__(self,parent = None):
        self.app = app

class PDBScan():

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
        mvars.pdbname = os.path.splitext(os.path.basename(mvars.pdbfile))[0]
        
        return

    def generate_full_report(self, report_lines):
        """
        @type  report_lines:
        @param report_lines:
        @return:
        """

        mkd_file = self.mvars.pdbname + '.mkd'
        mkd_filepath = os.path.join(self.runpath, mkd_file)

        out_file = open(mkd_filepath,'w')

        for line in report_lines:
            print(line,file=out_file)

        return mkd_filepath

    def report_to_user(self, report_lines):

        for line in report_lines:
            self.run_utils.print_gui(line)

        return


    def run_scan(self):

        mvars = self.mvars
        pgui = self.run_utils.print_gui
        log = self.log

        log.debug('in run_scan')
    
        pgui('-'*50)
        pgui('PDB Scan')
        pgui('-'*50)

        mol = pdbscan.SasMolScan()
        mol.read_pdb(mvars.pdbfile)

        pgui('Initiating scan')
        mol.run_scan()
        mol.copy_biomt_segments()

        short_report, long_report = report.generate_reports(mol)

        self.report_to_user(short_report)

        mkd_filepath = self.generate_full_report(long_report)
        html_filepath = pdbscan_utils.convert_mkd_html(mkd_filepath)


        
        
