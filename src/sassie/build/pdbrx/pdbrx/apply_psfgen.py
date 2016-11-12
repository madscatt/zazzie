# -*- coding: utf-8 -*-
"""
Create a input script to psfgen and execute in order to create a CHARMM
formatted PSF/PDB pair from an input SasMolScan object.

Based in part on the FIX_PDB script coded by Joseph E. Curtis - 05/10/2013

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

"""

import logging
import os
import subprocess

import sasmol.sasmol as sasmol
from sassie.util import sasconfig
import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils


class PsfgenDriver():
    '''
    Class to create psfgen script from an input SasMolScan object and run
    using psfgen installed in SASSIE sasconfig.__binpath__. This should
    produce a CHARMM formated PSF/PDB pair.
    '''

    def __init__(self, mol, segname_info, top_file_path, out_path, out_prefix):
        '''
        @type  mol:            SasMol
        @param mol:            Molecule to be processed
        @type  segname_info:   pdbscan.data_struct.Info
        @param segname_info:   Information about contents of each segment
        @type  top_file_path:  string
        @param top_file_path:  Full path to the topology file to use in PSFGEN
        @type  out_path:       string
        @param out_path:       Path to store output PDB/PSF
        @type  out_prefix:     string
        @param out_prefix:     Prefix to use for all output files
        '''

        self.logger = logging.getLogger(__name__)

        self.mol = mol

        self.segname_info = segname_info

        self.top_file_path = top_file_path

        self.out_path = out_path

        self.out_prefix = out_prefix

        self.script_lines = []

        self.psfgen_exe = os.path.join(sasconfig.__bin_path__, "psfgen")

        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        return

    def run_psfgen(self):
        '''
        Run PSFGEN using the script held in the self.script_lines.
        '''

        script_filename = os.path.join(self.out_path, 'psfgen_pdbrx.in')
        log_filename = os.path.join(self.out_path, 'psfgen_pdbrx.log')

        self.create_psfgen_script()

        outfile = open(script_filename, 'w')

        for line in self.script_lines:
            outfile.write('{0:s}\n'.format(line))

        outfile.close()

        script_file = open(script_filename)

        log_file = open(log_filename, 'w')

        # Run PSFGEN
        failure = subprocess.call(
            [self.psfgen_exe], stdin=script_file, stdout=log_file, stderr=subprocess.STDOUT)

        log_file.close()

        if failure:
            self.logger.critical(
                'ERROR: Unsuccessful PSFGEN run, check log for details')

        return

    def create_psfgen_script(self):
        '''
        Create lines of the PSFGEN script to parameterize the molecule stored
        in self.mol, output to self.script_lines.
        '''

        segnames = self.mol.segnames()
        out_path = self.out_path
        out_prefix = self.out_prefix

        out_lines = self.script_lines

        out_lines.append('')
        out_lines.append('topology {0:s}'.format(self.top_file_path))
        out_lines.append('')

        self.add_standard_aliases_to_script()
        self.add_nucleic_aliases_to_script()

        for segname in segnames:
            self.create_segname_script_section(segname)

        self.add_disulphides_script()

        out_lines.append('')
        out_lines.append('guesscoord')
        out_lines.append('')

        out_location_stub = os.path.join(out_path, out_prefix)

        out_lines.append('writepsf {0:s}_xplor.psf'.format(out_location_stub))
        out_lines.append('writepdb {0:s}.pdb'.format(out_location_stub))

        out_lines.append('writepsf charmm {0:s}.psf'.format(out_location_stub))

        out_lines.append('')

        return

    def create_segname_script_section(self, segname):
        '''
        Prepare segname to be read into PSFGEN. This involves both saving a
        PDB of the segment structure and adding lines to the self.script_lines
        list.

        @type  segname:  string
        @param segname:  Name of the segment to be prepared for PSFGEN
        '''

        mol = self.mol
        info = self.segname_info
        logger = self.logger
        out_lines = self.script_lines

        # Create segment PDB file for PSFGEN run
        out_filename = '{0:s}.pdb'.format(segname)
        output_pdb = os.path.join(self.out_path, out_filename)

        seg_mol = sasmol.SasMol(0)

        basis_filter = 'segname[i] == "{0:s}"'.format(segname)
        error, mask = mol.get_subset_mask(basis_filter)

        for line in error:
            logger.warning(line)

        error = mol.copy_molecule_using_mask(seg_mol, mask, 0)

        seg_mol.write_pdb(output_pdb, 0, 'w')

        for line in error:
            logger.warning(line)

        moltype = seg_mol.moltype()[0]
        nterm_resname = seg_mol.resname()[0]

        # Create psfgen input section for segname
        out_lines.append('segment {0:s} {{'.format(segname))

        out_lines.append('\tpdb {0:s}'.format(output_pdb))

        if moltype == 'protein':
            self.add_nterm_patch(nterm_resname, True)
            out_lines.append('\tlast CTER')
        if moltype in ['nucleic', 'rna', 'dna']:
            out_lines.append('\tfirst 5TER')
            out_lines.append('\tlast 3TER')

        out_lines.append('}')
        out_lines.append('')

        out_lines.append('coordpdb {0:s} {1:s}'.format(output_pdb, segname))

        if moltype == 'dna':

            out_lines.append('')

            res_desc = utils.uniquify_list(
                zip(seg_mol.resid(), seg_mol.resname()))
            for resid, resname in res_desc:
                self.apply_dna_patch(segname, resid, resname)

            out_lines.append('')

        return

    def add_disulphides_script(self):
        '''
        Create lines in the psfgen script to add patches for disulphide bonds
        identified in the input SasMol.

        @return:
        '''

        out_lines = self.script_lines

        chosen = self.mol.segnames()

        disulphides = self.segname_info.disulphides

        for bond in disulphides:

            bound = bond.bound
            residue1 = bound[0]
            residue2 = bound[1]

            if (residue1['subdiv'] in chosen) and (residue2['subdiv'] in chosen):
                out_lines.append('patch DISU {0:s}:{1:d} {2:s}:{3:d}'.format(residue1['subdiv'], residue1['resid'],
                                                                             residue2['subdiv'], residue2['resid']))

        return

    def add_standard_aliases_to_script(self):
        '''
        Many residues and atoms have differences in naming between the PDB
        standard and those in CHARMM. Aliases are used to convert these for
        standard residues. Add lines to do this to the psfgen input script.

        @return:
        '''

        out_lines = self.script_lines

        out_lines.append('')
        out_lines.append('alias residue HIS HSE')
        out_lines.append('alias atom ILE CD1 CD')
        out_lines.append('alias atom SER HG HG1')
        out_lines.append('alias atom CYS HG HG1')
        out_lines.append('alias residue HOH TIP3')

        out_lines.append('alias residue UNK GLY')
        out_lines.append('')

        return

    def add_nucleic_aliases_to_script(self):
        '''
        Many DNA/RNA residues and atoms have differences in naming between the
        PDB standard and those in CHARMM. Aliases are used to convert these for
        standard residues. Add lines to do this to the psfgen input script.

        @return:
        '''

        out_lines = self.script_lines

        out_lines.append('')

        out_lines.append('alias residue G GUA')
        out_lines.append('alias residue C CYT')
        out_lines.append('alias residue A ADE')
        out_lines.append('alias residue T THY')
        out_lines.append('alias residue U URA')

        out_lines.append('alias residue DG GUA')
        out_lines.append('alias residue DC CYT')
        out_lines.append('alias residue DA ADE')
        out_lines.append('alias residue DT THY')

        for base in ['GUA', 'CYT', 'ADE', 'THY', 'URA']:
            out_lines.append("alias atom {0:s} O5* O5'".format(base))
            out_lines.append("alias atom {0:s} C5* C5'".format(base))
            out_lines.append("alias atom {0:s} O4* O4'".format(base))
            out_lines.append("alias atom {0:s} C4* C4'".format(base))
            out_lines.append("alias atom {0:s} C3* C3'".format(base))
            out_lines.append("alias atom {0:s} O3* O3'".format(base))
            out_lines.append("alias atom {0:s} C2* C2'".format(base))
            out_lines.append("alias atom {0:s} O2* O2'".format(base))
            out_lines.append("alias atom {0:s} C1* C1'".format(base))
            out_lines.append("alias atom {0:s} O2* O2'".format(base))
            out_lines.append("alias atom {0:s} OP1 O1P".format(base))
            out_lines.append("alias atom {0:s} OP2 O2P".format(base))

        out_lines.append('alias atom THY C7 C5M')
        out_lines.append('alias atom THY C5A C5M')

        out_lines.append('')

    def add_nterm_patch(self, nterm_resname, flag):
        '''
        CHARMM uses patches to provide changes to n-terminal residues. Create
        appropriate lines for the selected terminal residues in the psfgen
        input script.

        @type nterm_resname :  str
        @param nterm_resname:  Residue name for terminus under consideration
        @type flag :  bool
        @param flag:  Should GLY and PRO be handled as special cases
        @return:
        '''

        out_lines = self.script_lines

        if not flag:

            if nterm_resname == 'GLY' or nterm_resname == 'PRO':

                out_lines.append('\tfirst NONE')

            else:

                out_lines.append('\tfirst NTER')

        elif flag:

            if nterm_resname == 'GLY':

                out_lines.append('\tfirst GLYP')

            elif nterm_resname == 'PRO':

                out_lines.append('\tfirst PROP')

            else:

                out_lines.append('\tfirst NTER')

        return

    def apply_dna_patch(self, segname, resid, resname):
        '''
        CHARMM treats DNA as modified RNA. Need to tell it to apply the correct
        patches for residues in segments of DNA.

        @type segname :  str
        @param segname:  Name of segment residue in question belongs to
        @type resid :  int
        @param resid:  Residue number for relevant residue
        @type resname :  str
        @param resname:  Name of the residue under consideration
        @return:
        '''

        out_lines = self.script_lines

        if resname in ['NUSA', 'NUSG', 'DA', 'DG', 'ADE', 'GUA']:

            # Purine bases
            out_lines.append('patch DEO2 {0:s}:{1:d}'.format(segname, resid))

        else:

            # Pyramidine bases
            out_lines.append('patch DEO1 {0:s}:{1:d}'.format(segname, resid))

        return
