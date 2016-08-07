#       FIX_PDB
#
#       05/10/2013      --      initial coding                  :       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        FIX_PDB is a script to read in a PDB file that has an official header
	and it itemizes chains, patches (disulfide and amino-terminus), writes
	a psfgen input script, runs psfgen on the file to create a cleaned-up
	pdb file and associated psf file.  Hydrogens and some missing atoms are
	added as well.  Atom names are changed to comply with CHARMM naming
	conventions.
	
	This program will process all files in a specified directory sequentially.

	Input values are listed at the bottom of this file.

	TODO:

	Methods to peform similar analyses with files without standard PDB headers
	are needed.  For example, if CHAINS are not listed in a header the program
	does nothing.

	Also, a systematic analysis to identify and replace missing residues and
	atoms not covered in this script needs to be incorporated.  Do not assume
	that the MISSING residue section is correct ... need to validate.


'''

import sasmol.sasmol as sasmol
from sassie.util import sasconfig
import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils

import logging
import os, sys, string, glob, numpy
import subprocess


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

        log_file = open(log_filename,'w')

        # Run PSFGEN
        failure = subprocess.call([self.psfgen_exe], stdin=script_file, stdout=log_file, stderr=subprocess.STDOUT)

        log_file.close()

        if failure:
            self.logger.critical('ERROR: Unsuccessful PSFGEN run, check log for details')

        return

    def create_psfgen_script(self):
        '''
        Create lines of the PSFGEN script to parameterize the molecule stored
        in self.mol, output to self.script_lines.
        '''

        segnames = self.mol.segnames()
        out_path = self.out_path
        out_prefix = self.out_prefix

        st = self.script_lines

        st.append('')
        st.append('topology {0:s}'.format(self.top_file_path))
        st.append('')

        self.add_standard_aliases_to_script()
        self.add_nucleic_aliases_to_script()

        for segname in segnames:
            self.create_segname_script_section(segname)

        self.add_disulphides_script()

        st.append('')
        st.append('guesscoord')
        st.append('')

        out_location_stub = os.path.join(out_path,out_prefix)

        st.append('writepsf {0:s}_xplor.psf'.format(out_location_stub))
        st.append('writepdb {0:s}.pdb'.format(out_location_stub))


        st.append('writepsf charmm {0:s}.psf'.format(out_location_stub))

        st.append('')

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
        st = self.script_lines

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
        st.append('segment {0:s} {{'.format(segname))

        st.append('\tpdb {0:s}'.format(output_pdb))

        if moltype == 'protein':
            self.add_nterm_patch(nterm_resname, True)
            st.append('\tlast CTER')
        if moltype in ['nucleic','rna','dna']:
            st.append('\tfirst 5TER')
            st.append('\tlast 3TER')

        st.append('}')
        st.append('')

        st.append('coordpdb {0:s} {1:s}'.format(output_pdb, segname))

        if moltype == 'dna':

            st.append('')

            res_desc = utils.uniquify_list(zip(seg_mol.resid(), seg_mol.resname()))
            for resid, resname in res_desc:
                self.apply_dna_patch(segname, resid, resname)

            st.append('')

        return

    def add_disulphides_script(self):

        st = self.script_lines

        chosen = self.mol.segnames()

        disulphides = self.segname_info.disulphides

        for bond in disulphides:

            bound = bond.bound
            residue1 = bound[0]
            residue2 = bound[1]

            if (residue1['subdiv'] in chosen) and (residue2['subdiv'] in chosen):
                st.append('patch DISU {0:s}:{1:d} {2:s}:{3:d}'.format(residue1['subdiv'], residue1['resid'],
                                                                      residue2['subdiv'], residue2['resid']))

        return

    def add_standard_aliases_to_script(self):

        st = self.script_lines

        st.append('')
        st.append('alias residue HIS HSE')
        st.append('alias atom ILE CD1 CD')
        st.append('alias atom SER HG HG1')
        st.append('alias atom CYS HG HG1')
        st.append('alias residue HOH TIP3')

        st.append('alias residue UNK GLY')
        st.append('')

        return

    def add_nucleic_aliases_to_script(self):

        st = self.script_lines

        st.append('')

        st.append('alias residue G GUA')
        st.append('alias residue C CYT')
        st.append('alias residue A ADE')
        st.append('alias residue T THY')
        st.append('alias residue U URA')

        st.append('alias residue DG GUA')
        st.append('alias residue DC CYT')
        st.append('alias residue DA ADE')
        st.append('alias residue DT THY')

        for base in ['GUA','CYT', 'ADE', 'THY', 'URA']:
            st.append("alias atom {0:s} O5* O5'".format(base))
            st.append("alias atom {0:s} C5* C5'".format(base))
            st.append("alias atom {0:s} O4* O4'".format(base))
            st.append("alias atom {0:s} C4* C4'".format(base))
            st.append("alias atom {0:s} C3* C3'".format(base))
            st.append("alias atom {0:s} O3* O3'".format(base))
            st.append("alias atom {0:s} C2* C2'".format(base))
            st.append("alias atom {0:s} O2* O2'".format(base))
            st.append("alias atom {0:s} C1* C1'".format(base))
            st.append("alias atom {0:s} O2* O2'".format(base))
            st.append("alias atom {0:s} OP1 O1P".format(base))
            st.append("alias atom {0:s} OP2 O2P".format(base))

        st.append('alias atom THY C7 C5M')
        st.append('alias atom THY C5A C5M')

        st.append('')

    def add_nterm_patch(self, nterm_resname, flag):

        st = self.script_lines

        if not flag:

            if (nterm_resname == 'GLY' or nterm_resname == 'PRO'):

                st.append('\tfirst NONE')

            else:

                st.append('\tfirst NTER')

        elif flag:

            if nterm_resname == 'GLY':

                st.append('\tfirst GLYP')

            elif nterm_resname == 'PRO':

                st.append('\tfirst PROP')

            else:

                st.append('\tfirst NTER')

        return

    def apply_dna_patch(self, segname, resid, resname):

        st = self.script_lines

        if resname in ['NUSA','NUSG','DA','DG','ADE','GUA']:
            st.append('patch DEO2 {0:s}:{1:d}'.format(segname,resid))
        else:
            st.append('patch DEO1 {0:s}:{1:d}'.format(segname,resid))

        return
