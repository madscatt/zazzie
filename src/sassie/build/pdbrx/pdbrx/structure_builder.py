# -*- coding: utf-8 -*-
"""

"""

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

import structure_builder_pyrosetta as build_pyrosetta
import sassie.build.pdbscan.pdbscan as pdbscan
import sassie.build.pdbscan.pdbscan.pdbscan_utils as utils

import sasmol.sasmol as sasmol

import logging
import os
import numpy as np
import copy

class StructureBuilder():
    '''
    Drive the building of structures.

    Currently limited to using PyRosetta for protein segments.
    '''

    def __init__(self, mol, tmp_dir, builder='pyrosetta'):
        '''

        @type  mol:      sassie.pdbscan.SasMolScan
        @param mol:      Molecule containing coordinates for basis of complete
                         model and information on regions to be completed
        @type  tmp_dir:  string
        @param tmp_dir:  Directory in which temporary files will be stored
        '''

        self.logger = logging.getLogger(__name__)

        self.input_mol = mol

        self.completed_segnames = []

        self.completed_segment_structures = {}

        self.tmp_dir = tmp_dir

        self.builder = builder

        return

    def complete_structure(self, *args, **kwargs):
        '''
        Wrapper to allow selection of different methods of model completion.

        Currently limited to PyRosetta. If self.builder is not set to
        'pyrosetta' then None is returned

        @rtype : sasmol.SasMol
        @return: Molecule containing models of all segments, including regions
                 modelled from gap_descriptions
        '''

        if self.builder == 'pyrosetta':

            out_mol = self.complete_structure_pyrosetta(self)

        else:

            out_mol = None

        return out_mol


    def complete_structure_pyrosetta(self, *args, **kwargs):
        '''
        Produce complete models for all segments in the input structure,
        self.mol using PyRosetta. Where segments are incomplete and
        descriptions of gaps are provided in (self.segname_info) missing
        sections are modelled in.

        Currently limited to completion of protein segments.

        @rtype : sasmol.SasMol
        @return: Molecule containing models of all segments, including regions
                 modelled from gap_descriptions
        '''

        mol = self.input_mol

        moltype_list = mol.moltypes()

        if 'protein' in moltype_list and len(moltype_list) == 1:

            protein_only = True

        else:

            protein_only = False

        segname_info = mol.segname_info

        completed_segment_structures = self.completed_segment_structures

        self.logger.info('Preparing to complete model')

        for segname in mol.segnames():

            # Obtain list containing a list for each missing region of the form:
            # [pre_anchor, post_anchor, pre_flank, gap, post_flank]
            gap_descriptions = segname_info.seqs_for_completion(segname)

            if gap_descriptions:

                self.logger.info('Dealing with gaps in segment {0:s}'.format(segname))

                scaffold_pdb = os.path.join(self.tmp_dir, 'tmp.pdb')

                # Scaffold PDB created in which the segment be edited is
                # labelled chain E, others F
                self.create_scaffold_pdb(segname,scaffold_pdb)

                # Build missing regions using PyRosetta
                # Note: pass teh fact that we need to edit chain E
                builder = build_pyrosetta.StructureBuilderPyRosetta(scaffold_pdb, gap_descriptions, 'E', protein_only=protein_only)

                builder.complete_structure(self.tmp_dir, 'completed.pdb')

                # Load PDB created by PyRosetta modelling
                completed_pdb = os.path.join(self.tmp_dir, 'completed.pdb')

                completed_model = sasmol.SasMol(0)

                completed_model.read_pdb(completed_pdb)

                # Filter out only the newly completed chain
                # Use to create a correctly labelled molecule of just selected segname
                basis_filter = 'chain[i] == "E" and element[i] != "H"'
                err, mask = completed_model.get_subset_mask(basis_filter)

                for line in err:
                    self.logger.critical(line)

                seg_mol = sasmol.SasMol(0)

                err = completed_model.copy_molecule_using_mask(seg_mol,mask,0)

                for line in err:
                    self.logger.critical(line)

                natoms = seg_mol.natoms()
                seg_list = np.array([segname] * natoms)

                seg_mol.setSegname(seg_list)
                seg_mol.setSegnames([segname])

                chain_label = segname[0]
                chain_list = np.array([chain_label] * natoms)
                seg_mol.setChain(chain_list)
                seg_mol.setChains([chain_label])

                # Store completed segments in a dictionary
                completed_segment_structures[segname] = seg_mol

                # Store segment name in order of reading
                self.completed_segnames.append(segname)

        # Combine all completed chains into a single model
        updated_mol = self.combine_component_models()

        return updated_mol

    def create_scaffold_pdb(self, segname, out_filename, frame = 0):
        '''
        Create a scaffold for missing region completion in segname segment
        of self.mol. The segment to be edited is renamed chain E (for edited)
        and the other segments all renamed F (for fixed) this allows processing
        of files in programs that don't understand segments (such as PyRosetta).

        @type  segname:       string
        @param segname:       Name of segment being processed
        @type  out_filename:  string
        @param out_filename:  Name of the PDB file to be created
        @type  frame:         integer
        @param frame:         Frame number to copy from input mol (self.mol)
        '''

        input_mol = self.input_mol
        natoms = input_mol.natoms()

        mask = np.array([1] * natoms)

        # Copy molecule to new version to allow chain re-labelling without
        # fouling the original data
        new_mol = sasmol.SasMol(0)

        input_mol.copy_molecule_using_mask(new_mol,mask,frame)

        chains = new_mol.chain()
        segnames = new_mol.segname()

        # Re-label chains: - E for editing, F for fixed
        for ndx in range(natoms):

            if segnames[ndx] == segname:
                chains[ndx] = 'E'
            else:
                chains[ndx] = 'F'

        new_mol.setChain(chains)

        new_mol.write_pdb(out_filename,0,'w')

        return

    def combine_component_models(self):
        '''
        Combine all non edited segments in the original model with those
        in self.completed_segment_structures which have been completed to
        create a full model of the system.

        @rtype :  sasmol.SasMol
        @return:  Molecule containing all completed segment models
        '''

        mol = self.input_mol
        completed_segment_structures = self.completed_segment_structures

        updated_mol = None

        # Get list of segnames in order foung in input model
        ordered_segnames = utils.uniquify_list(mol.segname())

        self.logger.info('Combining component models')

        for segname in ordered_segnames:

            if segname in completed_segment_structures:
                # Deal with segments where we needed to fill in missing regions

                seg_mol =  completed_segment_structures[segname]

            else:
                # Deal with segments which were already complete

                basis_filter = 'segname[i] == "{0:s}"'.format(segname)
                err, mask = mol.get_subset_mask(basis_filter)

                for line in err:
                    self.logger.critical(line)

                seg_mol = sasmol.SasMol(0)

                err = mol.copy_molecule_using_mask(seg_mol,mask,0)

                for line in err:
                    self.logger.critical(line)

            if updated_mol:
                tmp = sasmol.SasMol(0)

                err = tmp.merge_two_molecules(updated_mol, seg_mol)

                for line in err:
                    self.logger.critical(line)

                updated_mol = copy.deepcopy(tmp)

            else:

                updated_mol = copy.deepcopy(seg_mol)


            # TODO: Need to fix original_index
            # TODO: Need to copy conect records

        return updated_mol
