# -*- coding: utf-8 -*-
"""
Build the scaffold PDB to be used as the basis of the final parameterized
model. Also need to store data about regions to be completed.

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

from __future__ import division  # You don't need this in Python3

import logging
import sassie.build.pdbscan.pdbscan as pdbscan
import sassie.build.pdbscan.pdbscan.report as report
from . import segment_choice
from . import biomt_choice
from . import altloc_choice
import sasmol.sasmol as sasmol
import numpy as np


class ScaffoldBuilder():
    '''
    Build structure to form scaffold for building. Scaffold contains only
    regions selected by user and apply appropriate BIOMTs.
    '''

    def __init__(self, *args, **kwargs):

        self.logger = logging.getLogger(__name__)

        self.prep_report = []
        self.biomt_report = []

        self.selected_segnames = []
        self.selected_altlocs = {}
        self.selected_biomt = []

        self.segname_map = {}

        self.scaffold_model = None

        self.selected_mol = None

        if 'mol' in kwargs:

            self.mol = kwargs['mol']

            self.prep_report = report.generate_simulation_prep_report(self.mol)

            if self.mol.segname_info.biomt:
                self.biomt_report = report.create_biomt_summary(
                    self.mol.segname_info.biomt)

        else:
            # TODO: Make a better decision on what to do here
            self.mol = None

        if 'ui' in kwargs:
            self.ui = kwargs['ui']
        else:
            self.ui = 'terminal'

        return

    def check_segment_status(self, segnames_to_check):
        """
        Determine if the segments listed in self.selected_segnames
        are passed as ready for CHARMM parameterization

        @rtype :  boolean
        @return:  Are all segnames composed of CHARMM compatible residues only
        """

        sim_ready_checks = self.mol.sim_ready

        accepted_segnames = []

        for segname in segnames_to_check:

            checks = sim_ready_checks[segname]
            if checks['charmm']:
                accepted_segnames.append(segname)

        return accepted_segnames == segnames_to_check

    def process_non_ff(self):
        """
        Stub which will eventually process non-standard CHARMM residues
        using CGENFF or CHARM-GUI glycan reader for example

        """

        return

    def create_excess_mask(self, model_no=1):
        """
        Create a SasMol style mask to select all excess atoms

        @type model_no :  int
        @param model_no:  Model number
        @rtype : numpy.array
        @return: Integer mask for all atoms (1 = excess, 0 = not)
        """

        excess = self.selected_mol.segname_info.excess_atoms[model_no]
        mol = self.mol

        mask = np.zeroes(mol.natoms())

        for segname, res_info in excess.iteritems():

            for resid, excess_info in res_info.iteritems():

                for name in excess_info['atoms']:
                    sel_txt = ('segname[i] == "{0:s}" and resid[i] == {1:d}'
                               ' and name[i] == "{2:s}"'.format(
                        segname, resid, name))

                    err, tmp_mask = mol.get_subset_mask(sel_txt)

                    mask = np.logical_or(tmp_mask, mask)

        return mask.astype(int)

    def select_specified_regions(self):
        """
        Create self.selected_mol as a sasmol.SasMol object containing only the
        segments and residue conformations selected by the user to be made in the final model

        """

        selected_segnames = ['"{0:s}"'.format(x) for x in self.selected_segnames]
        selected_altlocs = self.selected_altlocs

        natoms = self.mol.natoms()

        err, alt_loc_mask = self.mol.get_subset_mask('loc[i] == " "')

        for line in err:
            self.logger.warning(line)

        for segname, resids in selected_altlocs.iteritems():

            for resid, loc in resids.iteritems():
                sel_txt = 'segname[i] == "{0:s}" and resid[i] == {1:d} and loc[i] == "{2:s}"'.format(
                    segname, resid, loc)
                err, tmp_mask = self.mol.get_subset_mask(sel_txt)

                alt_loc_mask = np.logical_or(
                    alt_loc_mask, tmp_mask).astype(int)

        sel_txt = 'segname[i] in [{0:s}]'.format(','.join(selected_segnames))
        err, segname_mask = self.mol.get_subset_mask(sel_txt)

        for line in err:
            self.logger.warning(line)

        # Combine segname and altloc masks
        # Note: sasmol uses 0/1 based masks not boolean
        mask = np.logical_and(segname_mask, alt_loc_mask).astype(int)

        self.selected_mol = pdbscan.SasMolScan()

        frame = self.mol.model_no - 1
        self.mol.copy_molecule_using_mask(self.selected_mol, mask, frame)

        self.selected_mol.segname_info = self.mol.segname_info

        self.fix_residue_numbering()

        # Purge non-selected segments from segname_info
        for segname in self.selected_mol.segname_info.subdivs:

            if segname not in selected_segnames:
                self.selected_mol.segname_info.purge_subdiv(segname)

        return

    def fix_residue_numbering(self):
        """
        Renumber selected segments to start at residue 1.

        @return:
        """

        segname_info = self.selected_mol.segname_info
        selected_segnames = self.selected_segnames

        for segname in selected_segnames:
            segname_info.subdiv_renumber_from_one(segname)

        return

    def apply_biomt_to_model(self):
        """
        Apply the BIOMTs in self.selected_biomt to self.selected_mol. Each
        model generated by application of a transformation is created as a
        new SasMol object and returned in a list. No new model is created
        for the identity transform.

        @rtype :  List
        @return:  List of sasmol.SasMol objects containing BIOMT transformed structures

        @todo:  Consider moving some of this code into sasmol -
                create a copy and transform method perhaps
        """

        mol = self.mol
        selected_mol = self.selected_mol
        biomts = mol.segname_info.biomt
        selected_biomt = self.selected_biomt
        selected_segnames = self.selected_segnames

        frame = mol.model_no - 1

        new_mols = []

        for biomt_no in selected_biomt:

            biomt = biomts[biomt_no]

            segs_to_transform = list(
                set(biomt['subdivs']).intersection(selected_segnames))

            if segs_to_transform:

                seg_list = ['"{0:s}"'.format(x) for x in segs_to_transform]
                sel_txt = 'segname[i] in [{0:s}]'.format(','.join(seg_list))

                err, mask = selected_mol.get_subset_mask(sel_txt)

                for line in err:
                    self.logger.warning(line)

                for i in range(len(biomt['trans'])):

                    u = biomt['rot'][i]
                    m = biomt['trans'][i]

                    if not ((u == np.identity(3)).all() and (m == np.array([0.0, 0.0, 0.0])).all()):
                        new_mols.append(sasmol.SasMol(0))

                        selected_mol.copy_molecule_using_mask(
                            new_mols[-1], mask, frame)

                        new_mols[-1].apply_biomt(frame, sel_txt, u, m)

        return new_mols

    def update_models_segname_original_index(self, models):
        """
        Give all SasMol models passed in segment names which do not conflict
        with those selected by the user from the original model.

        @type  models:  List
        @param models:  List of SasMol objects
        @rtype :        Dictionary
        @return:        Dictionary mapping new segment names to those in
                        the original structure
        """

        mol = self.mol

        last_original_index = mol.original_index()[-1]
        existing_segnames = list(self.selected_segnames)

        segname_map = {}

        for model in models:

            tmp_ndx = model.original_index() + last_original_index
            model.setOriginal_index(tmp_ndx)

            last_original_index = model.original_index()[-1]

            model_segnames = model.segnames()

            tmp_segnames = model.segname()

            for segname in model_segnames:
                new_segname = mol.next_segname_generator(
                    existing=existing_segnames)

                tmp_segnames = [new_segname if x ==
                                               segname else x for x in tmp_segnames]

                existing_segnames.append(new_segname)

                segname_map[new_segname] = segname

            model.setSegname(tmp_segnames)

        return segname_map

    def merge_models(self, models):
        """
        Create a merged model from self.selected_mol and the input list of
        SasMol models. The merged model is stored as self.scaffold_model as
        a SasMolScan object with segname_info copied from the original
        input model with an addition of self.segname_map.

        @type  models:  List
        @param models:  List of SasMol objects
        """

        segname_info = self.mol.segname_info

        mol = self.selected_mol

        combined_model = None

        if models:

            for model in models:

                if combined_model is None:
                    
                    combined_model = model

                else:

                    tmp_model = sasmol.SasMol(0)

                    tmp_model.merge_two_molecules(combined_model, model)

                    combined_model = tmp_model

            self.scaffold_model = pdbscan.SasMolScan()
            self.scaffold_model.merge_two_molecules(mol, combined_model)
            self.scaffold_model.segname_info = segname_info
            self.scaffold_model.segname_info.subdiv_map = self.segname_map

        else:

            self.scaffold_model = mol
            self.scaffold_model.segname_info = segname_info

        return

    def terminal_segment_selection(self):
        '''
        Get user to select segments from the terminal

        @return:
        '''

        mol = self.mol
        segname_list = mol.segnames()
        prep_report = self.prep_report

        choice_made = False

        while not choice_made:

            if len(segname_list) > 1:

                selected_segnames = segment_choice.select_segnames(
                    segname_list, prep_report)

            else:

                selected_segnames = [segname_list[0]]

            if selected_segnames:

                if not self.check_segment_status(selected_segnames):

                    # TODO: This is just a stub - needs filling out
                    self.process_non_ff()

                    if self.check_segment_status(selected_segnames):
                        choice_made = True
                else:
                    choice_made = True

            else:

                print("No segnames selected")

            self.selected_segnames = selected_segnames

        return

    def get_segnames_altlocs_selected(self):
        '''
        Return list of segnames fro segments containing altlocs from those
        selected for final model.

        @return:
        '''

        altlocs = self.mol.segname_info.altloc
        alt_segnames = []

        for segname in self.selected_segnames:

            if segname in altlocs and len(altlocs[segname]) > 0:
                alt_segnames.append(segname)

        return alt_segnames

    def _defaut_altloc(self):
        '''
        Chose the first altloc for all residues in which one is present

        @return:
        '''

        altlocs = self.mol.segname_info.altloc

        for segname in self.selected_segnames:

            self.selected_altlocs[segname] = {}

            if segname in altlocs:

                seg_altlocs = altlocs[segname]

                for resid, locs in seg_altlocs.iteritems():
                    # Select first non-blank loc label
                    chosable_locs = [x for x in locs if x != ' ']
                    chosen_loc = chosable_locs[0]

                    self.selected_altlocs[segname][resid] = chosen_loc

    def create_scaffold_model(self):
        '''
        Create the scaffold model - extract the selected regions and apply
        BIOMT transforms. Application of the BIOMT transforms creates
        moved duplicates of the selected segments - these are initially
        created as separate models and then merged.

        @return:
        '''

        self.select_specified_regions()

        if self.selected_biomt:

            models = self.apply_biomt_to_model()

            self.segname_map = self.update_models_segname_original_index(
                models)

            self.merge_models(models)

        else:

            self.scaffold_model = self.selected_mol

        # self.scaffold_model.write_pdb('combined_test.pdb', 0, 'w')

        return

    def create_default_scaffold(self):
        '''
        Create a scaffold without user input. Choose all segments which
        contain only standard CHARMM residues, apply author suggested BIOMT
        and the first AltLoc available.

        @return:
        '''

        sim_ready_checks = self.mol.sim_ready
        biomt = self.mol.segname_info.biomt

        # Select all valid segments
        selected_segnames = []

        for segname in sim_ready_checks.keys():

            if sim_ready_checks[segname]['charmm']:
                selected_segnames.append(segname)

        self.selected_segnames = selected_segnames

        # Select first AltLoc where necessary
        self._defaut_altloc()

        # Apply all BIOMT records suggested by PDB author
        if biomt:

            auth_biomt = []

            for biomt_ndx, info in biomt.iteritems():
                if info['auth_bio_unit']:
                    auth_biomt.append(biomt_ndx)

            self.selected_biomt = auth_biomt

        self.create_scaffold_model()

        return

    def user_system_selection(self):
        '''
        Get user to select the regions of the protein to be included in the
        final model. Also select which AltLoc to use if multiple available.

        @return:
        '''

        mol = self.mol
        altlocs = self.mol.segname_info.altloc
        biomol_report = self.biomt_report

        if self.ui == 'terminal':

            self.terminal_segment_selection()

            alt_segnames = self.get_segnames_altlocs_selected()

            if len(alt_segnames) > 0:

                print("Multiple conformations were found for residues "
                      "in the following segments:")
                print(str(alt_segnames))
                print("Are you happy using the first conformation (altloc)"
                      " for all residues? (answer [y]es/[n]o)")

                choice = ''

                while choice.lower() not in ['y', 'n', 'yes', 'no']:
                    choice = raw_input().lower()

                if choice.lower() in ['n', 'no']:

                    for segname in self.selected_segnames:

                        self.selected_altlocs[segname] = {}

                        if segname in altlocs:

                            seg_altlocs = altlocs[segname]

                            for resid, locs in seg_altlocs.iteritems():

                                chosable_locs = [x for x in locs if x != ' ']

                                chosen_altloc = altloc_choice.select_altloc(
                                    segname, resid, chosable_locs)

                                if chosen_altloc:
                                    self.selected_altlocs[segname][
                                        resid] = chosen_altloc[0]
                                else:
                                    self.selected_altlocs[segname][
                                        resid] = list(locs).remove(' ')[0]

                else:

                    self._defaut_altloc()

            if self.mol.segname_info.biomt:

                biomt_list = self.mol.segname_info.biomt.keys()

                if len(biomt_list) > 1:

                    sel_biomt = biomt_choice.select_biomt(
                        biomt_list, biomol_report)
                    self.selected_biomt = [int(x) for x in sel_biomt]

                else:

                    self.selected_biomt = [biomt_list[0]]

        else:
            pass

        self.create_scaffold_model()

        return
