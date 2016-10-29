# -*- coding: utf-8 -*-
"""
Scanner

Module containing classes to scan the coordinates read from a PDB
(SasMolScan) and extract data by chains (for reconciliation with PDB headers)
and segnames to be evaluate the system as a candidate for simulation in the
CHARMM forcefield.
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

import os
import numpy as np
import logging
import copy

from itertools import groupby
import itertools

import sasmol.sasmol as sasmol
from sassie.util import sasconfig

from . import data_struct
from . import pdbscan_utils as utils
from . import header_reader
from . import reconcile


class SasMolScan(sasmol.SasMol):

    def __init__(self, top_file=None, top_path=None, model_no=1):
        """
        Initialise the sasmol_scan object. If a topology is provided this is
        used instead of the default CHARMM 27 as provided by NAMD/VMD. Creates
        dictionaries of atoms which are expected for different residue types.

        @type  top_file:  string
        @param top_file:  Path to a CHARMM format topology file
        """

        super(SasMolScan, self).__init__(0)

        self.logger = logging.getLogger(__name__)

        if not top_file:
            top_path = os.path.join(sasconfig.__bin_path__, "toppar")
            top_file = "top_all27_prot_na.inp"
            self.top_file = os.path.join(top_path, top_file)
        elif not top_path:
            top_path = '.'

        self.topology = sasmol.saspdbrx.Topology()
        self.topology.read_charmm_topology(
            topology_file_path=top_path,
            topology_file_name=top_file)

        # model number to be used in analysis
        self.model_no = 1

        self.chain_info = data_struct.Info(scan_type='chain')
        self.segname_info = data_struct.Info(scan_type='segname')

        self.pdbname = ''

        self.header_reconciled_chains = {}
        self.header_reconciliation_status = {}
        self.chains_not_in_header = []

        self.hets_classified = False

        self.sim_ready = {}

        self.anchors = self._build_anchor_dict()


        return

    def _build_anchor_dict(self):
        """
        Create dictionary detailing the start and end of connected residues
        in different moltypes and the max/min distances for connecting bonds.

        @rtype   :  dictionary
        @return  :  Dictionary containing information on bonds connecting
                    residues for different moltypes
        """

        bond_anchors = {}
        bond_anchors['protein'] = {
            'pre': 'C',
            'post': 'N',
            'min_dist': 1.1,
            'max_dist': 1.5
        }
        bond_anchors['nucleic'] = {
            'pre': "O3'",
            'post': 'P',
            'min_dist': 1.1,
            'max_dist': 1.5
        }

        return bond_anchors

    def read_pdb(self, filename, **kwargs):
        """
        Read coordinates and associated information from a PDB. Sequence and
        heterogen information is extracted. Also checks to ensure that some
        atoms are present and that chain names have been provided.

        @type  filename:  string
        @param filename:  Path to PDB file
        """

        self.pdbname = os.path.basename(filename)

        kwargs['pdbscan'] = True
        super(SasMolScan, self).read_pdb(filename, **kwargs)

        self.header_data = header_reader.PdbHeader(
            text=self.header(), parse=True)

        self.charmm = [False] * self.natoms()
        self.md_ready = [False] * self.natoms()

        self.basic_checks()

        return

    def basic_checks(self):
        """
        Conduct basic checks on any PDB that has been read in. Currently just 
        that we have some atoms and that all contiguous atoms within a
        residue have agreeing chain/segnames.
        """

        if self.pdbname:

            if self.natoms() == 0:
                raise IOError('Invalid PDB ({0:s}) - No Atoms!'.format(
                    self.pdbname))

            self.resid_subdivision_match()

        return

    def resid_subdivision_match(self):
        """
        Check that all atoms within each residue have matching chain and
        segname labels.

        Scans through atoms capturing lists of segname and chain. When resid
        changes the lists are checked to ensure they contain only one
        consistent label. Where possible a fix is attempted.

        @todo: Check this works in the case of a single residue as first chain
               followed by a second starting with the same resid.
        """

        resids = self.resid()
        segnames = self.segname()
        chains = self.chain()

        natoms = self.natoms()

        # Initialise scan_resid to ensure capture of first atoms data
        scan_resid = resids[0]

        # Lists to capture atomic data for each residue
        scan_ndxs = []
        scan_chains = []
        scan_segnames = []

        # Lists to hold lists of atoms in residues with an issue
        broken_chain_ndxs = []
        broken_segname_ndxs = []

        for i in range(natoms):

            # Check lists on change of residue or end of molecule
            if resids[i] != scan_resid or i == natoms - 1:

                # All atoms in a residue should have the same chain and segname
                if len(set(scan_chains)) > 1:
                    broken_chain_ndxs.append(scan_ndxs)

                if len(set(scan_segnames)) > 1:
                    broken_segname_ndxs.append(scan_ndxs)

                # Initialize scan of next residue
                scan_resid = resids[i]

                scan_chains = [chains[i]]
                scan_segnames = [segnames[i]]
                scan_ndxs = [i]

            else:

                scan_chains.append(chains[i])
                scan_segnames.append(segnames[i])
                scan_ndxs.append(i)

        # Try to fix any residues which have multiple subdivision labels
        # Update chain/segname lists which have been altered (i.e. when
        # subdivisions are removed).

        if broken_chain_ndxs:
            self.correct_multiple_resid_subdivision(broken_chain_ndxs, 'chain')
            self.setChains(list(set(self.chain())))
            self.setNumber_of_chains(len(self.chains()))

        if broken_segname_ndxs:
            self.correct_multiple_resid_subdivision(
                broken_segname_ndxs, 'segname')
            self.setSegnames(list(set(self.segname())))
            self.setNumber_of_segnames(len(self.segnames()))

        return

    def correct_multiple_resid_subdivision(self, broken, subdiv_type):
        """
        Take a list containing lists of the atom indices of residues with
        multiple subdivision (segname or chain) labels and try to fix the
        labelling so there is only one label for each residue.

        @type broken :       list
        @param broken:       List of lists. Each of the sublists containing
                             atom indices for a single residue
        @type subdiv_type :  string
        @param subdiv_type:  Type of subdivision to correct - 'chain' or
                             'segname'
        """

        n_broken = len(broken)

        for i in range(n_broken):

            # Do not attempt to fix if the following residue is also broken
            # as we rely upon surrounding residues for fixing
            if (i == n_broken -
                    1) or (not broken[i][-1] + 1 == broken[i + 1][0]):
                self.correct_resid_subdivision(broken[i], subdiv_type)
            else:
                resid1 = self.resid()[n_broken[i][0]]
                resid2 = self.resid()[n_broken[i + 1][0]]
                raise Exception(
                    'Residues around residues {0:d} and {1:d} have multiply defined {2:s}s '.format(
                        resid1, resid2, subdiv_type))

        return

    def correct_resid_subdivision(self, ndxs, subdiv_type):
        """
        If possible we want to assign the same subdivision (segname or chain)
        label to each atom of the residue whose atom indices are listed in
        ndxs. The label is matched to that of the flanking residues.

        @type ndxs        :  list
        @param ndxs       :  List containing atom indices for a single residue.
        @type subdiv_type :  string
        @param subdiv_type:  Type of subdivision to correct - 'chain' or
                             'segname'.
        """

        resids = self.resid()

        # Get the resid for the residue being corrected
        resid = resids[ndxs[0]]

        subdivs = self.subdiv(subdiv_type)

        try:
            prev_subdiv = subdivs[ndxs[0] - 1]
            prev_resid = resids[ndxs[0] - 1]
        except:
            prev_subdiv = None
            prev_resid = None

        try:
            next_subdiv = subdivs[ndxs[-1] + 1]
            next_resid = resids[ndxs[-1] + 1]
        except:
            next_subdiv = None
            next_resid = None

        continuous = utils.is_sequential_numbers(
            [prev_resid, resid, next_resid])

        if continuous and prev_subdiv == next_subdiv:
            # Naming error happens in the middle of otherwise well ordered
            # section of the subdivision
            for ndx in ndxs:
                subdivs[ndx] = prev_subdiv
        else:
            # Other cases can get complicated - lets just throw up our hands
            # (at least for now)
            subdiv_list = '/'.join(set(subdivs[ndxs[0]:ndxs[-1]]))
            raise Exception(
                'Multiple {0:s}s ({1:s}) assigned within residue {2:d}'.format(
                    subdiv_type, subdiv_list, resid))

        return

    def subdiv(self, subdiv_type):
        """
        Return list of subdiv(isions) for each atom according to chosen type
        (chain or segname)

        @type subdiv_type  :  string
        @param subdiv_type:  Choice of subdivision type - chain or segname
        @rtype             :  list
        @return:           :  List of subdivision of chain labels for each atom
        """

        if subdiv_type == 'segname':
            subdivs = self.segname()
        elif subdiv_type == 'chain':
            subdivs = self.chain()
        else:
            raise ValueError(
                'An invalid subdivision type was selected - must be either chain or segname')

        return subdivs

    def select_resid_subdiv(self, resid, subdiv, subdiv_type='chain'):
        """
        Create a selection text to use to obtain information about the atoms
        within a residues specified by residue number within a subdivision
        (chain or segname).

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param subdiv      :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @rtype             :  string
        @return:           :  Selection text to use to select atoms in the
                              residue specified
        """

        sel_text = '{0:s}[i] == "{1:s}" and resid[i] == {2:d}'.format(
            subdiv_type, subdiv, resid)

        return sel_text

    def set_property_residue(
            self, resid, subdiv, attribute, value, subdiv_type='chain'):
        """
        Set selected property for a chosen residue to value.

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param  subdiv     :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type attribute    :  string
        @param  attribute  :  Property to be set for each residue
        @type value        :  ?
        @param  value      :  Value to give to each atom in the residue
        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname

        """

        target = getattr(self, attribute)

        ndxs = self.ndxs_of_residue(resid, subdiv, subdiv_type=subdiv_type)

        for ndx in ndxs:
            target()[ndx] = value

        return

    def ndxs_residue_from_ndx(self, ndx):
        """
        Get list of indexes in the same residue as the input index.

        @type ndx  :  integer
        @param  ndx:  Index of selected atom in SasMol object
        @rtype     :  numpy.array
        @return    :  Indices of the atoms in the same residue as input
        """

        if ndx < self.natoms():
            resid = self.resid()[ndx]
            chain = self.chain()[ndx]
            segname = self.segname()[ndx]

            sel_text = 'chain[i] == "{0:s}" and segname[i] == "{1:s}" and resid[i] == {2:d}'.format(
                chain, segname, resid)
            ndxs = self.ndxs_from_selection(sel_text)

        else:

            ndxs = []

        return ndxs

    def ndxs_from_selection(self, selection):
        """
        Get indices for the given selection

        @type selection  :  string
        @param  selection:  Text to select atoms from SasMol object
        @rtype           :  numpy.array
        @return:         :  Indices of the selected atoms
        """

        err, mask = self.get_subset_mask(selection)
        if not err:
            ndxs = self.get_indices_from_mask(mask)
        elif err[0].startswith('found no atoms using filter selection'):
            ndxs = np.array([])
        else:
            raise Exception('\n'.join(err))

        return ndxs

    def ndxs_of_residue(self, resid, subdiv, subdiv_type='chain'):
        """
        Get indices the atoms within a residues specified by residue number
        within a subdivision (chain or segname).

        @type resid        :  integer
        @param  resid      :  Residue number to use in selection
        @type subdiv       :  string
        @param  subdiv     :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname
        @rtype             :  numpy.array
        @return:           :  Indices of the selected atoms
        """

        selection = self.select_resid_subdiv(
            resid, subdiv, subdiv_type=subdiv_type)

        return self.ndxs_from_selection(selection)

    def get_initial_chain_info(self):
        """
        Populates the self.chain_info object with data based on the resids
        found in each chain. Also classifies HETATMs depending on their
        connectivity to protein/nucleic acid residues.

          -  if in chain with protein/nucleic moltype residues and have
             appropriate bonds the change moltype to match rest of chain and set
             residue_flag
          -  if have a CONECT record to a protein nucleic acid but no backbone
             bond then just set residue_flag
        """

        info = self.chain_info

        info.n_models = self.number_of_frames()
        info.subdivs = self.chains()

        # Populates info.heterogens and reclassifies heterogens that are
        # bonded as part of protein/nucleic chains
        self.process_hets_by_subdiv('chain')

        # Evaluates following CONECT record information:
        # Numbering gaps, non-chain links HET to protein/nucleic, disulphide
        # bonds
        self.basic_conect_analysis()

        self.extract_sequence_info(subdiv_type='chain')

        return

    def extract_sequence_info(self, subdiv_type='chain', selected_subdivs=[]):
        """
        Scan atoms to get sequences for all protein/nucleic acid subdiv(ision)s
        (i.e. chain or sgnames as selected). As chains are built up missing
        residues detected from gaps in resid numbering are also recorded.

        @type subdiv_type      : string
        @param subdiv_type     : Choice of subdiv(ision) type - chain or segname
        @type selected_subdivs : list
        @param selected_subdivs: List of subdiv(isions) for which sequence 
                                 needs to be extracted
        """

        subdivs = self.subdiv(subdiv_type)
        resids = self.resid()
        resnames = self.resname()
        moltypes = self.moltype()

        if subdiv_type == 'segname':
            info = self.segname_info
        else:
            info = self.chain_info

        # Get list of tuples (subdiv, resid, resnames, moltypes) for all
        # protein and nucleic residues so we can build sequences
        resid_description = [
            x for x in zip(
                subdivs,
                resids,
                resnames,
                moltypes) if x[3] in [
                'protein',
                'rna',
                'dna',
                'nucleic']]
        resid_description = utils.uniquify_list(resid_description)

        for subdiv, residues in groupby(resid_description, lambda x: x[0]):

            if not selected_subdivs or subdiv in selected_subdivs:

                info.initialize_missing_resids_subdiv(subdiv)

                last_seen = None

                for residue in residues:

                    current_resid = residue[1]
                    current_resname = residue[2]

                    if last_seen:

                        expected_resid = last_seen + 1

                        # Record missing residues when resids non-seqential
                        if current_resid != expected_resid:

                            # Don't record as missing if previously found numbering
                            # is non-sequential but a bond exists
                            if not info.known_number_gap(
                                    subdiv, expected_resid, current_resid):

                                resid_miss = range(
                                    expected_resid, current_resid)
                                resname_miss = [''] * len(resid_miss)

                                info.add_missing_resids(
                                    subdiv, resid_miss, resname_miss)

                                info.add_residues_to_sequence(
                                    subdiv, resname_miss, resid_miss)

                    # Record residue in sequence, format:
                    # {subdiv:[(resid, resname),..]}
                    info.add_residue_to_sequence(
                        subdiv, current_resid, current_resname)

                    last_seen = current_resid

        return

    def check_atoms_disulphide(self, ndx1, ndx2):
        """
        Check if the atoms at the two input indices are of type SG and in
        cysteine.

        @type ndx1  :  integer
        @param  ndx1:  Atom index believed to be part of disulphide bond
        @type ndx2  :  integer
        @param  ndx2:  Atom index believed to be part of disulphide bond
        """

        resnames = self.resname()
        names = self.name()

        bond1 = (resnames[ndx1] == 'CYS') and (names[ndx1] == 'SG')
        bond2 = (resnames[ndx2] == 'CYS') and (names[ndx2] == 'SG')

        return (bond1 and bond2)

    def sasmol_index_from_original_index(self, original_ndx):
        """
        Take an original atom index and return the index for the same atom in
        the form used by sasmol.

        @type original_ndx :  integer
        @param original_ndx:  Original atom index as read from PDB
        @rtype             :  integer
        @return            :  Index of atom in sasmol
        """

        original_ndxs = self.original_index()

        # TODO: Need to check robustness of this
        match = np.where(original_ndxs == original_ndx)[0]

        if match:
            ndx = match[0]
        else:
            ndx = None
            self.logger.warning(
                'No atom located for original index {0:d} found in CONECT record'.format(original_ndx))

        return ndx

    def basic_conect_analysis(self, subdiv_type='chain'):
        """
        Run through CONECT records to determine if any HET atoms are connected
        to proteins and flag appropriately, disulphide bonds and gaps in chain
        numbering.
        If CONECT found between atoms in residues of moltype 'other' and
        'protein' then change residue_flag of HET to True.

        @type subdiv_type :  string
        @param subdiv_type:  Choice of subdiv(ision) type - chain or segname
        """

        conect_recs = self.conect()
        moltypes = self.moltype()
        resids = self.resid()
        subdivs = self.subdiv(subdiv_type)
        residue_flags = self.residue_flag()

        if subdiv_type == 'segname':
            info = self.segname_info
        else:
            info = self.chain_info

        # CONECT record formatted as dictionary:
        # {base_atom_index: [linked_index1, linked_index2,...]}
        # Indices are original indices as read from the source PDB
        for base, linked in conect_recs.iteritems():

            base_ndx = self.sasmol_index_from_original_index(base)

            if base_ndx:

                base_moltype = moltypes[base_ndx]
                base_resid = resids[base_ndx]
                base_subdiv = subdivs[base_ndx]

                # Check for HET links to protein/nucleic residues
                if base_moltype == 'other' and not residue_flags[base_ndx]:

                    for old_ndx in linked:

                        new_ndx = self.sasmol_index_from_original_index(
                            old_ndx)

                        if new_ndx:

                            if moltypes[new_ndx] in ['protein', 'rna', 'dna', 'nucleic']:

                                # Set residue_flag for the HET residue containing base
                                # atom
                                self.set_property_residue(base_resid, base_subdiv,
                                                          'residue_flag', True,
                                                          subdiv_type=subdiv_type)

                                break

                elif base_moltype in ['protein', 'rna', 'dna', 'nucleic']:

                    for old_ndx in linked:

                        new_ndx = self.sasmol_index_from_original_index(
                            old_ndx)

                        if new_ndx:

                            new_moltype = moltypes[new_ndx]
                            new_resid = resids[new_ndx]
                            new_subdiv = subdivs[new_ndx]

                            # Process HET link to protein/nucleic
                            if new_moltype == 'other' and not residue_flags[new_ndx]:

                                # Set residue_flag for HET residue containing
                                # linked atom
                                self.set_property_residue(new_resid, new_subdiv,
                                                          'residue_flag', True,
                                                          subdiv_type=subdiv_type)

                            elif (new_moltype == base_moltype
                                  ) and (
                                    base_moltype in ['protein', 'rna', 'dna', 'nucleic']):

                                # Check for number gaps
                                if (new_subdiv == base_subdiv):

                                    if abs(base_resid - new_resid) != 1:

                                        if self.check_chain_bond(
                                                base_resid, new_resid, base_subdiv,
                                                base_moltype, subdiv_type=subdiv_type):

                                            info.add_number_gap(
                                                base_subdiv, base_resid, new_resid)

                                # Check for disulphides
                                if self.check_atoms_disulphide(base_ndx, new_ndx):

                                    info.add_disulphide(
                                        base_subdiv, base_resid,
                                        new_subdiv, new_resid,
                                        subdiv_type=subdiv_type)

        return

    def process_hets_by_subdiv(self, subdiv_type='chain'):
        """
        Record list of heterogens in self.chain_info.heterogens and reassign
        moltype from 'other' if connected to a recognized chain type. Residues
        with altered type are flagged in residue_flag.

        @type subdiv_type : string
        @param subdiv_type: Choice of subdiv(ision) type - chain or segname
        """

        subdivs = self.subdiv(subdiv_type)
        moltypes = self.moltype()
        resids = self.resid()
        resnames = self.resname()

        if subdiv_type == 'segname':
            info = self.segname_info
        else:
            info = self.chain_info

        hets = [
            x for x in zip(
                subdivs,
                resids,
                resnames,
                moltypes) if x[3] == 'other']

        hets = utils.uniquify_list(hets)

        for het in hets:
            info.add_heterogen(het[0], het[1], het[2])

        self.check_subdivs_other_moltype(subdiv_type)

        return

    def check_subdivs_other_moltype(self, subdiv_type):
        """
        Check subdiv(isions) that contain recognised and 'other' residues to
        see if the other residues are part of the chain. If so, moltype is
        altered to be consistent with the recognized chain type and the
        residue_flag is set to True for all atoms in the relevant residues.

        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname
        """

        moltypes = self.moltype()
        subdivs = self.subdiv(subdiv_type)
        resids = self.resid()
        residue_flags = self.residue_flag()

        get_next_residue_ndxs = self.get_next_residue_ndxs

        # Get array of indices of atoms before changes in moltype
        tmp_moltypes = np.array(moltypes)
        ndxs_moltype_switch = np.where(
            tmp_moltypes[:-1] != tmp_moltypes[1:])[0]

        for ndx in ndxs_moltype_switch:

            # Need index of first atom after the moltype changes
            switch_ndx = ndx + 1

            # Only looking to check moltype changes within subdivisions
            if subdivs[ndx] == subdivs[switch_ndx]:

                subdiv = subdivs[ndx]

                this_moltype = moltypes[ndx]
                next_moltype = moltypes[switch_ndx]

                if this_moltype == 'other':

                    # If 'other' before moltype switch we need to scan lower
                    # resids in chain for residues to which HET might bond
                    direction = -1
                    # Select first residue to scan (before switch)
                    selection = self.select_resid_subdiv(
                        resids[ndx], subdiv, subdiv_type=subdiv_type)
                    current_ndxs = self.ndxs_from_selection(selection)

                elif next_moltype == 'other':

                    # If 'other' after moltype switch we need check higher
                    # resids in the chain for residues to which HET might bond
                    direction = 1
                    # Select first residue to scan (after switch)
                    selection = self.select_resid_subdiv(
                        resids[switch_ndx], subdiv, subdiv_type=subdiv_type)
                    current_ndxs = self.ndxs_from_selection(selection)

                else:
                    # Only interested if moltype change involves 'other'
                    # moltype - i.e. HET residues
                    continue

                # Scan along chains checking if HET atoms are connected to
                # protein/nucleic residues using standard bonds (i.e. N-C)
                while True:

                    if not current_ndxs.any():
                        break

                    if not moltypes[current_ndxs[0]] == 'other':
                        break

                    current_resid = resids[current_ndxs[0]]

                    # Check to see if current resid is linked to a recognised
                    # type by a standard backbone bond - if so get new_moltype
                    # (otherwise None)
                    if direction == 1:
                        new_moltype = self.check_moltype_other_residue(
                            current_resid, subdiv, subdiv_type=subdiv_type, check_dir='pre')

                    elif direction == -1:
                        new_moltype = self.check_moltype_other_residue(
                            current_resid, subdiv, subdiv_type=subdiv_type, check_dir='post')

                    # If bond to chain found update residue moltype
                    if new_moltype:

                        for i in current_ndxs:

                            moltypes[i] = new_moltype
                            residue_flags[i] = True

                        current_ndxs = get_next_residue_ndxs(
                            current_resid, subdiv, direction=direction, subdiv_type=subdiv_type)

                    else:
                        break

        if 'other' not in moltypes:
            self.setMoltypes(sorted(list(set(moltypes))))

        self.hets_classified = True

        return

    def get_next_residue_ndxs(
            self, resid, subdiv, direction=1, subdiv_type='chain'):
        """
        Return the indexes of the atoms in the next residues within the
        subdivison (chain or segname). Next can mean the next in either
        direction along the subdivision - direction = 1 means the direction of
        increasing resid/index,  -1 means the opposite direction.

        @type resid        :  integer
        @param  resid      :  Residue number
        @type subdiv       :  string
        @param  subdiv     :  Subdiv(ision), either chain or segname, label
        @type direction    :  integer
        @param  direction  :  1 = direction of increasing resid/index, -1 =
                              direction of decreasing resid/index
        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname
        @rtype             :  np.array
        @return            :  Indices of atoms in the next residue
        """

        subdivs = self.subdiv(subdiv_type)
        resids = self.resid()
        natoms = self.natoms()

        selection = self.select_resid_subdiv(
            resid, subdiv, subdiv_type=subdiv_type)
        #selection = '{0:s} {1:s} and resid {2:d}'.format(subdiv_type, subdiv, resid)
        current_ndxs = self.ndxs_from_selection(selection)

        found_ndxs = np.array([])

        if direction == 1:
            new_ndx = current_ndxs[-1] + 1
        else:
            new_ndx = current_ndxs[0] - 1

        if 0 <= new_ndx < natoms:

            if subdivs[new_ndx] == subdiv:

                new_resid = resids[new_ndx]
                selection = self.select_resid_subdiv(
                    new_resid, subdiv, subdiv_type=subdiv_type)
                #selection = '{0:s} {1:s} and resid {2:d}'.format(subdiv_type, subdiv, new_resid)
                found_ndxs = self.ndxs_from_selection(selection)

        return found_ndxs

    def check_moltype_other_residue(
            self, resid, subdiv, subdiv_type='chain', check_dir='both'):
        """
        Return a moltype for the current residue based on those of the residues
        connected directly to it within the same subdivision (chain/segname).

        @type resid        :  integer
        @param  resid      :  Residue number to use in selection
        @type subdiv       :  string
        @param  subdiv     :  Subdiv(ision), either chain or segname, label
        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname
        @type check_dir    :  string
        @param  check_dir  :  Direction ofto check for bonds. 'post' = check
                              for a bond with a higher resid, 'pre' = check for
                              a bond with a lower resid, 'both' = check both
                              directions
        """

        resids = self.resid()
        natoms = self.natoms()

        new_moltype = None

        get_resid_moltype = self.get_resid_moltype

        selection = self.select_resid_subdiv(
            resid, subdiv, subdiv_type=subdiv_type)
        #selection = '{0:s} {1:s} and resid {2:d}'.format(subdiv_type, subdiv, resid)
        ndxs = self.ndxs_from_selection(selection)

        # Initialize flags for detected bonds
        pre_bond = False
        post_bond = False

        if check_dir in ['pre', 'both']:

            if ndxs[0] - 1 >= 0:
                pre_resid = resids[ndxs[0] - 1]

                # Moltype of previous residue
                pre_moltype = get_resid_moltype(pre_resid, subdiv, subdiv_type)
                pre_bond = self.check_chain_bond(pre_resid, resid, subdiv,
                                                 pre_moltype, subdiv_type)

        if check_dir in ['post', 'both']:

            if ndxs[-1] + 1 < natoms:
                post_resid = resids[ndxs[-1] + 1]
                # Moltype of following residue
                post_moltype = get_resid_moltype(
                    post_resid, subdiv, subdiv_type)
                post_bond = self.check_chain_bond(resid, post_resid, subdiv,
                                                  post_moltype, subdiv_type)

        if pre_bond and post_bond and pre_moltype == post_moltype:
            new_moltype = pre_moltype
        elif pre_bond:
            new_moltype = pre_moltype
        elif post_bond:
            new_moltype = post_moltype
        else:
            new_moltype = None

        return new_moltype

    def get_resid_moltype(self, resid, subdiv, subdiv_type='chain'):
        """
        Get the moltype for the specified residue

        @param  resid      :  Residue number to use in selection
        @type subdiv       :  string
        @param  subdiv     :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param  subdiv_type:  Choice of subdivision type - chain or segname
        @rtype             :  string
        @return            :  moltype of selected residue
        """

        moltypes = self.moltype()
        sel_moltype = None

        selection = self.select_resid_subdiv(
            resid, subdiv, subdiv_type=subdiv_type)
        #selection = '{0:s} {1:s} and resid {2:d}'.format(subdiv_type, subdiv, resid)
        ndxs = self.ndxs_from_selection(selection)

        if ndxs.any():
            sel_moltype = moltypes[ndxs[0]]

        return sel_moltype

    def check_chain_bond(self, resid1, resid2, subdiv,
                         moltype, subdiv_type='chain'):
        """
        Check whether the geometry of the input residues is consistent with
        being bound as part of a polymer of moltype.

        @type resid1       : string
        @param resid1      : Previous residue ID
        @type resid2       : string
        @param resid2      : Current residue ID
        @type subdiv       : string
        @param subdiv      : Fragment name
        @type moltype      : string
        @param moltype     : moltype of bond to check
        @type subdiv_type  : string
        @param  subdiv_type: Choice of subdivision type - chain or segname
        @rtype             : boolean
        @return            : Does the geometry indicate residues are bound
                             together in a chain?
        """

        logger = self.logger

        bond_exists = False

        if moltype in self.anchors:

            bond_anchors = self.anchors[moltype]

            coors = self.coor()[0]

            # Find anchor atoms for bonds between residues
            selection1 = '{0:s}[i] == "{1:s}" and resid[i] == {2:d} and name[i] == "{3:s}"'.format(
                subdiv_type, subdiv, resid1, bond_anchors['pre'])
            #selection1 = '{0:s} {1:s} and resid {2:d} and name {3:s}'.format(subdiv_type, subdiv, resid1, bond_anchors['pre'])
            pre_ndx = self.ndxs_from_selection(selection1)

            selection2 = '{0:s}[i] == "{1:s}" and resid[i] == {2:d} and name[i] == "{3:s}"'.format(
                subdiv_type, subdiv, resid2, bond_anchors['post'])
            #selection2 = '{0:s} {1:s} and resid {2:d} and name {3:s}'.format(subdiv_type, subdiv, resid2, bond_anchors['post'])
            post_ndx = self.ndxs_from_selection(selection2)

            if len(pre_ndx) and len(post_ndx):

                pre_ndx = pre_ndx[0]
                post_ndx = post_ndx[0]

                coor1 = coors[pre_ndx]
                coor2 = coors[post_ndx]

                bond_length = np.linalg.norm(coor1 - coor2)

                bond_exists = bond_anchors[
                    "min_dist"] < bond_length < bond_anchors["max_dist"]

            else:
                if not pre_ndx.any():
                    logger.warning(
                        'Checked resid {0:d} in {1:s} for {2:s} bond but found no {3:s} atom'.format(
                            resid1, subdiv, moltype, bond_anchors['pre']))

                if not post_ndx.any():
                    logger.warning(
                        'Checked resid {0:d} in {1:s} for {2:s} bond but found no {3:s} atom'.format(
                            resid2, subdiv, moltype, bond_anchors['post']))

        return bond_exists

    def get_duplicate_resids(self, subdiv_type='chain'):
        """
        Check for duplicate resids within subdivisions of specified type
        (segname or chain).

        @type subdiv_type :  string
        @param subdiv_type:  Type of subdivision to correct - 'chain' or
                             'segname'.
        @rtype            :  list
        @return           :  List of tuples (subdiv_id, resid)

        @todo: Decide if I need this now we have changed how segnames and
               altlocs are read ans assigned
        """

        duplicates = []

        subdivs = self.subdiv(subdiv_type)

        resids = self.resid()
        natoms = self.natoms()

        # Create list of tuples of the form: ((subdiv_id, resid), index)
        # Sort brings together atoms with same (subdiv_id, resid) if separated
        atom_ids = sorted(zip(zip(subdivs, resids), range(natoms)))

        # Group by (subdiv_id, resid)
        for subdiv_res, grpd_atoms in itertools.groupby(atom_ids,
                                                        key=lambda x: x[0]):

            # Create a list of the atom indices for each (subdiv_id, resid)
            # pair
            subdiv_res_ndxs = [atom[1] for atom in grpd_atoms]

            if not self.is_sequential_numbers(subdiv_res_ndxs):
                duplicates.append(subdiv_res)

        return duplicates

    def check_chains(self):

        # TODO: Do I need this anymore?
        # TODO: convert print statements to exceptions

        if not self.chains_exist():
            pass

        else:

            chains = self.chain()

            if ' ' in chains:
                print("Empty chain IDs are a real pain!")
            else:
                dup_resids = self.get_duplicate_resids('chain')
                print(dup_resids)

                # Check for numbering gaps & Check missing residues

        return

    def next_segname_generator(self, existing=None):
        """
        Return the next unused name for segments. A list of existing names can
        be provided, otherwise the default is to use self.segnames(). The names
        are simple combinations of alphanumeric characters (A-Z then 0-9)
        growing in length from one to four characters.

        type existing :  list
        param existing:  List of the segment names which have already been used
        rtype         :  string
        return        :  Name for new segment
        """

        # Use only alphanumeric characters in segnames
        alpha_chars = list(map(chr, range(ord('A'), ord('Z') + 1)))
        valid_chars = alpha_chars + [str(x) for x in range(10)]

        if not existing:
            segnames = self.segnames()
        else:
            segnames = existing

        next_segname = ''
        found_new = False

        for i in range(1, 5):

            if not found_new:

                for x in itertools.product(valid_chars, repeat=i):

                    name = ''.join(x)

                    if name not in segnames:
                        next_segname = name
                        found_new = True
                        break
            else:
                break

        return next_segname

    def name_split_segment(self, old_name, existing=None, **kwargs):
        """
        Create names for segments generated through splitting existing
        segments. If possible keep some encoding based on the original
        segment name. When not possible or if a blank old name is provided then
        the names are simple combinations of alphanumeric characters
        (A-Z then 0-9) growing in length from one to four characters. In all
        cases existing provides a list of names to avoid, if no l;ist is passed
        then self.segnames() is used to provide a list of existing names.

        type old_name :  string
        param old_name:  Name of the existing segment (the one being split)
        type existing :  list
        param existing:  List of the segment names which have already been used
        rtype         :  string
        return        :  Name for new segment

        """

        if not existing:
            existing = self.segnames()

        if 'moltype' in kwargs:
            moltype = kwargs['moltype']
        else:
            moltype = None

        # Use only alphanumeric characters in segnames
        alpha_chars = list(map(chr, range(ord('A'), ord('Z') + 1)))
        valid_chars = alpha_chars + [str(x) for x in range(10)]

        new_name = ''
        found_new = False

        old_name = old_name.strip()

        # If the old segname is not blank we prefer a related name
        if old_name:

            if moltype:
                # Format: 1st letter old segname + letter indicating moltype +
                # ??
                prefix = old_name[0] + moltype[0].upper()
                maxval = 99
                free_chars = 2

            elif (len(old_name) > 1) and (
                    old_name[0] in alpha_chars) and (
                    old_name[1] in alpha_chars):
                # We give 1st and 2nd character significance - so preserve if found
                # Format: 1st letter old segname + 2nd letter old segname + ??

                prefix = old_name[0:2].upper()
                maxval = 99
                free_chars = 2

            else:
                # Format: 1st letter old segname + ???
                prefix = old_name[0]
                maxval = 999
                free_chars = 3

            # Evaluate trial segnames:
            #    - add number + spaces to fill ???
            #    - accept if doesn't already exist
            for x in range(1, maxval):

                tmp_name = prefix + str(x).zfill(free_chars)

                if tmp_name not in existing:
                    new_name = tmp_name
                    found_new = True
                    break

        # If we don't have a new name yet - try combinations of 1, 2, 3 and 4
        # characters selected from valid_chars
        if not found_new:
            new_name = self.next_segname_generator(existing=existing)

        if not new_name:
            raise Exception('Ran out of segnames')

        return new_name

    def update_segnames(self, new_segnames):
        """
        Use the new_segnames list to update the segnames for all atoms. Also
        update the list of unique segnames.

        @type  new_segnames:  list
        @param new_segnames:  List of segnames for each atom in the structure
        """

        self.setSegname(new_segnames)
        self.setSegnames(sorted(list(set(new_segnames))))

        return

    def copy_chains_segnames(self):
        """
        Copy the chain names for each atom into the segname.
        """

        chains = self.chain()
        self.update_segnames(chains)

        return

    def suggest_subdivisions(self, subdiv_type, skip=[]):
        """
        Use the new_segnames list to update the segnames for all atoms. Also
        update the list of unique segnames.

        @type  subdiv_type:  string
        @param subdiv_type:  Which subdivision type ('segname' or 'chain') is
                             to be updated
        @type  skip       :  list
        @param skip       :  List of subdivisions (chains/segnames) that
                             should not be altered
        @rtype            :  list
        @return           :  List of updated subdivisions (chains/segnames) for
                             each atom
        """

        logger = self.logger
        logger.debug('in suggest_subdivisions')

        # TODO make moltype checks
        # TODO improve moltype checks to deal with MODRES

        subdivs = self.subdiv(subdiv_type)
        start_subdivs = list(set(subdivs))
        resids = self.resid()

        # Record sections - runs of residues where resid goes up
        sections = []

        last_seen = -1000

        for r in resids:
            if r >= last_seen and sections:
                sections.append(sections[-1])
            elif not sections:
                sections.append(0)
            else:
                sections.append(sections[-1] + 1)
            last_seen = r

        # subsections = matched sections with subdivs (chains/segnames)
        sub_secs = zip(subdivs, sections)

        unq_sub_secs = set(sub_secs)
        n_unq_sub_secs = len(unq_sub_secs)

        # If there are more unique subsections that subdivs then
        # we need to break the subdivs down
        # if n_unq_sub_secs == len(set(subdivs)):
        if n_unq_sub_secs == len(start_subdivs):
            if start_subdivs not in [[''], [' ']]:
                new_subdiv_names = subdivs
            else:
                subdivs = ['A'] * len(subdivs)
        else:
            if n_unq_sub_secs > len(sections):

                sections = []

                last_seen = sections[0]
                curr_section = 0

                for sub_sec in sub_secs:
                    if sub_sec != last_seen:
                        curr_section += 1
                    sections.append(curr_section)

            sub_secs = zip(subdivs, sections)

            new_subdiv_names = []

            seen = {}

            for subdiv in skip:
                seen[subdiv] = subdiv

            for sub_sec in sub_secs:

                if sub_sec[0] in seen:
                    # Preserve skipped
                    new_subdiv_names.append(sub_sec[0])

                elif sub_sec in seen:
                    new_subdiv_names.append(seen[sub_sec])

                else:
                    new_name = self.name_split_segment(
                        sub_sec[0], existing=seen.values())

                    seen[sub_sec] = new_name
                    new_subdiv_names.append(new_name)

        return new_subdiv_names

    def split_segments_moltype(self):
        """
        Update segnames such that there is only one moltype in each. If HET
        atoms have not been scanned to check if they are part of another chain
        this is performed prior to segment splitting.
        """

        logger = self.logger
        logger.debug('in split_segments_moltype')

        if not self.hets_classified:
            self.check_subdivs_other_moltype('segname')

        segnames = self.segname()
        moltypes = self.moltype()
        natoms = self.natoms()

        new_segnames = []

        added_segnames = []
        edit = False

        last_segname = segnames[0]
        last_moltype = moltypes[0]

        segtype = {}

        for ndx in range(natoms):

            segname = segnames[ndx]
            moltype = moltypes[ndx]

            if segname not in segtype:
                segtype[segname] = moltype

            if segtype[segname] != moltype:

                if segname == last_segname and moltype == last_moltype:

                    new_segnames.append(new_segnames[-1])

                else:

                    added_segnames.append(
                        self.name_split_segment(
                            segname,
                            existing=self.segnames() +
                            new_segnames))

                    new_segnames.append(added_segnames[-1])
            else:
                new_segnames.append(segname)

            last_segname = segname
            last_moltype = moltype

        self.update_segnames(new_segnames)

        return

    def subdivisions_exist(self):
        """
        Check that either segname or chain labels exist in molecule.

        rtype :  boolean
        return:  Does the loaded data contain entries for either chain or
                 segname for any atoms?
        """

        return self.chains_exist() or self.segments_exist()

    def chains_exist(self):
        """
        Check that chain labels exist in molecules.

        rtype :  boolean
        return:  Does the loaded data contain chain entries for any atoms?
        """

        return not self.chains() == [' ']

    def segments_exist(self):
        """
        Check that segname labels exist in molecules.

        rtype :  boolean
        return:  Does the loaded data contain segname entries for any atoms?
        """

        return not (self.segnames() in [[''], [' ']])

    def segname_chain_same(self):
        """
        Return True if the entries for segname and chain align for all atoms
        and False if not.

        @rtype :  boolean
        @return:  Are entries for segname and chain aligned for all atoms?
        """

        same = False

        if len(self.chains()) == len(self.segnames()):

            chains = self.chain()
            segnames = self.segnames()

            chain_segname_pairs = utils.uniquify_list(zip(chains, segnames))

            if len(chain_segname_pairs) == len(chains):
                same = True

        return same

    def all_resids_have_subdivs(self, subdiv_type='chain'):
        """
        Return True if all atoms have a non-blank entry for the chosen
        subdiv_type (chain or segname) and False otherwise.

        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @rtype             :  boolean
        @return            :  Do all atoms have a atoms have a non-blank entry
                              for the chosen subdiv_type (chain or segname)
        """

        all_assigned = True

        subdivs = list(self.subdiv(subdiv_type))

        if '' in subdivs or ' ' in subdivs:
            all_assigned = False

        return all_assigned

    def copy_chain_segname(self, chain):
        """
        Copy the values from self.chain to self.segname for the chosen chain.

        @type chain :  string
        @param chain:  Chain ID to be copied to segname
        """

        segnames = self.segname()
        chains = self.chain()

        ndxs = self.ndxs_from_selection('chain[i] == "{0:s}"'.format(chain))

        for ndx in ndxs:
            segnames[ndx] = chains[ndx]

        self.setSegnames(sorted(list(set(segnames))))

        return

    def blank_segment_copy_chain(self):
        """
        Copy the values from self.chain to self.segname for all blank segnames.

        """

        segnames = self.segname()
        chains = self.chain()

        ndxs = self.ndxs_from_selection('segname[i] == ""')

        for ndx in ndxs:
            segnames[ndx] = chains[ndx]

        self.setSegnames(sorted(list(set(segnames))))

        return

    def segname_chain_one2one(self, segname):
        """
        Do the segname and chain values line up one to one per atom for the 
        selected segname.

        @type segname :  string
        @param segname:  Segname selected to be checked for match with chain.

        """

        one2one = ''

        segnames = self.segname()
        chains = self.chain()

        chain_segname_map = set(zip(chains, segnames))
        chain_from_seg = [x for x in chain_segname_map if x[1] == segname]
        seg_from_chain = [
            x for x in chain_segname_map if x[0] in chain_from_seg]

        if len(chain_from_seg) == 1 and len(seg_from_chain) == 1:
            one2one = chain_from_seg[0]

        return one2one

    def segname_moltype_subsection_chain(self, segname):
        """
        Is the segname assigned to a subset of a chain where all the atoms 
        have the same moltype.

        @type segname :  string
        @param segname:  Segname selected to be checked for match with chain.

        """
        subsection = ''

        segnames = self.segname()
        chains = self.chain()
        moltypes = self.moltype()

        chain_segname_map = set(zip(chains, segnames))
        segname_moltype = [x for x in set(
            zip(segnames, moltypes)) if x[0] == segname]
        chain_from_seg = [x for x in chain_segname_map if x[1] == segname]

        if len(chain_from_seg) == 1 and len(segname_moltype) == 1:
            subsection = chain_from_seg[0]

        return subsection

    def initial_segmentation(self):
        """
        Creates an initial guessed segmentation of the protein if none is
        present. If no segnames are present and chains are available they are
        copied into the segname. If neither segnames or chains are available
        segnames guessed on assumption that runs of increasing resids are
        indicative of segments. If incomplete segnames are found then first
        chains are used to fill and then resid runs.

        The data behind self.segname and self.segnames are both modifield.
        """

        logger = self.logger
        logger.debug('in initial_segmentation')

        header_reconciled_chains = self.header_reconciled_chains

        chains_exist = self.chains_exist()
        segments_exist = self.segments_exist()

        #segments_all = self.all_resids_have_subdivs(subdiv_type='segname')

        segments_to_save = []
        matched_chains = set([])
        segments_to_edit = []

        if segments_exist and header_reconciled_chains:

            for segname in self.segnames():

                matched_chain = self.segname_chain_one2one(segname)

                if not matched_chain:
                    matched_chain = self.segname_moltype_subsection_chain(
                        segname)

                if matched_chain in header_reconciled_chains:
                    segments_to_save.append(segname)
                    matched_chains.add(matched_chain)
                else:
                    segments_to_edit.append(segname)

        elif chains_exist and not segments_exist:
            self.copy_chains_segnames()

        elif chains_exist and segments_exist:

            if not self.segname_chain_same():

                if '' in self.segnames():

                    self.blank_segment_copy_chain()

                else:

                    logger.warning(
                        'Chain/Segname mismatch found in initial structure')

        #self.update_segnames(self.suggest_subdivisions('segname',skip = segments_to_save))

        if not self.all_resids_have_subdivs('segname'):

            skip_segnames = self.segnames()
            if '' in skip_segnames:
                skip_segnames.remove('')
            if ' ' in skip_segnames:
                skip_segnames.remove(' ')

            new_segnames = self.suggest_subdivisions('segname',
                                                     skip=segments_to_save)

            self.update_segnames(new_segnames)

        return

    def strip_hydrogens_list(self, atom_list):
        """
        Remove hydrogen atoms from list of atom names. Hydrogen atoms are
        identified as having names where the first alphabetical character is H.

        type atom_list :  list
        param atom_list:  List of atom names
        rtype          :  list
        return         :  List of all of the non-hydrogen atoms in input list
        """

        return [x for x in atom_list if x.strip('0123456789')[0] != 'H']

    def strip_heavy_list(self, atom_list):
        """
        Remove non-hydrogen atoms from list of atom names. Hydrogen atoms are
        identified as having names where the first alphabetical character is H.

        type atom_list :  list
        param atom_list:  List of atom names
        rtype          :  list
        return         :  List of all of the hydrogen atoms in input list
        """

        return [x for x in atom_list if x.strip('0123456789')[0] == 'H']

    def resid_altloc_label(self, resid, subdiv, subdiv_type='chain'):
        """
        Assign loc labels to copies of coordinates with a residue (starting
        from 'A' as in PDB AltLoc). Return list of the locations labelled.

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param subdiv      :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @type loc          :  string
        @param loc         :  Select a altloc for which to get the atom names
        @rtype             :  list
        @return:           :  List of the locations labelled

        @todo: Check if this is robust if first altloc in a properly specified
               PDB has less atoms than later ones
        """

        cnt = {}

        ndxs = self.ndxs_of_residue(resid, subdiv, subdiv_type=subdiv_type)

        locs = self.loc()
        names = self.name()

        chr_start = ord('A') - 1

        seen_locs = set([])

        for ndx in ndxs:

            cur_name = names[ndx]

            if cur_name in cnt:

                cnt[cur_name] += 1

                loc_label = chr(cnt[cur_name] + chr_start)

            else:

                cnt[cur_name] = 1

                cur_label = locs[ndx]

                if cur_label in [' ', 'A']:
                    loc_label = cur_label
                else:
                    loc_label = chr(cnt[cur_name] + chr_start)

            locs[ndx] = loc_label

            seen_locs.add(loc_label)

        return list(seen_locs)

    def get_atom_names_residue(
            self,
            resid,
            subdiv,
            subdiv_type='chain',
            loc=None):
        """
        Return list of all atom names within the selected residue (if chosen
        give for the specified AltLoc).

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param subdiv      :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @type loc          :  string
        @param loc         :  Select a altloc for which to get the atom names
        @rtype             :  list
        @return:           :  List of atom names for selected residue
        """

        names = self.name()
        locs = self.loc()

        ndxs = self.ndxs_of_residue(resid, subdiv, subdiv_type=subdiv_type)

        checked_names = []

        for ndx in ndxs:
            if not loc or locs[ndx] in [loc, ' ']:
                checked_names.append(names[ndx])

        return checked_names

    def get_resid_resname(self, resid, subdiv, subdiv_type='chain'):
        """
        Return the resname of the selected residue.

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param subdiv      :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @rtype             :  string
        @return:           :  Residue name for selected residue
        """

        resnames = self.resname()
        ndxs = self.ndxs_of_residue(resid, subdiv, subdiv_type=subdiv_type)

        return resnames[ndxs[0]]

    def residue_atom_check(self, resid, subdiv, subdiv_type='chain'):
        """
        Check residue to find missing atoms (heavy and hydrogen) and determine
        if it is available in CHARMM forcefield and complete and ready for
        simulation. Missing hydrogens are only counted as naming is so variable
        in possible inputs.

        Two substitutions are made to correct for CHARMM naming oddities:
           1. name CD1 -> CD in ILE residues
           2. resname HIS -> HSE - standard choice for histidine protonation

        @type resid        :  integer
        @param resid       :  Residue number to use in selection
        @type subdiv       :  string
        @param subdiv      :  Subdiv(ision), either chain or segname, label to
                              use in selection
        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @rtype             :  list, list, integer, integer, boolean, boolean
        @return:           :    1.  List of heavy atoms missing in coordinates
                                2.  List of extra heavy atoms in coordinates
                                3.  Number of missing hydrigen atoms
                                4.  Number of excess hydrogen atoms
                                5.  Is the residue found in the CHARMM topology
                                    file
                                6.  Is the residues complete - all heavy and
                                    hydrogen atoms and the latter correctly
                                    named for use in CHARMM forcefield


        @todo: Check that terminal residues are handled correctly
        """

        resname = self.get_resid_resname(
            resid, subdiv, subdiv_type=subdiv_type)

        coor_atoms = self.get_atom_names_residue(
            resid, subdiv, subdiv_type=subdiv_type)

        if resname in utils.nucleic_res_charmm_dict:
            resname = utils.nucleic_res_charmm_dict[resname]
            for i in range(len(coor_atoms)):
                if coor_atoms[i][-1] == '*':
                    coor_atoms[i] = coor_atoms[i][:-1] + "'"
                elif coor_atoms[i] == 'OP1':
                    coor_atoms[i] = 'O1P'
                elif coor_atoms[i] == 'OP2':
                    coor_atoms[i] = 'O2P'
                elif resname == 'THY' and coor_atoms[i] == 'C7':
                    coor_atoms[i] = 'C5M'

        heavy_missing = []
        heavy_excess = []

        n_hyd_missing = 0
        n_hyd_excess = 0

        altlocs = []

        # Standard HIS in CHARMM is HSE (epsilon protonated)
        if resname == 'HIS':
            charmm = False
            resname = 'HSE'
        else:
            charmm = True

        hyd_naming_correct = False

        t = self.topology

        # Extract the atoms for each residue from CHARMM forcefield topology
        # If residue not found in topology then it is not available in CHARMM
        try:
            ff_atoms = set([x[0] for x in t.topology_info[resname]['ATOM']])

        except KeyError:

            charmm = False

            return heavy_missing, heavy_excess, n_hyd_missing, n_hyd_excess, hyd_naming_correct, charmm, altlocs

        # Get list of heavy and hydrogen atoms from CHARMM forcefield
        ff_heavy = set(self.strip_hydrogens_list(ff_atoms))
        ff_hyd = set(self.strip_heavy_list(ff_atoms))

        coor_heavy = self.strip_hydrogens_list(coor_atoms)

        # There may be multiple conformers (AltLocs) found in the coordinates
        # Select only the first
        if len(coor_heavy) != len(set(coor_heavy)):

            altlocs = self.resid_altloc_label(
                resid, subdiv, subdiv_type=subdiv_type)

            first_altloc = next(s for s in altlocs if s.strip())

            coor_atoms = self.get_atom_names_residue(
                resid, subdiv, subdiv_type=subdiv_type, loc=first_altloc)

        unq_coor_names = set(coor_atoms)

        if not ff_atoms.symmetric_difference(unq_coor_names):

            # TODO: Do we want to pass this residue as ready for simulation?
            hyd_naming_correct = True

        else:

            coor_heavy = set(self.strip_hydrogens_list(unq_coor_names))
            coor_hyd = set(self.strip_heavy_list(unq_coor_names))

            # Compare forcefield and coordinate heavy atoms - record
            # discrepancies
            if ff_heavy != coor_heavy:

                heavy_missing = list(ff_heavy - coor_heavy)
                heavy_excess = list(coor_heavy - ff_heavy)

            # Compare forcefield and hydrogens - record discrepancies or full
            # match
            if coor_hyd == ff_hyd:

                hyd_naming_correct = True

            else:

                n_ff_hyd = len(ff_hyd)
                n_coor_hyd = len(coor_hyd)

                if n_ff_hyd == n_coor_hyd:
                    pass

                elif n_ff_hyd - n_coor_hyd < 0:
                    n_hyd_excess = n_coor_hyd - n_ff_hyd
                else:
                    n_hyd_missing = n_ff_hyd - n_coor_hyd

        # Deal with CHARMM/PDB inconsistencies with naming of atoms
        if resname == 'ILE' and 'CD' in heavy_missing and 'CD1' in heavy_excess:
            heavy_missing.remove('CD')
            heavy_excess.remove('CD1')

        # Do the atoms in the coordinates perfectly match the forcefield
        if hyd_naming_correct and not (
                n_hyd_missing or n_hyd_excess) and not (
                heavy_missing or heavy_excess):
            complete = True
        else:
            complete = False

        return heavy_missing, heavy_excess, n_hyd_missing, n_hyd_excess, charmm, complete, altlocs

    def check_all_residues(self, subdiv_type='segname', skip=[]):
        """
        Check all residues to find missing atoms (heavy and hydrogen) and 
        determine if they are available in CHARMM forcefield and complete and 
        ready for simulation. Missing hydrogens are only counted as naming is 
        so variable in possible inputs.

        @type subdiv_type  :  string
        @param subdiv_type :  Choice of subdivision type - chain or segname
        @type skip         :  list
        @param skip        :  List of subdiv(isions) not to be scanned
        """

        logger = self.logger

        logger.debug('in check_all_residues')

        resids = self.resid()
        resnames = self.resname()
        subdivs = self.subdiv(subdiv_type)

        # TODO: how do we handle non-standard chains/segments?

        if subdiv_type == 'segname':
            info = self.segname_info
        else:
            info = self.chain_info

        # Loop through each residues of each subdivision to check it has the
        # correct atoms
        for subdiv in set(subdivs):

            if subdiv in skip:
                continue

            logger.debug('Checking residues in {0:s} {1:s}'.format(
                subdiv_type, subdiv))

            first_ndx = subdivs.index(subdiv)
            resid = resids[first_ndx]
            ndxs = self.ndxs_of_residue(resid, subdiv, subdiv_type=subdiv_type)

            while len(ndxs):

                resid = resids[ndxs[0]]
                resname = resnames[ndxs[0]]

                # Check to see if the resname is in the CHARMM forcefield, if
                # so compare the coordinate atoms to those described in the
                # forcefield topology
                # AltLocs also processed to ensure only one used for comparison
                (heavy_missing,
                 heavy_excess,
                 n_hyd_missing,
                 n_hyd_excess,
                 charmm,
                 complete,
                 altlocs) = self.residue_atom_check(resid,
                                                    subdiv,
                                                    subdiv_type=subdiv_type)

                # If the resid is in the CHARMM forcefield topology record
                # check if any atoms are missing or named incorrectly
                if charmm:

                    if heavy_missing:
                        info.add_missing_atoms(
                            subdiv, resid, resname, heavy_missing)

                    if n_hyd_missing:
                        info.add_missing_h(
                            subdiv, resid, resname, n_hyd_missing)

                    if heavy_excess:
                        info.add_excess_atoms(
                            subdiv, resid, resname, heavy_excess)

                    if n_hyd_excess:
                        info.add_excess_h(subdiv, resid, resname, n_hyd_excess)

                    for ndx in ndxs:
                        self.charmm[ndx] = True

                        if complete:
                            self.md_ready[ndx] = True

                if len(altlocs) > 1:
                    info.add_altlocs(subdiv, resid, altlocs)

                # Get the indices for the next resid (or an empty list if
                # none)
                ndxs = self.get_next_residue_ndxs(
                    resid, subdiv, subdiv_type=subdiv_type)

        return

    def get_segment_sequence_info(self):
        """
        Extract sequence from coordinates and then fill in the segname_info with
        resnames taken from chain_info (if this has been filled from header
        information).

        """

        # Get coordinate sequence and gap information
        self.extract_sequence_info(subdiv_type='segname')

        model_no = self.model_no

        resids = self.resid()
        segnames = self.segname()
        chains = self.chain()
        moltypes = self.moltype()
        natoms = self.natoms()

        chain_info = self.chain_info
        segname_info = self.segname_info

        seg_chain_map = np.array(zip(segnames, chains))

        # Get indexes where segname-chain combination changes
        seg_chain_ends = np.where(seg_chain_map[:-1] != seg_chain_map[1:])[0]

        # Get indices for the start of every contiguous segname-chain run
        seg_chain_starts = seg_chain_ends + 1
        seg_chain_starts = np.append([0], seg_chain_starts)

        # Add final index for completeness
        seg_chain_ends = np.append(seg_chain_ends, [natoms - 1])

        # Need the start and end of chains to know when to copy terminal regions
        # Same logic applied as for segname-chain combination
        chain_map = np.array(zip(segnames, moltypes))
        chain_ends = np.where(chain_map[:-1, 1] != chain_map[1:, 1])[0]

        chain_starts = chain_ends + 1
        chain_starts = np.append([0], chain_starts)

        chain_ends = np.append(chain_ends, [natoms - 1])

        # Need to copy all information from chains that map to segnames so
        # loop through all segname-chain combination regions
        for i in range(len(seg_chain_ends)):

            start_ndx = seg_chain_starts[i]
            end_ndx = seg_chain_ends[i]

            if moltypes[start_ndx] in ['nucleic', 'protein']:

                segname = segnames[start_ndx]
                chain = chains[start_ndx]

                chain_missing_res = chain_info.missing_resids[model_no][chain]

                if chain in chain_info.number_gaps:
                    chain_num_gaps = chain_info.number_gaps[chain]
                else:
                    chain_num_gaps = {}

                    # If we are at a chain terminus then chain_info may contain
                # information of preceding/following residues - copy this in
                if start_ndx in chain_starts:

                    start_resid = resids[start_ndx]

                    for resid, resname in chain_missing_res.iteritems():

                        if resid < start_resid:

                            segname_info.add_missing_resid(
                                segname, resid, resname, model_no=model_no)

                if end_ndx in chain_ends:

                    end_resid = resids[end_ndx]

                    for resid, resname in chain_missing_res.iteritems():

                        if resid > end_resid:

                            segname_info.add_missing_resid(
                                segname, resid, resname, model_no=model_no)

                # Copy in information on gaps in the coordinate sequence

                start_resid = resids[start_ndx]
                end_resid = resids[end_ndx]

                for resid in range(start_resid, end_resid + 1):

                    if resid in chain_missing_res.keys():

                        if ((resid not in segname_info.sequence[segname]) or
                                (segname_info.sequence[segname][resid] == '')):

                            segname_info.add_missing_resid(segname,
                                                           resid,
                                                           chain_missing_res[resid],
                                                           model_no=model_no)

                    if resid in chain_num_gaps.keys():

                        segname_info.add_number_gap(segname, resid, chain_num_gaps[resid])

        # Complete segment sequences with resnames from missing residues
        seg_missing = segname_info.missing_resids[model_no]

        for segname in seg_missing.keys():

            seq = segname_info.sequence[segname]
            resid_list = [x[0] for x in seq]

            pre_seq = []

            for resid, resname in seg_missing[segname].iteritems():

                if (resid, '') in seq:

                    ndx = seq.index((resid, ''))
                    seq[ndx] = (resid, resname)

                elif resid < resid_list[0]:

                    pre_seq.append((resid, resname))

                elif resid > resid_list[-1]:

                    segname_info.add_residue_to_sequence(
                        segname, resid, resname)

            segname_info.prepend_residues_to_sequence(segname, pre_seq)

        return

    def advanced_het_recognition(self):
        # This will one day scan coordinates and CONECT records to see if
        # we can identify classes of HET such as glycans.

        return

    def copy_biomt_segments(self):
        """
        Copy BIOMT data from self.chain_info to self.segname_info. All 
        segments associated with each chain is included separately in the 
        segname_info version of the information.
        """

        biomt = self.chain_info.biomt

        # If there is BIOMT information from the chain scan copy to segname
        # Replacing the chains to be transformed by the equivalent segnames
        if biomt:

            self.segname_info.biomt = copy.deepcopy(biomt)
            seg_biomt = self.segname_info.biomt

            segnames = self.segname()
            chains = self.chain()

            # Set up chain -> segname mapping
            ch_seg = zip(chains, segnames)

            ch_seg = utils.uniquify_list(ch_seg)

            chain_seg_map = {}

            for chain, segment in ch_seg:

                if chain not in chain_seg_map:
                    chain_seg_map[chain] = [segment]
                else:
                    chain_seg_map[chain].append(segment)

            # Replace the list of chains to which each BIOMT record applies
            # with relevant segment names
            for rec_no, record in seg_biomt.items():

                new_subdivs = []

                for chain in record['subdivs']:
                    new_subdivs += chain_seg_map[chain]

                record['subdivs'] = new_subdivs

        return

    def check_residues_charmm_ready(self):
        """
        Check that the residues in each segment are all CHARMM compatible.
        Returns a dictionary giving the status of each segment.

        @rtype :  dictionary
        @return:  Boolean for each segment stating if all residues are CHARMM
                  compatible.
        """

        segnames = self.segname()
        charmm_ready = self.charmm

        seg_charmm_valid = {}

        seg_charmm = zip(segnames, charmm_ready)

        for segment, charmm in groupby(seg_charmm, lambda x: x[0]):

            charmm_checks = list(charmm)

            if sum([x[1] for x in charmm_checks]) == len(charmm_checks):
                seg_charmm_valid[segment] = True
            else:
                seg_charmm_valid[segment] = False

        return seg_charmm_valid

    def check_residues_md_ready(self):
        """
        Check that the residues in each segment are all ready for MD
        simulation, this requires both CHARMM compatible residue and
        a match of atoms to the topology. Returns a dictionary
        giving the status of each segment.

        @rtype :  dictionary
        @return:  Boolean for each segment stating if all residues are CHARMM
                  compatible.
        """

        segnames = self.segname()
        md_ready = self.md_ready

        seg_md_valid = {}

        seg_md = zip(segnames, md_ready)

        for segment, md in groupby(seg_md, lambda x: x[0]):

            if sum([x[1] for x in md]):
                seg_md_valid[segment] = True
            else:
                seg_md_valid[segment] = False

        return seg_md_valid

    def check_segname_simulation_preparedness(self):
        """
        Creates a dictionary summarizing the structures preparedness in several
        categories (which are used keys) for each segname:
        chain - True if there are no missing residues (i.e. chains complete)
        charmm - True if all residues are found in the CHARMM topology
        md - True if all residues in the CHARMM topology and all atoms present with no additions
        single conformer - True if charmm criteria met and no altloc flags are found

        @attention: Alters self.sim_ready
        """

        seg_charmm_valid = self.check_residues_charmm_ready()
        seg_md_valid = self.check_residues_md_ready()

        info = self.segname_info

        altloc = self.segname_info.altloc

        sim_ready = self.sim_ready

        for segname in self.segnames():

            sim_ready[segname] = {}

            sim_ready[segname]['charmm'] = seg_charmm_valid[segname]

            sim_ready[segname]['single_conformer'] = segname not in altloc

            if info.no_missing_resids(segname) == 0:

                mc_ready = self.check_segment_links(segname)

                sim_ready[segname]['chain'] = mc_ready

            else:

                sim_ready[segname]['chain'] = False

            sim_ready[segname]['md'] = (seg_md_valid[segname] and
                                        sim_ready[segname]['chain'] and
                                        sim_ready[segname]['single_conformer'])

            # We should warn people if start residue is not 1
            if segname in self.segname_info.sequence:
                first_resid = self.segname_info.sequence[segname][0][0]
            else:
                # TODO: This is a default - can something less problematic
                first_resid = 0

            sim_ready[segname]['start'] = (first_resid == 1)

        return

    def check_protein_residue_dihed_atoms(self, ndxs):
        """
        Test if the atoms represented by the input residue list contain
        the types necessary for protein backbone dihedral Monte Carlo moves in
        SASSIE, namely 'CA','C' and 'N'.

        @type  ndxs:  list
        @param ndxs:  List of indices corresponding to atoms in resdiue being checked
        @rtype     :  bool
        @return    :  Are the backbone atoms needed for dihedral rotation found in the structure
        """

        names = self.name()

        needed = set(['C', 'CA', 'N'])

        found = set([])

        for ndx in ndxs:
            found.add(names[ndx])

        return needed.issubset(found)

    def check_rna_residue_dihed_atoms(self, ndxs):
        """
        Test if the atoms represented by the input residue list contain
        the types necessary for dihedral Monte Carlo moves in SASSIE,
        namely either:
            "O3'","P","O5'","C5'","C4'","C3'"
        or  "O3*","P","O5*","C5*","C4*","C3*"

        @type  ndxs:  list
        @param ndxs:  List of indices corresponding to atoms in resdiue being checked
        @rtype     :  bool
        @return    :  Are the backbone atoms needed for dihedral rotation found in the structure
        """

        names = self.name()

        needed = set(["O3'", "P", "O5'", "C5'", "C4'", "C3'"])

        found = set([])

        for ndx in ndxs:
            found.add(names[ndx])

        correct = needed.issubset(found)

        if not correct:

            needed = set(["O3*", "P", "O5*", "C5*", "C4*", "C3*"])
            correct = needed.issubset(found)

        return correct

    def check_segment_links(self, segname):
        """
        @type  segname:  string
        @param segname:  Selected segment name
        @rtype:          bool
        @return:         Are all residues in the segment connected in a way
                         suitable for MC simulation in SASSIE
        """

        segnames = self.segname()
        resids = self.resid()
        moltypes = self.moltype()

        bound = []
        gaps = []
        mc_ready = []

        first_ndx = segnames.index(segname)
        resid = resids[first_ndx]
        ndxs = self.ndxs_of_residue(resid, segname, subdiv_type='segname')

        moltype = moltypes[ndxs[0]]

        if moltype not in ['protein', 'rna', 'dna', 'nucleic']:
            # return bound, gaps
            return False
        elif moltype == 'protein':
            check_atoms = self.check_protein_residue_dihed_atoms
            mc_atoms_present = check_atoms(ndxs)
        elif moltype == 'rna':
            check_atoms = self.check_rna_residue_dihed_atoms
            mc_atoms_present = check_atoms(ndxs)
        else:
            mc_atoms_present = False

        bound.append(resid)

        n_resid = 1

        if mc_atoms_present:
            mc_ready.append(resid)

        while len(ndxs):

            resid = resids[ndxs[0]]

            if moltype in ['protein', 'rna']:
                mc_atoms_present = check_atoms(ndxs)
            else:
                mc_atoms_present = False

            next_ndxs = self.get_next_residue_ndxs(
                resid, segname, direction=1, subdiv_type='segname')

            if len(next_ndxs) > 0:
                next_resid = resids[next_ndxs[0]]

                bonded = self.check_chain_bond(
                    resid, next_resid, segname, moltype, subdiv_type='segname')

                if bonded:
                    bound.append(next_resid)

                    if mc_atoms_present:
                        mc_ready.append(resid)

                else:
                    gaps.append(range(resid + 1, next_resid + 1))

                # Get the indices for the next resid (or an empty list if
                # none)
                ndxs = self.get_next_residue_ndxs(
                    resid, segname, subdiv_type='segname')

                n_resid += 1

            else:
                ndxs = []

        if len(mc_ready) == n_resid:
            all_mc_ready = True
        else:
            all_mc_ready = False

        return all_mc_ready

    def chain_segment_map(self):
        """
        Generate a dictionary mapping each chain to the corresponding
        segnames.

        @rtype :  dictionary
        @return:  Mapping of chains to segnames
        """

        segnames = self.segname()
        chains = self.chain()

        chain_segment_map = {}

        tmp_map = set(zip(chains, segnames))

        for chain, segname in tmp_map:

            if chain in chain_segment_map:
                chain_segment_map[chain].append(segname)
            else:
                chain_segment_map[chain] = [segname]

        return chain_segment_map

    def chain_segment_map2(self):
        """
        Generate a dictionary mapping each chain to the corresponding
        segnames.

        @rtype :  dictionary
        @return:  Mapping of chains to segnames
        """

        segnames = self.segname()
        chains = self.chain()
        moltypes = self.moltype()
        original_indices = self.original_index()
        resids = self.resid()

        chain_segment_map = {}

        tmp_map = zip(chains, segnames, resids, original_indices, moltypes)

        for chain, chain_atom_info in groupby(tmp_map, lambda x: x[0]):

            if chain not in chain_segment_map:
                chain_segment_map[chain] = {}

            chain_atom_info = list(chain_atom_info)
            seg_list = set(map(lambda x: x[1], chain_atom_info))

            for segname in seg_list:

                chain_segment_map[chain][segname] = {}

                atom_list = [(x[2], x[3], x[4])
                             for x in chain_atom_info if x[1] == segname]

                chain_segment_map[chain][segname]['resid'] = (
                    atom_list[0][0], atom_list[-1][0])
                chain_segment_map[chain][segname]['original_index'] = (
                    atom_list[0][1], atom_list[-1][1])
                chain_segment_map[chain][segname]['moltype'] = atom_list[0][2]

        return chain_segment_map

    def segment_chain_map(self):
        """
        Generate a dictionary mapping each segname to the corresponding
        chain(s).

        @rtype :  dictionary
        @return:  Mapping of segnames to chains
        """
        segnames = self.segname()
        chains = self.chain()

        segment_chain_map = {}

        tmp_map = set(zip(segnames, chains))

        for segname, chain in tmp_map:

            if chain in segment_chain_map:
                segment_chain_map[segname].append(chain)
            else:
                segment_chain_map[segname] = [chain]

        return segment_chain_map

    def calc_terminal_distance(self, filter_txt):
        """
        Calculate the distance between the termini of residues matching
        the input filter text.

        @type  filter_txt:   string
        @param filter_txt:   Filter describing atoms to be considered.
        """

        anchors = self.anchors

        error, mask = self.get_subset_mask(filter_txt)

        if mask.sum() == 0:
            raise ValueError('No atoms selected by filter text: ' + filter_txt)

        mol_sub = sasmol.SasMol(1)
        error = self.copy_molecule_using_mask(mol_sub, mask, 0)

        moltypes = set(mol_sub.moltype())

        if len(moltypes) > 1:
            raise NotImplementedError(
                "Cannot calculate end-to-end distances of selections with mixed moltypes")
        else:
            moltype = moltypes.pop()

        if moltype not in anchors:
            raise NotImplementedError(
                "Unsupported moltype " +
                moltype +
                " in get_terminal_distance")

        if len(error):
            raise Exception(
                "Something wrong when copying molecule for moltype " +
                moltype)

        anchors_moltype = anchors[moltype]

        resids = mol_sub.resid()
        names = mol_sub.name()
        coors = mol_sub.coor()[0]

        for i in range(mol_sub.natoms()):
            if resids[i] == resids[0] and names[i] == anchors_moltype["pre"]:
                coor_1 = coors[i]
                continue
            if resids[i] == resids[-1] and names[i] == anchors_moltype["post"]:
                coor_2 = coors[i]
                continue

        if "coor_1" not in locals():
            raise NotImplementedError("There are missing " +
                                      anchors_moltype["pre"] +
                                      " atom in resid " +
                                      str(resids[0]) +
                                      ", and we can't determine missing residues before it!")

        if "coor_2" not in locals():
            raise NotImplementedError("There are missing " +
                                      anchors_moltype["post"] +
                                      " atom in resid " +
                                      str(resids[-1]) +
                                      ", and we can't determine missing residues after it!")

        distance = np.linalg.norm(coor_1 - coor_2)

        return distance

    def calc_protein_subdiv_terminal_distances(self, subdiv, subdiv_type='chain'):
        """
        Calculate the distance between N and C termini for protein resdues in
        given chain.

        @type  subdiv:        string
        @param subdiv:        Choice of subdivision to be measured
        @type  subdiv_type:   string
        @param subdiv_type:   Which type of subdivision to use in selection -
                              chain or segname.
        @rtype:               float
        @return:              Distance between the termini of any protein
                              sections of the selected subdiv
        """

        filter_txt = '{0:s}[i]=="{1:s}" and moltype[i] == "protein"'.format(
            subdiv_type, subdiv)

        try:
            dist = self.calc_terminal_distance(filter_txt)
        except Exception as e:
            txt = str(e)
            self.logger.warning(txt)
            dist = 0.0

        return dist

    def build_molecule_dictionary(self, filter_text=''):
        """
        Create a dictionary of statistics calculated from the structure. Only
        the atoms which are selected by filter_text are included in the
        calculations.

        @type  filter_txt:   string
        @param filter_txt:   Filter describing atoms to be considered.
        """

        from collections import OrderedDict

        if filter_text:

            m1 = sasmol.SasMol(1)
            error, mask = self.get_subset_mask(filter_text)
            error = self.copy_molecule_using_mask(m1, mask, 0)

        else:
            m1 = self

        m1.initialize_children()

        molecule_dictionary = {}

        try:
            wt = m1.calcmass()
            molecule_dictionary[
                'Molecular weight'] = '%.0f Da : %.2f kDa' % (wt, wt / 1000)
        except:
            molecule_dictionary[
                'Molecular weight'] = 'FAILED TO CALCULATE MASS'

        try:
            com = m1.calccom(0)
            com = ['%.2f' % item for item in com]
            molecule_dictionary[
                'Center of mass (x,y,z)'] = '(' + ', '.join(com) + ') Angstrom'
        except:
            molecule_dictionary[
                'Center of mass (x,y,z)'] = 'FAILED TO CALCULATE CENTER OF MASS'

        try:

            formula = m1.calcformula()
            name_order = ['C', 'H', 'O', 'N', 'S', 'P']
            standard_formula = dict(
                (item,
                 value) for item,
                value in formula.items() if item in name_order)
            non_standard_formula = dict(
                (item,
                 value) for item,
                value in formula.items() if item not in name_order)
            standard_formula = OrderedDict(
                sorted(
                    standard_formula.items(),
                    key=lambda x: name_order.index(
                        x[0])))
            formula = [
                str(k) +
                '<sub>' +
                str(v) +
                '</sub>' for k,
                v in standard_formula.items()]
            non_standard_formula = OrderedDict(non_standard_formula.items())
            formula += [str(k) +
                        '<sub>' +
                        str(v) +
                        '</sub>' for k, v in non_standard_formula.items()]
            molecule_dictionary['Formula'] = ''.join(formula)
        except:
            molecule_dictionary[
                'Formula'] = 'FAILED TO CALCULATE CHEMICAL FORMULA'

        try:
            rg = m1.calcrg(0)
            molecule_dictionary['Radius of gyration'] = '%.2f Angstrom' % rg
        except:
            molecule_dictionary[
                'Radius of gyration'] = 'FAILED TO CALCULATE RADIUS OF GYRATION'

        try:
            range = m1.calcminmax()
            range1 = ['%.2f' % item for item in range[0]]
            range1 = '(' + ', '.join(range1) + ')'
            range2 = ['%.2f' % item for item in range[1]]
            range2 = '(' + ', '.join(range2) + ')'
            molecule_dictionary[
                'Coordinate range (x,y,z)'] = 'Low: ' + range1 + ', High: ' + range2 + ' Angstrom'
        except:
            molecule_dictionary[
                'Coordinate range (x,y,z)'] = 'FAILED TO CALCULATE MIN & MAX DIMENSIONS'

        name_order = [
            'Formula',
            'Non-redundant residue names',
            'Non-redundant atom names',
            'Molecular weight',
            'Center of mass (x,y,z)',
            'Coordinate range (x,y,z)',
            'Radius of gyration',
            'PMI']

        if len(m1.chains()) == 1:
            chid = m1.chains()[0]
            dist = self.calc_protein_subdiv_terminal_distances(
                chid, subdiv_type='chain')
            description = 'Distance between N- & C-terminals'
            name_order.append(description)
            molecule_dictionary[description] = '{0:.2f} Angstrom'.format(dist)

        ordered_molecule_dictionary = OrderedDict(
            sorted(
                molecule_dictionary.items(),
                key=lambda x: name_order.index(
                    x[0])))

        return ordered_molecule_dictionary

    def segment_scan(self, initialize=True):
        """
        Perform a scan to determine the contents of the files using 
        organization by segname, rather than chain. In many cases segnames 
        will not be provided, here they are created inline with how segments 
        are used in CHARMM - i.e. each segment must be a continuous chain of a 
        particular molecule type. 

        Information is copied from self.chain_info to self.segname_info (with 
        relevant mappings from chains to segnames).

        @type    initialize:  boolean
        @keyword initialize:  Is it necessary to create an initial
                              segmentation
        """

        info = self.segname_info
        logger = self.logger

        logger.debug('in segment_scan')
        logger.info('Beginning segment based scan')

        info.n_models = self.number_of_frames()

        if initialize:
            self.initial_segmentation()
            self.split_segments_moltype()

        info.subdivs = self.segnames()

        # Populates info.heterogens and reclassifies heterogens that are
        # bonded as part of protein/nucleic chains
        self.process_hets_by_subdiv(subdiv_type='segname')

        # Evaluates following CONECT record information:
        # Numbering gaps, non-chain links HET to protein/nucleic, disulphide
        # bonds
        self.basic_conect_analysis(subdiv_type='segname')

        if 'other' in self.moltypes():
            # Hook for future recognition of glycans, etc.
            self.advanced_het_recognition()

        self.check_all_residues(subdiv_type='segname')

        # self.get_segment_sequence()
        #
        # if self.header_data.read_valid_header:
        #     self.get_segment_missing_resid_info()

        self.get_segment_sequence_info()

        return

    def any_charmm_ready_segments(self):

        ready = False

        for segname in self.segnames():

            if self.sim_ready[segname]['charmm']:
                ready = True
                break

        return ready

    def run_scan(self):
        """
        A full scan involves the comparison of information on PDB contents from 
        the header and coordinates. The first scan checks that any header data
        corresponds to the coordinates in the PDB. A second scan works on the 
        basis of trying to form CHARMM style segments from the coordinates. The
        second scan also checks at the atomic detail that residue contents 
        agree with what is expected in the CHARMM topology read into the 
        scanner.
        """

        logger = self.logger
        header_data = self.header_data

        self.get_initial_chain_info()

        if not header_data.read_valid_header:
            logger.info(
                'No valid header detected so proceeding to segname based scan')

        elif not header_data.has_seq_info:
            logger.info(
                'No sequence information detected in header so proceeding to segname based scan')

        else:

            chain_reconciler = reconcile.Chain_Info_Reconciler(
                header_data.chain_info, self.chain_info)
            chain_reconciler.reconcile_chain_info()

            self.header_reconciled_chains = chain_reconciler.reconciled
            self.header_reconciliation_status = chain_reconciler.issues
            self.chains_not_in_header = chain_reconciler.chains_not_in_head

            if True in chain_reconciler.reconciled.values():
                chain_reconciler.copy_biomt()

        self.segment_scan()

        self.check_segname_simulation_preparedness()

        return

# Structure preparation methods - more for PDB Rx than Scan

    def offset_original_index(self, offset):

        original_indices = self.original_index()

        for index in original_indices:

            index += offset

        return
