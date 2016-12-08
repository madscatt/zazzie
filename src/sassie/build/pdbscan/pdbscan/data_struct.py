# -*- coding: utf-8 -*-
"""
Common data structures to hold information extracted from PDB coordinates and
headers

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

from itertools import groupby
from operator import itemgetter

from . import pdbscan_utils as utils


class Disulphide():
    """
    Class to hold information about disulphide bond from a PDB. Provides
    methods to access.
    """

    def __init__(self, subdiv1, resid1, subdiv2, resid2,
                 subdiv_type='chain'):
        """
        Initialize a disulphide records - recording the subdiv(ision), i.e
        chain or segname, and resid of the two residues involved. Also record
        subdiv_type.

        @type subdiv1      :  string
        @param subdiv1     :  Subdivision (chain or segname) label for first
                              residue
        @type resid1       :  integer
        @param resid1      :  Resid of first residue
        @type subdiv2      :  string
        @param subdiv2     :  Subdivision (chain or segname) label for second
                              residue
        @type resid2       :  string
        @param resid2      :  Resid of second residue
        @type subdiv_type  :  string
        @param subdiv_type :  Are the subdivisions 'chain' or 'segname'
        """

        if subdiv_type not in ['chain', 'segname']:
            raise ValueError('subdiv_type must be "chain" or "segname"')

        self.bound = [{'subdiv': subdiv1,
                       'resid': resid1},
                      {'subdiv': subdiv2,
                       'resid': resid2}]

        self.subdiv_type = subdiv_type

        return

    def subdivs(self):
        """
        Return list of subdivisions involved in the bond
        """
        return list({self.bound[0]['subdiv'], self.bound[1]['subdiv']})

    def intra_subdiv(self):
        """
        Return True if both residues are within the same subdivision
        """
        return self.bound[0]['subdiv'] == self.bound[1]['subdiv']

    def __str__(self):
        """
        Report bond as string: {subdiv2}{resid2}-{subdiv2}{resid2}.
        Note: curly brackets characters are just seperators, not part of
        output string
        """

        return "{0:s}{1:d}-{2:s}{3:d}".format(self.bound[0]['subdiv'],
                                              self.bound[0]['bond'],
                                              self.bound[1]['subdiv'],
                                              self.bound[1]['bond'])


class Info():
    """
    Holds information from scanning PDB subdivisions (either chain or segname)
    and provides access to information and summary counts etc.
    """

    def __init__(self, scan_type='chain'):
        """
        Initialize object to hold information from scanning PDB subdivisions
        (chains or segments/segnames).
        """

        self.logger = logging.getLogger(__name__)

        self.n_models = 1
        self.subdivs = []

        # key = subdiv, value = list of (resid,resname) - 3 character resnames
        self.sequence = {}

        # Missing residues dictionay of the form:
        # missing_resids = {model_no:{subdiv: {resid:resname}, ..}}
        self.missing_resids = {}

        # TODO: Need to decide on format to summarize breaks in chain
        self.gaps = {}

        # Many PDBs have gaps in numbering without gaps in the sequence
        # number_gaps[subdiv] = {resid before gap : resid after gap}
        self.number_gaps = {}

        # Missing atom dictionary of the form:
        # missing_atoms = {model_no:{subdiv: {resid:
        #                                          {'atoms':[atm1, atm2, ...],
        #                                           'resname': resname}
        #                 }}}
        self.missing_atoms = {}

        # Excess atom dictionary of the form:
        # excess_atoms = {model_no:{subdiv: {resid:
        #                                          {'atoms':[atm1, atm2, ...],
        #                                           'resname': resname}
        #                 }}}
        self.excess_atoms = {}

        # Number of missing hydrogens dictionary of the form:
        # n = {model_no:{subdiv: {resid: {'n': number, 'resname': resname}}}}
        self.n_missing_hydrogens = {}

        # Number of hydrogens foundin excess of CHARMM FF, dictionary of the form:
        # n = {model_no:{subdiv: {resid: {'n': number, 'resname': resname}}}}
        self.n_excess_hydrogens = {}

        # heterogens = {subdiv: {resid:resname}}
        self.heterogens = {}

        # alternative locations = {subdiv: {resid:[label1, label2, ...]}}
        self.altloc = {}

        # List of Disulphides
        self.disulphides = []

        # Dictionary of form:
        # biomt[biomol_no] = {
        #   'subdivs' : [],
        #   'auth_bio_unit' : '',
        #   'soft_bio_unit' : '',
        #   'rot' : [],
        #   'trans' : []
        # }
        # rot = list of np.array((3,3))
        # trans = list of np.array(3)
        self.biomt = {}

        # Map of subdivisions created by transformations from
        # those detailed in other attributes
        self.subdiv_map = {}

        # chain or segname based scan?
        if scan_type in ['chain', 'segname']:
            self.scan_type = scan_type
        else:
            raise ValueError('scan_type must be "chain" or "segname"')

        return

    def no_missing_atoms(self, subdiv=None, model_no=1):
        """
        Return the number of missing atoms for the molecule or individual
        subidiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) identifier
        @type model_no :  integer
        @param model_no:  Model number for which to get data
        @rtype         :  integer
        @return        :  Number of missing atoms in chosen model/subdivision
        """

        if model_no in self.missing_atoms:
            model_missing = self.missing_atoms[model_no]

            if subdiv:

                n_missing = sum([len(x['atoms']) for x in model_missing[subdiv]])
            else:
                missing = []
                for sub in model_missing:
                    missing += [len(x['atoms']) for x in model_missing[sub]]
                n_missing = sum(missing)

        else:
            n_missing = 0

        return n_missing

    def no_missing_h(self, subdiv=None, model_no=1):
        """
        Return the number of missing hydrogens for the molecule or individual
        subidiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) identifier
        @type model_no :  integer
        @param model_no:  Model number for which to get data
        @rtype         :  integer
        @return        :  Number of missing atoms in chosen model/subdivision
        """

        if model_no in self.n_missing_hydrogens:
            model_missing = self.n_missing_hydrogens[model_no]

            if subdiv:

                n_missing = sum([x['n']
                                 for resid, x in model_missing[subdiv].items()])
            else:
                n_missing = 0
                for sub in model_missing:
                    self.logger.info(
                        'Missing hydrogens in ' + model_missing[sub])
                    n_missing += sum([x['n']
                                      for resid, x in model_missing[sub].items()])

        else:
            n_missing = 0

        return n_missing

    def no_missing_resids(self, subdiv=None, model_no=1):
        """
        Return the number of missing residues for the molecule or individual
        subdiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) identifier
        @type model_no :  integer
        @param model_no:  Model number for which to get data
        @rtype         :  integer
        @return        :  Number of missing residues in chosen model/subdivision
        """

        if model_no in self.missing_resids:
            model_missing = self.missing_resids[model_no]

            if subdiv:

                if subdiv in model_missing:
                    n_missing = len(model_missing[subdiv])
                else:
                    n_missing = 0

            else:
                n_missing = sum([len(model_missing[x]) for x in model_missing])

        else:
            n_missing = 0

        return n_missing

    def no_inter_disulphides(self, subdiv=None):
        """
        Return the number of inter-chain/segname disulphides for the molecule
        or individual subdiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) identifier
        @rtype         :  integer
        @return        :  Number of inter-subdivision disulphide bonds in
                         model/subdivision
        """

        if subdiv:
            is_inter = [not x.intra_subdiv()
                        for x in self.disulphides if subdiv in x.subdivs()]
        else:
            is_inter = [not x.intra_subdiv() for x in self.disulphides]

        return sum(is_inter)

    def no_intra_disulphides(self, subdiv=None):
        """
        Return the number of intra-chain/segname disulphides for the molecule
        or individual subdiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) identifier
        @rtype         :  integer
        @return        :  Number of intra-subdivision disulphide bonds in
                         model/subdivision
        """

        if subdiv:
            is_intra = [x.intra_subdiv()
                        for x in self.disulphides if subdiv in x.subdivs()]
        else:
            is_intra = [x.intra_subdiv() for x in self.disulphides]

        return sum(is_intra)

    def create_biomol(self, biomol_no, subdivs=[], auth='', soft=''):
        """
        Create record to hold information on biological unit transformation

        @type biomol_no :  integer
        @param biomol_no:  Numerical biomolecule identifier
        @type subdivs   :  list
        @param subdivs  :  List of chains to which transformation applies
        @type auth      :  string
        @param auth     :  Author provided unit description
        @type soft      :  string
        @param soft     :  Software predicted unit description
        """

        self.biomt[biomol_no] = {
            'subdivs': subdivs,
            'auth_bio_unit': auth,
            'soft_bio_unit': soft,
            'rot': [],
            'trans': []
        }

        return

    def add_biomt(self, biomol_no, rot, trans):
        """
        Add a transformation (rotation + translation matrices) to biomolecule
        biomol_no.

        @type biomol_no :  integer
        @param biomol_no:  Numerical biomolecule identifier
        @type rot       :  numpy.array
        @param  rot     :  3 x 3 matrix describing molecular rotation
        @type trans     :  numpy.array
        @param  trans   :  vector describing molecular translation
        """

        # TODO add check for array types
        self.biomt[biomol_no]['rot'].append(rot)
        self.biomt[biomol_no]['trans'].append(trans)

        return

    def add_missing_atoms(self, subdiv, resid, resname, atoms, model_no=1):
        """
        Add missing atom list for a specific residue to the structure holding
        the information for the entire model:
        missing_atoms[model_no][subdiv][resid] = atoms

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) label
        @type resid    :  integer
        @param resid   :  Residue number
        @type resname  :  string
        @param resname :  Residue name
        @type atoms    :  list
        @param atoms   :  List of missing atom names
        @type model_no :  integer
        @param model_no:  Numeric model identifier
        """

        missing_atoms = self.missing_atoms

        if model_no not in missing_atoms:
            missing_atoms[model_no] = {}

        if subdiv not in missing_atoms[model_no]:
            missing_atoms[model_no][subdiv] = {}

        missing_atoms[model_no][subdiv][resid] = {
            'resname': resname,
            'atoms': atoms,
        }

        return

    def add_excess_atoms(self, subdiv, resid, resname, atoms, model_no=1):
        """
        Add excess atom list for a specific residue to the structure holding
        the information for the entire model:
        missing_atoms[model_no][subdiv][resid] = atoms

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) label
        @type resid    :  integer
        @param resid   :  Residue number
        @type resname  :  string
        @param resname :  Residue name
        @type atoms    :  list
        @param atoms   :  List of missing atom names
        @type model_no :  integer
        @param model_no:  Numeric model identifier
        """

        excess_atoms = self.excess_atoms

        if model_no not in excess_atoms:
            excess_atoms[model_no] = {}

        if subdiv not in excess_atoms[model_no]:
            excess_atoms[model_no][subdiv] = {}

        excess_atoms[model_no][subdiv][resid] = {
            'resname': resname,
            'atoms': atoms,
        }

        return

    def add_missing_h(self, subdiv, resid, resname, n_missing, model_no=1):
        """
        Add count of missing hydrogens for a specific residue to the structure
        holding the information for the entire model:
        missing_atoms[model_no][subdiv][resid] = atoms

        @type subdiv    :  string
        @param subdiv   :  Subdivision (chain/segname) label
        @type resid     :  integer
        @param resid    :  Residue number
        @type resname   :  string
        @param resname  :  Residue name
        @type n_missing :  integer
        @param n_missing:  Number of missing hydrogens
        @type model_no  :  integer
        @param model_no :  Numeric model identifier
        """

        n_missing_hydrogens = self.n_missing_hydrogens

        if model_no not in n_missing_hydrogens:
            n_missing_hydrogens[model_no] = {}

        if subdiv not in n_missing_hydrogens[model_no]:
            n_missing_hydrogens[model_no][subdiv] = {}

        n_missing_hydrogens[model_no][subdiv][resid] = {'resname': resname,
                                                        'n': n_missing}

        return

    def add_excess_h(self, subdiv, resid, resname, n_missing, model_no=1):
        """
        Add count of excess hydrogens for a specific residue to the structure
        holding the information for the entire model:
        missing_atoms[model_no][subdiv][resid] = atoms

        @type subdiv    :  string
        @param subdiv   :  Subdivision (chain/segname) label
        @type resid     :  integer
        @param resid    :  Residue number
        @type resname   :  string
        @param resname  :  Residue name
        @type n_missing :  integer
        @param n_missing:  Number of missing hydrogens
        @type model_no  :  integer
        @param model_no :  Numeric model identifier
        """

        n_excess_hydrogens = self.n_excess_hydrogens

        if model_no not in n_excess_hydrogens:
            n_excess_hydrogens[model_no] = {}

        if subdiv not in n_excess_hydrogens[model_no]:
            n_excess_hydrogens[model_no][subdiv] = {}

        n_excess_hydrogens[model_no][subdiv][resid] = {'resname': resname,
                                                       'n': n_missing}

        return

    def add_missing_resid(self, subdiv, resid, resname, model_no=1):
        """
        Add missing residue to the structure holding the information for the
        entire model:
        missing_resids[model_no][subdiv][resid] = resname

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) label
        @type resid    :  integer
        @param resid   :  Residue number
        @type resname  :  string
        @param resname :  Residue name corresponding to the resid
        @type model_no :  integer
        @param model_no:  Numeric model identifier
        """

        missing_resids = self.missing_resids

        if model_no not in missing_resids:
            missing_resids[model_no] = {}

        if subdiv not in missing_resids[model_no]:
            missing_resids[model_no][subdiv] = {}

        missing_resids[model_no][subdiv][resid] = resname

        return

    def get_first_coor_resid(self, segname, model_no=1):

        missing_resids = self.missing_resids[model_no]

        if segname in self.sequence:

            seq = self.sequence[segname]

            ndx = 0
            guess = seq[ndx][0]

            if segname in missing_resids:

                while guess in missing_resids[segname].keys():
                    ndx += 1
                    guess = seq[ndx][0]

        return seq[ndx]

    def add_missing_resids(self, subdiv, resids, resnames, model_no=1):
        """
        Add missing list of missing residues to the structure holding the
        information for the entire model:
        missing_resids[model_no][subdiv][resid] = resname

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) label
        @type resids   :  list
        @param resids  :  List of residue numbers
        @type resnames :  list
        @param resnames:  Residue names corresponding to resids
        @type model_no :  integer
        @param model_no:  Numeric model identifier
        """

        missing_resids = self.missing_resids

        if model_no not in missing_resids:
            missing_resids[model_no] = {}

        if subdiv not in missing_resids[model_no]:
            missing_resids[model_no][subdiv] = {}

        missing = zip(resids, resnames)

        for resid, resname in missing:

            missing_resids[model_no][subdiv][resid] = resname

        return

    def add_heterogen(self, subdiv, resid, resname):
        """
        Add heterogen information to the structure holding the information for
        the entire model:
        heterogens[subdiv][resid] = resname

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain/segname) label
        @type resid    :  integer
        @param resid   :  Residue number
        @type resname  :  string
        @param resname :  Residue name corresponding to the resid
        """

        heterogens = self.heterogens

        if subdiv not in heterogens:
            heterogens[subdiv] = {}

        heterogens[subdiv][resid] = resname

        return

    def add_altloc(self, subdiv, resid, loc_id):
        """
        Add alternative location information to the dictionary holding the
        information for the entire model:
        altloc[subdiv][resid] = [loc_id1, ....]

        @type subdiv :  string
        @param subdiv:  Subdivision (chain/segname) label
        @type resid  :  integer
        @param resid :  Residue number
        @type loc_id :  string
        @param loc_id:  Alternative location label
        """

        altloc = self.altloc

        if subdiv not in altloc:
            altloc[subdiv] = {}

        if resid in altloc[subdiv]:
            altloc[subdiv][resid].append(loc_id)
        else:
            altloc[subdiv][resid] = [loc_id]

        return

    def add_altlocs(self, subdiv, resid, loc_ids):
        """
        Add several alternative locations to the dictionary holding the
        information for the entire model:
        altloc[subdiv][resid] = [loc_id1, ....]

        @type subdiv  :  string
        @param subdiv :  Subdivision (chain/segname) label
        @type resid   :  integer
        @param resid  :  Residue number
        @type loc_ids :  list
        @param loc_ids:  Alternative location labels
        """

        for loc in sorted(loc_ids):
            self.add_altloc(subdiv, resid, loc)

        return

    def add_disulphide(self, subdiv1, resid1, subdiv2,
                       resid2, subdiv_type='chain'):
        """
        Add a disulphide record to list for entire model

        @type subdiv1      :  string
        @param subdiv1     :  Subdivision (chain or segname) label for first
                              residue
        @type resid1       :  integer
        @param resid1      :  Resid of first residue
        @type subdiv2      :  string
        @param subdiv2     :  Subdivision (chain or segname) label for second
                              residue
        @type resid2       :  string
        @param resid2      :  Resid of second residue
        @type subdiv_type  :  string
        @param subdiv_type :  Are the subdivisions 'chain' or 'segname'
        """

        if not self.disulphide_exists(
                subdiv1, resid1, subdiv2, resid2, subdiv_type='chain'):

            dis = Disulphide(subdiv1, resid1, subdiv2, resid2, subdiv_type)
            self.disulphides.append(dis)

        return

    def disulphide_exists(self, subdiv1, resid1, subdiv2,
                          resid2, subdiv_type='chain'):
        """
        Check to see if a disulphide record exists for the two specified
        residues

        @type subdiv1      :  string
        @param subdiv1     :  Subdivision (chain or segname) label for first
                              residue
        @type resid1       :  integer
        @param resid1      :  Resid of first residue
        @type subdiv2      :  string
        @param subdiv2     :  Subdivision (chain or segname) label for second
                              residue
        @type resid2       :  string
        @param resid2      :  Resid of second residue
        @type subdiv_type  :  string
        @param subdiv_type :  Are the subdivisions 'chain' or 'segname'
        """

        exists = False

        bound1 = {'subdiv': subdiv1,
                  'resid': resid1}

        bound2 = {'subdiv': subdiv2,
                  'resid': resid2}

        for dis in self.disulphides:
            if (bound1 in dis.bound) and (bound2 in dis.bound):
                exists = True
                break

        return exists

    def add_number_gap(self, subdiv, resid1, resid2):
        """
        Add item to dictionary of number gaps (self.number_gaps[subdiv]):
        key = resid1, value = resid2. resid1 = residue before gap, resid2 that
        after.

        @type subdiv  :  string
        @param subdiv :  Subdivision (chain or segname) label for first residue
        @type resid1  :  integer
        @param resid1 :  Resid of first residue
        @type resid2  :  string
        @param resid2 :  Resid of second residue
        """
        number_gaps = self.number_gaps

        # Missing resids are between the two input found resids
        first_resid = min([resid1, resid2])
        last_resid = max([resid1, resid2])

        if subdiv not in number_gaps:
            number_gaps[subdiv] = {}

        number_gaps[subdiv][first_resid] = last_resid

        return

    def known_number_gap(self, subdiv, resid1, resid2):
        """
        Check if a number gap is known to exist between the two input resids in
        subdiv.

        @type subdiv  :  string
        @param subdiv :  Subdivision (chain or segname) label for first residue
        @type resid1  :  integer
        @param resid1 :  Resid of first residue
        @type resid2  :  string
        @param resid2 :  Resid of second residue
        """

        known = False

        if subdiv in self.number_gaps:

            if resid1 in self.number_gaps[subdiv]:

                if self.number_gaps[subdiv][resid1] == resid2:

                    known = True

        return known

    def prepend_residue_to_sequence(self, subdiv, resid, resname):
        """
        Add a single residue to the start of the sequence of the selected
        subdiv(ision).

        @type subdiv  :  string
        @param subdiv :  Subdivision (chain or segname) label for residue of
                         interest
        @type resid   :  integer
        @param resid  :  Resid of the residue of interest
        @type resname :  string
        @param resname:  Residue name
        """

        sequence = self.sequence

        residue = (resid, resname)

        if subdiv not in sequence:
            sequence[subdiv] = [residue]
        else:
            sequence[subdiv].insert(0, residue)

        return

    def prepend_residues_to_sequence(self, subdiv, residues):
        """
        Add list of residues to the start of the sequence of the selected
        subdiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain or segname) label for residue of
                          interest
        @type residues:   list
        @param residues:  list of (resid, resname) pairs
        """

        for residue in reversed(residues):
            self.prepend_residue_to_sequence(subdiv, residue[0], residue[1])

        return

    def initialize_missing_resids_subdiv(self, subdiv, model_no=1):

        missing_resids = self.missing_resids

        if model_no not in missing_resids:
            missing_resids[model_no] = {}

        if subdiv not in missing_resids[model_no]:
            missing_resids[model_no][subdiv] = {}

    def add_residue_to_sequence(self, subdiv, resid, resname):
        """
        Add a single residue to the sequence of the selected subdiv(ision).

        @type subdiv  :  string
        @param subdiv :  Subdivision (chain or segname) label for residue of 
                         interest
        @type resid   :  integer
        @param resid  :  Resid of the residue of interest
        @type resname :  string
        @param resname:  Residue name
        """

        sequence = self.sequence

        residue = (resid, resname)

        if subdiv not in sequence:
            sequence[subdiv] = [residue]
        else:
            sequence[subdiv].append(residue)

        return

    def add_residues_to_sequence(self, subdiv, resnames, resids):
        """
        Add list of residues to the sequence of the selected subdiv(ision).

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain or segname) label for residue of 
                         interest
        @type resnames :  string
        @param resnames:  Residue names for all residues
        @type resids   :  integer
        @param resids  :  Resid of the residues of interest
        """

        residues = zip(resids, resnames)

        for residue in residues:
            self.add_residue_to_sequence(subdiv, residue[0], residue[1])

        return

    def add_subdiv_sequence(self, subdiv, resnames, resids=[]):
        """
        Add a sequence for the selected subdiv(ision) containing the selected 
        residue names. If no list of resids is provided (as is typical from 
        SEQRES records) then None is entered as the resid in the sequence.

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain or segname) label for residue of 
                          interest
        @type resnames :  string
        @param resnames:  Residue names for all residues
        @type resids   :  integer
        @param resids  :  Resid of the residues of interest
        """

        sequence = self.sequence

        if list(resids):

            sequence[subdiv] = zip(resids, resnames)

        else:

            sequence[subdiv] = [(None, x) for x in resnames]

        return

    def sequence_to_fasta(self, subdiv, model_no=1, missing_lower=False, for_matching=False):
        """
        Create a FASTA format sequence string form the sequence helf for the 
        selected subdiv(ision).

        @type model_no :  int
        @param model_no:  Model number
        @type subdiv  :  string
        @param subdiv :  Subdivision (chain or segname) label for first residue
        @type missing_lower   :  boolean
        @keyword missing_lower:  Should missing residues be presented in
                                 lowercase?
        @type for_matching    :  boolean
        @keyword for_matching :  Only the coordinate residues to be output
                                 (i.e. filter out missing residues - replace
                                 gaps with single '.')
        @rtype        :  string
        @return       :  FASTA formatted sequence of selected subdiv(ision)
        """

        fasta = ''

        # sequence as list - [(resid, resname), ...]
        sequence = self.sequence[subdiv]

        if (missing_lower and subdiv in self.missing_resids[model_no]) or for_matching:

            missing = self.missing_resids[model_no][subdiv].keys()

        else:
            missing = []

        for residue in sequence:

            aa = utils.conv_aa3to1(residue[1])

            if residue[0] in missing:
                if for_matching:
                    fasta += '.'
                else:
                    fasta += aa.lower()
            else:
                fasta += aa

        return fasta

    def is_sequence_complete(self, subdiv, model_no=1):
        """
        Check if the sequence of the selected subdiv(ision) is complete - i.e. 
        all residues have both resid and resname.

        @type subdiv   :  string
        @param subdiv  :  Subdivision (chain or segname) label for first
                          residue
        @type model_no :  integer
        @param model_no:  Number of the model from which the sequence should
                          be derived
        @rtype         :  boolean
        @return        :  Do all residues in selected subdiv(ision) sequence
                          have both resid and resname?
        """

        complete = True

        seqs = self.sequence

        if model_no in seqs and subdiv in seqs[model_no]:
            seq = seqs[model_no][subdiv]
            for res in seq:
                if not res[0] or not res[1]:
                    complete = False
        else:
            complete = False

        return complete

    def seqs_for_completion(self, subdiv, model_no=1):
        """
        Return the sequence elements needed to model gaps - the sequence to
        fill in and the two residues flanking this at either end
        (if available). Makes use of the subdiv_map if the subdiv entered
        is present.

        @type  subdiv:    string
        @param subdiv:    Subdiv(ision) for which to get gap sequences
        @type model_no :  integer
        @param model_no:  Numeric model identifier
        @rtype:           list
        @return:          List containing a list for each gap;
                          [pre_anchor, post_anchor, pre_flank, gap, post_flank] sequences.
        """

        # Use subdiv map for cases where segments created using BIOMT records
        if subdiv in self.subdiv_map:
            chosen_subdiv = self.subdiv_map[subdiv]
        else:
            chosen_subdiv = subdiv

        seqs = self.sequence
        missing_resids = self.missing_resids

        if chosen_subdiv in seqs:
            seq = dict(seqs[chosen_subdiv])
        else:
            seq = {}

        if model_no in missing_resids and chosen_subdiv in missing_resids[model_no]:
            missing_seq = missing_resids[model_no][chosen_subdiv]
        else:
            missing_seq = {}

        seq_gaps = utils.group_sequential_nums(sorted(missing_seq.keys()))

        out_fragments = []

        if chosen_subdiv in self.number_gaps:

            num_gap_starts = self.number_gaps[chosen_subdiv].keys()

        else:

            num_gap_starts = []

        for gap in seq_gaps:

            pre_flank = ''
            post_flank = ''

            start = gap[0]
            anchor = start - 1

            if anchor not in num_gap_starts:

                if start - 1 in seq:

                    pre_anchor = start - 1

                    for i in range(-2, 0):
                        resid = start + i
                        if start + i in seq:
                            pre_flank += utils.conv_aa3to1(seq[resid])
                else:
                    pre_anchor = 0

                end = gap[-1]
                if end + 1 in seq:

                    post_anchor = end + 1

                    for i in range(1, 3):
                        resid = end + i
                        if end + i in seq:
                            post_flank += utils.conv_aa3to1(seq[resid])

                else:
                    post_anchor = 0

                gap_seq = ''

                for resid in gap:
                    gap_seq += utils.conv_aa3to1(seq[resid])

                out_fragments.append(
                    [pre_anchor, post_anchor, pre_flank, gap_seq, post_flank])

        return out_fragments

    def subdiv_renumber_from_one(self, subdiv):

        new_sequence = []

        old_to_new_resid = {}

        new_resid = 0

        for old_resid, resname in self.sequence[subdiv]:

            if resname:

                new_resid += 1

                new_sequence.append((new_resid, resname))

                old_to_new_resid[old_resid] = new_resid

        self.sequence[subdiv] = new_sequence

        self.update_after_renumber(subdiv, old_to_new_resid)

        return old_to_new_resid

    def update_after_renumber(self, subdiv, resid_mapping):

        self.number_gaps[subdiv] = {}

        for model_no in self.missing_resids.keys():

            if subdiv in self.missing_resids[model_no]:

                new_missing = {}

                for old_resid, resname in self.missing_resids[model_no][subdiv].iteritems():

                    if resname:
                        new_missing[resid_mapping[old_resid]] = resname

                self.missing_resids[model_no][subdiv] = new_missing

        for model_no in self.missing_atoms:

            if subdiv in self.missing_atoms[model_no]:

                new_missing = {}

                for old_resid, res_info in self.missing_atoms[model_no][subdiv].iteritems():

                    new_missing[resid_mapping[old_resid]] = res_info

                self.missing_atoms[model_no][subdiv] = new_missing

        for model_no in self.excess_atoms:

            if subdiv in self.excess_atoms[model_no]:

                new_excess = {}

                for old_resid, res_info in self.excess_atoms[model_no][subdiv].iteritems():

                    new_excess[resid_mapping[old_resid]] = res_info

                self.excess_atoms[model_no][subdiv] = new_excess

        for model_no in self.n_missing_hydrogens:

            if subdiv in self.n_missing_hydrogens[model_no]:

                new_missing = {}

                for old_resid, res_info in self.n_missing_hydrogens[model_no][subdiv].iteritems():

                    new_missing[resid_mapping[old_resid]] = res_info

                self.n_missing_hydrogens[model_no][subdiv] = new_missing

        for model_no in self.n_excess_hydrogens:

            if subdiv in self.n_excess_hydrogens[model_no]:

                new_excess = {}

                for old_resid, res_info in self.n_excess_hydrogens[model_no][subdiv].iteritems():

                    new_excess[resid_mapping[old_resid]] = res_info

                self.n_excess_hydrogens[model_no][subdiv] = new_excess

        for ndx in range(len(self.disulphides)):

            if subdiv in self.disulphides[ndx].subdivs():

                if self.disulphides[ndx].bound[0]['subdiv'] == subdiv:

                    old_resid = self.disulphides[ndx].bound[0]['resid']

                    self.disulphides[ndx].bound[0]['resid'] = resid_mapping[old_resid]

                if self.disulphides[ndx].bound[1]['subdiv'] == subdiv:

                    old_resid = self.disulphides[ndx].bound[1]['resid']

                    self.disulphides[ndx].bound[1]['resid'] = resid_mapping[old_resid]

        return

    def purge_subdiv(self, subdiv):
        """
        Remove all references to specified subdiv

        @type  subdiv:    string
        @param subdiv:    Subdiv(ision) to remove
        """

        if subdiv in self.subdivs:

            self.subdivs.remove(subdiv)

        elif subdiv in self.subdiv_map:

            del self.subdiv_map[subdiv]

        if subdiv in self.sequence:
            del self.sequence[subdiv]

        for model_no in self.missing_resids.keys():

            if subdiv in self.missing_resids[model_no]:

                del self.missing_resids[model_no][subdiv]

        if subdiv in self.number_gaps:

            del self.number_gaps[subdiv]

        for model_no in self.missing_atoms:

            if subdiv in self.missing_atoms[model_no]:

                del self.missing_atoms[model_no][subdiv]

        for model_no in self.excess_atoms:

            if subdiv in self.excess_atoms[model_no]:

                del self.excess_atoms[model_no][subdiv]

        for model_no in self.n_missing_hydrogens:

            if subdiv in self.n_missing_hydrogens[model_no]:

                del self.n_missing_hydrogens[model_no][subdiv]

        for model_no in self.n_excess_hydrogens:

            if subdiv in self.n_excess_hydrogens[model_no]:

                del self.n_excess_hydrogens[model_no][subdiv]

        if subdiv in self.heterogens:

            del self.heterogens[subdiv]

        if subdiv in self.altloc:
            del self.altloc[subdiv]

        for ndx in range(len(self.disulphides) - 1, -1, -1):

            if subdiv in self.disulphides[ndx].subdivs():

                del self.disulphides[ndx]

        for biomol_no in self.biomt.keys():

            if subdiv in self.biomt[biomol_no]['subdivs']:

                self.biomt[biomol_no]['subdivs'].remove(subdiv)

        return

# Here be Monsters

    def missing_resids_to_gaps(self, model_no=1):
        # TODO convert list of missing residues to gaps for output/consistency

        gaps = self.gaps
        missing_resids = self.missing_resids

        if model_no in missing_resids.keys():

            for subdiv in missing_resids[model_no].keys():

                resids = missing_resids[model_no][subdiv]
                ranges = []
                for k, g in groupby(enumerate(resids),
                                    lambda i_x: i_x[0] - i_x[1]):
                    group = map(itemgetter(1), g)
                    ranges.append((group[0], group[-1]))
                gaps[subdiv] = ranges

        return

    def missing_gaps_to_resids(self, model_no=1):
        # TODO convert list of gaps to missing residues output/consistency

        missing = {model_no: {}}

        for subdiv in self.gaps.keys():
            ranges = self.gaps[subdiv]
            for gap in ranges:
                missing[model_no][subdiv] = {}
                for resid in range(gap[0], gap[-1] + 1):
                    missing[model_no][subdiv][resid] = 'UNK'

        self.missing_resids = missing

        return

