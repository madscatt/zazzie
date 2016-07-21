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

import copy
import logging

from . import pdbscan_utils as utils
from .exceptions import NoCommonChainError

# TODO: Check robustness to missing chains in header or coordinates

class Chain_Info_Reconciler():
    """
    Reconcile data_struct.Info objects generated from scanning PDB header and
    coordinates. Note: As headers only refer to chains no segnames are 
    involved.
    Information derived from coordinates is updated to reflect additional 
    data taken from the header.
    """

    def __init__(self, head_info, coor_info):
        """
        Initialize with the two Info object to be compared and variables to 
        store results.
        
        @type head_info :  data_struct.Info
        @param head_info:  Chain information from PDB header
        @type coor_info :  data_struct.Info
        @param coor_info:  Chain information from PDB coordinates
        """
        
        self.head_info = head_info
        self.coor_info = coor_info
        
        self.coord_updated = False
        
        self.issues = {}
        
        self.logger = logging.getLogger(__name__)
    
        self.reconciled = {}

        self.common_chains = []
        self.chains_not_in_coor = []
        self.chains_not_in_head = []

        pass

    def reconcile_chain_info(self, model_no=1):
        """
        Run the reconciliation of the header and coordinate information. Check 
        that we have common chains and if so compare the sequences (including 
        gaps for missing resids and numerical gaps) and HETATM contents.
        Note: sequences only available for protein/nucleic chains, so HET check 
        is not the same thing as sequence check in this context.
        
        @type  model_no: integer
        @param model_no: number of the model whose data we are interested in.
        """        
        
        self.get_common_chains()
        
        self.setup_issue_tracker()
        
        if self.common_chains:
            self.compare_sequences(model_no=model_no)            
            self.check_hets()
        else:
            txt = 'No matching chains found in header/coordinates.'
            self.logger.error(txt)
            raise NoCommonChainError(txt)

        self.check_chain_reconciliation()
               
        return
        
    def check_chain_reconciliation(self):
        """
        Summarize which chains have been successfully reconciled between header 
        and coordinates in self.reconciled dictionary (key = chain, value = 
        boolean).
        
        """
        
        for chain, report in self.issues.items():
            
            if (not report['seq_mismatch'] and not report['het_mismatch'] and not report['missing']):
                self.reconciled[chain] = True
            else:
                self.reconciled[chain] = False
                
        return

    def get_reconciled_chains(self):
        """
        List the chains which have been successfully reconciled

        @rtype  :  list
        @return :  List of chain labels which have been reconciled
        """
        
        all_good = [x for x,y in self.reconciled.items() if y]        
        
        return sorted(all_good)
        
    def get_common_chains(self):
        """
        Determine the chains which are common to the header and coordinates 
        (hopefully all of them). Put there in self.common_chains. Lists of the 
        chains found in only one of the other also stored (in 
        chains_not_in_head & chains_not_in_coor respectively).
        """

        head_info = self.head_info

        coor_chains = set(self.coor_info.subdivs)
        head_chains = set(head_info.subdivs)
        head_seq_chains = set(head_info.sequence.keys() + head_info.heterogens.keys())

        if head_chains.symmetric_difference(head_seq_chains):
            head_chains = head_chains.union(head_seq_chains)

        self.chains_not_in_head = list(coor_chains.difference(head_chains))    
        self.chains_not_in_coor = list(head_chains.difference(coor_chains))    
        self.common_chains = list(head_chains.intersection(coor_chains))

        return

    def setup_issue_tracker(self):
        """
        Create dictinaries for each chain to hold outcome of scans.
        """
        
        issues = self.issues
        
        for chain in self.common_chains:
            issues[chain] = {'seq_mismatch' : False,
                             'het_mismatch' : [],
                             'missing'  : []
                            }
                            
        for chain in self.chains_not_in_head + self.chains_not_in_coor:
            issues[chain] = {'seq_mismatch' : True,
                             'het_mismatch' : [],
                             'missing'  : []
                            }
        
        return

    def compare_sequences(self, model_no = 1):
        """
        Compare sequences        
        
        @type  model_no: integer
        @param model_no: number of the model whose data we are interested in.        
        """

        logger = self.logger

        # Header info provides resnames for missing residues in coordinates
        # These are added to sequence and missing_resids
        self.update_coord_missing(model_no=model_no)

        
        for chain in self.common_chains:
            
            agree, removed_resids = self.compare_sequence(chain, model_no=model_no)
                        
            if agree and removed_resids:
                self.update_number_gaps(chain, removed_resids)
            elif not agree:
                self.issues[chain]['seq_mismatch'] = True
                logger.warning('Failure to match header and coordinate sequences in chain {0:s}'.format(chain))
                
        return

    def check_missing(self, model_no=1):
        """
        Check if the same residues are listed as missing in the header and 
        coordinates for each chain.
        
        @type  model_no: integer
        @param model_no: number of the model whose data we are interested in. 
        """
        
        head_info = self.head_info
        coor_info = self.coor_info        
   
        # missing_resids = {model_no:{subdiv: {resid:resname}, ..}}
        if model_no in head_info.missing_resids:
            head_missing = head_info.missing_resids[model_no]
        else:
            head_missing = {}
            
        if model_no in coor_info.missing_resids:
            coord_missing = coor_info.missing_resids[model_no]
        else:
            coord_missing = {}
        
        agree = {}
    
        head_chains = head_missing.keys()
        coord_chains = coord_missing.keys()
        
        if head_chains != coord_chains:
            
            set_head_chains = set(head_chains)
            set_coord_chains = set(coord_chains)
            
            non_common_chains = set_head_chains.symmetric_difference(set_coord_chains)
            common_chains = set_head_chains.intersection(set_coord_chains)
            
            for chain in non_common_chains:
                agree[chain] = False                            
        else:
            
            common_chains = head_chains
                
        for chain in common_chains:
            
            if head_missing[chain].keys() == coord_missing[chain].keys():
                agree[chain] = True
            else:
                agree[chain] = False
        
        return agree

    def update_coord_missing(self, model_no=1):
        """
        Reconcile the missing residues from header and coordinates. Where 
        possible add resname information to the coordinate missing residues 
        and sequence.
        
        @type  model_no: integer
        @param model_no: number of the model whose data we are interested in.
        
        @todo:  Make sure this is robust when chains between header and 
                coordinates don't match
        """
        
        coor_info = self.coor_info
        head_info = self.head_info
                
        # missing_resids = {model_no:{subdiv: {resid:resname, ..}}
        if model_no in head_info.missing_resids:
            head_missing = head_info.missing_resids[model_no]
        else:
            head_missing = {}
            
        if model_no in coor_info.missing_resids:
            coord_missing = coor_info.missing_resids[model_no]
        else:
            coord_missing = {}
        
        # NOTE info sequence format is:
        # key = subdiv, value = list of (resid,resname) - 3 character resnames
    
        # Check that the same residues are reported missing from coordinates and
        # header scans
        agree = self.check_missing(model_no=model_no)
        
        for chain, residues_same in agree.items():
            
    
            if residues_same:
            # If the same residues are missing then we copy the resnames from the 
            # header into the coordinate base information            
            
                coord_missing[chain] = head_missing[chain]
                
                self.resname_missing_resids_coor(head_missing[chain].keys(), head_missing[chain], coor_info.sequence[chain])
                
                for resid, resname in head_missing[chain].iteritems():
                    
                    coord_seq = coor_info.sequence[chain]
                    ndx = coord_seq.index((resid,''))
                    coord_seq[ndx] = (resid,resname)
                    self.coord_updated = True
                
            else:
        
                if chain in head_missing:
                    head_miss_resids = set(head_missing[chain].keys())
                else:
                    head_miss_resids = set([])
                    
                if chain in coord_missing:
                    coord_miss_resids = set(coord_missing[chain].keys())
                else:
                    coord_miss_resids = set([])
                
                only_missing_head = head_miss_resids.difference(coord_miss_resids)
                only_missing_coord = coord_miss_resids.difference(head_miss_resids)
                missing_both = head_miss_resids.intersection(coord_miss_resids)

                coord_seq = coor_info.sequence[chain]

                for resid in sorted(missing_both):

                    resname = head_missing[chain][resid]                                      
                    coord_missing[chain][resid] = resname
                                        
                    ndx = coord_seq.index((resid,''))
                    coord_seq[ndx] = (resid,resname)

                if missing_both:
                    self.coord_updated = True                                        
    
                if only_missing_head:                 
                    
                    first_coord_resid = coord_seq[0][0]
                    last_coord_resid = coord_seq[-1][0]
                    
                    missing_after = [x for x in sorted(only_missing_head) if x > last_coord_resid]
                    missing_before = [x for x in sorted(only_missing_head, reverse = True) if x < first_coord_resid]
                    missing_mid = [x for x in sorted(only_missing_head) if first_coord_resid < x < last_coord_resid]
                    
                    # Add header residues listed as missing after the final 
                    # coordinate resid to end of coordinate sequence (and to 
                    # coordinate missing list)
                    for resid in missing_after:
                        
                        resname = head_missing[chain][resid]
                        
                        coor_info.add_missing_resid(chain, resid, resname, model_no=model_no)
                        coord_seq.append((resid,resname))

                    # Add header residues listed as missing before the first 
                    # coordinate resid to start of coordinate sequence (and to 
                    # coordinate missing list)
                    for resid in missing_before:
                        
                        resname = head_missing[chain][resid]
                        
                        coor_info.add_missing_resid(chain, resid, resname, model_no=model_no)
                        coord_seq.insert(0,(resid,resname))
                        
                    if missing_after or missing_before:
                        self.coord_updated = True 

                    # Reported missing in the header but not coordinates is a 
                    # bit odd but not something we need to fix - log issue                        
                    for resid in missing_mid:
                    
                        resname = head_missing[chain][resid]
                        txt = 'Header claims resid {0:d} is missing in chain {1:s} but residue {2:s} found in coordinates'.format(resid, chain, resname)
                        self.logger.warning(txt)
                        
                for resid in sorted(only_missing_coord):
                                               
                    self.issues[chain]['missing'].append(resid)
               
        return

    def check_hets(self):
        """
        Check if the heterogens listed in teh coordinate and header information
        agree (residue by residue).
        """
        
        # heterogens = {subdiv: {resid:resname}}
        head_hets = self.head_info.heterogens        
        coor_hets = self.coor_info.heterogens
        logger = self.logger

        for chain, hets in head_hets.iteritems():
    
            if chain in coor_hets:
                
                for resid, resname in hets.iteritems():
                    
                    if (resid in coor_hets[chain]) and (resname != coor_hets[chain][resid]):
                        
                        coor_resname = coor_hets[chain][resid]
                                                    
                        issue = {
                                'resid': resid,
                                'head_resname': resname,
                                'coor_resname': coor_resname
                                }
                                
                        self.issues[chain]['het_mismatch'].append(issue)
                        
                        txt = 'Residue name mismatch between header ({0:s}) and coordinates ({1:s}) for chain {2:s} resid {3:d}.'
                        logger.warning(txt.format(resname, coor_resname, chain, resid))
               
        return

    def remove_unknown_from_sequence(self, sequence):
        """
        Remove residues from sequence if the resname is unknown. Return both 
        the updated sequence and the resids removed.
        
        Note: sequence is a list of tuples (resid, resname)

        @type sequence :  list
        @param sequence:  List of tuples (resid, resname)
        @rtype         :  list, list
        @return        :  
                          1.  List of tuples (resid, resname) with blank 
                              resnames removed.
                          2.  List of resids (integers)
        """
        
        edited_sequence = list(sequence)

        gaps = []

        blanks = [x for x in edited_sequence if x[1] == '']
           
        for blank in blanks:
                
            resid = blank[0]                
            gaps.append(resid)
            edited_sequence.remove(blank)

        return edited_sequence, gaps

    def compare_sequence(self, chain, model_no = 1):
        """
        Remove residues from sequence if the resname is unknown. Return both 
        the updated sequence and the resids removed. In this comparison all 
        HET atoms found in the sequence are treated as the same (converted to 
        X in fasta format in which the comparison is conducted).
        
        Note: sequence is a list of tuples (resid, resname)

        @type chain :  string
        @param chain:  Chain identifier
        @rtype      :  boolean, list
        @return     :  
                       1.  Do the sequences for the selected chains agree in 
                           coordinate and header reports
                       2.  List of missing resids removed from the coordinate 
                           sequence (integers).
        """
        
        agree = False
        removed_resids = []

        head_info = self.head_info
        coor_info = self.coor_info

        if chain in head_info.sequence:
            
            # Get target header sequence in fasta format     
            head_fasta = self.head_info.sequence_to_fasta(chain)
            
            # As there is no way of knowing if missing resids in coordinates 
            # are real or just numbering gaps prior to applying a sequence 
            # there may be rogue missing residues in sequence, we'll try 
            # removing them
            coor_seq = coor_info.sequence[chain]
    
            # keep old sequence in case of non-match
            bak_coor_seq = list(coor_seq)
    
            edited_seq, removed_resids = self.remove_unknown_from_sequence(coor_seq)
            coor_info.sequence[chain] = edited_seq
    
            # Get coordinate sequence in fasta format for comparison
            coor_fasta = coor_info.sequence_to_fasta(chain)
  
            # If they match then we can edit missing residues and flag agreement,
            # otherwise revert to saved sequence.
            if coor_fasta == head_fasta:
                
                agree = True

                if (model_no in coor_info.missing_resids and chain in coor_info.missing_resids[model_no]):
                
                    missing = coor_info.missing_resids[model_no][chain]
                    for resid in removed_resids:
                        del missing[resid]
    
            else:
                coor_info.sequence[chain] = bak_coor_seq

        return agree, removed_resids

    def update_number_gaps(self, chain, removed_resids):
        """
        Residues removed in sequence comparison for each chain represent number 
        gaps if the chain comparison is successful. Convert the list
        removed_resids into the end points either side of runs of resids which 
        do not exist in the sequence and add to the coordinate information 
        number_gaps.

        @type chain          :  string
        @param chain         :  Chain identifier
        @type removed_resids :  list
        @param removed_resids:  List of resids removed from the coordinate 
                                sequence which we now believe to be number gaps
        """
        
        for resid in removed_resids:
            self.issues[chain]['missing'].remove(resid)
            
        grouped_gaps = utils.group_ranges(removed_resids)
        
        for gap in grouped_gaps:
            # number_gap records the residues either side of the gap
            self.coor_info.add_number_gap(chain, gap[0] - 1 , gap[1] + 1)
            
        return

    def copy_biomt(self):
        """
        Copy biological unit information from the header to the coordinate 
        information.
        """
        
        self.coor_info.biomt = copy.deepcopy(self.head_info.biomt)
        
        return


