# -*- coding: utf-8 -*-
"""
header_reader

This module will read and parse PDB header records

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

import logging
import itertools
import numpy as np

from . import data_struct as data_struct
from .pdb_record_reader import (parse_line, process_runon_line, rec_schemas)
from .exceptions import (NoHeaderReadError, IncompleteBiomtError)


class ReferenceInfo:
    """
    Contains the information extracted from the PDB header about the structure.
    """

    def __init__(self):
        """
        Initialise all of the variables to store header information.
        Naming is either self explanatory (disulphides) or based on the PDB 
        records the information is derived from.
        """

        self.title = ''
        self.citation = {}
        self.compnd = {}
        self.formuls = {}
        self.hets = {}
        self.metrics = {}


class PdbHeader:
    """
    Reads and then stores PDB header information.
    """

    def __init__(self, pdbfile=None, parse=True, sasmol=None, text=None):
        """
        Initialise all of the variables to store header information, warnings 
        and errors.
        If provided read in a pdbfile and if the parse flag is True then 
        process to populate the variables.

        @type  pdbfile:  string
        @param pdbfile:  PDB filename
        @type  parse:    boolean
        @param pdbfile:  PDB filename        
        """

        # All of the record types for which pdb_record_reader has schemas
        self.head_types = rec_schemas.keys()

        self.coord_types = ['ATOM',
                            'HETATM',
                            'TER',
                            'END',
                            'MODEL',
                            'ENDMDL',
                            'CONECT',
                            'ANISOU']

        self.set_blank_values()

        self.logger = logging.getLogger(__name__)

        if pdbfile:
            self.read_pdb(pdbfile)
        elif sasmol:
            self.process_sasmol_header(sasmol)
        elif text:
            self.process_header_text(text)

        if (pdbfile or sasmol or text) and parse:
            try:
                self.parse_header()
            except NoHeaderReadError:
                pass
            except:
                if pdbfile:
                    raise ValueError(
                        'Unable to parse PDB header: {0:s}'.format(pdbfile))
                else:
                    raise ValueError(
                        'Unable to parse PDB header: {0:s}'.format(sasmol.pdbname))

        return

    def set_blank_values(self):
        """
        Blank values to be filled by processed header
        """

        self.pdb_recs = {}
        for rec_type in self.head_types:
            self.pdb_recs[rec_type] = []

        self.reference_info = ReferenceInfo()

        self.chain_info = data_struct.Info(scan_type='chain')

        self.read_valid_header = False
        self.has_seq_info = False

        return

    def process_header_line(self, line):
        """
        Process single line from a PDB header. Processed data is appended to 
        the appropriate list in self.pdb_recs dictionary (key is the record 
        type).

        @type  line:  string
        @param line:  Line from a PDB header
        """

        err = None
        rec_type = line[0:6].strip()

        if rec_type not in self.coord_types:

            # Adjust records where naming discrepancy has occured
            if rec_type[0:5] in ['MTRIX', 'ORIGX', 'SCALE']:
                rec_type = rec_type[0:5] + 'n'

            try:
                # Parse line in to variables using schema
                vals, err = parse_line(line, rec_schemas[rec_type])
                self.pdb_recs[rec_type].append(vals)

            except:
                # If no schema available or parsing fails use NONSTD
                # schema which just reads in the line as a text field
                vals, err = parse_line(line, rec_schemas['NONSTD'])
                self.pdb_recs['NONSTD'].append(vals)

        if err:
            err_txt = '\n'.join(err)
            raise IOError(
                'Unable to parse PDB header line:\n{0:s}\nError:\n{1:s}'.format(line, err_txt))

        return

    def process_header_text(self, header_txt):
        """
        Parse list of header lines lines and place by record type in
        self.pdb_recs dictionary

        @type  header_txt:  list
        @param header_txt:  List of stings containing PDB records
        """

        self.set_blank_values()

        try:

            for line in header_txt:
                self.process_header_line(line)
        except Exception as err:
            py_err = str(err)
            raise IOError(
                'Unable to read header line from SasMol object: {0:s}\n{1:s}'.format(line, py_err))

        return

    def process_sasmol_header(self, sasmol):
        """
        Read PDB file, parse lines and place by record type in self.pdb_recs
        dictionary

        @type  sasmol:  sasmol.SasMol
        @param sasmol:  SasMol containing data read in from a PDB
        """

        header_txt = sasmol.header()

        self.process_header_text(header_txt)

        return

    def read_pdb(self, pdbfile):
        """
        Read PDB file, parse lines and place by card type in self.pdb_recs
        dictionary

        @type  pdbfile:  string
        @param pdbfile:  PDB filename
        """

        self.set_blank_values()

        # line count used for reporting purposes
        line_no = 1

        try:

            with open(pdbfile, 'r') as f:

                for line in f:
                    self.process_header_line(line)
                    line_no += 1

        except IOError as err:
            raise IOError('Unable to read PDB: {0:s}, line {1:d}\n{2:s}'.format(
                pdbfile, line_no, err))

        return

    def is_obsolete(self):
        """
        Is the PDB obsolete according to the header?

        @rtype:   boolean
        @return:  Flag to say PDB is noted as being obsolete in header
        """

        obsolete = False

        if 'OBSLTE' in self.pdb_recs:
            if self.pdb_recs['OBSLTE']:
                obsolete = True

        return obsolete

    def is_split(self):
        """
        Is the PDB part of a structure split across multiple files?

        @rtype:   boolean
        @return:  Flag to say PDB is noted as being part of a split 
                  structure in the header
        """

        split = False

        if 'SPLIT' in self.pdb_recs:
            if self.pdb_recs['SPLIT']:
                split = True

        return split

    def parse_header(self):
        """
        Parse all header records in self.pdb_recs.
        Populate self.reference_info and the following attributes of 
        self.chain_info: sequence, missing_resids, missing_atoms, heterogens, 
        disulphides and n_models.
        """

        pdb_recs = self.pdb_recs

        has_header_info = sum(len(v) for v in pdb_recs.itervalues())

        if has_header_info:

            if pdb_recs['NUMMDL']:
                self.logger.info('Multiple models (' +
                                 str(pdb_recs['NUMMDL'][0]['no_models']) +
                                 ') detected in the file.')
                self.chain_info.n_models = int(
                    pdb_recs['NUMMDL'][0]['no_models'])

            if pdb_recs['TITLE']:
                self.reference_info.title = process_runon_line(
                    pdb_recs['TITLE'], 'text')

            # SEQRES records contain full sequence information
            self.process_seqres()

            # Remarks contain quality metrics, BIOMT and missing residue/atom
            # information
            self.parse_remarks()

            # Citation information
            self.parse_jrnl()

            if pdb_recs['SSBOND']:
                self.process_disulphides()

            # Information about non-standard residues
            self.process_header_het()

            self.process_compnd()

            self.read_valid_header = True
            self.has_seq_info = self.check_header_seq()

        else:
            raise NoHeaderReadError("No header has been read into object")

        return

    def check_header_seq(self):
        """
        Check to see if parsed header information included a sequence for one 
        or more chains

        @rtype  : boolean
        @return : Flag if sequence information has been read in

        @todo: Need to ensure this is an adequate check, could probably check 
               if missing_resids is consistent with seqres
        """

        has_seq = False

        if self.chain_info.sequence:
            has_seq = True

        return has_seq

    def parse_remarks(self):
        """
        Parse all REMARK records.
        Extracts missing residue (self.chain_info.missing_resids)
        missing atom (self.chain_info.missing_atoms), 
        BIOMT (self.chain_info.biomt) and experimental quality metrics 
        (self.reference_info.metrics).

        """

        if self.pdb_recs['REMARK']:
            self.get_quality_metrics()
            self.process_biomolecule()
            self.process_missing_res()
            self.process_missing_atoms()

        else:
            self.logger.info('No REMARK lines found in header:')
            self.logger.info(
                'BIOMT and missing residues cannot be evaluated')
        return

    def process_missing_res(self):
        """
        Parse REMARK 465 records from a PDB to obtain missing residues.
        Populate self.chain_info.missing_resids with a dictionary of the form:

        {model_no:{chain: [{resid:resname}, ..]}}
        """

        # missing_resids = self.chain_info.missing_resids
        chain_info = self.chain_info

        remarks465 = [x for x in self.pdb_recs['REMARK'] if x['num'] == 465]

        # Extract missing residue data from text field of REMARK 465 records
        missing_rec = self._remark465_missing_residues(remarks465)

        # Create a dictionay for the form:
        # missing_resids = {model_no:{chain: {resid:resname}, ..}}
        # from the parsed REMARK lines

        for model, grpd in itertools.groupby(
                missing_rec, key=lambda x: x['model']):

            for chain, residues in itertools.groupby(
                    grpd, key=lambda y: y['chain']):
                residues = list(residues)
                resids = [z['resid'] for z in residues]
                resnames = [z['resname'] for z in residues]

                chain_info.add_missing_resids(chain, resids, resnames, model)

                n_missing = chain_info.no_missing_resids(subdiv=chain,
                                                         model_no=model)

                self.logger.info(
                    str(n_missing) +
                    ' missing residues in chain ' +
                    chain)

        return

    def _remark465_missing_residues(self, remarks465):
        """
        Extract information from PDB REMARK 465 records for further processing 
        using the schema:

        missing_schema = (
        ('model', 0, 3, None),       # 0,  model
        ('resname', 4, 7, None),     # 1,  residue name
        ('chain', 8, None, None),    # 2,  chain ID
        ('resid', 9, 16, int),       # 3,  residue number
        ('insert', 16, None, None),  # 4,  insertion code
        )

        @type remarks465 :  list
        @param remarks465:  List containing PDB REMARK records after basic 
                            parsing
        @rtype           :  list
        @return          :  List of dictionaries extracting data from the 
                            text field of the original PDB REMARK 465 records
        """

        # Schema is for parsing the text section extracted from the original
        # REMARK record not the original PDB record
        missing_schema = (
            ('model', 0, 3, None),  # 0,  model
            ('resname', 4, 7, None),  # 1,  residue name
            ('chain', 8, None, None),  # 2,  chain ID
            ('resid', 9, 16, int),  # 3,  residue number
            ('insert', 16, None, None),  # 4,  insertion code
        )

        missing_rec = []

        # Useful records are preceeded by a header but also possibly other
        # nonsense - only start reading after get to usable lines
        start = False

        for remark in remarks465:

            if start:

                rec, err = parse_line(remark['text'], missing_schema)

                if err:
                    self.logger.warning(
                        'Possible malformed missing residue remark: ' +
                        remark['text'])

                else:
                    try:
                        rec['model'] = int(rec['model'])
                    except:
                        rec['model'] = 1

                    missing_rec.append(rec)

            elif remark['text'].startswith('  MODELS'):
                self.logger.warning(
                    'Missing report for NMR' +
                    remark['text'].lower() +
                    '\n')

            # Check for header of the residue list
            # Note: M column not used as not in new NMR table format
            elif 'RES C SSSEQI' in remark['text']:
                start = True

        return missing_rec

    def process_missing_atoms(self):
        """
        Parse REMARK 470 records from a PDB to obtain missing atoms.
        Populate self.chain_info.missing_atoms with a dictionary of the form:

        missing_atms = {model_no:{chain: {resid: {'atoms':[atm1, atm2, ...],'resname': resname}}}}

        """

        missing_atoms = self.chain_info.missing_atoms
        remarks470 = [x for x in self.pdb_recs['REMARK'] if x['num'] == 470]

        missing_rec = self._remark470_missing_atoms(remarks470)

        # Create a dictionay for the form:
        # missing_atoms = {model_no:{chain: {resid:
        #                               {'atoms':[atm1, atm2, ...],
        #                                'resname': resname}
        #                }}}
        # from the parsed REMARK lines

        for model, grpd in itertools.groupby(
                missing_rec, key=lambda z: z['model']):

            missing_atoms[model] = {}

            for chain, resids in itertools.groupby(
                    grpd, key=lambda y: y['chain']):

                missing_atoms[model][chain] = {}

                for res in resids:
                    missing_atoms[model][chain][res['resid']] = {
                        'resname': res['resname'],
                        'atoms': res['atoms']
                    }

        return

    def _remark470_missing_atoms(self, remarks470):
        """
        Extract information from PDB REMARK 470 records for further processing 
        using schema:

        missing_schema = (
        ('model', 0, 3, None),       # 0,  model
        ('resname', 4, 7, None),     # 1,  residue name
        ('chain', 8, None, None),    # 2,  chain ID
        ('resid', 9, 16, int),       # 3,  residue number
        ('insert', 16, 16, None),    # 4,  insertion code
        ('atoms', 17, 80, None),     # 5,  atom names
        )

        @type remarks470 :  list
        @param remarks470:  List containing PDB REMARK records after basic 
                            parsing
        @rtype           :  list
        @return          :  List of dictionaries extracting data from the 
                            text field of the original PDB REMARK 470 records
        """

        # Schema is for parsing the text section extracted from the original REMARK
        # record not the original PDB record
        missing_schema = (
            ('model', 0, 3, None),  # 0,  model
            ('resname', 4, 7, None),  # 1,  residue name
            ('chain', 8, None, None),  # 2,  chain ID
            ('resid', 9, 16, int),  # 3,  residue number
            ('insert', 16, 16, None),  # 4,  insertion code
            ('atoms', 17, 80, None),  # 5,  atom names
        )

        missing_rec = []
        # Useful records are preceeded by a header but also possibly other
        # nonsense
        start = False

        for remark in remarks470:

            if start:
                rec, err = parse_line(remark['text'], missing_schema)
                rec['atoms'] = remark['text'][17:-1].split()

                if err:
                    self.logger.warning(
                        'Possible malformed missing atom remark: ' +
                        remark['text'])
                else:
                    try:
                        rec['model'] = int(rec['model'])
                    except:
                        rec['model'] = 1

                missing_rec.append(rec)

            elif 'RES CSSEQI' in remark['text']:
                start = True

        return missing_rec

    def process_seqres(self):
        """
        Parse SEQRES records to provide a sequence for each chain
        Populate self.chain_info.sequence as a dictionary of the form:

        {chain: [resname, ..]}

        resnames are three letter codes taken from SEQRES records
        """

        seqres_recs = self.pdb_recs['SEQRES']

        # Group SEQRES records by chain
        for chain, grpd in itertools.groupby(
                seqres_recs, key=lambda x: x['chain']):

            chain_seq = []

            for entry in grpd:
                chain_seq += entry['resnames']

            self.chain_info.add_subdiv_sequence(chain, chain_seq)

        return

    def get_quality_metrics(self):
        """
        Parse REMARK 2 and 3 records from PDB to obtain R values and resolution
        Store in self.reference_info.metrics as a dictionary containing:
        'resolution', 'r' and 'r free'
        """

        remarks = self.pdb_recs['REMARK']

        metrics = {}

        for remark in remarks:

            text = remark['text']

            if remark['num'] == 2:

                if len(text.split()) != 0:
                    try:
                        metrics['resolution'] = float(text[15:19])
                    except:
                        metrics['resolution'] = None

            if remark['num'] == 3:

                if 'R VALUE            (WORKING SET) :' in remark['text']:

                    try:
                        metrics['r'] = float(text.split()[5])
                    except:
                        metrics['r'] = None

                elif 'FREE R VALUE                     :' in remark['text']:

                    try:
                        metrics['r_free'] = float(text.split()[4])
                    except:
                        metrics['r_free'] = None

        self.reference_info.metrics = metrics

        return

    def process_biomolecule(self):
        """
        Parse REMARK 300 and 350 records from PDB to obtain biological unit
        specification.
        Populate self.chain_info.biomt with the following format:

        biomt[biomol_no] = {
        'subdivs' : [],
        'auth_bio_unit' : '',
        'soft_bio_unit' : '',
        'rot' : [],
        'trans' : []
        }

        rot = list of np.array((3,3))
        trans = list of np.array(3)
        """

        biomt = self.chain_info.biomt

        # REMARK 300 section is a free text description of any biomolecules
        # We extract a list of the numberical biomolecule labels described
        biomol_300 = self._parse_biomol_300()

        # REMARK 350 section contains the transforms for each biomolecule and
        # information on how the description was arrived at.
        # Transformation read in from BIOMT records
        self._parse_biomol_350()

        if len(biomol_300) != len(biomt.keys()):
            self.logger.warning(
                'No. biomolecules suggested in REMARK 300 and supplied in REMARK 350 records are inconsistent!')

        if biomt:
            self.logger.warning(
                'BIOMT present - PDB may not represent the biological unit')

        return

    def _parse_biomol_300(self):
        """
        Parse REMARK 300 records from PDB to obtain list of biological unit 
        labels which we expect to find specified in the REMARK 350 records.
        """

        biomol_300 = []
        remarks300 = [x for x in self.pdb_recs['REMARK'] if x['num'] == 300]

        # REMARK 300 section is a free text description of any biomolecules
        # described in REMARK 350 records
        # We just want the biomolecule ID numbers to be described
        for remark in remarks300:

            if remark['text'].startswith('BIOMOLECULE:'):

                for biomol_no in remark['text'][13:].split(','):
                    biomol_300.append(int(biomol_no))

        return biomol_300

    def _parse_biomol_350(self):
        """
        Parse REMARK 350 records from PDB to obtain biological unit 
        specification.

        Populate self.chain_info.biomt with the following format:

        biomt[biomol_no] = {
        'subdivs' : [],
        'auth_bio_unit' : '',
        'soft_bio_unit' : '',
        'rot' : [],
        'trans' : []
        }

        rot = list of np.array((3,3))
        trans = list of np.array(3)
        """

        logger = self.logger 

        chain_info = self.chain_info
        biomt = self.chain_info.biomt

        remarks350 = [x for x in self.pdb_recs['REMARK'] if x['num'] == 350]

        in_biomol = False

        # REMARK 350 records contain data on origin of the biomolecule
        # description along with the transformations needed to create it from
        # the ATOM/HETATM coordinates
        for remark in remarks350:

            content = remark['text']

            if content.startswith('BIOMOLECULE:'):

                bm_no = int(content[13:])
                chain_info.create_biomol(bm_no, subdivs=[])

                in_biomol = True
                last_seen_row = None

            elif in_biomol and content.strip():

                if content.startswith('AUTHOR DETERMINED'):
                    biomt[bm_no]['auth_bio_unit'] = content.split(':')[
                        1].strip()

                elif content.startswith('SOFTWARE DETERMINED'):
                    biomt[bm_no]['soft_bio_unit'] = content.split(':')[
                        1].strip()

                elif content[0:31] in ['APPLY THE FOLLOWING TO CHAINS: ',
                                       '                   AND CHAINS: ']:

                	content = content.split(":")[1]
                        content = content.replace(' ', '')
                        for chain in content.split(','):
                            if not chain.isspace() and len(chain) > 0:
                            	biomt[bm_no]['subdivs'].append(chain.strip())
                                logging.debug('biomt[bm_no] ' +
                                          biomt[bm_no]['subdivs'][-1])
                            else:
                                logging.debug('chain == empty_space = ' +
                                                chain)

                elif content.startswith('  BIOMT'):

                    # If we have not yet read any rows from a BIOMT matrix line
                    # initialize r(otation) and t(ranslation) arrays
                    if last_seen_row is None:
                        r = np.identity(3)
                        t = np.zeros(3)

                    # Read BIOMT record and parse to add information to r & t
                    last_seen_row = self._add_biomt_row(
                        content, last_seen_row, r, t)

                    # If we have finished reading a transformation add it to
                    # the chain_info object
                    if last_seen_row is None:
                        chain_info.add_biomt(bm_no, r, t)

        return

    def _add_biomt_row(self, biomt_text, last_seen, rot, trans):
        """
        Parse the BIOMT information from a PDB REMARK 350 record and add to 
        rotation matrix and translation vector.
        Example record:
        BIOMT1   1  1.000000  0.000000  0.000000        0.00000

        Column 0 = Record name, number indicated transform row
        Column 1 = Transform number within this biomolecule
        Columns 2,3,4 = rotation matrix row
        Column 5 = translation vector row

        @type biomt_text :  string
        @param biomt_text:  Text containing content of the REMARK 350 
                            containing a BIOMT line (record type and number 
                            pre-stripped)
        @type last_seen  :  integer
        @param last_seen :  Last row of the rotation/translation information 
                            read or None if this should be the first line
        """

        cols = biomt_text.split()

        # Convert PDB row number to python zero based index for rot and trans
        ndx = int(cols[0][-1]) - 1

        # Should be first row or the one after the last read row
        if (last_seen is None and ndx == 0) or (last_seen == ndx - 1):

            rot[ndx] = np.array([float(cols[2]),
                                 float(cols[3]),
                                 float(cols[4])])

            trans[ndx] = float(cols[5])

            if ndx == 2:
                # Only three rows - set last_seen for a new matric on next read
                last_seen = None
            else:
                last_seen = ndx

        else:
            raise IncompleteBiomtError('Incomplete BIOMT matrix encountered')

        return last_seen

    def parse_jrnl(self):
        """
        Parse an JRNL record line from a PDB. Extract the authors, title, 
        citation data, Pubmed IS and DOI for the primary citation if present.
        Populates: self.reference_info.citation as a dictionary with keys based 
        on sub-record names:

        {'authors': '', 'title': '', 'citation': '', 'pmid': '', 'doi': ''}
        """

        jrnl_recs = self.pdb_recs['JRNL']

        if jrnl_recs:

            ref_data = {
                'authors': '',
                'title': '',
                'citation': '',
                'pmid': '',
                'doi': '',
            }

            for rec in jrnl_recs:

                rec_type = rec['text'][0:4].strip()
                rec_val = rec['text'][7:]

                if rec_type == 'AUTH':

                    if rec['text'][5] != ' ':
                        ref_data['authors'] += ',' + rec_val
                    else:
                        ref_data['authors'] = rec_val

                elif rec_type == 'TITL':

                    if rec['text'][5] != ' ':
                        ref_data['title'] += ' ' + rec_val
                    else:
                        ref_data['title'] = rec_val

                elif rec_type == 'REF':
                    ref_data['citation'] += ' '.join(rec_val.split())

                elif rec_type == 'PMID':
                    ref_data['pmid'] = rec_val.strip()

                elif rec_type == 'DOI':
                    ref_data['doi'] = rec_val.strip()

            self.reference_info.citation = ref_data

        return

    def process_disulphides(self):
        """
        Extract disulphide bond information from SSBOND records
        Populates self.chain_info.disulphides as a list of 
        data_struct.Disulphide objects.
        """

        ssbond_recs = self.pdb_recs['SSBOND']

        if ssbond_recs:

            for rec in ssbond_recs:
                bond = data_struct.Disulphide(rec['chain1'], rec['resid1'],
                                              rec['chain2'], rec['resid2'],
                                              'chain')

                self.chain_info.disulphides.append(bond)

        return

    def process_compnd(self):
        """
        Parse COMPND records to get description of the macromolecular contents of
        the PDB.
        Populates self.reference_info.compnd dictionary with key/value pairs:

        {mol_no: 'chains': [], 'name': '', 'fragment': '', 'type': 'protein'}
        """

        compnds = {}

        # Booleans used to determine whether we have finished reading a multi-line
        # entry
        open_name = False
        open_fragment = False

        for rec in self.pdb_recs['COMPND']:

            content = rec['text']

            # Title lines separate different molecule IDs
            if 'MOL_ID:' in content:

                mol_no = int(content.split(':')[1].strip(';'))

                compnds[mol_no] = {
                    'chains': [],
                    'name': '',
                    'fragment': '',
                    'type': 'protein',
                }

            elif 'MOLECULE:' in content:

                text = content.split(':')[1].strip()

                if text:

                    if 'PROTEIN' in text:

                        compnds[mol_no]['type'] = 'protein'

                    elif text.startswith('RNA'):

                        compnds[mol_no]['type'] = 'rna'

                    elif text.startswith('DNA'):

                        compnds[mol_no]['type'] = 'dna'

                    # Complete entries end with a ';' run on lines do not
                    if text[-1] != ';':

                        open_name = True

                    else:

                        open_name = False
                else:

                    open_name = True

                compnds[mol_no]['name'] += text.strip(';')

            elif 'CHAIN:' in content:

                # CHAIN entry contains a comma separated list of chain IDs
                chain_list = [c.strip() for c in content.strip(
                    ';').strip(',')[6:].split(',')]
                self.chain_info.subdivs = chain_list
                compnds[mol_no]['chains'] = chain_list

            elif 'FRAGMENT:' in content:

                text = content.split(':')[1].strip()

                if text[-1] != ';':

                    open_fragment = True

                else:

                    open_fragment = False

                compnds[mol_no]['fragment'] += text.strip(';')

            elif open_name:

                compnds[mol_no]['name'] += content.strip(';')

                if content[-1] == ';':
                    open_name = False

            elif open_fragment:

                compnds[mol_no]['fragment'] += content.strip(';')

                if content[-1] == ';':
                    open_fragment = False

        self.reference_info.compnd = compnds

        if 1 in compnds:
            self.chain_info.chains = compnds[1]['chains']

        return

    def create_resid_dict(self, records):
        """
        Create a dictionary mapping resnames to descriptions from a list of
        dictionaries containing 'resname' and 'text' keys. Used to describe
        HETATMs.

        @type  records:   list
        @param records:   List of parsed PDB lines (dictionaries)
        @rtype:           dictionary
        @return:          Dictionary with residue name as key and text description
                          values
        """

        dic = {}

        for rec in records:
            if rec['resname'] in dic:
                dic[rec['resname']] += rec['text']
            else:
                dic[rec['resname']] = rec['text']

        return dic

    def process_header_het(self):
        """
        Parse HETATM related records to get descriptions of the heterogens.
        Populates self.chain_info.heterogens with a dictionary of the form:

        {chain: {resid: resname}}

        and self.reference_info.hets with:

        {resname: text description}

        and self.reference_info.formuls with:

        {resname: text formula}        

        """

        pdb_recs = self.pdb_recs

        formuls = self.reference_info.formuls

        heterogens = self.chain_info.heterogens

        self.reference_info.hets = self.create_resid_dict(pdb_recs['HETNAM'])
        hets = self.reference_info.hets

        formuls = self.create_resid_dict(pdb_recs['FORMUL'])

        for het in pdb_recs['HET']:

            if het['resname'] not in hets:
                self.logger.info(
                    'HETATM ' +
                    het['resname'] +
                    ' not given description in header')

            if het['resname'] not in formuls:
                self.logger.info(
                    'HETATM ' +
                    het['resname'] +
                    ' formula not given in header')

            chain = het['chain']

            if chain not in heterogens:
                heterogens[chain] = {het['resid']: het['resname']}
            else:
                heterogens[chain][het['resid']] = het['resname']

        return
