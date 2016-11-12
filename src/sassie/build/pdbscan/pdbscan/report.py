# -*- coding: utf-8 -*-
"""
Report
Generates a report detailing the contents and simulation preparedness of input
SasMolScan object.
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

from . import pandoctable as pdt
from . import pdbscan_utils as utils
from textwrap import TextWrapper


def generate_reports(mol):
    """
    Create reports on the contents of the structure. A long report to be saved
    to file detailing all of the contents and a shorter summary focused on
    which, if any, parts of the system are ready for simulation.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @rtype:           list, list
    @return:            1.  List of strings containing summary report focused
                            on the suitability of the system for simulation
                        2.  List of strings containing full report on the
                            contents of the input structure.
    """

    short_report = []
    long_report = []

    short_report.append('# PDB Scan report for ' + mol.pdbname + '\n')

    reconciliation_report = generate_reconciliation_report(mol)

    if reconciliation_report:

        short_report.append(
            '##Warning header information does not match coordinates\n')
        short_report += generate_reconciliation_report(mol)

    short_report += generate_summary_report(mol, model_no=1)

    short_report += generate_simulation_prep_report(mol)

    long_report += short_report

    long_report += create_chain_report(mol, model_no=1)

    return short_report, long_report


def generate_reconciliation_report(mol):
    """
    Create a report detailing any discrepancies between the header and
    coordinate description of the molecule.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @rtype:           list
    @return:          List of strings containing report detailing any
                      discrepancies between the header and coordinate
                      description of the molecule
    """

    reconcile_errors = []

    for chain, reconciled in mol.header_reconciled_chains.iteritems():

        if not reconciled:

            issues = mol.header_reconciliation_status[chain]

            if issues['seq_mismatch']:

                if issues['missing']:

                    missing_ranges = ",".join(
                        utils.get_grouped_ranges(issues['missing']))

                    reconcile_errors.append(
                        'Chain {0:s}: Sequences from header SEQRES records and coordinates do not match (mismatch in missing residue(s): {1:s}).'.format(chain, missing_ranges))

                else:

                    reconcile_errors.append(
                        'Chain {0:s}: Sequences from header SEQRES records and coordinates do not match.'.format(chain))

            if issues['het_mismatch']:

                reconcile_errors.append(
                    'Chain {0:s}: HETATM information in header coordinates do not match.'.format(chain))

    return reconcile_errors


def boolean_to_yn(bool):
    """
    Create a report detailing any discrepancies between the header and
    coordinate description of the molecule.

    @type  bool:      boolean
    @param bool:      Boolean value to convert to a text 'Y'/'N'
    @rtype:           string
    @return:          'Y' if input is True, otherwise 'N'
    """

    if bool:
        text = 'Y'
    else:
        text = 'N'

    return text


def generate_simulation_prep_report(mol):
    """
    Create a report detailing any discrepancies between the header and
    coordinate description of the molecule.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @rtype:           list
    @return:          List of strings containing report detailing the
                      simulation preparedness of the system
    """

    sim_ready_checks = mol.sim_ready

    chain_segment_map = mol.chain_segment_map2()

    rep = []

    rep.append('### Proposed segment to chain mapping\n')

    rep.append('Segments rather than chains are used to create CHARMM ')
    rep.append('parameterized structures. PDB Scan suggests the following ')
    rep.append('mapping of residues from the input chains to segments.\n')

    widths = [6, 7, 4, 13, 16]
    just = ['l', 'l', 'l', 'l', 'l']
    contents = []

    chain_list = sorted(chain_segment_map.keys())
    segname_list = []

    for chain in chain_list:

        segment_info = chain_segment_map[chain]

        tmp_contents = []

        for segname, atom_info in segment_info.iteritems():

            resid_range = str(atom_info['resid'][0]) + \
                '-' + str(atom_info['resid'][1])
            index_range = str(atom_info['original_index'][
                              0]) + '-' + str(atom_info['original_index'][1])
            seg_moltype = atom_info['moltype'][:3]

            tmp_contents.append(
                [chain, segname, seg_moltype, resid_range, index_range])

        contents += sorted(tmp_contents, key=lambda x: x[3])

    for line in contents:
        segname_list.append(line[1])

    header = ['Chain', 'Segname', 'Type', 'Resids', 'Indices']

    rep += pdt.create_pandoc_table(header, contents, widths, just)

    rep.append('### Proposed segment simulation readiness\n')

    rep.append('Each proposed segment has been checked to see if it is ready for ')
    rep.append('simulation through MD (in CHARMM) or SASSIE dihedral MC.\n')

    widths = [7, 11, 6, 6, 6]
    header = ['Segname', 'Single Conformer', 'CHARMM', 'MC', 'MD']
    just = ['c', 'c', 'c', 'c', 'c']
    contents = []

    start_warnings = []

    for segname in segname_list:

        checks = sim_ready_checks[segname]

        conf = boolean_to_yn(checks['single_conformer'])
        charmm = boolean_to_yn(checks['charmm'])
        mc = boolean_to_yn(checks['chain'])
        md = boolean_to_yn(checks['md'])

        line = [segname, conf, charmm, mc, md]

        contents.append(line)

        if not checks['start']:
            start_warnings.append(segname)

    rep += pdt.create_pandoc_table(header, contents, widths, just)

    if start_warnings:
        rep.append('\n')
        warn_txt = ('WARNING: Segments {:s} do not start with resid 1, check '
                    'sequence is correct\n'.format(','.join(start_warnings)))

        rep.append(warn_txt)

    return rep


def generate_summary_report(mol, model_no=1):
    """
    Create report summarizing the content of the molecule - outputs description
    from header, table detailing polymer (nucleic or proteinatious at this
    point) chains, any hetergens and the biological unit(s) given in the
    header.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @type model_no :  integer
    @param model_no:  Model number to get information from and to use in
                      storage dictionaries.
    @rtype:           list
    @return:          List of strings containing summary data on the whole
                      molecule
    """

    header_data = mol.header_data
    chain_info = header_data.chain_info
    ref_info = header_data.reference_info

    summary_report = []

    summary_report.append('## Structure contents summary\n\n')

    # If there is a header we will believe the heterogen information from there
    if mol.header_data.read_valid_header:

        het_res = chain_info.heterogens
        het_dict = header_data.reference_info.hets
        biomt = chain_info.biomt

        summary_report += create_header_summary(header_data)

    else:

        het_res = mol.chain_info.heterogens
        # If there is no header we have no information about what the HET
        # resids represent
        het_dict = {}

        biomt = {}

    summary_report += create_polymer_table(mol, model_no)

    summary_report += create_heterogen_table(het_res, het_dict)

    summary_report += create_biomt_summary(biomt)

    return summary_report


def get_expt_string(metrics, space_group):
    """
    Convert information about experimental outcome from PDB header to a string 
    for output.

    @type  metrics:       dictionary
    @param metrics:       Dictionary containing information parsed from header
    @type  space_group:   string
    @param space_group:   Crystallographic space group
    @rtype:               string
    @return:              String containing experimental quality report
    """

    expt_str = 'Resolution '

    try:
        expt_str += '{0:.2f} '.format(metrics['resolution'])
    except:
        expt_str += '* '

    expt_str += 'Angstrom; R-factor '

    try:
        expt_str += '{0:.2f}'.format(metrics['r'] * 100)
    except:
        expt_str += '* '

    expt_str += '; Free R-factor '

    try:
        expt_str += '{0:.2f}'.format(metrics['r_free'] * 100)
    except:
        expt_str += '* '

    expt_str += '; Space Group {0:s}'.format(space_group)

    return expt_str


def create_header_summary(header):
    """
    Create a table containing a summary header data. The table is output
    as a list of lines in the Pandoc Markdown multi-line table format.

    @type  header:  header_reader.PdbHeader
    @param header:  Information extracted from PDB header records
    @rtype:         list
    @return:        list of strings containing report summarizing the 
                    general information gained about the system contained in 
                    the PDB from its header
    """

    ref_info = header.reference_info
    pdb_recs = header.pdb_recs

    rep = []

    title = ref_info.title

    wrapper = TextWrapper(width=78)

    for line in wrapper.wrap(title):
        rep.append('*' + line + '*')
    rep.append('\n\n')

    # metrics contains resolution, r, and r free values
    metrics = ref_info.metrics

    # primary citation information
    cite_data = ref_info.citation

    if pdb_recs['HEADER']:
        class_txt = pdb_recs['HEADER'][0]['class']
    else:
        class_txt = 'Unknown'

    table_header = ['Class', class_txt]

    # contents will contain line of the table
    contents = []

    if cite_data:

        cite_text = cite_data['authors'] + ' "' + \
            cite_data['title'] + '" ' + cite_data['citation']
        cite_text += ' [Link](http://dx.doi.org/' + cite_data['doi'] + ')'

        contents.append(['Citation', cite_text])

    # EXPDTA contains information on the type of experiment
    if pdb_recs['EXPDTA']:

        expt_type = pdb_recs['EXPDTA'][0]['text']
        contents.append(['Method', expt_type.title()])

        if expt_type.startswith('X-RAY'):

            try:

                space_group = pdb_recs['CRYST1'][0]['space_group']

            except:

                space_group = 'Unknown'

        else:
            space_group = 'N/A'

        if metrics:

            expt_str = get_expt_string(metrics, space_group)

            contents.append(['Experiment', expt_str])

    widths = [20, 59]
    just = ['l', 'l']
    rep += pdt.create_pandoc_table(table_header, contents, widths, just)

    return rep


def create_polymer_table(mol, model_no=1):
    """
    Create a table containing a summary of each chain taken from header data.
    Returns header and list of lists (one for each row) containing column
    values.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @type model_no :  integer
    @param model_no:  Model number to get information from and to use in 
                      storage dictionaries.
    @rtype:           list
    @return:          List of strings containing table of basic information 
                      about polymer chains for a PDB
    """

    contents = []

    if mol.header_data.read_valid_header:

        header = mol.header_data

        # Calculate % of each chain observed in the experiment

        observed = {}

        for chain, seq in header.chain_info.sequence.iteritems():

            tot = len(seq)

            if ((model_no in header.chain_info.missing_resids) and
                (chain in header.chain_info.missing_resids[model_no])):

                # if chain in
                # header.chain_info.missing_resids[model_no].keys():

                miss = len(header.chain_info.missing_resids[model_no][chain])
                observed[chain] = 1.0 - float(miss) / tot

            else:

                observed[chain] = 1.0

        tab_header = [
            'Chain ID',
            'Name',
            'Type',
            'No. Sequence Residues',
            '% Residues Observed']

        for mol_id, mol_info in header.reference_info.compnd.iteritems():

            for ch_id in mol_info['chains']:

                no_res = str(len(header.chain_info.sequence[ch_id]))

                observed_pc = '{0:.0f}'.format(observed[ch_id] * 100)

                if mol_info['type'] in ['rna', 'dna']:
                    mol_type = mol_info['type'].upper()
                else:
                    mol_type = mol_info['type'].title()

                row = [ch_id, mol_info['name'], mol_type, no_res, observed_pc]

                contents.append(row)

        widths = [10, 36, 10, 10, 10]
        just = ['l', 'l', 'l', 'l', 'l']

    else:

        chain_info = mol.chain_info

        tab_header = [
            'Chain ID',
            'Type',
            'No. Residues Found',
            'Gaps Detected']

        chains = mol.chains()
        atom_chains = mol.chain()
        moltypes = mol.moltype()
        resids = mol.resid()

        chain_resids = set(zip(atom_chains, resids))

        chain_types = report_moltypes_chain(chains, atom_chains, moltypes)

        for chid in chains:

            if chain_types[chid] != 'Others' and chain_types[chid] != 'Water':

                n_res = str(len([y for y in chain_resids if y[0] == chid]))

                chain_type = chain_types[chid]

                if chid in chain_info.missing_resids[model_no]:
                    gaps = 'Yes'
                else:
                    gaps = 'No'

                row = [chid, chain_type, n_res, gaps]

                contents.append(row)

        widths = [10, 10, 10, 10]
        just = ['l', 'l', 'l', 'l']

    if contents:

        table = pdt.create_pandoc_table(tab_header, contents, widths, just)

    else:

        table = []

    return table


def report_moltypes_chain(chains, atom_chains, atom_moltypes):
    """
    Ascertain the moltypes of each chain. Mixed type chains reported with a
    '/' separating the different types.

    @type chains         :  list
    @param chains        :  List of all unique chains in the molecule
    @type atom_chains    :  list
    @param atom_chains   :  List of chain ID for each atom
    @type moltype_chains :  list
    @param moltype_chains:  List of chain ID for each atom
    @rtype               :  dictionary
    @return              :  Dictionary with key as chain and value as moltype
    """

    chain_moltype_pairs = set(zip(atom_chains, atom_moltypes))

    chain_types = {}

    for chid, moltype in chain_moltype_pairs:
        if chid in chain_types:
            chain_types[chid].append(moltype)
        else:
            chain_types[chid] = [moltype]

    for chid in chains:
        chain_type = ' / '.join([x.title() for x in chain_types[chid]])
        chain_types[chid] = chain_type

    return chain_types


def create_biomt_summary(biomols):
    """
    Report on the biological units found when parsing the header. The table is
    output as a list of lines in the Pandoc Markdown multi-line table format.

    @type  head_info:   dictionary
    @param head_info:   Dictionary containing information parsed from header
    @rtype:             list
    @return:            list of strings containing report for all biomolecules

    """

    rep = []

    need = True
    if len(biomols) == 1:
        for n, entry in biomols.items():
            if (entry['auth_bio_unit'] == 'MONOMERIC' or
                    entry['soft_bio_unit'] == 'MONOMERIC'):
                need = False

    if need:

        rep.append('\n###Biological Unit\n\n')

        if len(biomols) == 0:

            rep.append('No information detected\n')

        else:

            for biomol_id in sorted(biomols.keys()):
                headline, table = create_biomol_entry(biomol_id, biomols)
                rep.append(headline + '\n')
                if table:
                    rep += table

    return rep


def create_biomol_entry(biomol_id, biomols):
    """
    Report for a biological unit found when parsing the header. Table is
    output as a list of lines in the Pandoc Markdown multiline table format.

    @type  biomol_id:   integer
    @param biomol_id:   Number of the biomolecule to report on
    @type  biomols:     dictionary
    @param biomols:     Dictionary containing BIOMT information parsed from a
                        PDB header
    @rtype:             string, list
    @return:            1. String containing biomolecule report
                        2. list of strings containing table detailing the 
                           transformation suggested to create the biological
                           unit in the PDB header
    """

    entry = ''
    table = []

    biomol = biomols[biomol_id]

    entry += '* Biomolecule: ' + str(biomol_id) + ' - '

    nonmonomer = True
    if biomol['auth_bio_unit'] == 'MONOMERIC' or biomol[
            'soft_bio_unit'] == 'MONOMERIC':
        nonmonomer = False

    # What are the transformations to be applied to?
    if len(biomol['subdivs']) == 1:

        ch_text = 'chain ' + ' '.join(biomol['subdivs'])

    else:

        ch_text = 'chains ' + ' '.join(biomol['subdivs'])

    # Origin of the suggested biological unit

    if biomol['auth_bio_unit'] and biomol['soft_bio_unit']:

        entry += 'Biological units agree for ' + ch_text + \
            ' from both author and software: ' + biomol['auth_bio_unit'] + '\n'

    elif biomol['auth_bio_unit']:

        entry += 'Author asserted biological unit for ' + \
            ch_text + ': ' + biomol['auth_bio_unit'].title() + '\n'

    elif biomol['soft_bio_unit']:

        entry += 'Software suggested biological unit for ' + \
            ch_text + ': ' + biomol['soft_bio_unit'].title() + '\n'

    if nonmonomer:
        # Create table detailing transformations
        head_cols = ['#', 'Rotation', 'Translation']

        contents = []

        for biomt_no in range(len(biomol['rot'])):

            rot = biomol['rot'][biomt_no]
            trans = biomol['trans'][biomt_no]

            for line_no in range(len(rot)):
                rot_txt = '{0:8.3f} {1:8.3f} {2:8.3f}'.format(
                    rot[line_no, 0], rot[line_no, 1], rot[line_no, 2])
                trans_txt = '{0:8.3f}'.format(trans[line_no])
                contents.append([str(biomt_no + 1), rot_txt, trans_txt])

        widths = [5, 47, 26]
        just = ['l', 'l', 'l']

        table = pdt.create_pandoc_table(head_cols, contents, widths, just)

    return entry, table


def create_chain_report(mol, model_no=1):
    """
    Create a report on each individual chain which has a recognized sequence
    in the input molecule. Report contains a FASTA sequence table (with any
    missing residues shown in lower case), a list of the gaps in the sequence,
    a table detailing the location of any heterogens, any missing atoms and
    the disulphide bonds detected.

    @type  mol:       SasMolScan
    @param mol:       Molecular coordinates and associated information
    @type model_no :  integer
    @param model_no:  Model number to get information from and to use in
                      storage dictionaries.
    @rtype:           list
    @return:          List of strings containing
    """

    chain_report = []

    chain_info = mol.chain_info
    ref_info = mol.header_data.reference_info

    for chain in mol.chain_info.subdivs:

        if chain in chain_info.sequence:

            sequence = chain_info.sequence[chain]

            fasta = chain_info.sequence_to_fasta(chain, missing_lower=True)
            start_resid = sequence[0][0]

            chain_report.append("## Sequence for chain {0:s}\n".format(chain))

            chain_report += create_sequence_table(fasta, start_resid)

            chain_report += create_gap_summary(chain,
                                               chain_info.missing_resids,
                                               chain_info.number_gaps)

            chain_report += heterogen_report_chain(chain,
                                                   chain_info.heterogens,
                                                   ref_info.hets)

            chain_report += missing_atoms_chain(chain,
                                                chain_info.missing_atoms)

            chain_report += create_disulphide_table(
                chain, chain_info.disulphides)

            sel_text = 'chain[i] == "{0:s}"'.format(chain)
            chain_report += report_pdb_statistics(mol,
                                                  selection=sel_text, extra={})

    return chain_report


def create_sequence_table(seq, start_resid):
    """
    Produce table containing a single letter (amino acid) sequence. Heterogens 
    are reported as 'X'. 

    @type  seq:  string
    @param seq:  FASTA sequence
    @rtype:      list
    @return:     List of strings making up a table formatting a FASTA style 
                 sequence
    """

    rep = ['### Sequence:\n\n']

    seq_length = 50
    widths = [10, seq_length + 2]

    padded_fasta, row_start = create_padded_fasta(seq, start_resid, widths[1])

    fasta_table = create_fasta_table(padded_fasta, widths, row_start)

    seq_line = '`' + \
        ('123456789|' * ((seq_length / 10) + 1))[:seq_length] + '`'
    header = ['', seq_line]
    rep += pdt.create_pandoc_table(header, fasta_table, widths, ['r', 'l'])

    return rep


def create_padded_fasta(seq, start_resid, width):
    """
    Create a FASTA style sequence from a list of (resid, resname, aa) items,
    where resname are three letters and aa one letter residue names. The start
    of the sequence is padded with - characters to make printing start at a
    round number.

    @type  seq:  string
    @param seq:  FASTA sequence
    @type  width:      integer
    @param width:      Line width (characters) at which to wrap the input
    @rtype:            list, integer
    @return:
                         1.  List of strings of length width containing the wrapped text
                         2.  Starting residue number (including padding)

    """

    offset = start_resid % 10
    start_pad_length = offset - 1
    row_start = start_resid - offset

    padded_seq = '-' * start_pad_length + seq

    return padded_seq, row_start


def create_fasta_table(fasta, col_widths, row_start):
    """
    Create a table containing a FASTA style sequence and a preceding column
    indicating the starting residue number for each line. The table is output
    as a list of lines in the Pandoc Markdown multiline table format.

    @type  seq:          list
    @param seq:          List of (resid, resname, aa) items, where resname are
                         three letters and aa one letter residue names.
    @type  col_widths:   list
    @param col_widths:   List of integer column widths (in characters)
    @rtype:              list
    @return:             Lines forming a Pandoc Markdown multiline table
                         containing a formatted version of the input sequence

    """

    fasta_width = col_widths[1] - 2

    formatted_fasta = format_fasta(fasta, fasta_width)

    table = []

    for line in formatted_fasta:

        # Add `'s to indicate fixed width to Pandoc
        table.append([str(row_start), '`' + line + '`'])

        row_start += fasta_width

    return table


def format_fasta(sequence, width):
    """
    Wrap the of a FASTA type sequence to fit into a column defined by the
    input width.

    @type  sequence:   string
    @param sequence:   Sequence as a continuous string of single letter codes
    @type  width:      integer
    @param width:      Line width (characters) at which to wrap the input
    @rtype:            list
    @return:           List of strings of length width containing the wrapped
                       text
    """

    joined_sequence = ''.join(sequence)

    wrapper = TextWrapper(width=width)

    formatted_sequence = wrapper.wrap(joined_sequence)

    formatted_sequence[-1] += " " * (width - len(formatted_sequence[-1]))

    return formatted_sequence


def create_gap_summary(chain, missing_res_all, numbering_gaps_all, model_no=1):
    """
    Creates a report on the disulphides contained in the specified PDB chain. 

    @type  chain:           string  
    @param chain:           Chain ID  
    @type  missing_res:     dictionary
    @param missing_res:     Missing residues by chain
    @type  numbering_gaps:  dictionary
    @param numbering_gaps:  Gaps in sequential residue numbering by chain
    @type model_no :        integer
    @param model_no:        Model number to get information from and to use in 
                            storage dictionaries.
    @rtype:                 list
    @return:                List of strings containing grouped ranges of 
                            missing residues and any regions which appear to 
                            have non-sequential residue numbering

    """

    if model_no in missing_res_all:
        missing_res = missing_res_all[model_no]
    else:
        missing_res = {}

    if model_no in numbering_gaps_all:
        numbering_gaps = numbering_gaps_all[model_no]
    else:
        numbering_gaps = {}

    rep = ['### Gaps:\n\n']

    if chain in missing_res and missing_res[chain]:
        gaps = utils.get_grouped_ranges(sorted(missing_res[chain].keys()))
        rep.append(', '.join(gaps))
        rep.append('\n')
    else:
        rep.append("None\n")

    if chain in numbering_gaps:
        gaps = utils.get_grouped_ranges(sorted(numbering_gaps[chain].keys()))
        rep.append('### Potential numbering issues:\n\n')
        rep.append(', '.join(gaps))
        rep.append('\n')

    return rep


def heterogen_report_chain(chain, hets, het_dict):
    """
    Creates a report on the heterogens found in the specified PDB chain. 

    @type  chain:     string
    @param chain:     Chain ID
    @type  hets:      dictionary
    @param hets:      Hetrogens by chain
    @type  het_dict:  dictionary
    @param het_dict:  Description for heterogens by resid
    @rtype:           list
    @return:          List of strings containing table detailing heterogens in 
                      the selected chain
    """

    rep = ['\n### Heterogens:\n\n']

    het_head = ['Residue No.', 'Resdiue ID', 'Description']

    content = []

    water_res = []

    if chain in hets:

        het_resids = sorted(hets[chain].keys())

    else:

        het_resids = []
        rep.append('\nNone Found')

    for resid in het_resids:

        resname = hets[chain][resid]

        if resname == 'HOH':
            # Water is a special case - very common and usually multiple
            # molecules together
            # Add residue number to list of water
            water_res.append(resid)

        else:

            if resname in het_dict:

                desc = het_dict[resname]

            else:

                desc = 'No description available'

            content.append([str(resid), resname, desc])

    if het_resids:

        # Group water molecules ranges for output
        if water_res:
            wat_list = ''
            wat_ranges = utils.get_grouped_ranges(water_res)
            wat_list += ', '.join(wat_ranges)
            content.append([wat_list, 'HOH', 'Water'])

        het_widths = [20, 10, 48]
        het_just = ['l', 'l', 'l']

        rep += pdt.create_pandoc_table(het_head, content, het_widths, het_just)

    rep.append('\n')

    return rep


def missing_atoms_chain(chain, missing_atms, model_no=1):
    """
    Creates a report on the missing atoms for the specified PDB chain. 

    @type  chain:          string
    @param chain:          Chain ID
    @type  missing_atms:   dictionary
    @param missing_atms:   Missing atoms by chain
    @type model_no :       integer
    @param model_no:       Model number to get information from and to use in 
                           storage dictionaries.
    @rtype:                list
    @return:               List of strings containing table detailing which 
                           residues are missing which atoms in the selected 
                           chain
    """

    rep = []

    rep.append("\n### Missing heavy atoms:\n\n")

    if model_no in missing_atms and chain in missing_atms[model_no]:

        if missing_atms and chain in missing_atms.keys():

            resids = sorted(missing_atms[chain].keys())
            contents = []

            for resid in resids:

                res = missing_atms[chain][resid]
                contents.append(
                    [str(resid), res['resname'], ' '.join(res['atoms'])])

            rep += pdt.create_pandoc_table(['Resid',
                                            'Resname',
                                            'Missing'],
                                           contents,
                                           [8,
                                            8,
                                            50],
                                           ['l',
                                            'l',
                                            'l'])

    else:
        rep.append('No missing heavy atoms detected in recognized residues')

    rep.append('\n')

    return rep


def create_disulphide_table(chain, dis):
    """
    Creates a report on the disulphides contained in the specified PDB chain. 

    @type  chain:     string
    @param chain:     Chain ID
    @type  dis:       list
    @param dis:       List of disulphide bonds
    @rtype:           list
    @return:          List of strings containing table detailing the disulphide 
                      bonds in the chain
    """

    rep = []

    if dis:

        contents = []

        for bond in dis:

            res0 = bond.bound[0]
            res1 = bond.bound[1]

            if res0['subdiv'] == chain or res1['subdiv'] == chain:
                contents.append([res0['subdiv'],
                                 str(res0['resid']),
                                 res1['subdiv'],
                                 str(res1['resid'])])

        if contents:
            rep.append("\n*Disulphide bonds:*\n\n")
            headers = ['Chain 1', 'Resid 1', 'Chain 2', 'Resid 2']
            col_widths = [8, 8, 8, 8]
            col_align = ['l', 'l', 'l', 'l']
            rep += pdt.create_pandoc_table(headers,
                                           contents,
                                           col_widths,
                                           col_align)
            rep.append('\n')

    return rep


def create_heterogen_table(hets, het_dict):
    """
    Create a table describing all heterogens, listing which chains they occur 
    in.

    @type  hets:       dictionary
    @param hets:       Hetrogens by chain
    @type  het_dict:   dictionary
    @param het_dict:   Dictionary containing heterogen descriptions from header
    @rtype:            list
    @return:           List of strings containing a table of the heterogens

    """

    contents = []

    het_list = {}

    table = ['\n### Heterogens \n']

    for chain in sorted(hets.keys()):

        het_resids = hets[chain]

        for resid in het_resids:

            resname = hets[chain][resid]

            if resname != 'HOH':

                if resname not in het_list.keys():

                    het_list[resname] = [chain]

                else:

                    if chain not in het_list[resname]:
                        het_list[resname] += chain

    if het_list:

        for res, chains in het_list.iteritems():

            if res in het_dict.keys():

                het_line = [res, het_dict[res], ','.join(chains)]

            else:

                het_line = [res, 'No header description', ','.join(chains)]

            contents.append(het_line)

        header = ['Residue', 'Name', 'Chain IDs']

        widths = [10, 58, 10]
        just = ['l', 'l', 'l']

        table += pdt.create_pandoc_table(header, contents, widths, just)

    else:

        table += ["None"]

    return table


def report_pdb_statistics(mol, selection='', extra={}):
    """
    Creates a summary structural statistics for the selection provided within
    the molecule given.

    @type  mol:         SasMolScan
    @param mol:         Molecular coordinates and associated information
    @type  selection:   string
    @param selection:   SASSIE style selection text
    @type  extra:       dictionary
    @param extra:       Extra lines for the summary table provided - key =
                        metric name, value = entry for metric
    @rtype:             list
    @return:            List of strings containing a table detailing structural
                        metrics
    """

    molecule_dictionary = mol.build_molecule_dictionary(selection)

    widths = [50, 100]
    just = ['l', 'l']

    contents = []

    # Table header, not PDB header
    header = None

    for key, value in molecule_dictionary.iteritems():
        str_property = str(key)
        str_value = str(value)
        contents.append([str_property, str_value])

    if len(extra):
        for key, value in extra.items():
            contents.append([key, value])

    rep = ['### Statistics:\n\n']
    rep += pdt.create_pandoc_table(header, contents, widths, just)

    return rep
