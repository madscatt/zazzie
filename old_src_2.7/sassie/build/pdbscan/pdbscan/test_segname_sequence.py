import sassie.build.pdbscan.pdbscan as pdbscan
import numpy as np

pdbfile = '../testing/data/s.pdb'

mol = pdbscan.SasMolScan()
mol.read_pdb(pdbfile)
mol.run_scan()

print 'chain info'
print 'type(mol.chain_info) = ' + str(type(mol.chain_info))
print 'type(mol.chain_info.sequence) = ' + str(type(mol.chain_info.sequence))
for k,v in mol.chain_info.sequence.iteritems():
    chain_sequence =  v
    print 'chain : ' + k + '\n\tsequence = %s\n' % (chain_sequence,)

print 'segname info'
print 'type(mol.segname_info) = ' + str(type(mol.segname_info))
print 'type(mol.segname_info.sequence) = ' + str(type(mol.segname_info.sequence))

for k,v in mol.segname_info.sequence.iteritems():
    segname_sequence =  v
    print 'segname : ' + k + '\n\tsequence = %s\n' % (segname_sequence,)


class dum():

    def __init__(self):
        pass

    def get_segment_sequence_info(self, mol):
        """
        Extract sequence from coordinates and then fill in the segname_info with
        resnames taken from chain_info (if this has been filled from header
        information).

        """

        # Get coordinate sequence and gap information
        mol.extract_sequence_info(subdiv_type='segname')

        model_no = mol.model_no

        resids = mol.resid()
        segnames = mol.segname()
        chains = mol.chain()
        moltypes = mol.moltype()
        natoms = mol.natoms()

        chain_info = mol.chain_info
        segname_info = mol.segname_info

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

            print 'start_ndx = ', start_ndx

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

                        if resid < start_resid and resname != '':

                            segname_info.add_missing_resid(
                                segname, resid, resname, model_no=model_no)

                if end_ndx in chain_ends:

                    end_resid = resids[end_ndx]

                    for resid, resname in chain_missing_res.iteritems():

                        if resid > end_resid and resname != '':

                            segname_info.add_missing_resid(
                                segname, resid, resname, model_no=model_no)

                # Copy in information on gaps in the coordinate sequence

                start_resid = resids[start_ndx]
                end_resid = resids[end_ndx]

                for resid in range(start_resid, end_resid + 1):

                    if resid in chain_missing_res.keys():

                        # if ((resid not in segname_info.sequence[segname]) or
                        #         (segname_info.sequence[segname][resid] == '')):
                        if (resid not in segname_info.sequence[segname]):

                            segname_info.add_missing_resid(segname,
                                                           resid,
                                                           chain_missing_res[resid],
                                                           model_no=model_no)

                    if resid in chain_num_gaps.keys():

                        segname_info.add_number_gap(segname, resid, chain_num_gaps[resid])

        # Complete segment sequences with resnames from missing residues
        seg_missing = segname_info.missing_resids[model_no]

        for segname in sorted(seg_missing.keys()):

            seq = segname_info.sequence[segname]
            resid_list = [x[0] for x in seq]

            pre_seq = []


            for resid in sorted(seg_missing[segname].keys()):

                resname = seg_missing[segname][resid]

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

d = dum()
d.get_segment_sequence_info(mol)

print 'final segname info'
print 'type(mol.segname_info) = ' + str(type(mol.segname_info))
print 'type(mol.segname_info.sequence) = ' + str(type(mol.segname_info.sequence))

for k,v in mol.segname_info.sequence.iteritems():
    segname_sequence =  v
    print 'segname : ' + k + '\n\tsequence = %s\n' % (segname_sequence,)

