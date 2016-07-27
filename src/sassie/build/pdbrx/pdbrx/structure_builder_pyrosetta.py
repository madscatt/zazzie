import rosetta

#args = '-chemical:exclude_patches VirtualDNAPhosphate'
#args = '-add_orbitals'
#rosetta.init(extra_options=args)

rosetta.init()

from rosetta.protocols.grafting import CCDEndsGraftMover
from rosetta.protocols.loops.loop_closure.ccd import CCDLoopClosureMover
from rosetta.protocols.loops.loop_mover.refine import LoopMover_Refine_CCD

import salign

import os
import logging

class StructureBuilderPyRosetta():

    def __init__(self, scaffold_pdb, gap_descriptions, chain, protein_only = True):

        self.logger = logging.getLogger(__name__)

        # Scoring function to determine probability of moves
        self.scorefxn = rosetta.core.scoring.get_score_function()

        self.res_type = 'fa_standard'

        # Load pose from input PDB file
        self.scaffold_pose =  rosetta.pose_from_pdb(scaffold_pdb)

        # Gap descriptions as list of lists:
        # [[pre_anchor, post_anchor, pre_flank, gap, post_flank], ...]
        self.gap_descriptions = gap_descriptions
        self.gap_descriptions.sort(key=lambda x:x[0],reverse=True)

        self.Loops = rosetta.Loops()

        self.loop_list = []

        self.chain = chain

        self.protein_only = protein_only

        return

    def get_last_residue_id_chain(self, pose, chain):
        """
        Get the PyRosetta resdiue ID of the last residue of the selected chain
        in the given pose.

        @type  pose:   Pose
        @param pose:   PyRosetta pose from which to determine the residue ID
        @type  chain:  string
        @param chain:  Chain letter of interest (NOT the Pyrosetta numerical ID)
        @rtype :       integer
        @return:       PyRosetta resdiue ID of the end of chain
        """

        # Rosetta uses numbers for chains
        chid = rosetta.core.pose.get_chain_id_from_chain(chain,pose)

        # PyRosetta residue numbering rather than PDB
        last_res_no = pose.conformation().chain_end(chid)

        return last_res_no


    def get_first_residue_id_chain(self, pose, chain):
        """
        Get the PyRosetta resdiue ID of the first residue of the selected chain
        in the given pose.

        @type  pose:   Pose
        @param pose:   PyRosetta pose from which to determine the residue ID
        @type  chain:  string
        @param chain:  Chain letter of interest (NOT the Pyrosetta numerical ID)
        @rtype :       integer
        @return:       PyRosetta resdiue ID of the start of chain
        """

        # Rosetta uses numbers for chains
        chid = rosetta.core.pose.get_chain_id_from_chain(chain,pose)

        # PyRosetta residue numbering rather than PDB
        first_res_no = pose.conformation().chain_begin(chid)

        return first_res_no

    def insert_single_internal_residue(self, chain, anchors, seq, pre_flank_seq):

        info = self.scaffold_pose.pdb_info()

        start_res = info.pdb2pose(chain,anchors[0])

        seq_overhang = pre_flank_seq + seq

        overhang_len = len(pre_flank_seq)

        if overhang_len == 2:

            scaffold_range = [start_res]

        else:

            scaffold_range = [start_res -1, start_res]

        frag_range = range(1,overhang_len + 1)

        frag_pose = rosetta.pose_from_sequence(seq_overhang,self.res_type, auto_termini=False)

        aligned_frag = salign.kabsch_alignment(self.scaffold_pose, frag_pose, scaffold_range, frag_range)

        for i in range(2):
            aligned_frag.delete_residue_range_slow(1,1)

        self.scaffold_pose.append_polymer_residue_after_seqpos(aligned_frag.residue(1),start_res, False)

        self.scaffold_pose.pdb_info.set_resinfo(start_res + 1, chain, anchors[0] + 1)

        self.scaffold_pose.pdb_info().obsolete(False)

        loop = rosetta.Loop(start_res, start_res+2, start_res+1)

        return loop

    def graft_loop_double_end(self, chain, anchors, seq, pre_flank_seq, post_flank_seq):
        """
        Graft a loop fragment into a gap in the self.scaffold_structure.
        Requires flanking residues to be present at both ends.

        @type  chain:    string
        @param chain:    Chain into which loop will be grafted
        @type  anchors:  list
        @param anchors:  Resids before and after region to be inserted
        @type  seq:      string
        @param seq:      Sequence of the inserfed fragment
        @type  pre_flank_seq:  string
        @param pre_flank_seq:  Sequence of the residues preceding the inserted fragment
        @type  post_flank_seq:  string
        @param post_flank_seq:  Sequence of the residues following the inserted fragment
        @rtype:          Pyrosetta Loop
        @return:         Loop describing the region added/remodelled
        """

        loop_length = len(seq)

        if loop_length == 1:

            loop = self.insert_single_internal_residue(chain, anchors, seq, pre_flank_seq)

        else:

            info = self.scaffold_pose.pdb_info()

            start_res = info.pdb2pose(chain,anchors[0])
            end_res = info.pdb2pose(chain,anchors[1])

            #nter_overlap_length = len(pre_flank_seq)
            #cter_overlap_length = len(post_flank_seq)

            if loop_length < 4:

                seq = pre_flank_seq[-1] + seq + pre_flank_seq[0]

                loop_length = len(seq)

                pre_flank_seq = pre_flank_seq[:-1]
                post_flank_seq = post_flank_seq[1:]

                self.scaffold_pose.delete_polymer_residue(end_res)
                self.scaffold_pose.delete_polymer_residue(start_res)

                self.scaffold_pose.pdb_info().obsolete(False)

                anchors[0] = anchors[0] - 1
                anchors[1] = anchors[1] + 1

                start_res = info.pdb2pose(chain,anchors[0])
                end_res = info.pdb2pose(chain,anchors[1])

                flex_length_nter = 1
                flex_length_cter = 1

            elif loop_length == 4:

                flex_length_nter = 1
                flex_length_cter = 1

            elif loop_length > 12:
                flex_len = int((loop_length/2) - 4)
                flex_length_nter = flex_len
                flex_length_cter = flex_len

            else:
                flex_length_nter = 2
                flex_length_cter = 2


            seq_overhang = pre_flank_seq + seq + post_flank_seq
            nter_overlap_length = len(pre_flank_seq)
            cter_overlap_length = len(post_flank_seq)

            frag_pose = rosetta.pose_from_sequence(seq_overhang,self.res_type, auto_termini=False)

            mover = CCDEndsGraftMover(start_res, end_res, frag_pose, nter_overlap_length, cter_overlap_length)

            mover.set_insert_flexibility(flex_length_nter,flex_length_cter)

            mover.apply(self.scaffold_pose)

            # Give grafted residues correct chain and numbering info

            for ii in range(loop_length):
                res = start_res + 1 + ii
                pdb_res = anchors[0] + 1 + ii
                self.scaffold_pose.pdb_info().set_resinfo(res,chain,pdb_res)

            # Needed for updated information to be accessible
            self.scaffold_pose.pdb_info().obsolete(False)

            # Setup loop object for later refinement
            loop_begin = self.scaffold_pose.pdb_info().pdb2pose(chain, anchors[0]) + 1
            loop_end = self.scaffold_pose.pdb_info().pdb2pose(chain,anchors[1]) - 1

            loop = rosetta.Loop(loop_begin,loop_end,(loop_begin + loop_end)/2)

        return loop

    def graft_c_terminal_loop(self, chain, anchors, seq, pre_flank_seq):
        """
        Graft a C-terminal loop to self.scaffold_structure.
        Requires at least one flanking residue to be given.

        @type  chain:    string
        @param chain:    Chain into which loop will be grafted
        @type  anchors:  list
        @param anchors:  Resids before and after region to be inserted
        @type  seq:      string
        @param seq:      Sequence of the inserfed fragment
        @type  pre_flank_seq:  string
        @param pre_flank_seq:  Sequence of the residues preceding the inserted fragment
        @rtype:          Pyrosetta Loop
        @return:         Loop describing the region added/remodelled
        """

        frag_seq = pre_flank_seq + seq

        frag_pose = rosetta.pose_from_sequence(frag_seq, 'fa_standard', auto_termini=False)

        # PDB resid of anchor
        start_anchor = anchors[0]

        # PyRosetta indices of the pre_flank_seq
        frag_range = range(1,len(pre_flank_seq) + 1)

        # PyRosetta index of the final residue
        # The terminal residue of the scaffold structure will be deleted and
        # replaced from the fragment
        #original_end = self.scaffold_pose.total_residue()
        original_end = self.get_last_residue_id_chain(self.scaffold_pose,chain)

        # If we can it helps to use an idealized residue in the overlap
        # so delete the original residue and generate replacement
        # The chain isn't built this way as it produces overlapping conformations
        if len(pre_flank_seq) > 1:
            self.scaffold_pose.delete_polymer_residue(original_end)
            self.scaffold_pose.append_polymer_residue_after_seqpos(frag_pose.residue(1),original_end -1, True)

            self.scaffold_pose.pdb_info().set_resinfo(original_end,chain,start_anchor)

        # If two flanking residues available use for alignment range
        if len(frag_range) == 1:
            scaffold_range = [original_end]
        else:
            scaffold_range = [original_end -1, original_end]

        aligned_frag = salign.kabsch_alignment(self.scaffold_pose, frag_pose,
                                               scaffold_range, frag_range)

        self.scaffold_pose.delete_polymer_residue(original_end)

        loop_length = len(frag_seq)

        # In most cases we should have two flanking residues
        # First is used for alignment only - so delete now
        if len(frag_range) != 1:
            aligned_frag.delete_residue_range_slow(1,1)
            loop_length -= 1

        post_edit_end = self.get_last_residue_id_chain(self.scaffold_pose, chain)
        #self.scaffold_pose.append_pose_by_jump(aligned_frag, post_edit_end)

        for ii in range(loop_length):

            frag_res = ii + 1

            new_res = post_edit_end + 1
            pdb_res = start_anchor + ii

            self.scaffold_pose.append_polymer_residue_after_seqpos(aligned_frag.residue(frag_res),post_edit_end, False)

            post_edit_end += 1

            self.scaffold_pose.pdb_info().set_resinfo(new_res,chain,pdb_res)

        # for ii in range(loop_length):
        #
        #     res = original_end + ii
        #     pdb_res = start_anchor + ii
        #     self.scaffold_pose.pdb_info().set_resinfo(res,chain,pdb_res)

        self.scaffold_pose.pdb_info().obsolete(False)

        end_res = self.get_last_residue_id_chain(self.scaffold_pose, chain)

        if len(frag_range) == 1:
            loop = rosetta.Loop(original_end ,end_res-1, (original_end + end_res)/2)
        else:
            loop = rosetta.Loop(original_end-1 ,end_res-1, (original_end + end_res)/2)

        return loop

    def graft_n_terminal_loop(self, chain, anchors, seq, post_flank_seq):
        """
        Graft a N-terminal loop to self.scaffold_structure.
        Requires at least one flanking residue to be given.

        @todo: Check if we need the replacement shenanigans used for
        C-terminal additions

        @type  chain:    string
        @param chain:    Chain into which loop will be grafted
        @type  anchors:  list
        @param anchors:  Resids before and after region to be inserted
        @type  seq:      string
        @param seq:      Sequence of the inserted fragment
        @type  post_flank_seq:  string
        @param post_flank_seq:  Sequence of the residues following the inserted fragment
        @rtype:          Pyrosetta Loop
        @return:         Loop describing the region added/remodelled
        """

        logger = self.logger

        logger.info("N-terminal Loop modelling")

        frag_seq = seq + post_flank_seq

        frag_pose = rosetta.pose_from_sequence(frag_seq, 'fa_standard', auto_termini=False)
        logger.info("Fragment built")

        # PDB resid of anchor
        start_anchor = anchors[1]

        # PyRosetta indices of the pre_flank_seq
        #frag_range = range(1,len(post_flank_seq) + 1)
        last_frag_resid = len(frag_seq)
        frag_range = [last_frag_resid-1,last_frag_resid]


        # PyRosetta index of the first residue
        # The terminal residue of the scaffold structure will be deleted and
        # replaced from the fragment
        original_term = self.get_first_residue_id_chain(self.scaffold_pose,chain)

        # If two flanking residues available use for alignment range
        if len(frag_range) == 1:
            scaffold_range = [original_term]
        else:
            scaffold_range = [original_term, original_term + 1]

        logger.info("Anchor alignment setup complete")
        aligned_frag = salign.kabsch_alignment(self.scaffold_pose, frag_pose,
                                               scaffold_range, frag_range)
        logger.info("Fragment aligned")

        #self.scaffold_pose.delete_polymer_residue(original_term)
        self.scaffold_pose.delete_residue_range_slow(original_term,original_term)
        logger.info("Anchor deleted")

        loop_length = len(frag_seq)

        # Generated with default chain label 'A'
        frag_term = self.get_last_residue_id_chain(aligned_frag, 'A')

        # In most cases we should have two flanking residues
        # First is used for alignment only - so delete now
        if len(frag_range) != 1:
            aligned_frag.delete_residue_range_slow(frag_term,frag_term)
            loop_length -= 1

        logger.info("Delete alignment residue")

        count = 0

        for ii in range(loop_length, 0 , -1):

            print "add residue " + str(ii)

            self.scaffold_pose.prepend_polymer_residue_before_seqpos(frag_pose.residue(ii),1,False)
            pdb_res = start_anchor - count
            self.scaffold_pose.pdb_info().set_resinfo(1,chain,pdb_res)
            count += 1

        self.scaffold_pose.pdb_info().obsolete(False)

        if len(frag_range) == 1:
            loop = rosetta.Loop(1 ,original_term, (1 + original_term)/2)
        else:
            loop = rosetta.Loop(1 ,original_term + 1, (1 + original_term + 1)/2)

        logger.info("Loops generated")

        return loop

    def prepare_loop_refine(self, loop, loops):
        """
        @type  loop:  Loop
        @param loop:  Loop to be refined
        @type  loops:  Loops
        @param loops:  List of loops contained in scaffold_pose
        @rtype :  movemap
        @return:  Movemap showing possible moves in the loop
        """

        logger = self.logger

        logger.info("Prepare refinement")

        rosetta.add_single_cutpoint_variant(self.scaffold_pose, loop)
        rosetta.set_single_loop_fold_tree(self.scaffold_pose, loop)

        loops.add_loop(loop)

        move_map = rosetta.move_map_from_loop(self.scaffold_pose,loop,False)

        logger.info("Preparation complete")

        return move_map

    def refine_loop(self, loop, loops, move_map, inner_cycles = 10):

        logger = self.logger

        logger.info("Peform refinement")

        loop_refine = LoopMover_Refine_CCD( loops, self.scorefxn )
        loop_refine.max_inner_cycles(inner_cycles)
        loop_refine.move_map(move_map)

        logger.info("Refinement set up")

        loop_refine.apply(self.scaffold_pose)
        logger.info("Refinement applied")

        ccd = CCDLoopClosureMover(loop, move_map)
        ccd.apply(self.scaffold_pose)
        logger.info("Loop closed")

        return

    def model_all_loops(self):

        chain = self.chain

        loops = self.Loops

        #orig_conf = self.scaffold_pose.conformation()

        # Gap descriptions as list of lists:
        # [[pre_anchor, post_anchor, pre_flank, gap, post_flank], ...]

        for loop_desc in self.gap_descriptions:

            anchors = [loop_desc[0], loop_desc[1]]

            seq = loop_desc[3]
            pre_flank_seq = loop_desc[2]
            post_flank_seq = loop_desc[4]

            if pre_flank_seq and post_flank_seq:

                loop = self.graft_loop_double_end(chain, anchors, seq, pre_flank_seq, post_flank_seq)

            elif pre_flank_seq:

                loop = self.graft_c_terminal_loop(chain, anchors, seq, pre_flank_seq)

            elif post_flank_seq:

                loop = self.graft_n_terminal_loop(chain, anchors, seq, post_flank_seq)

            else:

                raise Exception('At least one flanking region must be provided for each fragment')

            if self.protein_only:

                self.loop_list.append(loop)
                movemap = self.prepare_loop_refine(loop, loops)

                self.refine_loop(loop,loops, movemap)

        return

    def complete_structure(self, output_path, filename):

        logger = self.logger

        out_filename = os.path.join(output_path, filename)

        logger.info("About to model loops")

        self.model_all_loops()

        logger.info("Writing output PDB")

        # Save output as PDB
        self.scaffold_pose.dump_pdb(out_filename)

        return out_filename



