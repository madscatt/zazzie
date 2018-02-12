Altens v.1
# You should change import for altens

##### input files
1. pdb file: Current altens works only with a pdb file to extract N-H unit vector 
             For consistency, pdb should be compatible with sasmol and CHARMM.
             It is recommended to run pdbrx before running altens if a user use a pdb file from database. 

2. rdc file: Current altens is applicable to N-H vector case. 
             Multiple vector types are not considered in this version.  
             Dimension: N(residue_exp) X 2 
3. residue list file: list of residues that are considered for analysis and selected by a user.
             Dimension: 1 X N(residue_user) 
             Note that N(residue_user) can be different from N(residue_exp).
  
4. Optional
   1. mcon: Type - boolean 
         If mcon == True, Perform Monte-Carlo analysis for the calculated orientation of protein

##### Output files

1. runname_frame.txt: includes summary of inputs and outputs.
2. runname_calc_rdc_frame.txt: includes calculated and experimental RDCs.
3. runname_mc_analysis_frame.txt: includes summary of MC analysis for calculated orientation.
4. runname_mc_trajectory_frame.txt: include the trajectories of orientational information from MC.  

###### Example
Input file: 1D3Z_mod1_new.pdb was obtained from 1D3Z_mod1.pdb (David Fushiman provided) by running pdbrx. 
            See 1D3Z_mod1.pdb has different naming convention and numbering scheme with CHARMM.
Output file: 
            run_0_00001.txt 
            run_0_calc_rdc_00001.txt
            run_0_mc_analysis_00001.txt
            run_0_mc_trajectory_00001.txt

