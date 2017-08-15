#
# Hailiang Zhang
# NIST & UTK
#

import sys, os, locale, glob
import random
import numpy as np
#import h5py

import sasmol.sasmol as sasmol
# import sasproperties

import sassie.util.basis_to_python as basis_to_python

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'periodictable-1.3.0'))

import periodictable as pt
import periodictable.xsf as xxx

ELECTRON_RADIUS = pt.constants.electron_radius*1.0e14   #converts from units of m to units of 10^-12 cm

def set_atomic_radii(mol, radii, scale):

#	print "VDW radii scale factor: ", scale

#    atomic=sasproperties.Atomic()
    atomic = sasmol.sasproperties.Atomic()
    natoms = mol.natoms()
#    element_names = [atom_name[0] for atom_name in mol.name()]


    R_VDW = atomic.van_der_Waals_radii()

    for i in range(natoms):
        for key, item in R_VDW.items():
            if key == mol.element()[i]:
#                print key, item
                radii[i] = item*scale
    return radii

def calculate_solvent_sl(fraction_D2O):
    '''
    Functionality: calculate the H2O/D2O solvent scattering length density based on the fraction of D2O
    Input: fraction_D2O--fraction of D2O
    Return: solvent scattering length (x 10^-12 cm)
    Note: room temperature is assumed
    '''
    fraction_H2O = 1.0-fraction_D2O
    bH = float(pt.H.neutron.b_c)*0.1
    bD = float(pt.D.neutron.b_c)*0.1
    bO = float(pt.O.neutron.b_c)*0.1
#    print bH,bD,bO
    b_H2O = bH*2 + bO
    b_D2O = bD*2 + bO
    b_solvent = fraction_H2O*b_H2O + fraction_D2O*b_D2O
#    print b_H2O,fraction_H2O,b_D2O, fraction_D2O
#    print 'b_solvent: ',b_solvent

    return b_solvent

def calculate_solvent_electron_sl(qlist):
    '''
    Return: solvent electron scattering length (x 10^-12 cm) as a function of q
    Note: room temperature is assumed
    '''

#    print 'electron_radius: ', ELECTRON_RADIUS

    bxH = np.zeros(qlist.size)
    bxO = np.zeros(qlist.size)
    electron_b_solvent = np.zeros(qlist.size)

    for j in range(qlist.size):

        bxH[j] = xxx.Xray(pt.H).f0(qlist[j])*ELECTRON_RADIUS
        bxO[j] = xxx.Xray(pt.O).f0(qlist[j])*ELECTRON_RADIUS
#        print j, qlist[j], bxH[j], bxO[j]
        electron_b_solvent[j] = 2*bxH[j] + bxO[j]

#    print 'electron b_solvent: ',electron_b_solvent

    return electron_b_solvent

def find_active_Hs(mol):
    '''
    Functionality: find the active Hs for a molecule
    Input: mol--sasmol object
    Return: list of masks for the active Hs (0--is NOT a active H; 1--is a active H)
    Note: it is assumed that Hs always follows its attached heavy atoms in the pdb file
    '''
    natoms = mol.natoms()
    mask_active_Hs = [0]*natoms
#   need to ignore the H5T at the N-terminus of nucleic acids that appears before any "nonH" atom
    for i in xrange(natoms):
        atom_full_name = mol.name()[i]
        if atom_full_name == 'H5T':
            continue
        if atom_full_name[0] == 'H':
            if previous_nonH_atom_short_name != 'C':
                mask_active_Hs[i] = 1
            else:
                mask_active_Hs[i] = 0
        else:
#            mask_active_Hs[i] = 0
            previous_nonH_atom_short_name = atom_full_name[0]
    # print "number of non_aliphatic Hs found: ", mask_active_Hs.count(1)
    return mask_active_Hs

def find_inert_Hs(mol):
    '''
    Functionality: find the inert Hs for a molecule
    Input: mol--sasmol object
    Return: list of masks for the inert Hs (0--is NOT a inert H; 1--is a inert H)
    Note: it is assumed that Hs always follows its attached heavy atoms in the pdb file
    '''
    natoms = mol.natoms()
    mask_inert_Hs = [0]*natoms
#   need to ignore the H5T at the N-terminus of nucleic acids that appears before any "nonH" atom
    for i in xrange(natoms):
        atom_full_name = mol.name()[i]
        if atom_full_name == 'H5T':
            continue
        if atom_full_name[0] == 'H':
            if previous_nonH_atom_short_name == 'C':
                mask_inert_Hs[i] = 1
            else:
                mask_inert_Hs[i] = 0
        else:
#            mask_active_Hs[i] = 0
            previous_nonH_atom_short_name = atom_full_name[0]
    # print "number of aliphatic Hs found: ", mask_inert_Hs.count(1)
    return mask_inert_Hs

def set_atomic_neutron_bs(mol, B):
    '''
    Functionality: get the neutron b values for each individual atom
    Input: mol--sasmol object
           B--list of scattering lengths for each atom
    '''
    natoms = mol.natoms()

#    element_names = [atom_name[0] for atom_name in mol.name()]

    for i in range(natoms):
        for el in pt.elements:
# The following line is needed because periodic table lists a 0th element, "n" that could be confused with "N" for nitrogen.
            if el.symbol != 'n':
# Only the first letter is capitalized in the periodic table but all letters are capitalized in sasproperties
                if el.symbol.upper() == mol.element()[i]:
                    B[i] = float(el.neutron.b_c)*0.1
#                    print i, el.symbol, B[i]

def set_atomic_xray_bs(mol, B, qlist):
    '''
    Functionality: get the X-ray f0(q) values for each individual atom (Cromer-Mann formula)
    Input: mol--sasmol object
           B -- a list of f0 at different q values for each atom
    '''
    natoms = mol.natoms()
#   print 'natoms: ',natoms

#    element_names = [atom_name[0] for atom_name in mol.name()]


    for j in range(qlist.size):
#        print j
        for i in range(natoms):
            for el in pt.elements:
# The following line is needed because periodic table lists a 0th element, "n" that could be confused with "N" for nitrogen.
                if el.symbol != 'n':
# Only the first letter is capitalized in the periodic table but all letters are capitalized in sasproperties
                    if el.symbol.upper() == mol.element()[i]:
                        B[j, i] = xxx.Xray(el).f0(qlist[j])*ELECTRON_RADIUS
#                        print i, el.symbol, B[j, i]


def deuterate_neutron_bs(mol, B_neutron_vacuum, mask_active_Hs, mask_inert_Hs, mask_exH_regions, mask_deuteration_regions, mvars, D2O_fraction, o_pdbfile, ofile, ofile1):
    '''
    Functionality: deuterate the exchangeable H's based on the HD exhcange ratio in various ways
    Input: mol--sasmol object
           HD_exchange_ratio--H/D exchange ratio
    Return: None
    Note: this method is named in consistent with the CPP code and may be misleading, because there is no atomic atom type changed from H to D in the sasmol object
    '''
    of = open(ofile, 'w')
    of1 = open(ofile1, 'w')
    natoms = mol.natoms()
    for i in xrange(natoms):
        if mol.element()[i] == 'D':
            mol.element()[i] = 'H'
    for i in xrange(mvars.number_exH_regions):
        fraction_exH = mvars.fraction_exH_array[i]
        #
        mask_active_Hs_regional = [x and y for (x, y) in zip(mask_active_Hs, mask_exH_regions[i])]
        count_target_exH_active_Hs = round(mask_active_Hs_regional.count(1)*fraction_exH*D2O_fraction)
        indices_active_Hs = np.where(np.array(mask_active_Hs_regional) == 1)[0]
        np.random.shuffle(indices_active_Hs)
        indices_exH_active_Hs = indices_active_Hs[:count_target_exH_active_Hs]
        of.write("Fraction D2O: %.2f"%(D2O_fraction)+"; exH region: %3d"%(i+1)+"; fraction exH: %.2f"%(fraction_exH)+"; number exH: %5d"%(count_target_exH_active_Hs)+"\n")
        #
        for i in xrange(natoms):
            if i in indices_exH_active_Hs:
                of.write("Deuterate atom #%5d"%(i+1)+" (atom name: %5s"%(mol.name()[i])+" in residue: %5s"%(mol.resname()[i])+")\n")
                B_neutron_vacuum[i] = float(pt.D.neutron.b_c)*0.1
                mol.element()[i] = 'D'
    for i in xrange(mvars.number_deuteration_regions):
        fraction_deuterated = mvars.fraction_deuterated_array[i]
        #
        mask_inert_Hs_regional = [x and y for (x, y) in zip(mask_inert_Hs, mask_deuteration_regions[i])]
        count_target_deuterated_inert_Hs = round(mask_inert_Hs_regional.count(1)*fraction_deuterated)
        indices_inert_Hs = np.where(np.array(mask_inert_Hs_regional) == 1)[0]
        np.random.shuffle(indices_inert_Hs)
        indices_deuterated_inert_Hs = indices_inert_Hs[:count_target_deuterated_inert_Hs]
        of1.write("Fraction D2O: %.2f"%(D2O_fraction)+"; deuteration region: %3d"%(i+1)+"; fraction deuterated: %.2f"%(fraction_deuterated)+"; number deuterated H: %5d"%(count_target_deuterated_inert_Hs)+"\n")
        for i in xrange(natoms):
            if i in indices_deuterated_inert_Hs:
                of1.write("Deuterate atom #%5d"%(i+1)+" (atom name: %5s"%(mol.name()[i])+" in residue: %5s"%(mol.resname()[i])+")\n")
                B_neutron_vacuum[i] = float(pt.D.neutron.b_c)*0.1
                mol.element()[i] = 'D'
    mol.write_pdb(o_pdbfile, 0, 'w')

def preprocess_neutron_bs(mol, B, b_solvent, Btot, solvent_volume, radii):
    '''
    Functionality: preprocess the B values
    Input: mol--sasmol object
           B--shape==(3,natoms)
           Btot--sum of Bvac - Bsolv(corrected) over all atoms; used for Rg calculation
    Return: Btot

    '''
    natoms = B.shape[1]
    for i in xrange(natoms):
        B[1, i] = b_solvent * ((4./3.*np.pi*radii[i]*radii[i]*radii[i])/solvent_volume)
        B[2, i] = B[0, i] - B[1, i]
        Btot = Btot + B[2, i]
#        print i, B[2,i], Btot

# Need to return Btot because it is a scalar and its modified value isn't passed back to the calling routine.
# The modified B matrix is returned to the calling routine.

    return Btot

def preprocess_xray_bs(B, Q, b_solvent, Btot, solvent_volume, radii):
    '''
    Functionality: preprocess the B values
    Input: mol--sasmol object
           B--shape==(3,Nq,natoms)
           Btot--sum of Bvac - Bsolv(corrected) over all atoms for q = 0; used for Rg calculation
    Return: Btot(q=0)

    '''
    num_Q = B.shape[1]
    natoms = B.shape[2]
    for i in xrange(num_Q):
        for k in xrange(natoms):
            # B[1,i,k] = b_solvent * ((4./3.*np.pi*radii[k]*radii[k]*radii[k])/solvent_volume)
            B[1, i, k] = b_solvent[i] * ((4./3.*np.pi*radii[k]*radii[k]*radii[k])/solvent_volume)
            B[2, i, k] = B[0, i, k] - B[1, i, k]
            if i == 0:
                Btot = Btot + B[2, i, k]
#               print i,k,B[2,i,k],Btot

# Need to return Btot because it is a scalar and its modified value isn't passed back to the calling routine.
# The modified B matrix is returned to the calling routine.

    return Btot

def calc_cosm(mol, frame, Btot, B):
    '''
    Functionality: calculate center of scattering mass (contrast-dependent)
    Input:  mol--sasmol object
            frame--dcd frame
            B,Btot:  for neutrons, B = B[2] (complete) and Btot = sum of the B[2] values
                     for x-rays, B = B[2,0] (complete) and Btot = sum of the B[2,0] values (q = 0)
    Return: cosm

    '''

    natoms = mol.natoms()
    coor = mol.coor()
#    print 'natoms = ',natoms
#    print coor

#    print Btot
#    print B

    cmx = 0.0
    cmy = 0.0
    cmz = 0.0

    for k in range(natoms):
        cmx = cmx + (B[k]/Btot) * coor[frame, k, 0]
        cmy = cmy + (B[k]/Btot) * coor[frame, k, 1]
        cmz = cmz + (B[k]/Btot) * coor[frame, k, 2]
#    print cmx,cmy,cmz

    cosm = np.array([cmx, cmy, cmz], np.float)

#    print 'Center of scattering mass:  ',cosm

    return cosm

def calc_rg_contrast(mol, frame, Btot, B):
    '''
    Functionality: calculate scattering (contrast-dependent) radius of gyration, Rg
    Input:  mol--sasmol object
            frame--dcd frame
            B,Btot:  for neutrons, B = B[2] (complete) and Btot = sum of the B[2] values
                     for x-rays, B = B[2,0] (complete) and Btot = sum of the B[2,0] values (q = 0)
    Return: Rg

    '''

    natoms = mol.natoms()
    coor = mol.coor()
#    print 'natoms = ',natoms
#    print coor

#    print Btot
#    print B

    xcoor = []
    ycoor = []
    zcoor = []

    for k in range(natoms):
        this_xcoor = coor[frame, k, 0]
        this_ycoor = coor[frame, k, 1]
        this_zcoor = coor[frame, k, 2]
        xcoor.append(this_xcoor)
        ycoor.append(this_ycoor)
        zcoor.append(this_zcoor)

    cosm = calc_cosm(mol, frame, Btot, B)


#    print 'Center of scattering mass: ', cosm
#    print xcoor, ycoor, zcoor

    Rg2 = 0.0

    for j in range(natoms):
        Rg2 = Rg2 + (B[j]/Btot) * ((xcoor[j]-cosm[0])**2 + (ycoor[j]-cosm[1])**2 + (zcoor[j]-cosm[2])**2)
#	print 'Rg2 = ',Rg2
    if Rg2 < 0:
        Rg = np.nan
#		print 'Rg2 was negative.'
    else:
        Rg = np.sqrt(Rg2)

#    print 'Rg(contrast): ', Rg

    return Rg


def prepare_sascalc_inputs(mol, mvars, scvars, output_folder):
    '''
    Functionality: deuterate the exchangeable H's based on the HD exhcange ratio in various ways
    Input: mol--sasmol object
           mvars--module variables
    Return: B: array of atomic scattering length for each atom with H's deutrated based on the input deuteration_option
            Q: array of Q-values for scattering intensity calculation
            R: array of distance values for distance distribution calculation
    '''
    # folder preparation
    xon = mvars.xon
    if xon in ['neutron', 'neutron_and_xray']:
        for idx_contrast in xrange(mvars.number_contrast_points):
            output_dir = os.path.join(output_folder, 'neutron_D2Op_%d'%mvars.D2O_percentage_array[idx_contrast])
            os.mkdir(output_dir)
    if xon in ['xray', 'neutron_and_xray']:
        for i, idx_contrast in enumerate(range(mvars.xray_number_contrast_points)): ## @NOTE to developers: the output folder name is currently in-consistent with neutron and xray
            # output_dir = os.path.join(output_folder,'xray_D2Op_%d'%mvars.xray_D2O_percentage_array[idx_contrast])
            # output_dir = os.path.join(output_folder,'xray_condition_%d'%(i+1))
            output_dir = os.path.join(output_folder, 'xray')
            os.mkdir(output_dir)

    # get Q,R
    number_q_values = mvars.number_q_values
    q_max = mvars.q_max
    Q = np.zeros(number_q_values)
    for i in xrange(number_q_values):
        Q[i] = i*(q_max/(number_q_values-1));
    number_r_values = mvars.number_r_values
    minmax = mol.calcminmax()
    r_max = np.linalg.norm(minmax[1]-minmax[0])
    R = np.zeros(number_r_values)
    for i in xrange(number_r_values):
        R[i] = i*(r_max/(number_r_values-1));
    scvars.Q = Q
    scvars.R = R

    # get initial Bs from sasmol object (parsing the PDB file)
    natoms = mol.natoms()
    if xon in ['neutron', 'neutron_and_xray']:
        Btot_neutron_array = np.zeros(mvars.number_contrast_points)
        B_neutron_vacuum = np.zeros(natoms)
        set_atomic_neutron_bs(mol, B_neutron_vacuum)
        B_neutron_array = np.zeros((3, mvars.number_contrast_points, natoms)) ## vacuum, solvent, and complete
        B_neutron_array[0] = np.tile(B_neutron_vacuum, (mvars.number_contrast_points, 1))
    else:
        B_neutron_array = np.zeros((3, 0, natoms)) ## vacuum, solvent, and complete
        Btot_neutron_array = np.zeros(0) ## set array to zeros if only xray

    if xon in ['xray', 'neutron_and_xray']:
        Btot_xray_array = np.zeros(mvars.xray_number_contrast_points) ##Btot will be calculated for q=0 only, so the array is one-dimensional
        B_xray_vacuum = np.zeros((number_q_values, natoms))
        set_atomic_xray_bs(mol, B_xray_vacuum, Q)
        B_xray_array = np.zeros((3, mvars.xray_number_contrast_points, number_q_values, natoms)) ## vacuum, solvent, and complete
        B_xray_array[0] = np.tile(B_xray_vacuum, (mvars.xray_number_contrast_points, 1, 1))
    else:
        B_xray_array = np.zeros((3, 0, number_q_values, natoms)) ## vacuum, solvent, and complete
        Btot_xray_array = np.zeros(0) ##set array to zeros if only neutron


    # deuterate Bs for each contrast point
    radii = np.zeros(natoms)
    set_atomic_radii(mol, radii, mvars.VDW_scaling_factor)
#    print 'radii: ',radii
    if xon in ['neutron', 'neutron_and_xray']:
        mask_active_Hs = find_active_Hs(mol)
        mask_inert_Hs = find_inert_Hs(mol)
        mask_exH_regions = []
        mask_deuteration_regions = []
        for basis_string in mvars.exH_basis_string_array_python:
            error, mask = mol.get_subset_mask(basis_string)
            mask_exH_regions.append(mask)
        for basis_string in mvars.deuterated_basis_string_array_python:
            error, mask = mol.get_subset_mask(basis_string)
            mask_deuteration_regions.append(mask)
        for idx_contrast in xrange(B_neutron_array.shape[1]):
            D2O_fraction = mvars.D2O_percentage_array[idx_contrast]*0.01
            output_dir = os.path.join(output_folder, 'neutron_D2Op_%d'%mvars.D2O_percentage_array[idx_contrast])
            o_pdbfile = os.path.join(output_dir, "D2Op_%d.pdb"%(int(D2O_fraction*100)))
            ofile = os.path.join(output_dir, "HDexchange_Info_D2Op_%d.txt"%(int(D2O_fraction*100)))
            ofile1 = os.path.join(output_dir, "Deuteration_Info_D2Op_%d.txt"%(int(D2O_fraction*100)))
            deuterate_neutron_bs(mol, B_neutron_array[0, idx_contrast], mask_active_Hs, mask_inert_Hs, mask_exH_regions, mask_deuteration_regions, mvars, D2O_fraction, o_pdbfile, ofile, ofile1)
            b_solvent_neutron = calculate_solvent_sl(D2O_fraction)
            Btot_neutron_array[idx_contrast] = preprocess_neutron_bs(mol, B_neutron_array[:, idx_contrast,:], b_solvent_neutron, Btot_neutron_array[idx_contrast], mvars.solvent_volume, radii)
#            print 'B: ',idx_contrast,B_neutron_array[:,idx_contrast,:]
#            print 'B(complete): ', idx_contrast, B_neutron_array[2, idx_contrast,:]
#            print 'Btot: ', idx_contrast, Btot_neutron_array[idx_contrast]

    if xon in ['xray', 'neutron_and_xray']:
        for idx_contrast in xrange(B_xray_array.shape[1]):
            b_solvent_xray = calculate_solvent_electron_sl(Q)
            Btot_xray_array[idx_contrast] = preprocess_xray_bs(B_xray_array[:, idx_contrast,:,:], Q, b_solvent_xray, Btot_xray_array[idx_contrast], mvars.solvent_volume, radii)
#            print 'Bx:',idx_contrast,B_xray_array[:,idx_contrast,:,:]
#            print 'Bx[q=0]: ',idx_contrast,B_xray_array[:,idx_contrast,0,:]
#            print 'Bx[complete,q=0]: ', idx_contrast, B_xray_array[2, idx_contrast, 0,:]
#            print 'Bxtot[q=0]:', idx_contrast, Btot_xray_array[idx_contrast]

    scvars.B_neutron_array = B_neutron_array
    scvars.Btot_neutron_array = Btot_neutron_array
    scvars.B_xray_array = B_xray_array
    scvars.Btot_xray_array = Btot_xray_array

    return

def write_sascalc_log_file(log_filename, idx_contrast, mvars, scvars, results, frame, app, flag_neutron=True):
    '''
    Functionality: write the sascalc log file in a consistent way with xtal2sas
    Input: mol--sasmol object
           mvars--module variables
           scvars--sascalc input variables
           app--module name
    Return: None
    '''
#    print 'writing logfile'
#   'xray'
#    print idx_contrast, scvars.Btot_xray_array[idx_contrast]
#    print scvars.B_xray_array[2,idx_contrast,0,:]
#   'neutron'
#    print idx_contrast, scvars.Btot_neutron_array[idx_contrast]
#    print scvars.B_neutron_array[2,idx_contrast]


    fout = open(log_filename, 'w')
    output_folder = os.path.join(mvars.runname, app)
    natoms = scvars.mol.natoms()
    mw = scvars.mol.calcmass()
    molrg = scvars.mol.calcrg(frame)
    molcom = scvars.mol.calccom(frame)
    minmax = scvars.mol.calcminmax_frame(frame)
    r_max = np.linalg.norm(minmax[1]-minmax[0])

    fout.write("#Structural Information:\n")
    fout.write("1. PDB input file = %s\n"%mvars.pdbfile)
    fout.write("2. Number of atoms = %d; Mw = %.2f\n"%(natoms, mw))
    fout.write("3. Dimensions = x: %.2f, %.2f, y: %.2f, %.2f, z: %.2f, %.2f\n"%(minmax[0][0], minmax[1][0], minmax[0][1], minmax[1][1], minmax[0][2], minmax[1][2]))
    fout.write("4. Maximum radial dimension Dmax = %.2f A (contrast-independent)\n"%r_max)
#   f.out.write("5. Box size = x: %.2f, y: %.2f, z: %.2f\n"%(minmax[1][0]-minmax[0][0],minmax[1][1]-minmax[0][1],minmax[1][2]-minmax[0][2]))
    fout.write("5. Molecular center of mass = x: %.2f, y: %.2f, z: %.2f (contrast-independent)\n"%(molcom[0], molcom[1], molcom[2]))
    fout.write("6. Molecular Rg = %f A (contrast-independent)\n"%molrg)
    fout.write("#Scattering Intensity:\n")
    if flag_neutron:
        I0 = mvars.I0_array[idx_contrast]
        D2O_fraction = mvars.D2O_percentage_array[idx_contrast]*0.01
        fout.write("7. Source = neutron\n")
        fout.write("8. D2O %% = %7.3f\n"%(mvars.D2O_percentage_array[idx_contrast]))
        output_dir = os.path.join(output_folder, 'neutron_D2Op_%d'%mvars.D2O_percentage_array[idx_contrast])
        fout.write("9. Non-exchangeable H deuteration file = %s\n"%(os.path.join(output_dir, "HDexchange_Info_D2Op_%d.txt"%(int(D2O_fraction*100)))))
        fout.write("10. Exchangeable H deuteration file = %s\n"%(os.path.join(output_dir, "Deuteration_Info_D2Op_%d.txt"%(int(D2O_fraction*100)))))
        fout.write("11. PDB output file with explicit D atoms = %s\n"%(os.path.join(output_dir, "D2Op_%d.pdb"%(int(D2O_fraction*100)))))
        fout.write("12. I(q) .vs. q file = %s\n"%(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.iq'))
        if mvars.golden_vector_method_option == "fixed":
            fout.write("13. Number of golden vectors used = %d\n"%mvars.number_golden_vectors)
        if mvars.golden_vector_method_option == "converge":
            fout.write("13: Convergence tolerance used = %f\n"%mvars.golden_vector_method_converge_tolerance)
            fout.write("    Converged number of golden vectors (for complete scattering profile) = %d\n"%results.Ngv_converged_neutron_array[2, idx_contrast])
        fout.write("14. Io = %10.3f\n"%I0)
        scatcom = calc_cosm(scvars.mol, frame, scvars.Btot_neutron_array[idx_contrast], scvars.B_neutron_array[2, idx_contrast,:])
        fout.write("15. Center of Mass = %.2f, y: %.2f, z: %.2f (contrast-dependent)\n"%(scatcom[0], scatcom[1], scatcom[2]))
        scatrg = calc_rg_contrast(scvars.mol, frame, scvars.Btot_neutron_array[idx_contrast], scvars.B_neutron_array[2, idx_contrast,:])
        fout.write("16. Rg = %f A (contrast-dependent)\n"%scatrg)
    else:
        I0 = mvars.xray_I0_array[idx_contrast]
        fout.write("7. Source = x-ray\n")
        output_dir = os.path.join(output_folder, 'xray')
        fout.write("8. I(q) .vs. q file = %s\n"%(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.iq'))
        if mvars.golden_vector_method_option == "fixed":
            fout.write("9. Number of golden vectors used = %d\n"%mvars.number_golden_vectors)
        if mvars.golden_vector_method_option == "converge":
            fout.write("9: Convergence tolerance used = %f\n"%mvars.golden_vector_method_converge_tolerance)
            fout.write("   Converged number of golden vectors (for complete scattering profile) = %d\n"%results.Ngv_converged_xray_array[2, idx_contrast])
        fout.write("10. Io = %10.3f\n"%I0)
        # Center of Scattering Mass and Scattering Rg are calculated for q=0 only.
        scatcom = calc_cosm(scvars.mol, frame, scvars.Btot_xray_array[idx_contrast], scvars.B_xray_array[2, idx_contrast, 0,:])
        fout.write("11. Center of Mass = x: %.2f, y: %.2f, z: %.2f (contrast-dependent)\n"%(scatcom[0], scatcom[1], scatcom[2]))
        scatrg = calc_rg_contrast(scvars.mol, frame, scvars.Btot_xray_array[idx_contrast], scvars.B_xray_array[2, idx_contrast, 0,:])
        fout.write("12. Rg = %f (contrast-dependent)\n"%scatrg)


#   fout.write("17. Output: P(R) function = %s\n"%(os.path.join(mvars.runname,app,mvars.runname+'_'+str(frame+1).zfill(5))+'.pr'))
#   fout.write("18. Scattering Dmax = value (contrast-dependent)")
#   fout.write("19. Rg from P(r) = value (contrast-dependent)")

    fout.close()

def save_sascalc_outputs(mvars, scvars, results, frame, app):
    '''
    Functionality: save sascalc results to the output file
    Input: mol--sasmol object
           mvars--module variables
           scvars--sascalc input variables
           app--module name
    Return: None
    '''
    # neutron
    Iq_neutron_array = results.Iq_neutron_array
    Ncontrast_neutron = Iq_neutron_array.shape[1]
    if Ncontrast_neutron:
        for idx_contrast in xrange(Ncontrast_neutron):
            output_dir = os.path.join(mvars.runname, app, 'neutron_D2Op_%d'%mvars.D2O_percentage_array[idx_contrast])
            fout = open(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.iq', 'w')
            i0 = Iq_neutron_array[2, idx_contrast, 0]
            I0 = mvars.I0_array[idx_contrast]
            for q, i_vac, i_sol, i_com in zip(scvars.Q, Iq_neutron_array[0, idx_contrast], Iq_neutron_array[1, idx_contrast], Iq_neutron_array[2, idx_contrast]):
                fout.write("%8.6f  %8.6f  %8.6f  %15.5f %15.5f %15.5f\n"%(q, i_com/i0*I0, 0.0, i_com, i_vac, i_sol))
            fout.close()
            if mvars.complex_amplitudes:
                fout = open(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'_complex_amplitude.ciq', 'w')
                fout.write("# Q   (Ireal_complete, Iimag_complete)   (Ireal_vacuum, Iimag_vacuum)   (Ireal_solvent, Iimag_solvent)\n")
                for q, ireal_vac, ireal_sol, ireal_com, iimag_vac, iimag_sol, iimag_com in zip(scvars.Q, Iq_neutron_array[3, idx_contrast], Iq_neutron_array[4, idx_contrast], Iq_neutron_array[5, idx_contrast], Iq_neutron_array[6, idx_contrast], Iq_neutron_array[7, idx_contrast], Iq_neutron_array[8, idx_contrast]):
                    fout.write("%8.6f  (%15.5f, %15.5f)   (%15.5f, %15.5f)   (%15.5f, %15.5f)\n"%(q, ireal_com, iimag_com, ireal_vac, iimag_vac, ireal_sol, iimag_sol))
                fout.close()
            log_filename = os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.log'
            write_sascalc_log_file(log_filename, idx_contrast, mvars, scvars, results, frame, app)
    # xray
    Iq_xray_array = results.Iq_xray_array
    Ncontrast_xray = Iq_xray_array.shape[1]
    if Ncontrast_xray:
        for i, idx_contrast in enumerate(range(Ncontrast_xray)): ## @NOTE to developers: xray ouput folder names currently in-conistent with those of neutrons
            output_dir = os.path.join(mvars.runname, app, 'xray')
            fout = open(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.iq', 'w')
            i0 = Iq_xray_array[2, idx_contrast, 0]
            I0 = mvars.xray_I0_array[idx_contrast]
            for q, i_vac, i_sol, i_com in zip(scvars.Q, Iq_xray_array[0, idx_contrast], Iq_xray_array[1, idx_contrast], Iq_xray_array[2, idx_contrast]):
                fout.write("%8.6f  %8.6f  %8.6f  %15.5f %15.5f %15.5f\n"%(q, i_com/i0*I0, 0.0, i_com, i_vac, i_sol))
            fout.close()
            if mvars.complex_amplitudes:
                fout = open(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'_complex_amplitude.ciq', 'w')
                fout.write("# Q   (Ireal_complete, Iimag_complete)   (Ireal_vacuum, Iimag_vacuum)   (Ireal_solvent, Iimag_solvent)\n")
                for q, ireal_vac, ireal_sol, ireal_com, iimag_vac, iimag_sol, iimag_com in zip(scvars.Q, Iq_xray_array[3, idx_contrast], Iq_xray_array[4, idx_contrast], Iq_xray_array[5, idx_contrast], Iq_xray_array[6, idx_contrast], Iq_xray_array[7, idx_contrast], Iq_xray_array[8, idx_contrast]):
                    fout.write("%8.6f  (%12.5f, %12.5f)   (%12.5f, %12.5f)   (%12.5f, %12.5f)\n"%(q, ireal_com, iimag_com, ireal_vac, iimag_vac, ireal_sol, iimag_sol))
                fout.close()
            log_filename = os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.log'
            write_sascalc_log_file(log_filename, idx_contrast, mvars, scvars, results, frame, app, 0)


def ite(f, v=[]):
    for item in f.items():
        if  isinstance(f[item[0]], h5py.Dataset):
            v.append(f[item[0]].parent.name)
        else:
            ite(f[item[0]], v)
    return v

def reoganize_h5(output_folder):
    fin = os.path.join(output_folder, "results_raw.h5")
    fin = os.path.join(output_folder, "results.h5")
    v = ite(fin)
    #import pprint
    #pprint.pprint(v)
    data = {}
    Nq = f['/Q'].shape[0]
    Q = f['/Q'].value
    s=set([])
    sframe=set([])
    ll = []
    for i in v:
        #print i
        #if i=='/': print f['/Q'].value
        #else: print f[i+'/data'].value
        d=i.split('/')
        if len(d)==6:
            ll.append(d)
            print d
            s.add(d[2]+', '+d[4])
            sframe.add(d[3])

    fo = h5py.File(fout,'w')
    for i in s:
        ii = i.split(', ')
        grp = fo.create_group(i)
        for frame in sframe:
            grp.create_group(frame)
            dset = fo.create_dataset('/'+i+'/'+frame+'/data', (Nq,6), dtype='float64')
            dset[:,0] = Q
            for l in ll:
                if l[2] == ii[0] and l[4]== ii[1]:
                    if l[5] == 'normalized':
                        dset[:,1]=f['/'.join(l)+'/data'].value
                    if l[5] == 'error':
                        dset[:,2]=f['/'.join(l)+'/data'].value
                    elif l[5] == 'complete':
                        dset[:,3]=f['/'.join(l)+'/data'].value
                    elif l[5] == 'vacuum':
                        dset[:,4]=f['/'.join(l)+'/data'].value
                    elif l[5] == 'solvent':
                        dset[:,5]=f['/'.join(l)+'/data'].value
            dset.attrs['col 1'] = 'Q'
            dset.attrs['col 2'] = 'Iq normalized'
            dset.attrs['col 3'] = 'Iq error'
            dset.attrs['col 4'] = 'Iq complete'
            dset.attrs['col 5'] = 'Iq vacuum'
            dset.attrs['col 6'] = 'Iq solvent'
    fo.close()
