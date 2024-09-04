# this should, eventually, be a sasmol getter and live elsewhere

def sequence(residlist, resnamelist):
    marker = 0
    seq = []
    for i in range(len(residlist)):
        if marker != residlist[i]:
            seq.append(resnamelist[i])
            marker = residlist[i]
    return seq


def FASTA_sequence(infile):
    seq = []
    o = open(infile)
    for line in o:
        if line[0] != ">":
            for j in line:
                j = j.upper()
                if j in dna_sl() or j in rna_sl() or j in protein_sl():
                    seq.append(j)
    return seq


# Scattering properties for elements, proteins and nucleotides; need to update in sasproperties!
# Scattering lengths are in units of 10^-12 cm

def element_sl():

    # element name : [MW, vol Ang^3, eSL Ang, nSL Ang]
    # approx volumes are per Whitten 2008

    atomic_scattering = {
        'D': [2.014, 5.15, 0.282, 0.6671],
        'H': [1.008, 5.15, 0.282, -0.3741],
        'C': [12.01, 16.44, 1.692, 0.6646],
        'N': [14.01, 2.49, 1.974, 0.9360],
        'O': [16.00, 9.130, 2.256, 0.5803],
        'Na': [22.99, 4.45, 3.102, 0.3630],
        'Mg': [24.31, 1.56, 3.384, 0.5375],
        'K': [39.10, 11.01, 5.358, 0.3670],
        'Ca': [40.08, 4.19, 5.640, 0.4700],
        'Cl': [35.45, 24.84, 4.794, 0.9577],
        'Br': [79.90, 31.54, 9.870, 0.6795],
        'I': [126.9, 44.6, 14.946, 0.5280],
        'P': [30.97, 3.37, 4.230, 0.5130],
        'S': [32.07, 26.09, 4.512, 0.2847],
        'Fe': [55.85, 7.99, 7.332, 0.9450],
        'Co': [58.93, 7.99, 7.614, 0.2490],
        'Ni': [58.69, 8.18, 7.896, 1.0300],
        'Cu': [63.55, 8.78, 8.178, 0.7718],
        'Zn': [65.39, 9.85, 8.460, 0.5680],
        'Au': [196.97, 14.75, 22.278, 0.7630]}

    return atomic_scattering


def nucleotide_sl():

    # nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
    # vol = MW*(0.586 cm^3/g)/NA

    # electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
    # molecule formulae are for nucleotides as embedded in a chain, per Jacrot
    # 1976

    residue_scattering = {
        'DA': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
        'DT': [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12],  # 303
        'DG': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],  # 328
        'DC': [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11],  # 288
        'U': [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
        'A': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
        'G': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
        'C': [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11]}

    return residue_scattering


def dna_sl():

    # nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
    # vol = MW*(0.56 cm^3/g)/NA, where partial specific volume of 0.56 is per Hearst 1962, JMB 4, 415-417
    # psv depends on pH and salt, so there is a range between ~0.55-0.59

    # electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
    # molecule formulae are for nucleotides as embedded in a chain, per Jacrot
    # 1976

    residue_scattering = {
        'DA': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
        'DT': [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12],  # 303
        'DG': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],  # 328
        'DC': [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11],
        'ADE': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
        'THY': [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12],  # 303
        'GUA': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],  # 328
        'CYT': [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11],
        'A': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
        'T': [304.0, 294.9, 44.2, 8.61,  9.65,  1, 12],  # 303
        'G': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],  # 328
        'C': [289.0, 280.3, 41.9, 8.68,  10.77, 2, 11]}  # 288

    return residue_scattering


def rna_sl():

    # nucleotide name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang, #exchngH, #totalH]
    # vol = MW*(0.55 cm^3/g)/NA, where the partial specific volume of 0.55 is from Chien, et al. 2004,
    # Biochemistry 43, 1950-1962.  psv depends on pH and salt and whether RNA is ss or ds.  This value
    # is at the high end of the range for ssRNA, ~0.47-0.55.  psv is larger
    # for dsRNA, ~0.57.

    # electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
    # molecule formulae are for nucleotides as embedded in a chain, per Jacrot
    # 1976

    residue_scattering = {
        'U': [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
        'A': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
        'G': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
        'C': [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11],
        'URA': [305.0, 296.9, 44.2, 9.28,  11.36, 2, 10],
        'ADE': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
        'GUA': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
        'CYT': [304.0, 295.9, 44.2, 8.68,  10.77, 2, 11]}

    return residue_scattering


def protein_sl():

    # aa residue name : [MW, vol Ang^3, eSL Ang, SLprot Ang, SLdeut Ang,
    # #exchngH, #totalH]

    # electron scattering lengths = Thompson radius (2.82 x 10^-13 cm) * atomic number
    # molecule formulae are for residues embedded in a chain, per Jacrot 1976

    residue_scattering = {
        'ALA': [71.1,  88.6,  10.7, 1.645, 2.686, 1, 5],
        'ARG': [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
        'ASP': [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
        'ASN': [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
        'CYS': [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
        'GLU': [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
        'GLN': [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
        'GLY': [57.1,   60.1,  8.5, 1.728, 2.769, 1, 3],
        'HSD': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
        'HIS': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
        'HSE': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
        'HSP': [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
        'ILE': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
        'LEU': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
        'LYS': [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
        'MET': [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
        'PHE': [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
        'PRO': [97.1,  112.7, 14.7, 2.227, 2.227, 0, 7],
        'SER': [87.1,  89.0,  13.0, 2.225, 4.308, 2, 5],
        'THR': [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
        'TRP': [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
        'TYR': [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
        'VAL': [99.1,  140.0, 15.3, 1.478, 2.520, 1, 9],
        'A': [71.1,  88.6,  10.7, 1.645, 2.686, 1, 5],
        'R': [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
        'D': [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
        'N': [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
        'C': [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
        'E': [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
        'Q': [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
        'G': [57.1,   60.1,  8.5, 1.728, 2.769, 1, 3],
        'H': [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
        'I': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
        'L': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
        'K': [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
        'M': [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
        'F': [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
        'P': [97.1,  112.7, 14.7, 2.227, 2.227, 0, 7],
        'S': [87.1,  89.0,  13.0, 2.225, 4.308, 2, 5],
        'T': [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
        'W': [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
        'Y': [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
        'V': [99.1,  140.0, 15.3, 1.478, 2.520, 1, 9]}

    return residue_scattering

if __name__ == "__main__":

    import sasmol.sasmol as sasmol
    
    protein_dict = protein_sl()
    nuc_dict = nucleotide_sl()
    dna_dict = dna_sl()
    rna_dict = rna_sl()
    atomic_dict = element_sl()

    print(atomic_dict['H'])
    print(atomic_dict['D'])
    print(atomic_dict['O'])
    print(protein_dict['A'])
    print(dna_dict['THY'])
    print(rna_dict['URA'])
    print(nuc_dict['DG'])

    seq = FASTA_sequence('protein_sequence.txt')
    print(seq)

    m1 = sasmol.SasMol(0)
    m1.read_pdb('hiv1_gag.pdb')

    resids = m1.resid()
    resnames = m1.resname()

    seq = sequence(resids, resnames)
    print(seq)


    
