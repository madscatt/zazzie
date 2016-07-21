import string

class psf_data():
    pass

def parse_psf_file(psf_file_data, psf_data):

#    psf_keywords = ['PSF', '!NTITLE', '!NATOM', '!NBOND:', '!NTHETA:', '!NPHI:', '!NIMPHI:', '!NDON:', '!NACC:', '!NNB', '!NGRP']

    psf_data.header = [] ; psf_data.natom = [] ; psf_data.nbond = [] ; psf_data.ntheta = [] ; psf_data.nphi = [] 
    psf_data.nimphi = [] ; psf_data.ndon = [] ; psf_data.nacc = [] ; psf_data.nnb = [] ; psf_data.ngrp = []
  
    header_flag = True 
    natom_flag = nbond_flag = ntheta_flag = nphi_flag = nimphi_flag = ndon_flag = nacc_flag = nnb_flag = ngrp_flag = False
     
    for line in psf_file_data:
        this_line = string.split(line)

        if "!NATOM" in this_line:
            header_flag = False
            natom_flag = True
            psf_data.natom.append(line)
            continue

        elif "!NBOND:" in this_line:
            natom_flag = False
            nbond_flag = True
            psf_data.nbond.append(line)
            continue

        elif "!NTHETA:" in this_line:
            nbond_flag = False
            ntheta_flag = True
            psf_data.ntheta.append(line)
            continue

        elif "!NPHI:" in this_line:
            ntheta_flag = False
            nphi_flag = True
            psf_data.nphi.append(line)
            continue

        elif "!NIMPHI:" in this_line:
            nphi_flag = False
            nimphi_flag = True
            psf_data.nimphi.append(line)
            continue

        elif "!NDON:" in this_line:
            nimphi_flag = False
            ndon_flag = True
            psf_data.ndon.append(line)
            continue

        elif "!NACC:" in this_line:
            ndon_flag = False
            nacc_flag = True
            psf_data.nacc.append(line)
            continue

        elif "!NNB" in this_line:
            nacc_flag = False
            nnb_flag = True
            psf_data.nnb.append(line)
            continue

        elif "!NGRP" in this_line:
            nnb_flag = False
            ngrp_flag = True
            psf_data.ngrp.append(line)
            continue
    
        if header_flag:
            psf_data.header.append(line)
        elif natom_flag: 
            psf_data.natom.append(line)
        elif nbond_flag: 
            psf_data.nbond.append(line)
        elif ntheta_flag: 
            psf_data.ntheta.append(line)
        elif nphi_flag: 
            psf_data.nphi.append(line)
        elif nimphi_flag: 
            psf_data.nimphi.append(line)
        elif ndon_flag: 
            psf_data.ndon.append(line)
        elif nacc_flag: 
            psf_data.nacc.append(line)
        elif nnb_flag: 
            psf_data.nnb.append(line)
        elif ngrp_flag: 
            psf_data.ngrp.append(line)

    return


def get_group_psf(other_self,new_psf_file_name, psf_data, group_mol):

    log = other_self.log 
    mvars = other_self.mvars


    new_file = open(new_psf_file_name,'w')

    new_number = 1

    index_map = {}

    for index in group_mol.index():
        index_map[index] = new_number
        new_number += 1

#    psf_keywords = ['PSF', '!NTITLE', '!NATOM', '!NBOND:', '!NTHETA:', '!NPHI:', '!NIMPHI:', '!NDON:', '!NACC:', '!NNB', '!NGRP']

#    psf_data.header = [] ; psf_data.natom = [] ; psf_data.nbond = [] ; psf_data.ntheta = [] ; psf_data.nphi = [] 
#    psf_data.nimphi = [] ; psf_data.ndon = [] ; psf_data.nacc = [] ; psf_data.nnb = [] ; psf_data.ngrp = []

    for line in psf_data.header:
        new_file.write(line)

    flag = True
    for line in psf_data.natom:
        if flag:
            flag = False
            st = '{:8d} {:7s}'.format(group_mol.natoms(),'!NATOM') + '\n'      
            new_file.write(st)
        else:
            try: 
                this_index = int(string.split(line)[0])
                if (this_index in index_map):
                    st = '{:8d}{:62s}'.format(index_map[this_index],line[8:])
                    new_file.write(st)
            except:
                pass

#
#         1         2        3         4         5         6         7
#1234567890123456789012345678012345678901234567890123456789012345678901
#   19668 !NATOM
#       1 FIXL 1    GLN  N    NH3   -0.300000       14.0070           0
#    5252 !NATOM 
#       1 FIXL 1    GLN  N    NH3   -0.300000       14.0070           0
#

    # count bonds
    number_of_bonds = 0
              
    pairs = []
    for line in psf_data.nbond:                     
        this_line = string.split(line)
        if "!NBOND:" not in this_line and len(this_line) > 1:
            try:
                local_number_of_pairs = len(this_line)/2
                j = 0
                for i in xrange(local_number_of_pairs):
                    if (int(this_line[j]) in index_map) and (int(this_line[j+1]) in index_map):
                        pairs.append([index_map[int(this_line[j])],index_map[int(this_line[j+1])]])
                        number_of_bonds += 1
                    j += 2
            except:
                pass          

    log.debug('found '+str(number_of_bonds)+' bonds')
    
    st = '\n{:8d} {:13s}'.format(number_of_bonds,'!NBOND: bonds') + '\n'      
    new_file.write(st)

    count = 0
    st = ''
    for pair in pairs:
        st += '{:8d}{:8d}'.format(pair[0],pair[1])  
        count += 1
        if count == 4:
            new_file.write(st+'\n')
            count = 0 ; st = ''
    
    if st!='': new_file.write(st+'\n') ## @NOTE to ZHL: this is a temporary fix

    # count angles
    number_of_angles = 0
    
    triplets = []
    for line in psf_data.ntheta:                     
        this_line = string.split(line)
        if "!NTHETA:" not in this_line:
            try:
                local_number_of_triplets = len(this_line)/3
                j = 0
                for i in xrange(local_number_of_triplets):
                    if (int(this_line[j]) in index_map)  and (int(this_line[j+1]) in index_map) and (int(this_line[j+2]) in index_map):
                        triplets.append([index_map[int(this_line[j])],index_map[int(this_line[j+1])],index_map[int(this_line[j+2])]])
                        number_of_angles += 1
                    j += 3
            except:
                pass          
    
    log.debug('found '+str(number_of_angles)+' angles')
    
    st = '\n{:8d} {:8s}'.format(number_of_angles,'!NTHETA:') + '\n'      
    new_file.write(st)

    count = 0
    st = ''
    for triplet in triplets:
        st += '{:8d}{:8d}{:8d}'.format(triplet[0],triplet[1],triplet[2])  
        count += 1
        if count == 3:
            new_file.write(st+'\n')
            count = 0 ; st = ''
    
    if st!='': new_file.write(st+'\n') ## @NOTE to ZHL: this is a temporary fix

    # count dihedrals
    number_of_dihedrals = 0
    
    quartets = []
    for line in psf_data.nphi:                     
        this_line = string.split(line)
        if "!NPHI:" not in this_line:
            try:
                local_number_of_quartets = len(this_line)/4
                j = 0
                for i in xrange(local_number_of_quartets):
                    if (int(this_line[j]) in index_map) and (int(this_line[j+1]) in index_map) and (int(this_line[j+2]) in index_map) and (int(this_line[j+3]) in index_map):
                        quartets.append([index_map[int(this_line[j])],index_map[int(this_line[j+1])],index_map[int(this_line[j+2])],index_map[int(this_line[j+3])]])
                        number_of_dihedrals += 1
                    j += 4
            except:
                pass          
    
    log.debug('found '+str(number_of_dihedrals)+' dihedrals')
    
    st = '\n{:8d} {:6s}'.format(number_of_dihedrals,'!NPHI:') + '\n'      
    new_file.write(st)

    count = 0
    st = ''
    for quartet in quartets:
        st += '{:8d}{:8d}{:8d}{:8d}'.format(quartet[0],quartet[1],quartet[2],quartet[3])  
        count += 1
        if count == 2:
            new_file.write(st+'\n')
            count = 0 ; st = ''
    
    if st!='': new_file.write(st+'\n') ## @NOTE to ZHL: this is a temporary fix

    # count impropers
    number_of_impropers = 0
    
    quartets = []
    for line in psf_data.nimphi:                     
        this_line = string.split(line)
        if "!NIMPHI:" not in this_line:
            try:
                local_number_of_quartets = len(this_line)/4
                j = 0
                for i in xrange(local_number_of_quartets):
                    if (int(this_line[j]) in index_map) and (int(this_line[j+1]) in index_map) and (int(this_line[j+2]) in index_map) and (int(this_line[j+3]) in index_map):
                        quartets.append([index_map[int(this_line[j])],index_map[int(this_line[j+1])],index_map[int(this_line[j+2])],index_map[int(this_line[j+3])]])
                        number_of_impropers += 1
                    j += 4
            except:
                pass          
    
    log.debug('found '+str(number_of_impropers)+' impropers')
    
    st = '\n{:8d} {:8s}'.format(number_of_impropers,'!NIMPHI:') + '\n'      
    new_file.write(st)

    count = 0
    st = ''
    for quartet in quartets:
        st += '{:8d}{:8d}{:8d}{:8d}'.format(quartet[0],quartet[1],quartet[2],quartet[3])  
        count += 1
        if count == 2:
            new_file.write(st+'\n')
            count = 0 ; st = ''
    
    if st!='': new_file.write(st+'\n') ## @NOTE to ZHL: this is a temporary fix

    new_file.write('\n')

#    psf_data.nimphi = [] ; psf_data.ndon = [] ; psf_data.nacc = [] ; psf_data.nnb = [] ; psf_data.ngrp = []
    for line in psf_data.ndon:
        new_file.write(line)

    for line in psf_data.nacc:
        new_file.write(line)
   
    count = 0
#12345678
#       0       0       0       0       0       0       0       0  

    for line in psf_data.nnb:
        try:
            this_line = string.split(line)
            if "!NNB:" not in this_line and this_line[3] == '0':
                sum = len(this_line)
                count += sum
                if value in this_line != '0':
                    log.error('ERROR: NON ZERO VALUE FOUND FOR NNB')
                    log.error('ERROR: NOT CONFIGURED TO HANDLE THIS CASE')
                    log.error('ERROR: QUITTING NOW')
                    sys.exit()
            #else:
            #    new_file.write(line)     
        except:
            pass
   
    new_file.write("\n       0 !NNB\n\n")  
  
    st = '' 
    count = 0
    total_count = 0
    for number in xrange(group_mol.natoms()):
        st += '{:8d}'.format(0)
        count += 1 ; total_count += 1
        if count == 8:
            new_file.write(st+'\n')
            count = 0 ; st = ''
   
    new_file.write("\n")
             
    log.debug('wrote '+str(total_count)+' NNB entries')             
    
    for line in psf_data.ngrp:
        new_file.write(line)
    
    new_file.write("\n")
    
    new_file.close()


    return



if __name__ == "__main__":

    import sasmol.sasmol as sasmol
   
    frame = 0

    pdb_file_name = 'asa2.pdb'
    psf_file_name = 'asa2.psf'

    psf_file_data = open(psf_file_name).readlines()
    
    psf_data = psf_data()

    parse_psf_file(psf_file_data, psf_data)

    m = sasmol.SasMol(0)
    m.read_pdb(pdb_file_name)

    basis = 'name[i] == "N" or name[i] == "CA" or name[i] == "C" or name[i] == "O"'

    basis = 'segname[i] == "FIXL" or segname[i] == "FIXH"'
    error, mask = m.get_subset_mask(basis) 
    group_mol = sasmol.SasMol(0)
    error = m.copy_molecule_using_mask(group_mol,mask,frame) 
    new_psf_file_name = 'mab_1.psf'
#    new_psf_file_name = 'backbone.psf'
    get_group_psf(new_psf_file_name, psf_data, group_mol)

    index = []
    for i in xrange(group_mol.natoms()):
        index.append(i+1)
    group_mol.setIndex(index)
    group_mol.write_pdb("mab_1.pdb",frame,"w")
#    group_mol.write_pdb("backbone.pdb",frame,"w")
    
    basis = 'segname[i] == "FIXM" or segname[i] == "FIXK"'
    error, mask = m.get_subset_mask(basis) 
    group_mol = sasmol.SasMol(0)
    error = m.copy_molecule_using_mask(group_mol,mask,frame) 
    new_psf_file_name = 'mab_2.psf'
    get_group_psf(new_psf_file_name, psf_data, group_mol)

    index = []
    for i in xrange(group_mol.natoms()):
        index.append(i+1)
    group_mol.setIndex(index)
    group_mol.write_pdb("mab_2.pdb",frame,"w")
            

