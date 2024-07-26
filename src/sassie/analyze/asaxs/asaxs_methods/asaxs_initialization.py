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
import os
import string
import locale
import numpy
import io

import sasmol.system as system

local_debug = True

def process_pdb_file(other_self):

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    log.debug('in read_structure_files')
    mvars = other_self.module_variables
    avars = other_self.asaxs_variables

    mol = avars.mol

    print("in process_pdb_file: natoms = ", mol.natoms())

    index = mol.index()
    name = mol.name()
    resname = mol.resname()
    segname = mol.segname()

    x = mol.coor()[0,:,0]
    y = mol.coor()[0,:,1]
    z = mol.coor()[0,:,2]

    # matlab variables
    
    if local_debug:
        output_filename = 'temp_pdb_output.txt'
        outfile = open(output_filename, 'w')
        outfile.write('Atom Data:\n') 

    labels = []
    atom_set = ['C', 'N', 'O', 'S', 'P', 'H', 'ZN', 'AU', 'AN', 'TB', 'PT', 'TT']

    for i in range(mol.natoms()):
        atom_type = 0
        for atom_set_item in atom_set:
            n_char = len(atom_set_item)
        
            # Manually remove digits from the name
            alt_name = ''.join([char for char in name[i][:n_char] if not char.isdigit()])
        
            if ''.join([char for char in name[i] if not char.isdigit()])[:n_char] == atom_set_item:
                atom_type = atom_set.index(atom_set_item) + 1  # +1 to match MATLAB's 1-based indexing
                break
        if atom_type == 0:
            raise ValueError("Error: Atom type not found.")

        #Replace with alt_name
        if segname[i] == "LBL":
            labels.append([alt_name, atom_type, resname[i], x[i], y[i], z[i]])
        elif atom_type != 6: # Skip hydrogen atoms
            if local_debug:
                outfile.write(f'{alt_name} {atom_type} {resname[i]} {x[i]} {y[i]} {z[i]}\n')

    if local_debug:
        outfile.write('Label Data:\n')
        for i in range(len(labels)):
            for j in range(len(labels[i])):
                outfile.write(f'{labels[i][j]} ')
            outfile.write('\n')
        outfile.close()

    print("in process_pdb_file: atom_set = ", atom_set)

    return

