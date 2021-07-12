import sasmol.sasmol as sasmol
import sys

pdbfile = "two_lysozyme.pdb"

system = sasmol.SasMol(0)

system.read_pdb(pdbfile)

system.initialize_children()

number_of_segments = len(system.segnames())

print 'number of segments = ', number_of_segments

print 'segnames = ',system.segnames()

all_masks = []
all_molecules = []

frame = 0

### create new molecules from each segname in the system

for this_segname in system.segnames():
    basis = 'segname[i] == "' + this_segname + '"'
    error, this_mask = system.get_subset_mask(basis)
    if len(error) > 0:
        print 'ERROR : ' + str(error)
        sys.exit()

    all_masks.append(this_mask)    
    this_molecule = sasmol.SasMol(0)        
    system.copy_molecule_using_mask(this_molecule,this_mask,frame)

    ### write out the individual molecules as a test

    this_molecule.write_pdb("test_"+this_segname+".pdb",frame,"w")
    all_molecules.append(this_molecule)



### OUR "simulation"
### move second molecule to x,y,z = [90,90,90]

all_molecules[1].moveto(frame,[90,90,90])

all_molecules[1].write_pdb("test_moved_0002.pdb",frame,"w")

###
### END of "simulation"


### now put the coordinates back into the full "system" and write to a pdb file

error = system.set_coor_using_mask(all_molecules[1],frame,all_masks[1])

system.write_pdb('new_'+pdbfile, frame, "w")


