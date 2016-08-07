#!/Library/Frameworks/Python.framework/Versions/7.1/bin/python

import prdf
import numpy
#from sassie.sasmol import sasmol
#from sassie.simulate.rigid_body.two_body_grid.construct import construct_lmn
import sasmol.sasmol as sasmol
from construct import construct_lmn
import pprint

def Error(msg, error=[]):
	print(msg)
	print(error)
	exit(1)

if __name__=="__main__":
	# paramters to be modified
	pdbFile = "./data/rna.pdb"
	dcdFile = "./data/rna.dcd"
	solute_filter = 'moltype[i]=="rna"'
	solvent_filter = 'moltype[i]=="water"'

	# sasmol setup
	mol = sasmol.SasMol(0)
	mol.read_pdb(pdbFile)
	mol.read_dcd(dcdFile)
	tot_frames = mol.number_of_frames()

	# get the solute atom names
	error, solute_mask = mol.get_subset_mask(solute_filter)
	if (len(error)>0): Error("solute get_subset_mask error", error)
	solute_mol = sasmol.SasMol(1)
	error = mol.copy_molecule_using_mask(solute_mol, solute_mask, 0)
	if (len(error)>0): Error("solute copy_molecule_using_mask error", error)
	#solute_atom_types = solute_mol.name()
	#solute_atom_types = [i[0] for i in solute_mol.name()]
	#solute_atom_types = [t[0]+"-"+t[1] for t in zip(solute_mol.resname(),solute_mol.name())]
	solute_atom_types = []
	Q_1stline = True
	for name in solute_mol.name():
		if Q_1stline:
			solute_atom_types.append('HO')
			Q_1stline  = False
			continue
		atom_type = name[0]
		if (atom_type!='H'):
			solute_atom_types.append(atom_type)
			heavy_atom_type = atom_type
		else:
			solute_atom_types.append(atom_type+heavy_atom_type)
	# get the solute atom coordinates
	solute_atom_coors = []
	error, coors = mol.get_coor_using_mask(1, solute_mask) # get coordinates for frame 1 only
	if (len(error)>0): Error("solute get_coor_using_mask error", error)
	solute_atom_coors = numpy.array(coors[0], dtype=numpy.double)

	# find the surface atoms for the solute
	r = 1.8 # parameter to derive almn and blmn
	d = 2.5 # thickness of the surface layer
	rou = -15 # molecular 1 interior parameter
	eta = 1.0 # grid step size, 1.0-1.2
	N = 90 # Number of grid points
	qsurf = 1 # flag for surf atoms to be passed to construct_lmn
	percsurf = 0.3 # percentage of coverage to decide a surf atom in construct_lmn
	solute_atom_surf_mask = construct_lmn(solute_atom_coors, rou, N, r, d, eta, qsurf, percsurf)
	#solute_atom_surf_mask = numpy.array([0]*len(solute_mol.name()),numpy.int32)

	# get the solvent atom names
	error, solvent_mask = mol.get_subset_mask(solvent_filter)
	if (len(error)>0): Error("solvent get_subset_mask error", error)
	solvent_mol = sasmol.SasMol(3)
	error = mol.copy_molecule_using_mask(solvent_mol, solvent_mask, 0)
	if (len(error)>0): Error("solvent copy_molecule_using_mask error", error)
	solvent_atom_types = solvent_mol.name()
	# get the solvent atom coordinates
	solvent_atom_coors = []
	for frame in range(tot_frames):
		error, coors = mol.get_coor_using_mask(frame, solvent_mask)
		if (len(error)>0): Error("solvent get_coor_using_mask error", error)
		solvent_atom_coors.append(coors[0])
	solvent_atom_coors = numpy.array(solvent_atom_coors, dtype=numpy.double)
	
	# get the prdf parameters from the solvent atoms
	coor_minmax  = solvent_mol.calcminmax()
	(xmin,ymin,zmin)=coor_minmax[0]
	(xmax,ymax,zmax)=coor_minmax[1]
	cube_length = 1.0;
	nx = int((xmax-xmin)/cube_length)+1;
	ny = int((ymax-ymin)/cube_length)+1;
	nz = int((zmax-zmin)/cube_length)+1;
	cutoff_distance = 10.0;

	# wrap the solute atom coordinates in the box
	length_x = xmax-xmin
	length_y = ymax-ymin
	length_z = zmax-zmin
	for coor in solute_atom_coors:
		while (coor[0] <  xmin): coor[0]+=length_x
		while (coor[0] >= xmax): coor[0]-=length_x
		while (coor[1] <  ymin): coor[1]+=length_y
		while (coor[1] >= ymax): coor[1]-=length_y
		while (coor[2] <  zmin): coor[2]+=length_z
		while (coor[2] >= zmax): coor[2]-=length_z
	# test by printing out the wrapped solute atoms
	solute_mol.setCoor(numpy.array([solute_atom_coors],dtype=numpy.float))
	solute_mol.write_pdb("solute_wrapped.pdb",0,"w")

	# write out the solute_surf_mol
	solute_surf_mol = sasmol.SasMol(2)
	error = solute_mol.copy_molecule_using_mask(solute_surf_mol, solute_atom_surf_mask, 0)
	#pprint.pprint(list(solute_atom_surf_mask))
	if (len(error)>0): Error("solute_surf copy_molecule_using_mask error", error)
	solute_surf_mol.write_pdb("solute_surf.pdb",0,"w")

	# print-outs
	print("ZHL")
	pprint.pprint(solute_atom_types)
	print("ZHL")
	pprint.pprint(solute_atom_coors)
	print("ZHL")
	pprint.pprint(solvent_atom_types)
	print("ZHL")
	pprint.pprint(solvent_atom_coors)
	print("ZHL")

	# run it
	prdf.predict_prdf(solute_atom_types, solute_atom_coors, solute_atom_surf_mask, xmin, ymin, zmin, cube_length, nx, ny, nz, cutoff_distance)
