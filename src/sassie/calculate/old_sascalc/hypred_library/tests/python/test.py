import prdf
import numpy

if __name__=="__main__":
	solute_atom_types = ["C","NH"]
	solute_atom_coors = numpy.array([[0.,0.,0.],[1.,1.,1.]], dtype=numpy.double)
	solvent_atom_types = ["OH2","H1"]
	solvent_atom_coors = numpy.array([[[2.,2.,2.],[3.,3.,3.]],[[4.,4.,4.],[5.,5.,5.]]], dtype=numpy.double)
	x0=0.0;
	y0=0.0;
	z0=0.0;
	cube_length=0.5;
	nx=50;
	ny=50;
	nz=50;
	prdf_spacing = 0.5;
	prdf_n=100;
	prdf.construct_prdf(solute_atom_types, solute_atom_coors, solvent_atom_types, solvent_atom_coors, x0, y0, z0, cube_length, nx, ny, nz, prdf_spacing, prdf_n)
