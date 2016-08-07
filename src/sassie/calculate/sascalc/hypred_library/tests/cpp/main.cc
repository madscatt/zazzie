#include "../../include/Atom.h"
#include "../../include/coor.h"
#include "../../include/box.h"
#include "../../include/pRDF.h"

#include <iostream>

using namespace std;

int main()
{
	// types
	typedef Atom::SoluteAtom::SoluteAtomType solute_atom_type_t;
	typedef Atom::SolventAtom::SolventAtomType solvent_atom_type_t;

	// box parameters
	double cube_length=0.5;
	size_t nx=10; // ZHL OPEN
	size_t ny=10; // ZHL OPEN
	size_t nz=10; // ZHL OPEN
	Coor coor_box_origin(0.0, 0.0, 0.0); // ZHL OPEN
	// set up the box
	Box box(coor_box_origin, nx, ny, nz, cube_length);

	// prdf parameters
	size_t prdf_spacing = 0.5;
	size_t prdf_n=10; // ZHL OPEN
	// set up prdf
	pRDF rdfs(prdf_n, prdf_spacing);

	// set up the cubes in box
	std::vector<solute_atom_type_t> v_solute_atom_types; //ZHL OPEN
	std::vector<Coor> v_solute_atom_coordinates; //ZHL OPEN
	box.setupCube(v_solute_atom_types, v_solute_atom_coordinates); // set up the _p_cube_type_distance from the passed solute atom set

	// accumulate the solvent density into cubes
	size_t tot_frames=10; // ZHL OPEN
	std::vector<solvent_atom_type_t> v_solvent_atom_types; //ZHL OPEN
	std::vector<Coor> v_solvent_atom_coordinates;
	for (size_t frame = 0; frame<tot_frames; ++frame)
	{
		// ZHL OPEN
		// get the solvent atom coorindates for the current frame
		// ...
		// accumulate the solvent density into cubes
		box.accumulateCubeDensity(v_solvent_atom_types, v_solvent_atom_coordinates);
	}

	// finally update the pRDF based on the box
	rdfs.update(box);

	// output pRDF
	std::ofstream fout("pRDF.txt");
	rdfs.print(fout);
	fout.close();

	// clean up
	v_solute_atom_types.clear();
	v_solute_atom_coordinates.clear();
	v_solvent_atom_types.clear();
	v_solvent_atom_coordinates.clear();

	return 0;
}
