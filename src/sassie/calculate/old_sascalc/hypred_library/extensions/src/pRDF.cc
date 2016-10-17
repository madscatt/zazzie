// Hailiang Zhang


#include "../include/pRDF.h"

// dependencies
//#include "box.h"

// update the pRDF data based on the cube denisty in the current box
const pRDF &
pRDF::
update(Box & box)
{
	// local variables
	size_t i_cube;
	size_t i_distance;
	size_t n_density;
	solute_atom_type_t solute_atom_type;
	double solute_atom_radius;
	double distance;
	double density;

	// get the box parameters
	const size_t nx = box.nx();
	const size_t ny = box.ny();
	const size_t nz = box.nz();

	// loop over cubes
	for (i_cube=0; i_cube< nx*ny*nz; i_cube++)
	{
		// get the solute atom type of the cube
		solute_atom_type = box.p_cube_type_distance()[i_cube].first;
		if (strcmp(solute_atom_type.c_str(), "core")==0 || strcmp(solute_atom_type.c_str(), "null")==0) continue;
		std::string atom_type_short = solute_atom_type.substr(0,1);
		solute_atom_radius = Atom::AtomVDWRadius[atom_type_short];
		// get the index of distance (where this distance stays in the "pRDF_t*" array)
		distance = box.p_cube_type_distance()[i_cube].second + solute_atom_radius;
		i_distance = size_t(distance/_d_step);
		if (i_distance>_n_step)
		{
			std::cout<<i_distance<<" "<<_n_step<<std::endl;
			prdf::utils::Error("the index of the distance exceeds the total number of distances!");
		}
		// get the accumlated density value and counts of this cube
		density = box.p_cube_density_number()[i_cube].first;
		n_density = box.p_cube_density_number()[i_cube].second;
		//std::cout<<"box.p_cube_density_number @"<<i_cube<<" "<<density<<std::endl;// ZHL test
		// deposit this density/count into pRDF for this solute atom type
		getpRDF(solute_atom_type)[i_distance].first += density;
		//getpRDF(solute_atom_type)[i_distance].second+= n_density;//ZHL
		getpRDF(solute_atom_type)[i_distance].second++;
		//std::cout<<"getpRDF for "<<solute_atom_type<<" @"<<i_distance<<" "<<getpRDF(solute_atom_type)[i_distance].first<<" "<<getpRDF(solute_atom_type)[i_distance].second<<std::endl;//ZHL test

	} // end of "for i_cube"

	// all done
	return *this;
}

// print out the pRDF data
const pRDF &
pRDF::
print(std::ostream & os, Box & box, const size_t tot_frames)
{
	// local variables
	size_t i;
	pRDF_t rdf;
	solute_atom_type_t solute_atom_type;
	const double cube_length = box.cube_length();

	// print out _d_step and _n_step
	os<<"dstep "<<_d_step<<std::endl;
	os<<"nstep "<<_n_step<<std::endl;

	// loop over solute atom types
	for (pRDF_map_t::const_iterator it = _pRDF_map.begin(); it != _pRDF_map.end(); ++it)
	{
		// print out the solute atom type
		os<<it->first<<" ";
		// print out its RDF
		for (i=0; i<_n_step; i++)
		{
			rdf = (it->second)[i];
			if (rdf.second) os<<rdf.first/rdf.second/tot_frames/(cube_length*cube_length*cube_length)<<" "; 
			else os<<0.0<<" ";
		}
		os<<std::endl;
	}
	// all done
	return * this;
}


// destructor 
pRDF::
~pRDF()
{
	pRDF_map_t::iterator it;
	for (it=_pRDF_map.begin(); it!=_pRDF_map.end(); ++it) delete [] it->second;
}
