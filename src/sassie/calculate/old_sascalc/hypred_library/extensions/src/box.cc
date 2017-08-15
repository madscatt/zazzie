// Hailiang Zhang

#include "../include/box.h"

// dependencies
#include "../include/coor.h"
#include "../include/pRDF.h"
#include "../include/util.h"

// set up the _p_cube_type_distance from the passed solute atom set
Box & 
Box::
setupCube(const size_t n_solute_atoms, const solute_atom_type_t * solute_atom_types, const double * solute_atom_coors, const int * solute_atom_surf_mask, const double cutoff_distance) 
{
	// local variables
	int i_solute_atom, ix_solute_atom, iy_solute_atom, iz_solute_atom;
   	int i_cube, ix_cube, iy_cube, iz_cube, ix_cube_wrapped, iy_cube_wrapped, iz_cube_wrapped;
	double distance_to_nucleic, distance_to_surface;
	double solute_atom_radius;
	solute_atom_type_t solute_atom_type;
	Coor coor_cube_center, coor_solute_atom;
	int idx_span = int(cutoff_distance/_cube_length);

	// for each atom
	for (i_solute_atom=0; i_solute_atom<n_solute_atoms; i_solute_atom++)
	{
		coor_solute_atom = Coor(solute_atom_coors[i_solute_atom*3], solute_atom_coors[i_solute_atom*3+1], solute_atom_coors[i_solute_atom*3+2]);
		solute_atom_type = solute_atom_types[i_solute_atom];
		// ZHL note: type hard-coded
		// ZHL note: for a python-c setup compiler issue in mac, I am using a tempory string (atom_type_short)
		std::string atom_type_short = solute_atom_type.substr(0,1);
		solute_atom_radius = Atom::AtomVDWRadius[atom_type_short];
		//std::cout<<solute_atom_type<<" "<<solute_atom_radius<<std::endl;
		ix_solute_atom  = int((solute_atom_coors[i_solute_atom*3]-_c0.x())/_cube_length);
		iy_solute_atom  = int((solute_atom_coors[i_solute_atom*3+1]-_c0.y())/_cube_length);
		iz_solute_atom  = int((solute_atom_coors[i_solute_atom*3+2]-_c0.z())/_cube_length);
		for (ix_cube = ix_solute_atom-idx_span; ix_cube < ix_solute_atom+idx_span; ix_cube++)
		{
			for (iy_cube = iy_solute_atom-idx_span; iy_cube < iy_solute_atom+idx_span; iy_cube++)
			{
				for (iz_cube = iz_solute_atom-idx_span; iz_cube < iz_solute_atom+idx_span; iz_cube++)
				{
					// get the cube index
					ix_cube_wrapped=ix_cube; while(ix_cube_wrapped<0) ix_cube_wrapped += _nx; while(ix_cube_wrapped>=_nx) ix_cube_wrapped -= _nx;
					iy_cube_wrapped=iy_cube; while(iy_cube_wrapped<0) iy_cube_wrapped += _ny; while(iy_cube_wrapped>=_ny) iy_cube_wrapped -= _ny;
					iz_cube_wrapped=iz_cube; while(iz_cube_wrapped<0) iz_cube_wrapped += _nz; while(iz_cube_wrapped>=_nz) iz_cube_wrapped -= _nz;
					i_cube = iz_cube_wrapped*(_nx*_ny) + iy_cube_wrapped*_nx + ix_cube_wrapped;
					// get the cube center coordinate
					coor_cube_center = _c0+Coor((ix_cube+0.5)*_cube_length, (iy_cube+0.5)*_cube_length, (iz_cube+0.5)*_cube_length);
					// get the distance between the cube center and the solute atom scaled vdw surface
					distance_to_nucleic = coor_cube_center.distance(coor_solute_atom);
					distance_to_surface = distance_to_nucleic - solute_atom_radius;
					// if this distance is within the atomic radius
					if (distance_to_surface <= 0.0)
					{
						// set the cube_type_distance to be ("core",0.0)
						_p_cube_type_distance[i_cube] = std::make_pair(std::string("core").c_str(), 0.0);
					}
					// else, if this distance is smaller than the existing distance value in this cube
					else if (distance_to_nucleic<=cutoff_distance && distance_to_surface<_p_cube_type_distance[i_cube].second)
					{
						// if this a solute surface atom
						if (solute_atom_surf_mask[i_solute_atom])
						{
							// assign this cube
							_p_cube_type_distance[i_cube] = std::make_pair(solute_atom_type.c_str(), distance_to_surface);
						}
						// else, if this is a solute core atom
						else
						{
							// set the cube_type_distance to be ("core", distance_to_surface)
							_p_cube_type_distance[i_cube] = std::make_pair(std::string("core").c_str(), distance_to_surface);
						}
						//_p_cube_type_distance[i_cube] = std::make_pair(solute_atom_type.c_str(), distance_to_surface);//ZHL note:  no surf atoms
					}
				} // end for iz_cube
			} //end for iy_cube
		} // end for ix_cube
	} // end for atom

	// all done
	return * this;
}


// accumulate the solvent density into cubes 
Box & 
Box::
accumulateCubeDensity(const size_t n_solvent_atoms, const solvent_atom_type_t * solvent_atom_types, const double * solvent_atom_coors)
{
	// short namespace
	namespace SolventAtom = Atom::SolventAtom;

	// local variables
	size_t i_solvent_atom, i_cube, ix_cube, iy_cube, iz_cube;
	solvent_atom_type_t solvent_atom_type;
	Coor solvent_atom_coor;

	// for each solvent atom
	for (i_solvent_atom = 0; i_solvent_atom < n_solvent_atoms; i_solvent_atom++)
	{
		// determine the cube it resides
		ix_cube = size_t((solvent_atom_coors[i_solvent_atom*3]-_c0.x())/_cube_length);
		iy_cube = size_t((solvent_atom_coors[i_solvent_atom*3+1]-_c0.y())/_cube_length);
		iz_cube = size_t((solvent_atom_coors[i_solvent_atom*3+2]-_c0.z())/_cube_length);
		i_cube = iz_cube*(_nx*_ny) + iy_cube*_nx + ix_cube;
		// update _p_cube_density_number of that cube
		solvent_atom_type = solvent_atom_types[i_solvent_atom];
		_p_cube_density_number[i_cube].first += SolventAtom::SolventAtomDensity[solvent_atom_type];
		_p_cube_density_number[i_cube].second ++;
		//std::cout<<"_p_cube_density_number for "<<solvent_atom_type<<" @"<<i_cube<<" "<<_p_cube_density_number[i_cube].first<<" "<<_p_cube_density_number[i_cube].second<<std::endl;//ZHL test
	}
	//std::cout<<"new frame"<<std::endl;//ZHL test

	// all done
	return *this;
}

// predict the solvent density in cubes 
Box & 
Box::
predictCubeDensity(pRDF & prdfs, const double cutoff_distance)
{
	// local variables
	pRDF::pRDF_t * prdf;
	size_t i_cube, ix_cube, iy_cube, iz_cube;
	solute_atom_type_t solute_atom_type;
	double solute_atom_radius;
	double distance;
	size_t i_distance;
	double d_step = prdfs.d_step();
std::cout<<__LINE__<<std::endl;

	// loop over cubes
	for (ix_cube = 0; ix_cube < _nx; ix_cube++)
	{
		for (iy_cube = 0; iy_cube < _ny; iy_cube++)
		{
			for (iz_cube = 0; iz_cube < _nz; iz_cube++)
			{
				// get the solute_atom_type of this cube
				i_cube = iz_cube*(_nx*_ny) + iy_cube*_nx + ix_cube;
				solute_atom_type = _p_cube_type_distance[i_cube].first;
				// skip core or null cubes
				if (strcmp(solute_atom_type.c_str(), "core")==0 || strcmp(solute_atom_type.c_str(), "null")==0) continue;
std::cout<<__LINE__<<std::endl;
				// get the distance and index of distance for this cube
				std::string atom_type_short = solute_atom_type.substr(0,1);// ZHL note: for a python-c setup compiler issue in mac, I am using a tempory string (atom_type_short)
				solute_atom_radius = Atom::AtomVDWRadius[atom_type_short];
				distance = _p_cube_type_distance[i_cube].second + solute_atom_radius;
				i_distance = size_t(distance/d_step);
std::cout<<__LINE__<<std::endl;
std::cout<<solute_atom_type<<std::endl;
				// get the prdf profile for this atom type
				prdf = prdfs.getpRDF(solute_atom_type);
std::cout<<__LINE__<<std::endl;
				// set the cube density and number
				_p_cube_density_number[i_cube].first = prdf[i_distance].first/prdf[i_distance].second;
std::cout<<__LINE__<<std::endl;
				_p_cube_density_number[i_cube].second = 1;
std::cout<<__LINE__<<std::endl;
			}
		}
	}
std::cout<<__LINE__<<std::endl;

	// all done
	return *this;
}

// print out the types and distances values of all the not-null cubes in pdb format
void
Box::
fprintfTypeDistancePDB(const std::string & foutName) const
{
	// local variables
	int i_cube, ix_cube, iy_cube, iz_cube;
	solute_atom_type_t solute_atom_type;
	double distance, x, y, z;
	double x0=_c0.x();
	double y0=_c0.y();
	double z0=_c0.z();
	FILE * fin = fopen(foutName.c_str(),"w"); // ZHL note: a hardwired pdb name

	// loop over cubes
	for (ix_cube = 0; ix_cube < _nx; ix_cube++)
	{
		for (iy_cube = 0; iy_cube < _ny; iy_cube++)
		{
			for (iz_cube = 0; iz_cube < _nz; iz_cube++)
			{
				// get the solute_atom_type and distance of this cube
				i_cube = iz_cube*(_nx*_ny) + iy_cube*_nx + ix_cube;
				solute_atom_type = _p_cube_type_distance[i_cube].first;
				distance = _p_cube_type_distance[i_cube].second;
				// skip core or null cubes
				//if (strcmp(solute_atom_type.c_str(), "core")==0 || strcmp(solute_atom_type.c_str(), "null")==0) continue;
				if (strcmp(solute_atom_type.c_str(), "null")==0) continue;
				// get the cube center coordinate
				x = x0+(ix_cube+0.5)*_cube_length;
				y = y0+(iy_cube+0.5)*_cube_length;
				z = z0+(iz_cube+0.5)*_cube_length;
				// print them out
				fprintf(fin, "ATOM    145 %4s VAL A  25    %8.3f%8.3f%8.3f  1.00%6.2f      A1   N\n", solute_atom_type.c_str(), x, y, z, distance);
			}
		}
	}

	// clean up
	fclose(fin);
}

// print out the types and densities values of all the solvent layer cubes in pdb format
void
Box::
fprintfTypeDensityPDB(const std::string & foutName) const
{
	// local variables
	int i_cube, ix_cube, iy_cube, iz_cube;
	solute_atom_type_t solute_atom_type;
	double density, x, y, z;
	double x0=_c0.x();
	double y0=_c0.y();
	double z0=_c0.z();
	FILE * fin = fopen(foutName.c_str(),"w"); // ZHL note: a hardwired pdb name

	// loop over cubes
	for (ix_cube = 0; ix_cube < _nx; ix_cube++)
	{
		for (iy_cube = 0; iy_cube < _ny; iy_cube++)
		{
			for (iz_cube = 0; iz_cube < _nz; iz_cube++)
			{
				// get the solute_atom_type of this cube
				i_cube = iz_cube*(_nx*_ny) + iy_cube*_nx + ix_cube;
				solute_atom_type = _p_cube_type_distance[i_cube].first;
				// skip core or null cubes
				if (strcmp(solute_atom_type.c_str(), "core")==0 || strcmp(solute_atom_type.c_str(), "null")==0) continue;
				//if (strcmp(solute_atom_type.c_str(), "null")==0) continue;
				// get the density of this cube
				//density = (_p_cube_density_number[i_cube].first)/(_p_cube_density_number[i_cube].second);
				density = (_p_cube_density_number[i_cube].first);//ZHL test
				// get the cube center coordinate
				x = x0+(ix_cube+0.5)*_cube_length;
				y = y0+(iy_cube+0.5)*_cube_length;
				z = z0+(iz_cube+0.5)*_cube_length;
				// print them out
				fprintf(fin, "ATOM    145 %4s VAL A  25    %8.3f%8.3f%8.3f  1.00%6.2f      A1   N\n", solute_atom_type.c_str(), x, y, z, density);
			}
		}
	}

	// clean up
	fclose(fin);
}

// destructor
Box::
~Box()
{
	delete [] _p_cube_type_distance;
	delete [] _p_cube_density_number;
}
