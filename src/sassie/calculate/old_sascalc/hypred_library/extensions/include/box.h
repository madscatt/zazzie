// Hailiang Zhang

#ifndef BOX_H
#define BOX_H

// dependencies
#include "common.h"
#include "Atom.h"
#include "coor.h"
#include "util.h"

// externals
#include <vector>
#include <algorithm>
#include <utility>
#include <stdlib.h>
#include <stdio.h>

// forward declaration
class pRDF;

/*
   The box containing the solvent and solute information for a single frame.
   It is composed of cubes; it is orthogonal in cartisian system.
   Each cube contains either the solvent information (i.e. total number of protons for X-ray, or total cross-sections for neutron),
     or the solute information
*/


class Box
{
	// types
	public:
		typedef Atom::AtomType atom_type_t;
		typedef Atom::SoluteAtom::SoluteAtomType solute_atom_type_t;
		typedef Atom::SolventAtom::SolventAtomType solvent_atom_type_t;
		typedef std::pair<solute_atom_type_t, double> cube_type_distance_t;
		typedef std::pair<double, size_t> cube_density_number_t;
	// data
	private:
		cube_type_distance_t * _p_cube_type_distance; // the atomic type the cube belonging to, and the distance from the cube center to the atom
		cube_density_number_t * _p_cube_density_number; // the accumulated solvent density and the number of accumulation for each cube
  		int _nx, _ny, _nz; // the number of cubes along each dimension of the box
  		Coor _c0; // the starting left/bottom coordinates of the first cube
  		double _cube_length; // the length of each dimension of the cube

	// interfaces
	public:
		cube_type_distance_t * & p_cube_type_distance() {return _p_cube_type_distance;}
		cube_density_number_t * & p_cube_density_number() {return _p_cube_density_number;}
		inline size_t nx() const {return _nx;}
		inline size_t ny() const {return _ny;}
		inline size_t nz() const {return _nz;}
		inline double cube_length() const {return _cube_length;}

	// methods 
	public:
		virtual Box & setupCube(const size_t n_atoms, const solute_atom_type_t * solute_atom_types, const double * solute_atom_coors, const int * solute_atom_surf_mask, const double cutoff_distance); // set up the _p_cube_type_distance from the passed solute atom set
		virtual Box & accumulateCubeDensity(const size_t n_atoms, const solvent_atom_type_t *, const double *); // accumulate the solvent density into cubes 
		virtual Box & predictCubeDensity(pRDF &, const double); // predict the solvent density in cubes 
		virtual void fprintfTypeDistancePDB(const std::string &) const; // print out the types and distances values of all the not-null cubes in pdb format
		virtual void fprintfTypeDensityPDB(const std::string &) const; // print out the types and denisties values of all the solvent layer cubes in pdb format

	// meta methods
	public:
		inline Box(); // default constructor
		inline Box(const Coor &, const size_t, const size_t, const size_t, const double); // constructor
		inline Box(const std::vector<Coor> & solute_coordinates, const double cube_length, const double delta_expand); // construct the box from solute atoms
  		virtual ~ Box ();
	private: // disallowed
		inline Box(const Box &);
		inline Box & operator=(const Box &);
};

#define BOX_ICC
#include "box.icc"
#undef BOX_ICC

#endif
