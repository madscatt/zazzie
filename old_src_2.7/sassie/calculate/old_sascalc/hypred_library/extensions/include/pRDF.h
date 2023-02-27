// Hailiang Zhang

#ifndef PRDF_H
#define PRDF_H

// dependencies
#include "box.h"
#include "Atom.h"
#include "coor.h"

// externels
#include <vector>
#include <utility>
#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>

/*
   the proximal radial distribution functions for different solute atom types
*/

// forward declaration
class Box;

class pRDF
{
	// friend
	friend class Box;

	// type
	public:
		typedef Atom::SoluteAtom::SoluteAtomType solute_atom_type_t;
		typedef std::pair<double, size_t> pRDF_t;
		typedef std::map<solute_atom_type_t, pRDF_t*> pRDF_map_t;
	// data
	private:
		double _d_step; // the distance step for pRDF
		size_t _n_step; // the number of distance steps
		pRDF_map_t _pRDF_map; // the pRDF map: the key is the atom type, the value is the accumulated density and the number of the accumulation for each grid

	// interface
	public:
		inline pRDF_t * & getpRDF(const solute_atom_type_t solute_atom_type); // get the pRDF for a given solute atom type; if key not exists, insert value
		inline const double d_step() const;
		inline const size_t n_step() const;

	// methods
	public:
		virtual const pRDF & update(Box &); // update the pRDF data
		virtual const pRDF & print(std::ostream &, Box &, const size_t); // print pRDF data

	// meta methods
	public:
		inline pRDF(); // default constructor
		inline pRDF(const size_t n, const double spacing); // constructor
		inline pRDF(const std::string & pRDFFileName); // constructor from the externel pRDF file (previously written by pRDF::print)
		virtual ~pRDF(); // destructor
};


#define PRDF_ICC
#include "pRDF.icc"
#undef PRDF_ICC

#endif
