// Hailiang Zhang

#ifndef ATOM_H
#define ATOM_H

#include <map>
#include <string>

namespace Atom
{
	typedef std::string AtomType;
	extern std::map<AtomType, double> AtomVDWRadius;

	namespace SoluteAtom
	{
		typedef std::string SoluteAtomType;
	}// of namespace SoluteAtom

	namespace SolventAtom
	{
		typedef std::string SolventAtomType;
		extern std::map<SolventAtomType, double> SolventAtomDensity;
	}// of namespace SolventAtom

}// of namespace Atom

#define ATOM_ICC
#include "Atom.icc"
#undef ATOM_ICC

#endif
