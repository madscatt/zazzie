// Hailiang Zhang

#include "../include/Atom.h"

#include <map>

#define SCALE 0.53

namespace Atom
{
	std::map<AtomType, double> AtomVDWRadius
		=
		{
			{"null", 0.0},
			{"H", 1.2*SCALE},
			{"C", 1.7*SCALE},
			{"N", 1.55*SCALE},
			{"O", 1.52*SCALE},
			{"S", 1.80*SCALE}
		};

	namespace SoluteAtom
	{

	}// of namespace SoluteAtom

	namespace SolventAtom
	{
		/*
		std::map<SolventAtomType, double> createMap()
		{
			std::map<SolventAtomType, double> m;
			m["null"]=0.0;
			m["O"]=8.0;
			m["HOH"]=1.0;
			return m;
		}

		std::map<SolventAtomType, double> SolventAtomDensity = createMap();
		*/
		std::map<SolventAtomType, double> SolventAtomDensity
			=
			{
				{"null",0.0},
				{"OH2",8.0},
				{"H1",1.0},
				{"H2",1.0}
			};

	}// of namespace SolventAtom

}// of namespace Atom

