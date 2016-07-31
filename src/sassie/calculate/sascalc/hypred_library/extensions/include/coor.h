// Hailiang Zhang

#ifndef COOR_H
#define COOR_H

// externels
#include <cmath>

class Coor
{
	// data
	private:
		double _x;
		double _y;
		double _z;

	// interface
	public:
		inline double distance(const Coor &) const; // distance
		inline Coor operator+(const Coor &) const; // add operator
		inline double x() const; // get the x coordinate
		inline double y() const; // get the y coordinate
		inline double z() const; // get the z coordinate

	// constructor/destructor
	public:
		inline Coor();
		inline Coor(const double x, const double y, const double z);
		virtual ~Coor();
};


#define COOR_ICC
#include "coor.icc"
#undef COOR_ICC 

#endif
