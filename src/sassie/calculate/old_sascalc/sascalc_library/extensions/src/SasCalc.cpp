//
// Hailiang Zhang
// NIST & UTK
//

#include "SasCalc.h"
#include <math.h>
#include <cstdlib>
#include <iostream>

//////////////////////////////////////////////////////////
/// constructor 
//////////////////////////////////////////////////////////
sascalc::SasCalc::
SasCalc(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q):
_Natoms(Natoms),
_Nframes(Nframes),
_Nq(Nq),
_Nr(Nr),
_dr(dr),
_coordinates(coordinates),
_q(q)
{
    _Iq = (TYPE*)calloc(_Nq,sizeof(TYPE));
    _Pr = (TYPE*)calloc(_Nr,sizeof(TYPE));
}

//////////////////////////////////////////////////////////
/// SasCalc deallocation
//////////////////////////////////////////////////////////
sascalc::SasCalc::
~SasCalc()
{
    free(_Iq);
    free(_Pr);
}
