//
// Hailiang Zhang
// NIST & UTK
//

#include "Debye.h"
#include <math.h>
#include <cstdlib>
#include <iostream>

//////////////////////////////////////////////////////////
/// constructor 
//////////////////////////////////////////////////////////
sascalc::Debye::
Debye(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q):
SasCalc(Natoms,Nframes,Nq,Nr,dr,coordinates,q)
{
}

//////////////////////////////////////////////////////////
/// batch load
//////////////////////////////////////////////////////////
void
sascalc::Debye::
batch_load(const int offset, const int extend)
{
    _coordinates += offset*_Natoms*3;
}

//////////////////////////////////////////////////////////
/// Pr calculator
//////////////////////////////////////////////////////////
void
sascalc::Debye::
calculate_pr(const int frame) const
{
    const TYPE * coor = &_coordinates[frame*_Natoms*3];

    for (int i=0; i<_Nr; ++i) _Pr[i] = 0.0;

    // get p(r)
    // this slows things quite a bit, because p(r) needs looping over atom pairs, which is not requires by Debye
    int iatom, jatom, idx;
    TYPE r;
    TYPE xi,yi,zi,xj,yj,zj;
    for (iatom=0; iatom<_Natoms; ++iatom)
    {
        xi = coor[iatom];
        yi = coor[_Natoms + iatom];
        zi = coor[_Natoms*2 + iatom];
        for (jatom=0; jatom<_Natoms; ++jatom)
        {
            xj = coor[jatom];
            yj = coor[_Natoms + jatom];
            zj = coor[_Natoms*2 + jatom];
            r = sqrt(pow((xi-xj),2.0)+pow((yi-yj),2.0)+pow((zi-zj),2.0));
            idx = int(r/_dr);
            if (idx>=_Nr) continue;
            else _Pr[idx] += 1.;
        }
    }
}

//////////////////////////////////////////////////////////
/// debye calculator
//////////////////////////////////////////////////////////
void
sascalc::Debye::
calculate(const TYPE *const b, const int frame, const int xray)
{
    const TYPE * coor = &_coordinates[frame*_Natoms*3];
    const TYPE * B = b;

    int iatom, jatom, idx;
    TYPE r;
    TYPE xi,yi,zi,xj,yj,zj;
    TYPE bi,bj;
    int iq;
    TYPE q;
    TYPE I;
    for (iq=0; iq<_Nq; ++iq)
    {
        q = _q[iq];
        I = 0;
        for (iatom=0; iatom<_Natoms; ++iatom)
        {
            xi = coor[iatom];
            yi = coor[_Natoms + iatom];
            zi = coor[_Natoms*2 + iatom];
            if (xray) bi = B[_Natoms*iq+iatom];
            else bi = B[iatom];
            for (jatom=0; jatom<_Natoms; ++jatom)
            {
                xj = coor[jatom];
                yj = coor[_Natoms + jatom];
                zj = coor[_Natoms*2 + jatom];
                if (xray) bj = B[_Natoms*iq+jatom];
                else bj = B[jatom];
                r = sqrt(pow((xi-xj),2.0)+pow((yi-yj),2.0)+pow((zi-zj),2.0));
                if (q*r==0) I += bi*bj;
                else I += bi*bj*sin(q*r)/(q*r);
            }
        }
        _Iq[iq] = I;
    }
}

//////////////////////////////////////////////////////////
/// Debye deallocation
//////////////////////////////////////////////////////////
sascalc::Debye::
~Debye()
{
    free(_Iq);
    free(_Pr);
}
