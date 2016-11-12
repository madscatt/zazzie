//
// Hailiang Zhang
// NIST & UTK
//

#include "GV.h"
#include <math.h>
#include <cstdlib>
#include <iostream>

//////////////////////////////////////////////////////////
/// constructor 
//////////////////////////////////////////////////////////
sascalc::GV::
GV(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q):
SasCalc(Natoms,Nframes,Nq,Nr,dr,coordinates,q)
{
    _Ngv = 0;
    _gv = NULL;
    _Iq_real = (TYPE*)calloc(_Nq,sizeof(TYPE));
    _Iq_imag = (TYPE*)calloc(_Nq,sizeof(TYPE));
}

//////////////////////////////////////////////////////////
/// batch load
//////////////////////////////////////////////////////////
void
sascalc::GV::
batch_load(const int offset, const int extend)
{
    _coordinates += offset*_Natoms*3;
}

//////////////////////////////////////////////////////////
/// Golden vector calculator
//////////////////////////////////////////////////////////
void
sascalc::GV::
updateGV(const int Ngv)
{
    // assign Ngv
    _Ngv = Ngv;

    // free exisiting GV
    if (_gv) free(_gv);

    // setup GV
    if (_Ngv%2==0)
    {
        std::cout<<"The number of golden vectors should be an odd integer, and so it will be reset to be: "<<_Ngv+1<<std::endl;
        ++_Ngv;
    }
    _gv = (TYPE*)malloc(3*_Ngv*sizeof(TYPE));
    
    const TYPE phi_inv = 2.0/(1+sqrt(5.)); // golden ratio
    TYPE cos_theta, sin_theta, phi;
    TYPE qx,qy,qz;
    int igv;
    const int rank = _Ngv/2;
    for (int i=-rank; i<=rank; i++)
    {   
        sin_theta = cos(asin(2.0*i/_Ngv));
        cos_theta = 2.0*i/_Ngv;
        phi = 2*M_PI*i*phi_inv;
        igv = i + rank;
        _gv[igv] = sin_theta*cos(phi);
        _gv[_Ngv+igv] = sin_theta*sin(phi);
        _gv[2*_Ngv+igv] = cos_theta;
    }   
}


//////////////////////////////////////////////////////////
/// SasGolden vector calculator
//////////////////////////////////////////////////////////
void
sascalc::GV::
calculate(const TYPE *const B, const int frame, const int xray)
{
    const TYPE * coor = &_coordinates[frame*_Natoms*3];

    // locals
    int iatom,iq,igv;
    TYPE x, y, z, b;
    TYPE qx,qy,qz;
    TYPE q_dot_r;
    TYPE sine, cose;
    TYPE Ireal = 0.0;
    TYPE Iimag = 0.0;
    TYPE qmag;

    // summation
    for (iq=0; iq<_Nq; ++iq)
    {
        _Iq[iq] = 0.0;
        _Iq_real[iq] = 0.0;
        _Iq_imag[iq] = 0.0;
        qmag = _q[iq];
        //std::cout<<"Q: "<<qmag<<std::endl;
        for (igv=0; igv<_Ngv; ++igv)
        {
            qx = qmag * _gv[igv];
            qy = qmag * _gv[_Ngv+igv];
            qz = qmag * _gv[2*_Ngv+igv];
            Ireal = 0.0;
            Iimag = 0.0;
            for (iatom=0; iatom<_Natoms; ++iatom)
            {
                x = coor[iatom];
                y = coor[_Natoms + iatom];
                z = coor[_Natoms*2 + iatom];
                if (xray) b = B[_Natoms*iq+iatom];
                else b = B[iatom];
                q_dot_r = qx*x+qy*y+qz*z;
                //sincos_TYPE(q_dot_r, &sine, &cose);
                sine = sin(q_dot_r);
                cose = cos(q_dot_r);
                Ireal += b*cose;
                Iimag += b*sine;
                //if (iq==1 && igv==1) printf("%8.3f %8.3f %8.3f  %8.3f  %8.3f %8.3f %8.3f\n",x,y,z,b,qx,qy,qz);
            }
            _Iq[iq] += Ireal*Ireal + Iimag*Iimag;
            _Iq_real[iq] += Ireal*Ireal;
            _Iq_imag[iq] += Iimag*Iimag;
        }
        _Iq[iq] /= _Ngv;
        //_Iq[iq] = (_Iq_real[iq]+_Iq_imag[iq])/_Ngv;
        _Iq_real[iq] = pow(_Iq_real[iq]/_Ngv, 0.5);
        _Iq_imag[iq] = pow(_Iq_imag[iq]/_Ngv, 0.5);
    }
    //printf("%16.8f %16.8f %16.8f\n",_Iq[0],_Iq_real[0],_Iq_imag[0]);
}

//////////////////////////////////////////////////////////
/// SasGolden vector calculator
//////////////////////////////////////////////////////////
void
sascalc::GV::
calculate_pr(const int frame) const
{
    const TYPE * coor = &_coordinates[frame*_Natoms*3];

    for (int i=0; i<_Nr; ++i) _Pr[i] = 0.0;

    // get p(r)
    // this slows things quite a bit, because p(r) needs looping over atom pairs, which is not requires by GV
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
/// GV deallocation
//////////////////////////////////////////////////////////
sascalc::GV::
~GV()
{
}
