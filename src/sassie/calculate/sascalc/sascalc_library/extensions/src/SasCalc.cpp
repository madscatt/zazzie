//
// Hailiang Zhang
// NIST & UTK
//

#include "SasCalc.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
/*
//////////////////////////////////////////////////////////
/// constructor 
/// currently not implemented due to overhead of B array preparation
//////////////////////////////////////////////////////////
sascalc::SasCalc::
SasCalc(const SasMol & mol, const ScVars & scvar)
{
}
*/

//////////////////////////////////////////////////////////
/// constructor 
/// patched to take the python-preprocess B array as the input
/// assume that the B-array memory layout follows what has been prepared in sascalc_util.py
//////////////////////////////////////////////////////////
sascalc::SasCalc::
SasCalc(const int Natoms, const ScVars & scvar, const double *const B_neutron_array, const double *const B_xray_array):
_Natoms(Natoms),
_Qmax(scvar._Qmax),
_Nq(scvar._Nq)
{
    // get B factors and its descriptors
    const std::string items[] = {"vacuum", "solvent", "complete"};
    const double *pB;
    const int Ncontrast_neutron = scvar._v_neutron_contrasts.size();
    for (int idx_contrast=0; idx_contrast!=Ncontrast_neutron; ++idx_contrast)
    {
        for (int idx=0; idx!=3; ++idx)
        {
            pB = &B_neutron_array[idx*(Ncontrast_neutron*_Natoms) + idx_contrast*_Natoms];
            std::stringstream B_descriptor;
            B_descriptor<<scvar._v_neutron_contrasts[idx_contrast].first<<"pD2O "<<items[idx];
            _B_neutron_vector.push_back(std::make_pair(B_descriptor.str(), pB));
            _neutron_I0.push_back(scvar._v_neutron_contrasts[idx_contrast].second);
        }
    }
    const int Ncontrast_xray = scvar._v_xray_contrasts.size();
    for (int idx_contrast=0; idx_contrast!=Ncontrast_xray; ++idx_contrast)
    {
        for (int idx=0; idx!=3; ++idx)
        {
            pB = &B_xray_array[idx*(Ncontrast_xray*_Nq*_Natoms) + idx_contrast*_Nq*_Natoms];
            std::stringstream B_descriptor;
            B_descriptor<<scvar._v_xray_contrasts[idx_contrast].first<<"pD2O "<<items[idx];
            _B_xray_vector.push_back(std::make_pair(B_descriptor.str(), pB));
            _xray_I0.push_back(scvar._v_xray_contrasts[idx_contrast].second);
        }
    }
    // _q
    _q = new double[_Nq]; 
    for (int iq=0; iq<_Nq; ++iq)
    {
        _q[iq] = iq*(_Qmax/(_Nq-1));
    }
}


//////////////////////////////////////////////////////////
/// SasCalc deallocation
//////////////////////////////////////////////////////////
sascalc::SasCalc::
~SasCalc()
{
    delete [] _q;
}
