//
// Hailiang Zhang
// NIST & UTK
//

#include "ScVars.h"
#include <iostream>

//////////////////////////////////////////////////////////
/// constructor 
//////////////////////////////////////////////////////////
sascalc::ScVars::
ScVars(const int Nq, const double Qmax, const std::vector<neutron_contrast_t> & v_neutron_contrasts, const std::vector<neutron_exH_t> & v_neutron_exH, const std::vector<neutron_deut_t> & v_neutron_deut, const std::vector<xray_contrast_t> & v_xray_contrasts, const std::string gv_method, const double gv_parameter):
_Nq(Nq),
_Qmax(Qmax),
_v_neutron_contrasts(v_neutron_contrasts),
_v_neutron_exH(v_neutron_exH),
_v_neutron_deut(v_neutron_deut),
_v_xray_contrasts(v_xray_contrasts),
_gv_method(gv_method), 
_gv_parameter(gv_parameter)
{
}

//////////////////////////////////////////////////////////
/// display
//////////////////////////////////////////////////////////
void
sascalc::ScVars::
show()
{
    std::cout<<"Nq: "<<_Nq<<std::endl;
    std::cout<<"Qmax: "<<_Qmax<<std::endl;
    std::cout<<"neutron contrast:"<<std::endl;
    for (int i=0; i<_v_neutron_contrasts.size(); ++i) std::cout<<"    "<<_v_neutron_contrasts[i].first<<", "<<_v_neutron_contrasts[i].second<<std::endl;
    std::cout<<"neutron exH regions:"<<std::endl;
    for (int i=0; i<_v_neutron_exH.size(); ++i) std::cout<<"    "<<_v_neutron_exH[i].first<<", "<<_v_neutron_exH[i].second<<std::endl;
    std::cout<<"neutron deuterated regions:"<<std::endl;
    for (int i=0; i<_v_neutron_deut.size(); ++i) std::cout<<"    "<<_v_neutron_deut[i].first<<", "<<_v_neutron_deut[i].second<<std::endl;
    std::cout<<"x-ray contrast data:"<<std::endl;
    for (int i=0; i<_v_xray_contrasts.size(); ++i) std::cout<<"    "<<_v_xray_contrasts[i].first<<", "<<_v_xray_contrasts[i].second<<std::endl;
    std::cout<<"gv_method: "<<_gv_method<<std::endl; 
    std::cout<<"gv_parameter: "<<_gv_parameter<<std::endl;
}

//////////////////////////////////////////////////////////
/// ScVars deallocation
//////////////////////////////////////////////////////////
sascalc::ScVars::
~ScVars()
{
}
