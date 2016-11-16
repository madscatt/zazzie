#include <iostream>
#include "ScVars.h"

using namespace sascalc;

int main()
{
    const int Nq = 10;
    const double Qmax = 0.5;

    ///////////////////////////////////////
    // neutron inputs
    ///////////////////////////////////////
    // constrast points data
    ScVars::neutron_contrast_t neutron_contrast;
    std::vector<ScVars::neutron_contrast_t> neutron_contrasts_vector;
    neutron_contrast = std::make_pair(0.0, 1.0);
    neutron_contrasts_vector.push_back(neutron_contrast);
    neutron_contrast = std::make_pair(100.0, 1.0);
    neutron_contrasts_vector.push_back(neutron_contrast);
    // exH region data
    ScVars::neutron_exH_t neutron_exH;
    std::vector<ScVars::neutron_exH_t> neutron_exH_vector;
    neutron_exH = std::make_pair(std::string("moltype protein"), 0.95);
    neutron_exH_vector.push_back(neutron_exH);
    // deuterated region data
    ScVars::neutron_deut_t neutron_deut;
    std::vector<ScVars::neutron_deut_t> neutron_deut_vector;
    neutron_deut = std::make_pair(std::string("moltype protein"), 0.0);
    neutron_deut_vector.push_back(neutron_deut);
    
    ///////////////////////////////////////
    // xray inputs
    ///////////////////////////////////////
    // constrast points data
    std::vector<ScVars::xray_contrast_t> xray_contrasts_vector;

    ///////////////////////////////////////
    // advanced inputs
    ///////////////////////////////////////
    const std::string gv_method = std::string("fixed");
    const double gv_parameter = 35;


    ///////////////////////////////////////
    // initialize the ScVar object
    ///////////////////////////////////////
    ScVars * p_scvars = new ScVars(Nq, Qmax, neutron_contrasts_vector, neutron_exH_vector, neutron_deut_vector, xray_contrasts_vector, gv_method, gv_parameter);

    ///////////////////////////////////////
    // show the ScVar object
    ///////////////////////////////////////
    p_scvars->show();


}
