#include <string>
#include <vector>

#include <SasMol.h>
#include <ScVars.h>
#include <SasCalc.h>
#include <ScResults.h>

int main()
{
    // initialize the sasmol object
    SasMol * p_sasmol = sasmol::SasMol(0);
    p_sasmol->read_pdb("test.pdb");

    // initialize the scvars object
    const std::string golden_vector_method_option = "fixed";
    const int number_of_golden_vectors = 35;
    const std::string xon = "neutron";
    const int number_of_contrasts = 2;
    ScVars * p_scvars = sascalc::ScVars(golden_vector_method_option, number_of_golden_vectors, xon, number_of_contrasts);

    // initialize the sascalc object
    SasCalc * p_sascalc = sascalc::SasCalc(p_sasmol, p_scvars);

    // scattering calculation
    for (int i=0; i<3; ++i)
    {
        // generate some coordinates
        std::vector <std::vector <double> > coordinates = ...;
        // kernel calculation
        sascalc::ScResults * p_scresults = o_sascalc.calculate(coordinates);
        // save the results
        p_scresults->write_to_file("results_"+std::string(i)+".txt");
    }

    // epilogue/clean up
    p_sascalc_epilogue();

    // return
    return 0;
}
