import SasMol
import ScVars
import SasCalc

if __name__=='__main__':

    # initialize the sasmol object
    o_sasmol = sasmol.SasMol(0);
    o_sasmol.read_pdb("test.pdb");

    # initialize the scvars object
    golden_vector_method_option = "fixed"
    number_of_golden_vectors = 35
    xon = "neutron"
    number_of_contrasts = 2
    o_scvars = ScVars(golden_vector_method_option, number_of_golden_vectors, xon, number_of_contrasts)

    # initialize the sascalc object
    o_sascalc = SasCalc(o_sasmol, o_scvars)

    # scattering calculation
    for i in range(3):
        # generate some coordinates
        coordinates = ...
        # kernel calculation
        sascalc_results = o_sascalc.calculate(coordinates)
        # save the results
        sascalc_results.write_to_file("results_%d.txt"%i)

    # epilogue/clean up
    o_sascalc_epilogue()
