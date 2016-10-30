//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_SASCALC
#define H_SASCALC

#include "ScVars.h"
#include "ScResults.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class ScVars;
    class ScResults;
    class SasCalc;
}

// definitions
class sascalc::SasCalc
{
    // type
    public:
        typedef std::pair<std::string, const double *> B_t;

    // data
    protected:
        const int _Natoms;
        const int _Nq;
        const double _Qmax;
        double *_q;
        std::vector<B_t> _B_neutron_vector;
        std::vector<B_t> _B_xray_vector;
        std::vector<double> _neutron_I0;
        std::vector<double> _xray_I0;

    // interface
    public:
        //inline double * Is() const {return _Is;}

    // methods
    public:
        virtual sascalc::ScResults * calculate(const double *const coor, const int Nframes) = 0;

    // meta
    public:
        //SasCalc(const SasMol &, const ScVars &); //< currently not implemented due to overhead of B array preparation
        SasCalc(const int, const ScVars &, const double *const, const double *const); //< patched constructor to take the python-preprocess B array as the input
        virtual ~SasCalc();

    // disallow copy and assign
    private:
        inline SasCalc(const SasCalc &);
        inline const SasCalc & operator=(const SasCalc &);
};

#endif
