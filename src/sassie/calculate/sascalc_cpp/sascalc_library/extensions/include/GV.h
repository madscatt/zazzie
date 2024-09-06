//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_GV
#define H_GV

#include "SasCalc.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class ScVars;
    class GV;
}

// definitions
class sascalc::GV :
    public sascalc::SasCalc
{
    // data
    protected:
        const std::string _gv_method; 
        const double _gv_parameter; 
        double * _Iq_real;
        double * _Iq_imag;

    // interface
    public:
        inline double * Iq_real() const {return _Iq_real;}
        inline double * Iq_imag() const {return _Iq_imag;}

    // methods
    protected:
        double * _getGV(int Ngv) const;
        int _tolerance_to_ngv(const double tolerance) const;
        void _calculate_singleFrame_fixed(const bool *const flag_skip, const int Nitems, const double *const coor, double *const Iq, const int Ngv) const;
        void _calculate_singleFrame_converge(const double *const coor, double *const Iq, const double tolerance) const;

    public:
        virtual sascalc::ScResults * calculate(const double *const coor, const int Nframes) const;

    // meta
    public:
        //GV(const SasMol &, const ScVars &); //< currently not implemented due to overhead of B array preparation
        GV(const int, const ScVars &, const double *const, const double *const); //< patched constructor to take the python-preprocess B array as the input
        virtual ~GV();

    // disallow copy and assign
    private:
        inline GV(const GV &);
        inline const GV & operator=(const GV &);
};

#endif
