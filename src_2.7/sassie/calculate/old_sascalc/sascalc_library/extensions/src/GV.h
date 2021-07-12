//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_GV
#define H_GV

#define TYPE double

#include "SasCalc.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class GV;
}

// definitions
class sascalc::GV :
    public sascalc::SasCalc
{
    // data
    protected:
        int _Ngv;
        TYPE * _gv;
        TYPE * _Iq_real;
        TYPE * _Iq_imag;

    // interface
    public:
        inline TYPE * Iq_real() const {return _Iq_real;}
        inline TYPE * Iq_imag() const {return _Iq_imag;}

    // methods
    public:
        void updateGV(const int Ngv);
        virtual void batch_load(const int offset, const int extend);
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0);
        virtual void calculate_pr(const int frame = 0) const;

    // meta
    public:
        GV(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q);
        virtual ~GV();

    // disallow copy and assign
    private:
        inline GV(const GV &);
        inline const GV & operator=(const GV &);
};

#endif
