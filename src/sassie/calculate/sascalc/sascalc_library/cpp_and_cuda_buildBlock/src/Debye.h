//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_DEBYE
#define H_DEBYE

#define TYPE double

#include "SasCalc.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Debye;
}

// definitions
class sascalc::Debye :
    public sascalc::SasCalc
{
    // data
    protected:

    // interface
    public:

    // methods
    public:
        virtual void batch_load(const int offset, const int extend);
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0);
        virtual void calculate_pr(const int frame = 0) const;

    // meta
    public:
        Debye(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q);
        virtual ~Debye();

    // disallow copy and assign
    private:
        inline Debye(const Debye &);
        inline const Debye & operator=(const Debye &);
};

#endif
