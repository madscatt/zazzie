//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_SASCALC
#define H_SASCALC

#define TYPE double

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class SasCalc;
}

// definitions
class sascalc::SasCalc
{
    // data
    protected:
        const int _Nq;
        const int _Natoms;
        const int _Nframes;
        const int _Nr;
        const TYPE _dr;

        const TYPE * _q;
        TYPE * _Iq;
        TYPE * _Pr;

        const TYPE * _coordinates;

    // interface
    public:
        inline TYPE * Iq() const {return _Iq;}
        inline TYPE * Pr() const {return _Pr;}

    // methods
    public:
        virtual void batch_load(const int offset, const int extend) = 0;
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0) = 0;
        virtual void calculate_pr(const int frame = 0) const = 0;

    // meta
    public:
        SasCalc(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q);
        virtual ~SasCalc();

    // disallow copy and assign
    private:
        inline SasCalc(const SasCalc &);
        inline const SasCalc & operator=(const SasCalc &);
};

#endif
