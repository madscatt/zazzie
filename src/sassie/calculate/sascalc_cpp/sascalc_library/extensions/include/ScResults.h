//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_SCRESULTS
#define H_SCRESULTS

#include <string>
#include <vector>
#include <dirent.h>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <sys/mman.h>
#include <cmath>

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class ScResults;
    class SasCalc;
}

// definitions
class sascalc::ScResults
{
    // type
    public:
        typedef std::pair<std::string, const double *> Iq_t;

    // data
    protected:
        const int _Nframes;
        std::vector<Iq_t> _Iq_vector;
        const int _Nq;
        const double _Qmax;
        double * _Q;

    // interface
    public:
        inline void push_result(Iq_t Iq) {_Iq_vector.push_back(Iq);}

    // methods
    public:
        const sascalc::ScResults & save(const std::string &, const std::string &) const;
        const sascalc::ScResults & epilogue(const std::string &) const;
        friend std::ostream & operator<<(std::ostream &os, const sascalc::ScResults & results);

    // meta
    public:
        ScResults(const int, const int, const double);
        virtual ~ScResults();

    // disallow copy and assign
    private:
        inline ScResults(const ScResults &);
        inline const ScResults & operator=(const ScResults &);
};

namespace sascalc
{
    std::ostream & operator<<(std::ostream &os, const sascalc::ScResults & results);
}

#endif
