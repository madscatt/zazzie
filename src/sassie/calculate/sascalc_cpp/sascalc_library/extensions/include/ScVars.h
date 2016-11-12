//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_SCVARS
#define H_SCVARS

#include <iostream>
#include <string>
#include <vector>
#include <utility>

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class ScVars;
}

// definitions
class sascalc::ScVars
{
    // type
    public:
        typedef std::pair<double, double> neutron_contrast_t;
        typedef std::pair<double, double> xray_contrast_t;
        typedef std::pair<std::string, double> neutron_exH_t;
        typedef std::pair<std::string, double> neutron_deut_t;

    // data
    public:
        const int _Nq;
        const double _Qmax;
        const std::vector<neutron_contrast_t> _v_neutron_contrasts;
        const std::vector<xray_contrast_t> _v_xray_contrasts;
        const std::vector<neutron_exH_t> _v_neutron_exH;
        const std::vector<neutron_deut_t> _v_neutron_deut;
        const std::string _gv_method; 
        const double _gv_parameter; 

    // interface
    public:
        void show();

    // methods
    public:

    // meta
    public:
        ScVars(const int, const double, const std::vector<neutron_contrast_t> &, const std::vector<neutron_exH_t> &, const std::vector<neutron_deut_t> &, const std::vector<xray_contrast_t> &, const std::string, const double);
        virtual ~ScVars();

    // disallow copy and assign
    private:
        inline ScVars(const ScVars &);
        inline const ScVars & operator=(const ScVars &);
};

#endif
