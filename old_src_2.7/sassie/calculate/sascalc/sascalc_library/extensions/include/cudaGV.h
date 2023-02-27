//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDAGV
#define H_CUDA_CUDAGV

#include <cuda.h>
#include <cuda_runtime.h>

#include "GV.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaGV;
}

// definitions
class sascalc::cudaGV :
    public sascalc::GV
{
    // data
    protected:
        double * _gpu_gv;
        double * _gpu_q;
        double * _Is;
        double * _gpu_Is;

        double * _gpu_coor;
        double * _gpu_B_neutron;
        double * _gpu_B_xray;

        int * _gpu_flag_skip;

        int Ncontrast_neutron;
        int Ncontrast_xray;

    // methods
    protected:
        virtual void _updateGPU(const int *const flag_skip, const int Ngv, const double *const coor);
        virtual void _collectGPU(const int Ngv, const int *const flag_skip, double *const Iq);
        virtual void _calculate_singleFrame_fixed(const int *const flag_skip, const int Nitems, const double *const coor, double *const Iq, const int Ngv);

    // meta
    public:
        cudaGV(const int Natoms, const ScVars & scvar, const double *const B_neutron_array, const double *const B_xray_array);
        virtual ~cudaGV();

    // disallow copy and assign
    private:
        inline cudaGV(const cudaGV &);
        inline const cudaGV & operator=(const cudaGV &);
};

#endif
