//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_WAPPERCUDAKERNEL
#define H_WAPPERCUDAKERNEL

#include "cudaUtil.h"
#include <curand_kernel.h>

namespace sascalc
{
    namespace cuda
    {
        void wrapper_cudaGV_calcIq(const double * const gpu_coor, const double * const gpu_B_neutron, const double * const gpu_B_xray, const int *const gpu_flag_skip, const double * const gpu_gv, const double * const gpu_q, double * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int Ncontrast_neutron, const int Ncontrast_xray);

    }
}

#endif
