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
        void wrapper_cudaGV_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int xray = 0);

        void wrapper_cudaGV_calcIq_complex(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is_real, TYPE * const gpu_Is_imag, const int Natoms, const int Nq, const int Ngv, const int xray = 0);

        void wrapper_cudaDebye_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int xray = 0);

        void wrapper_cudaDebyePBC_calcIq(const TYPE boxl, const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int xray = 0);

    }
}

#endif
