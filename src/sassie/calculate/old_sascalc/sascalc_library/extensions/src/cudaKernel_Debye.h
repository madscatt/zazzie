//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_KERNEL_DEBYE
#define H_CUDA_KERNEL_DEBYE

// macros
#include "cudaCommon.h"

#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


namespace sascalc
{
    namespace cuda
    {
        __global__ void cudaKernel_Debye_calc(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int Xray = 0);
    }
}

#endif
