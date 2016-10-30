//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_KERNEL_GV
#define H_CUDA_KERNEL_GV

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
        __global__ void cudaKernel_GV_calc(const double * const gpu_coor, const double * const gpu_B_neutron, const double * const gpu_B_xray, const int *const gpu_flag_skip, const double * const gpu_gv, const double * const gpu_q, double * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int Ncontrast_neutron, const int Ncontrast_xray);
    }
}

#endif
