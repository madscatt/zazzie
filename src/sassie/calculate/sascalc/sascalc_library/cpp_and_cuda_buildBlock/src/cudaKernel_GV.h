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
        __global__ void cudaKernel_GV_calc(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const int xray = 0);
        __global__ void cudaKernel_GV_calc_complex(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I_real, TYPE * const gpu_I_imag, const int Natoms, const int Nq, const int Ngv, const int xray = 0);

        /*
        __global__ void cudaKernel_GV_calc_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent);

        __global__ void cudaKernel_GV_calc_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent, const TYPE scale_volume);
        */
    }
}

#endif
