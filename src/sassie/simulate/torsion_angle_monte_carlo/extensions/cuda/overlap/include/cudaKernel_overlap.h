//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDAKERNEL_OVERLAP
#define H_CUDAKERNEL_OVERLAP

#include <math.h>
#include <stdio.h>
#include "cudaUtil.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

/// @par Main functionality
/// check for overlap
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is X->Y->Z and slow dimension is natoms
__global__ void cudaKernel_overlap(const double * const gpu_coor, const double * const gpu_radius, const int natoms, const int npairs, const double scale, int * const gpu_noverlap, const int * const gpu_bond_table);

#endif
