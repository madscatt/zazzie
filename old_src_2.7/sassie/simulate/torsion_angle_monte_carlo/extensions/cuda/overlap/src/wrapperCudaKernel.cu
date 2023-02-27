//
// Hailiang Zhang
// NIST & UTK
//

#include <wrapperCudaKernel.h>
#include <cudaKernel_overlap.h>

#include <cudaUtil.h>

#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one atom pair
void
wrapper_cudaOverlap(const double * const gpu_coor, const double * const gpu_r, const int natoms, const double scale, int * const gpu_noverlap, const int * const gpu_bond_table)
{
	// set the CUDA block dimenstions
    const int npairs = (natoms*(natoms-1))/2;
    dim3 dim_grid((npairs-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(cudaKernel_overlap, cudaFuncCachePreferL1));
	// launch the kernel
    cudaKernel_overlap<<<dim_grid,dim_block>>>(gpu_coor, gpu_r, natoms, npairs, scale, gpu_noverlap, gpu_bond_table);
    // get the error
	CALL_CUDA(cudaGetLastError());
}
