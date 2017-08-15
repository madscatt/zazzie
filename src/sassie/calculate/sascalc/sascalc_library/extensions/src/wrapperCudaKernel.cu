//
// Hailiang Zhang
// NIST & UTK
//

#include "wrapperCudaKernel.h"

#include "cudaKernel_GV.h"

#include "cudaUtil.h"

#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq(const double * const gpu_coor, const double * const gpu_B_neutron, const double * const gpu_B_xray, const int *const gpu_flag_skip, const double * const gpu_gv, const double * const gpu_q, double * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int Ncontrast_neutron, const int Ncontrast_xray)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_B_neutron, gpu_B_xray, gpu_flag_skip, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, Ngv, Ncontrast_neutron, Ncontrast_xray);
	CALL_CUDA(cudaGetLastError());
}
