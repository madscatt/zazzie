//
// Hailiang Zhang
// NIST & UTK
//

#include "wrapperCudaKernel.h"

#include "cudaKernel_GV.h"
#include "cudaKernel_Debye.h"
#include "cudaKernel_DebyePBC.h"

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
wrapper_cudaGV_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int xray)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_gv, gpu_q, gpu_Is, Natoms, Nq, Ngv, xray);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaGV_calcIq_complex(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_Is_real, TYPE * const gpu_Is_imag, const int Natoms, const int Nq, const int Ngv, const int xray)
{
	// set the CUDA block dimenstions
    dim3 dim_grid(Ngv), dim_block(Nq);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_GV_calc_complex, cudaFuncCachePreferL1));
	// launch the kernel
	sascalc::cuda::cudaKernel_GV_calc_complex<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_gv, gpu_q, gpu_Is_real, gpu_Is_imag, Natoms, Nq, Ngv, xray);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaDebye_calcIq(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int xray)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_Debye_calc, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_Debye_calc<<<dim_grid, dim_block>>>(gpu_coor, gpu_b, gpu_q, gpu_Ia, Natoms, Nq, xray);
	CALL_CUDA(cudaGetLastError());
}

/// @par Main functionality
/// wrap the cuda function by a C++ interface
/// @par CUDA threads layout
///- one thread corresponds to one sample
void
sascalc::cuda::
wrapper_cudaDebyePBC_calcIq(const TYPE boxl, const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int xray)
{
	// set the CUDA block dimenstions
    dim3 dim_grid((Natoms-1)/BLOCKDIM+1), dim_block(BLOCKDIM);
	// cudakernel launch configuration
	CALL_CUDA(cudaFuncSetCacheConfig(sascalc::cuda::cudaKernel_DebyePBC_calc, cudaFuncCachePreferShared));
	// launch the kernel
	sascalc::cuda::cudaKernel_DebyePBC_calc<<<dim_grid, dim_block>>>(boxl, gpu_coor, gpu_b, gpu_q, gpu_Ia, Natoms, Nq, xray);
	CALL_CUDA(cudaGetLastError());
}

