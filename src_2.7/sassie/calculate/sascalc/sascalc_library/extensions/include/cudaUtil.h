//
// Hailiang Zhang
// NIST & UTK
//

// code guard
#ifndef H_CUDA_UTIL
#define H_CUDA_UTIL

// macros
#include "cudaCommon.h"

#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <iomanip>

//////////////////////////////////////////
// Error handlings
//////////////////////////////////////////
#define CALL_CUDA( err ) \
{ \
    cudaError_t cudaErr = err; \
    if (cudaErr != cudaSuccess) \
    { printf("cuda Error: \"%s\" in %s at line %d\n", cudaGetErrorString(cudaErr), __FILE__, __LINE__); exit(EXIT_FAILURE); } \
}


#define CALL_CUBLAS( err ) \
{ \
  cublasStatus_t cublasErr = err; \
  if (cublasErr != CUBLAS_STATUS_SUCCESS) \
  { printf("cublas Error code %d in %s at line %d\n", cublasErr, __FILE__, __LINE__); exit(EXIT_FAILURE); } \
}


#define CALL_CURAND( err ) \
{ \
  curandStatus_t curandErr = err; \
  if (curandErr != CURAND_STATUS_SUCCESS) \
  { printf("curand Error code %d in %s at line %d\n", curandErr, __FILE__, __LINE__); exit(EXIT_FAILURE); } \
}


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}

#endif
